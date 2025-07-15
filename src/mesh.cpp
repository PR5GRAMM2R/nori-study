/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { }

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

//void Mesh::activate() {
//    if (!m_bsdf) {
//        /* If no material was assigned, instantiate a diffuse BRDF */
//        m_bsdf = static_cast<BSDF *>(
//            NoriObjectFactory::createInstance("diffuse", PropertyList()));
//    }
//}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF*>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
    m_pdf = DiscretePDF(getTriangleCount());
    for (int i = 0; i < getTriangleCount(); i++) {
        m_pdf.append(surfaceArea(i));
    }
    m_pdf.normalize();
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}


// [ Sampling ]
// 
// 1. The sampled position p on the surface of the mesh.
// 2. The interpolated surface normal n at p computed from the per-vertex normals.
//      When the mesh does not provide per-vertex normals, compute and return the face normal instead.
// 3. The probability density of the sample.
//      This should be the reciprocal of the surface area of the entire mesh.
//
// Once a triangle is chosen, you can (uniformly) sample a barycentric coordinate
//  (alpha, beta, 1 - alpha - beta) using the mapping :
//      alpha = 1 - sqrt(1 - xi1)  
//      beta = xi2 * sqrt(1 - xi1)
//  where xi1 and xi2 are uniform variates in [0, 1).
float Mesh::samplePosition(Sampler* sampler, Point3f& position, Normal3f& normal) {
    Point2f s = sampler->next2D();
    uint32_t index = m_pdf.sample(s.x());

    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    float r = std::sqrt(1 - s.x());
    float a = 1 - r;
    float b = s.y() * r;

    normal = getFaceNormal(index, Point2f(a, b));
    float pd = 1.0f / allSurfaceArea();

    position = p0 * a + p1 * b + p2 * (1 - a - b);

    return pd;
    //return (i1 - i0) * sampler->next1D() + (i2 - i0) * sampler->next1D();
}

//float Mesh::samplePosition(const Point2f& sample, Point3f& position, Normal3f& normal)
//{
//    uint32_t index = m_pdf.sample(sample.x());
//
//    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
//    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
//
//    float r = std::sqrt(1 - sample.x());
//    float a = 1 - r;
//    float b = sample.y() * r;
//
//    normal = getFaceNormal(index, sample);
//    float pd = 1.0f / allSurfaceArea();
//}

Normal3f Mesh::getFaceNormal(uint32_t index, Point2f& i) {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    Normal3f normal;

    Point3f p0, p1, p2;
    if (m_N.size() == 0) {
        p0 = m_V.col(i0); p1 = m_V.col(i1); p2 = m_V.col(i2);
        normal = (p1 - p0).cross(p2 - p0);
        normal.normalize();
        return normal;
    }

    p0 = m_N.col(i0), p1 = m_N.col(i1), p2 = m_N.col(i2);

    normal = i.x() * p0 + i.y() * p1 + (1 - i.x() - i.y()) * p2;
    //normal.normalize();
    return normal;
}

float Mesh::allSurfaceArea() const
{
    float area = 0;

    for (int idx = 0; idx < getTriangleCount(); idx++)
        area += surfaceArea(idx);

    return area;
}

NORI_NAMESPACE_END
