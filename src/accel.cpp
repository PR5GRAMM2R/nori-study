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

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

//void Accel::addMesh(Mesh *mesh) {
//    if (m_mesh)
//        throw NoriException("Accel: only a single mesh is supported!");
//    m_mesh = mesh;
//    m_bbox = m_mesh->getBoundingBox();
//}

void Accel::build() {
    /* Nothing to do here for now */
}

BoundingBox3f getFullMeshBoundingBox(Mesh *m_mesh)
{
    BoundingBox3f meshBoundingBox;

    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx){
        BoundingBox3f temp = m_mesh->getBoundingBox(idx);
        meshBoundingBox.expandBy(temp);
    }

    return meshBoundingBox;
}

#define OCTREE_CHILDS 8
#define OCTREE_NODE_CAPACITY 10
#define OCTREE_MAX_DEPTH 7

Node::Node(){

}

Node::~Node(){
    for(int i = 0; i < OCTREE_CHILDS; i++){
        delete child[i];
    }

    delete this;
}

bool addNode(Mesh *m_mesh, Node* node, uint32_t index){
    if(node->depth > OCTREE_MAX_DEPTH)
         return false;

    node->totalSize++;

    if(node->depth != OCTREE_MAX_DEPTH && node->triangleIdxs.size() < OCTREE_NODE_CAPACITY && !node->hasChild){
        node->triangleIdxs.push_back(index);
        node->localSize++;

        return true;
    }
    else if (node->depth == OCTREE_MAX_DEPTH && !node->hasChild) {
        node->triangleIdxs.push_back(index);
        node->localSize++;
        
        return true;
    }
    else if(node->hasChild){
        Point3f idxCenter = m_mesh->getCentroid(index);

        for(int i = 0; i < OCTREE_CHILDS; i++){
            if(node->child[i]->box.contains(idxCenter, false)){
                addNode(m_mesh, node->child[i], index);

                break;
            }
        }

        return true;
    }
    else if(node->triangleIdxs.size() >= OCTREE_NODE_CAPACITY && !node->hasChild) {
        node->hasChild = true;
        node->localSize = 0;

        Node* newNodes = new Node[OCTREE_CHILDS];
        Point3f boxCenter = node->box.getCenter();
        std::vector<Point3f> boxCorners = node->box.get3DimCorners();

        /*cout << boxCenter.x() << ", " << boxCenter.y() << ", " << boxCenter.z() << endl << endl;
        for (int i = 0; i < 8; i++) {
            cout << boxCorners[i].x() << ", " << boxCorners[i].y() << ", " << boxCorners[i].z() << endl;
        }
        cout << endl;*/

        for(int i = 0; i < OCTREE_CHILDS; i++){
            newNodes[i].depth = node->depth + 1;
            newNodes[i].hasChild = false;
            newNodes[i].box = BoundingBox3f(boxCorners[i], boxCenter, true);
            node->child[i] = &newNodes[i];
        }

        node->triangleIdxs.push_back(index);
        for (int i = 0; i < node->triangleIdxs.size(); i++) {
            uint32_t idx = node->triangleIdxs[i];
            Point3f idxCenter = m_mesh->getCentroid(idx);

            for(int j = 0; j < OCTREE_CHILDS; j++){
                if(node->child[j]->box.contains(idxCenter, false)){
                    addNode(m_mesh, node->child[j], index);

                    break;
                }
            }
        }

        node->triangleIdxs.clear();

        return true;
    }
}

Node* buildOctree(Mesh* m_mesh)
{
    BoundingBox3f meshBoundingBox = getFullMeshBoundingBox(m_mesh);

    Node* octree = new Node;
    octree->box = meshBoundingBox;
    octree->hasChild = false;
    octree->depth = 0;
    
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx){
        addNode(m_mesh, octree, idx);
    }

    return octree;
}

void printOctree(Node* node)
{
    printf("%d Depth -> %lu Counts\n", node->depth, node->totalSize);

    for (int i = 0; i < node->depth; i++) {
        cout << "=====||";
    }
    cout << endl;

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            printOctree(node->child[i]);
        }
    }

    return;
}

//int temp = 0;
bool aaa = false;

void Accel::addMesh(Mesh* mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();

    {
        if (!aaa) {
            Node* octree = buildOctree(m_mesh);
            //printOctree(octree);
            aaa = true;
        }
    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    /*{
        if(!aaa){
            Node* octree = buildOctree(m_mesh);
            printOctree(octree);
            aaa = true;
        }
    }*/
    
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    //printf("RayIntersect %d, ", ++temp);
    
    /* Brute force search through all triangles */
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        float u, v, t;
        if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
            /* An intersection was found! Can terminate
               immediately if this is a shadow ray query */
            if (shadowRay)
                return true;
            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;
            foundIntersection = true;
        }
    }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

