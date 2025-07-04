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
#define OCTREE_NODE_MAX_CAPACITY 10
#define OCTREE_MAX_DEPTH 6

Node::Node(){

}

Node::~Node(){
    //for(int i = 0; i < OCTREE_CHILDS; i++){
        delete[] child;
    //}
}


bool nodeInsert(Mesh *m_mesh, Node* node, uint32_t index){
    assert(node != nullptr);
    assert(node->depth <= OCTREE_MAX_DEPTH);
    BoundingBox3f idxBox = m_mesh->getBoundingBox(index);
    //assert(node->box.contains(idxBox, false) == true);

    node->totalSize++;

    if (node->hasChild) {
        Point3f idxCenter = m_mesh->getCentroid(index);

        for (int i = 0; i < OCTREE_CHILDS; i++) {
            if (node->child[i]->boxForCenters.contains(idxCenter, false)) {
                return nodeInsert(m_mesh, node->child[i], index);
            }
        }

        cout << "ERROR" << endl;
        return false;
    }
    else if (node->triangleIdxs.size() >= OCTREE_NODE_MAX_CAPACITY && node->depth < OCTREE_MAX_DEPTH) {
        node->hasChild = true;
        node->localSize = 0;

        Node* newNodes = new Node[OCTREE_CHILDS];
        Point3f boxCenter = node->boxForCenters.getCenter();
        std::vector<Point3f> boxCorners = node->boxForCenters.get3DimCorners();

        //cout << boxCenter.x() << ", " << boxCenter.y() << ", " << boxCenter.z() << endl << endl;
        //for (int i = 0; i < 8; i++) {
        //    cout << boxCorners[i].x() << ", " << boxCorners[i].y() << ", " << boxCorners[i].z() << endl;
        //}
        //cout << endl;

        for (int i = 0; i < OCTREE_CHILDS; i++) {
            node->child[i] = &newNodes[i];
            node->child[i]->depth = node->depth + 1;
            node->child[i]->boxForCenters = BoundingBox3f(boxCorners[i], boxCenter, true);
            node->child[i]->parent = node;
        }

        node->triangleIdxs.push_back(index);

        for (int i = 0; i < node->triangleIdxs.size(); i++) {
            uint32_t idx = node->triangleIdxs[i];
            Point3f idxCenter = m_mesh->getCentroid(idx);

            for (int j = 0; j < OCTREE_CHILDS; j++) {
                if (node->child[j]->boxForCenters.contains(idxCenter, false)) {
                    nodeInsert(m_mesh, node->child[j], idx);

                    break;
                }
            }
        }

        node->triangleIdxs.clear();

        return true;
    }
    else{
        node->triangleIdxs.push_back(index);
        node->localSize++;

        return true;
    }

    cout << "ERROR" << endl;
    return false;
}

void makeActualOctreeBox(Mesh* m_mesh, Node* node)
{
    node->box.reset();

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            makeActualOctreeBox(m_mesh, node->child[i]);
            node->box.expandBy(node->child[i]->box);
        }
    }
    else {
        if (node->triangleIdxs.empty()) {
            node->box.expandBy(node->boxForCenters);
        }
        else {
            for (uint32_t idx = 0; idx < node->triangleIdxs.size(); ++idx) {
                BoundingBox3f temp = m_mesh->getBoundingBox(idx);
                node->box.expandBy(temp);
            }
        }

        /*Node* nodeTemp = node;
        for (int depth = node->depth; depth > 0; depth--) {
            nodeTemp = nodeTemp->parent;
            nodeTemp->box.expandBy(node->box);
        }*/
    }

    for (uint32_t idx = 0; idx < node->triangleIdxs.size(); ++idx) {
        BoundingBox3f temp = m_mesh->getBoundingBox(idx);
        assert(node->box.contains(temp, false) == true);
    }

    return;
}

Node* buildOctree(Mesh* m_mesh)
{
    BoundingBox3f meshBoundingBox = getFullMeshBoundingBox(m_mesh);

    Node* octree = new Node;
    octree->boxForCenters = meshBoundingBox;
    octree->hasChild = false;
    octree->depth = 0;
    
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx){
        nodeInsert(m_mesh, octree, idx);
    }

    return octree;
}

uint32_t nodeCount = 0;
uint32_t trianglesCount = 0;

uint32_t scanNodesOctree(Node* node)
{
    //static uint32_t nodeCount = 0;

    nodeCount++;

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            scanNodesOctree(node->child[i]);
        }
    }

    return nodeCount;
}

uint32_t scanTrianglesOctree(Node* node)
{
    //static uint32_t trianglesCount = 0;

    trianglesCount += node->localSize;

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            scanTrianglesOctree(node->child[i]);
        }
    }

    return trianglesCount;
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

void Accel::addMesh(Mesh* mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();

    {
        octree = buildOctree(m_mesh);
        makeActualOctreeBox(m_mesh, octree);

        std::cout << "root box min: " << octree->box.min.transpose()
            << "  max: " << octree->box.max.transpose() << std::endl;

        std::cout << "root centerBox min: " << octree->boxForCenters.min.transpose()
            << "  max: " << octree->boxForCenters.max.transpose() << std::endl;

        //printOctree(octree);
        scanNodesOctree(octree);
        cout << "Octree's node Count : " << nodeCount << endl;// scanNodesOctree(octree) << endl;
        scanTrianglesOctree(octree);
        cout << "Octree's triangle Count : " << trianglesCount << endl;// scanTrianglesOctree(octree) << endl;
    }
}

void sortFoundTriangles(std::vector<std::pair<float, IdxData>>& triangles)
{
    std::sort(triangles.begin(), triangles.end(),
        [](const auto& a, const auto& b) {
            return a.second.t < b.second.t;  // second(nearT) 기준 오름차순
        });

    return;
}

std::vector<std::pair<float, IdxData>> findIntersectedTriangles(Mesh* mesh, Node* node, const Ray3f& ray_)
{
    std::vector<std::pair<float, IdxData>> foundTrianglesTemp;
    float nearT, farT;

    if (node->box.rayIntersect(ray_, nearT, farT)) {
        if (node->hasChild) {
            std::vector<std::pair<float, IdxData>> temp;

            for (int i = 0; i < OCTREE_CHILDS; i++) {
                temp = findIntersectedTriangles(mesh, node->child[i], ray_);
                if (!temp.empty()) {
                    foundTrianglesTemp.insert(foundTrianglesTemp.end(), temp.begin(), temp.end());
                }
            }

            return foundTrianglesTemp;
        }
        else {
            std::vector<std::pair<float, IdxData>> foundTriangles;

            for (uint32_t idx : node->triangleIdxs) {
                float u, v, t;
                if (mesh->rayIntersect(idx, ray_, u, v, t)) {
                    IdxData i;
                    i.u = u;
                    i.v = v;
                    i.t = t;
                    foundTriangles.push_back(std::pair<float, IdxData>(idx, i));
                }
            }

            return foundTriangles;
        }
    }

    return {};
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {    
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)
    
    ///* Brute force search through all triangles */
    //{
    //    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //        float u, v, t;
    //        if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //            /* An intersection was found! Can terminate
    //               immediately if this is a shadow ray query */
    //            if (shadowRay)
    //                return true;
    //            ray.maxt = its.t = t;
    //            its.uv = Point2f(u, v);
    //            its.mesh = m_mesh;
    //            f = idx;
    //            foundIntersection = true;
    //        }
    //    }
    //}

    {
        std::vector<std::pair<float, IdxData>> foundTriangles = findIntersectedTriangles(m_mesh, octree, ray);

        sortFoundTriangles(foundTriangles);

        //cout << foundTriangles.size() << endl;
        int i = 0;

        for (int i = 0; i < foundTriangles.size(); ++i) {
            uint32_t idx = foundTriangles[i].first;

            if (shadowRay)
                return true;

            ray.maxt = its.t = foundTriangles[i].second.t;
            its.uv = Point2f(foundTriangles[i].second.u, foundTriangles[i].second.v);
            its.mesh = m_mesh;
            f = foundTriangles[i].first;
            foundIntersection = true;
        }

        /*
        for (uint32_t i = 0; i < foundTriangles.size(); ++i) {
            uint32_t idx = foundTriangles[i].first;

            if (shadowRay)
                return true;

            ray.maxt = its.t = foundTriangles[i].second.t;
            its.uv = Point2f(foundTriangles[i].second.u, foundTriangles[i].second.v);
            its.mesh = m_mesh;
            f = foundTriangles[i].first;
            foundIntersection = true;

            //float u, v, t;
            //if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
            //   //  An intersection was found! Can terminate
            //   //    immediately if this is a shadow ray query 
            //    if (shadowRay)
            //        return true;
            //    ray.maxt = its.t = t;
            //    its.uv = Point2f(u, v);
            //    its.mesh = m_mesh;
            //    f = idx;
            //    foundIntersection = true;
            //}
        }
        */
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

