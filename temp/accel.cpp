﻿/*
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
#include <chrono>

NORI_NAMESPACE_BEGIN

//BoundingBox3f getFullMeshBoundingBox(Mesh* m_mesh)
//{
//    BoundingBox3f meshBoundingBox = m_mesh->getBoundingBox(0);;
//
//    for (uint32_t idx = 1; idx < m_mesh->getTriangleCount(); ++idx) {
//        meshBoundingBox.expandBy(m_mesh->getBoundingBox(idx));
//    }
//
//    return meshBoundingBox;
//}

Node::Node() {

}

Node::~Node() {
    if(hasChild)
        for(int i = 0; i < OCTREE_CHILDS; i++)
            delete child[i];
}

Octree::Octree()
{

}

Octree::~Octree()
{
    delete rootNode;
}

bool octreeInsert(Mesh* m_mesh, Node* node, uint32_t index)
{
    assert(node != nullptr);
    assert(node->depth <= OCTREE_MAX_DEPTH);
    //BoundingBox3f idxCenterPointBox(m_mesh->getCentroid(index));
    //assert(node->tempBox.contains(idxCenterPointBox, false) == true);

    if (node->hasChild) {
        BoundingBox3f idxBox = m_mesh->getBoundingBox(index);

        node->box = node->tempBox;

        for (int i = 0; i < OCTREE_CHILDS; i++) {
			if (node->child[i]->tempBox.overlaps(idxBox, false)) {
				octreeInsert(m_mesh, node->child[i], index);
                node->box.expandBy(node->child[i]->box);
			}
		}
	}
    else {
        node->triangleIdxs[node->localSize] = index;
        node->localSize++;
        node->box.expandBy(m_mesh->getBoundingBox(index));

        assert(node->localSize != OCTREE_MAX_DEPTH_NODE_MAX_CAPACITY);

        Node* nodeTemp = node;
        for (int depth = node->depth; depth > 0; depth--) {
            nodeTemp->totalSize++;
            nodeTemp = nodeTemp->parent;
        }
        nodeTemp->totalSize++;

        if (node->depth == OCTREE_MAX_DEPTH) {
            return true;
        }
        else if (node->localSize > OCTREE_NODE_MAX_CAPACITY) {
            node->box.reset();

            Point3f boxCenter = node->tempBox.getCenter();
            std::vector<Point3f> boxCorners = node->tempBox.get3DimCorners();

            for (int i = 0; i < OCTREE_CHILDS; i++) {
                node->child[i] = new Node;
                node->child[i]->depth = node->depth + 1;
                node->child[i]->tempBox = BoundingBox3f(boxCorners[i], boxCenter, true);
                node->child[i]->parent = node;
                if (node->depth == OCTREE_MAX_DEPTH - 1) {
                    delete[] node->child[i]->triangleIdxs;
                    node->child[i]->triangleIdxs = new uint32_t[OCTREE_MAX_DEPTH_NODE_MAX_CAPACITY];
                }
                else
                    node->child[i]->triangleIdxs = new uint32_t[OCTREE_NODE_MAX_CAPACITY + 1];
            }

            for (int i = 0; i < node->localSize; i++) {
                uint32_t idx = node->triangleIdxs[i];
                BoundingBox3f idxBox = m_mesh->getBoundingBox(idx);

                for (int j = 0; j < OCTREE_CHILDS; j++) {
                    if (node->child[j]->tempBox.overlaps(idxBox, false)) {
                        octreeInsert(m_mesh, node->child[j], idx);
                        node->box.expandBy(node->child[j]->box);
                    }
                }
            }

            Node* nodeTemp = node;
            for (int depth = node->depth; depth > 0; depth--) {
                nodeTemp->totalSize -= node->localSize;
                nodeTemp = nodeTemp->parent;
            }
            nodeTemp->totalSize -= node->localSize;

            node->hasChild = true;
            node->localSize = 0;
        }
    }

    return true;
}

void makeActualOctreeBox(Mesh* mesh, Node* node)
{
    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            makeActualOctreeBox(mesh, node->child[i]);
            node->box.expandBy(node->child[i]->box);
        }
    }
    else {
        if (node->localSize == 0) {
            node->box = node->tempBox;
            return;
        }

        node->box = mesh->getBoundingBox(0);
        for (uint32_t idx = 1; idx < node->localSize; idx++) {
            BoundingBox3f temp = mesh->getBoundingBox(idx);
            node->box.expandBy(temp);
        }
    }

    /*for (uint32_t idx = 0; idx < node->localSize; ++idx) {
        BoundingBox3f temp = mesh->getBoundingBox(idx);
        assert(node->box.contains(temp, false) == true);
    }*/

    return;
}

bool Octree::buildOctree(Mesh* mesh)
{
    BoundingBox3f meshBoundingBox = mesh->getBoundingBox();

    rootNode = new Node;
    rootNode->tempBox = meshBoundingBox;
    rootNode->hasChild = false;
    rootNode->depth = 0;

    for (uint32_t idx = 0; idx < mesh->getTriangleCount(); ++idx) {
        octreeInsert(mesh, rootNode, idx);
    }

    //makeActualOctreeBox(mesh, rootNode);

    if (rootNode->box == meshBoundingBox)
        return true;

    return false;
}

void Octree::scanNodesOctree(Node* node)
{
    //static uint32_t nodeCount = 0;

    nodeCount++;

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            scanNodesOctree(node->child[i]);
        }
    }

    return;
}

void Octree::scanTrianglesOctree(Node* node)
{
    //static uint32_t trianglesCount = 0;

    trianglesCount += node->localSize;

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            scanTrianglesOctree(node->child[i]);
        }
    }

    return;
}

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();

    cout << "Loading mesh Completed." << endl;

    auto start = std::chrono::high_resolution_clock::now();

    bool isOctreeNormal = octree.buildOctree(m_mesh);
    assert(isOctreeNormal == true);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Making Octree took " << duration.count() << "ms." << std::endl;

    std::cout << "root box min: " << octree.rootNode->box.min.transpose()
        << "  max: " << octree.rootNode->box.max.transpose() << std::endl;

    std::cout << "root centerBox min: " << octree.rootNode->tempBox.min.transpose()
        << "  max: " << octree.rootNode->tempBox.max.transpose() << std::endl;

    octree.scanNodesOctree(octree.rootNode);
    cout << "Octree's node Count : " << octree.nodeCount << endl;// scanNodesOctree(octree) << endl;
    octree.scanTrianglesOctree(octree.rootNode);
    cout << "Octree's triangle Count : " << octree.trianglesCount << endl;// scanTrianglesOctree(octree) << endl;
}

void Accel::build() {
    /* Nothing to do here for now */
}

void findIntersectedBoxes
(Node* node, const Ray3f& ray_, std::vector<std::pair<Node*, float>>& foundBoxes)
{
    float nearT, farT;

    if (node->hasChild) {
        for (int i = 0; i < OCTREE_CHILDS; i++) {
            if (node->child[i]->box.rayIntersect(ray_))
                findIntersectedBoxes(node->child[i], ray_, foundBoxes);
        }
    }
    else {
        if (node->localSize != 0) {
            node->box.rayIntersect(ray_, nearT, farT);
            foundBoxes.push_back(std::pair<Node*, float>(node, farT));
        }
        return;
    }

    return;
}

//void findIntersectedBoxes
//    (Node* node, const Ray3f& ray_, std::vector<std::pair<Node*, float>>& foundBoxes, float nearT, float farT)
//{
//    if (node->hasChild) {
//        for (int i = 0; i < OCTREE_CHILDS; i++) {
//            if (node->child[i]->box.rayIntersect(ray_, nearT, farT))
//                findIntersectedBoxes(node->child[i], ray_, foundBoxes, nearT, farT);
//        }
//    }
//    else {
//        if (node->localSize != 0) {
//            foundBoxes.emplace_back(std::pair<Node*, float>(node, farT));
//        }
//        return;
//    }
//
//    return;
//}

//void findIntersectedBoxes
//    (Node* node, const Ray3f& ray_, std::pair<Node*, float>*& foundBoxes, int& foundBoxesSize, float nearT, float farT)
//{
//    if (foundBoxesSize > 100)
//        return;
//
//    if (node->hasChild) {
//        for (int i = 0; i < OCTREE_CHILDS; i++) {
//            if (node->child[i]->box.rayIntersect(ray_, nearT, farT))
//                findIntersectedBoxes(node->child[i], ray_, foundBoxes, foundBoxesSize, nearT, farT);
//        }
//    }
//    else {
//        if (node->localSize != 0) {
//            foundBoxes[foundBoxesSize] = std::pair<Node*, float>(node, farT);
//            foundBoxesSize++;
//        }
//        return;
//    }
//
//    return;
//}


void sortFoundBoxes(std::vector<std::pair<Node*, float>>& boxes)
{
    std::sort(boxes.begin(), boxes.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    return;
}

void sortFoundTriangles(std::vector<std::pair<uint32_t, float>>& triangles)
{
    std::sort(triangles.begin(), triangles.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    return;
}

//bool findIntersectedTriangles
//    (std::vector<std::pair<Node*, float>>& foundBoxes, const Ray3f& ray_, Mesh *mesh, 
//        Intersection &its, bool shadowRay)
//{
//    for (int i = 0; i < foundBoxes.size(); i++) {
//
//    }
//}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    ///* Brute force search through all triangles */
    //for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //    float u, v, t;
    //    if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //        /* An intersection was found! Can terminate
    //           immediately if this is a shadow ray query */
    //        if (shadowRay)
    //            return true;
    //        ray.maxt = its.t = t;
    //        its.uv = Point2f(u, v);
    //        its.mesh = m_mesh;
    //        f = idx;
    //        foundIntersection = true;
    //    }
    //}

    {
        //auto start = std::chrono::high_resolution_clock::now();                                    ////////////////////

        std::vector<std::pair<Node*, float>> foundBoxes;
        findIntersectedBoxes(octree.rootNode, ray, foundBoxes);
        float nearT, farT;
        //findIntersectedBoxes(octree.rootNode, ray, foundBoxes, nearT, farT);

        /*std::pair<Node*, float>* foundBoxes = new std::pair<Node*, float>[100];
        int foundBoxesSize = 0;
        findIntersectedBoxes(octree.rootNode, ray, foundBoxes, foundBoxesSize, nearT, farT);*/

        //auto end = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        //if (ray.o.x() - -65.605499 > 1e-3f)
            //printf("11111 findIntersectedBoxes %f => %lld us.\n", ray.o.x(), duration.count());

        //start = std::chrono::high_resolution_clock::now();                                    ////////////////////

        if (foundBoxes.size() > 0)
            sortFoundBoxes(foundBoxes);
        else
            return foundIntersection;

        //if (foundBoxesSize > 0) {
        //    std::sort(foundBoxes,                       // 시작 주소
        //        foundBoxes + foundBoxesSize,      // 끝   주소(한 칸 뒤)
        //        [](const std::pair<Node*, float>& a,
        //            const std::pair<Node*, float>& b) {
        //                return a.second < b.second; // 오름차순(nearT 가 작은 것부터)
        //        });
        //}
        //else
        //    return foundIntersection;

        //end = std::chrono::high_resolution_clock::now();
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        //if (ray.o.x() - -65.605499 > 1e-6f)
            //printf("22222 sortFoundBoxes %f => %lld us.\n", ray.o.x(), duration.count());

        uint32_t foundTriangle = -1;
        float foundTriangleDistance = INFINITY;

        //std::pair<Node*, float>* foundBox;
        Node* nodeTemp;

        //start = std::chrono::high_resolution_clock::now();                                    ////////////////////

        for (int i = 0; i < foundBoxes.size(); i++) {
        //for (int i = 0; i < foundBoxesSize; i++) {
            //foundBox = &foundBoxes[i];
            nodeTemp = foundBoxes[i].first;
            float u, v, t;

            uint32_t foundTriangleTemp;
            float foundTriangleDistanceTemp = foundTriangleDistance;

            for (uint32_t j = 0; j < nodeTemp->localSize; j++) {
                uint32_t idx = nodeTemp->triangleIdxs[j];

                if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
                    if (shadowRay)
                        return true;

                    if (t < foundTriangleDistance) {
                        foundTriangleTemp = idx;
                        foundTriangleDistanceTemp = t;
                        foundIntersection = true;
                    }
                }
            }

            if (foundIntersection) {
                if (foundTriangleDistanceTemp < foundTriangleDistance) {
                    foundTriangle = foundTriangleTemp;
                    foundTriangleDistance = foundTriangleDistanceTemp;
                }
                else {
                    break;
                }
            }
        }

        if (foundIntersection) {
            float u, v, t;
            uint32_t idx = foundTriangle;

            m_mesh->rayIntersect(idx, ray, u, v, t);

            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;

            //end = std::chrono::high_resolution_clock::now();
            //duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            //if (ray.o.x() - -65.605499 > 1e-6f)
                //printf("33333 foundTriangle %f => %lld us.\n", ray.o.x(), duration.count());
        }
    }

    //auto start = std::chrono::high_resolution_clock::now();                                    ////////////////////

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

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //if(ray.o.x() - -65.605499 > 1e-6f)
        //printf("44444 foundIntersection %f => %lld us.\n", ray.o.x(), duration.count());

    return foundIntersection;
}

NORI_NAMESPACE_END

