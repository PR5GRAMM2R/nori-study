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

#pragma once

#include <nori/mesh.h>
#include <nori/bbox.h>

#define OCTREE_CHILDS 8
#define OCTREE_NODE_MAX_CAPACITY 10
#define OCTREE_MAX_DEPTH_NODE_MAX_CAPACITY 1000
#define OCTREE_MAX_DEPTH 6

NORI_NAMESPACE_BEGIN

class Node {
public:
    Node();
    ~Node();

    BoundingBox3f box = BoundingBox3f();
    BoundingBox3f tempBox = BoundingBox3f();
    //std::vector<uint32_t> triangleIdxs;
    uint32_t* triangleIdxs = new uint32_t[OCTREE_NODE_MAX_CAPACITY + 1];
    Node* parent = nullptr;
    Node* child[8] = { nullptr, };
    bool hasChild = false;
    int depth;
    int localSize = 0;
    int totalSize = 0;

    int visitIdx;
};

class Octree {
public:
    Node* rootNode = new Node;
    uint32_t nodeCount = 0;
    uint32_t trianglesCount = 0;

    Octree();
    ~Octree();

    bool buildOctree(Mesh* mesh);
    void scanNodesOctree(Node* node);
    void scanTrianglesOctree(Node* node);
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    Octree octree;
};

NORI_NAMESPACE_END
