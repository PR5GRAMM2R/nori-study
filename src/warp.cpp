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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

float tentTransform(float x) {
    if (x < 0.5f)
        return std::sqrt(2.0f * x) - 1.0f;
    else
        return 1.0f - std::sqrt(2.0f - 2.0f * x);
}

Point2f Warp::squareToTent(const Point2f &sample) {
    //throw NoriException("Warp::squareToTent() is not yet implemented!");
    return Point2f(
        tentTransform(sample.x()),
        tentTransform(sample.y())
    );
}

float Warp::squareToTentPdf(const Point2f &p) {
    //throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
    if (p.x() < -1 || p.x() > 1 || p.y() < -1 || p.y() > 1)
        return 0.0f;

    return (1.0f - std::abs(p.x())) * (1.0f - std::abs(p.y()));
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    //throw NoriException("Warp::squareToUniformDisk() is not yet implemented!");
    if (std::pow(sample.x(), 2.0) + std::pow(sample.y(), 2.0) <= 1.0)               ////////////////////////
        return sample;
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    //throw NoriException("Warp::squareToUniformDiskPdf() is not yet implemented!");
    return (std::pow(p.x(), 2.0) + std::pow(p.y(), 2.0) <= 1.0) ? 1.0f : 0.0f;      ////////////////////////
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    throw NoriException("Warp::squareToUniformSphere() is not yet implemented!");
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    throw NoriException("Warp::squareToUniformHemisphere() is not yet implemented!");
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    throw NoriException("Warp::squareToUniformHemispherePdf() is not yet implemented!");
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    throw NoriException("Warp::squareToCosineHemisphere() is not yet implemented!");
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    throw NoriException("Warp::squareToCosineHemispherePdf() is not yet implemented!");
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END
