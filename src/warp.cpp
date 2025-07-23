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
    return Point2f(std::sqrt(sample.x()) * std::cos(2.0 * M_PI * sample.y()),
        std::sqrt(sample.x()) * std::sin(2.0 * M_PI * sample.y()));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    //throw NoriException("Warp::squareToUniformDiskPdf() is not yet implemented!");
    if (p.x() * p.x() + p.y() * p.y() <= 1.0f)
        return 1.0f / M_PI;
    else
        return 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    //throw NoriException("Warp::squareToUniformSphere() is not yet implemented!");
    float theta = 2.0f * M_PI * sample.x();
    float z = 1.0f - 2.0f * sample.y();
    float r = std::sqrt(std::max(0.0f, 1.0f - z * z));

    return Vector3f(
        r * std::cos(theta),
        r * std::sin(theta),
        z
    );
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    //throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
    if (v.x() * v.x() + v.y() * v.y() + v.z() * v.z() <= 1.0f)
        return 1.0f / (4.0 * M_PI);
    else
        return 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    //throw NoriException("Warp::squareToUniformHemisphere() is not yet implemented!");
    float theta = 2.0f * M_PI * sample.x();
    float z = 1.0f - 2.0f * sample.y();
    float r = std::sqrt(std::max(0.0f, 1.0f - z * z));

    Vector3f uniformSphere = Vector3f(
        r * std::cos(theta),
        r * std::sin(theta),
        z
    );

    return (uniformSphere.dot(Vector3f(0, 0, 1)) > 0) ? uniformSphere : -uniformSphere;
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    //throw NoriException("Warp::squareToUniformHemispherePdf() is not yet implemented!");
    if (v.x() * v.x() + v.y() * v.y() + v.z() * v.z() <= 1.0f && v.z() >= 0)
        return 1.0f / (2.0 * M_PI);
    else
        return 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    //throw NoriException("Warp::squareToCosineHemisphere() is not yet implemented!");
    Point2f d = Warp::squareToUniformDisk(sample);

    float x = d.x();
    float y = d.y();
    float z = std::sqrt(std::max(0.0f, 1.0f - x * x - y * y));

    return Vector3f(x, y, z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    //throw NoriException("Warp::squareToCosineHemispherePdf() is not yet implemented!");
    return v.z() > 0 ? v.z() / M_PI : 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    //throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
    float theta = std::atan(alpha * std::sqrt(-std::log(1.0 - sample.y())));
    float phi = 2.0 * M_PI * sample.x();

    return Vector3f(
        std::sin(theta) * std::cos(phi),
        std::sin(theta) * std::sin(phi),
        std::cos(theta)
    );
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    //throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
    return (0.5 * INV_PI) * (2 * std::exp(-std::pow(Frame::tanTheta(m), 2) / std::pow(alpha, 2)))
        / std::fmax(1e-6f, (std::pow(alpha, 2) * std::pow(Frame::cosTheta(m), 3)));
}

NORI_NAMESPACE_END
