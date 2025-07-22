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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

#include <nori/dpdf.h>
#include <random>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        Vector3f iRay = bRec.wi.normalized();
        Vector3f oRay = bRec.wo.normalized();
        Vector3f hRay = (iRay + oRay).normalized();
        Vector3f n = Vector3f(0, 0, 1);
        float inN = m_intIOR;
        float outN = m_extIOR;
        float alpha = m_alpha;

        Color3f kd = m_kd;
        float ks = m_ks;

        float thetaH = std::acos(hRay.dot(n));
        float D = (0.5 * INV_PI) * (2 * std::exp(-std::pow(std::tan(thetaH), 2)) / std::pow(alpha, 2)) 
            / (std::pow(alpha, 2) * std::pow(std::cos(thetaH), 3));

        float F = fresnel(iRay.dot(hRay), outN, inN);

        float thetaI = std::acos(iRay.dot(n));
        float thetaO = std::acos(oRay.dot(n));

        float bi = 1.0 / (alpha * std::tan(thetaI));
        float ci = iRay.dot(hRay) / iRay.dot(n);
        float Xi = (ci > 0) ? 1 : 0;
        float tempI = (bi < 1.6) ? ((3.535 * bi + 2.181 * std::pow(bi, 2)) / (1 + 2.276 * bi + 2.577 * std::pow(bi, 2))) : 1;
        float Gi = Xi * tempI;

        float bo = 1.0 / (alpha * std::tan(thetaO));
        float co = oRay.dot(hRay) / oRay.dot(n);
        float Xo = (co > 0) ? 1 : 0;
        float tempO = (bo < 1.6) ? ((3.535 * bo + 2.181 * std::pow(bo, 2)) / (1 + 2.276 * bo + 2.577 * std::pow(bo, 2))) : 1;
        float Go = Xo * tempO;

        float G = Gi * Go;

        Color3f f = (kd * INV_PI) + ks * (D * F * G / (4 * std::cos(thetaI) * std::cos(thetaO) * std::cos(thetaH))) * Color3f(1.0f);

        return f;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        Vector3f iRay = bRec.wi.normalized();
        Vector3f oRay = bRec.wo.normalized();
        Vector3f hRay = (iRay + oRay).normalized();
        //Vector3f n = Vector3f(0, 0, 1);
        float alpha = m_alpha;
        float ks = m_ks;

        float cosThetaO = Frame::cosTheta(oRay);

        if (cosThetaO <= 0.f)
            return 0.f;

        //float thetaH = std::acos(hRay.dot(n));
        float D = (0.5 * INV_PI) * (2 * std::exp(-std::pow(Frame::tanTheta(hRay), 2) / std::pow(alpha, 2)))
            / std::fmax(1e-6f, (std::pow(alpha, 2) * std::pow(Frame::cosTheta(hRay), 3)));

        float J = 1.0 / (4 * std::fmax(1e-6f, (hRay.dot(oRay))));

        float _pdf = ks * D * J + (1.0 - ks) * (cosThetaO / M_PI);

        return _pdf;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        float alpha = m_alpha;
        bool isSpecular = true;

        Point2f __sample;

        // Sample Reusing ÇÊ¿äÇÔ
        if (_sample[0] >= m_ks) {
            isSpecular = false;
            __sample.x() = (_sample.x() - m_ks) / (1 - m_ks);
        }
        else {
            __sample.x() = _sample.x() / m_ks;
        }
        __sample.y() = _sample.y();

        if (!isSpecular) {          // Diffuse
            if (Frame::cosTheta(bRec.wi) <= 0)
                return Color3f(0.0f);

            bRec.measure = ESolidAngle;
            bRec.wo = Warp::squareToCosineHemisphere(__sample);
            bRec.eta = 1.0f;
            return m_kd;
        }
        else {                      // Specular
            Vector3f n = Warp::squareToBeckmann(__sample, alpha);

            reflection(-(bRec.wi), n, bRec.wo);

            bRec.measure = ESolidAngle;
            bRec.eta = 1.f;

            float _pdf = pdf(bRec);
            if (_pdf <= 0)   return Color3f(0.f);

            return eval(bRec) * bRec.wo.dot(n) / pdf(bRec);
        }

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
