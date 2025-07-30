
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <pcg32.h>
#include <chrono>

NORI_NAMESPACE_BEGIN

#define MAX_SAMPLES 10

class PATHMATSIntegrator : public Integrator {
public:
    PATHMATSIntegrator(const PropertyList& props) {
    }

    // Recursive
    /*
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        float continueProb = 0.99;
        Color3f throughput = Color3f(1.f);

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        Point3f p = its.p;
        Color3f result;

        const BSDF* bsdf = its.mesh->getBSDF();

        if (n.dot(-ray.d) <= 0 && bsdf->isDiffuse()) {
            return Color3f(0.0f);
        }

        Vector2f sample = sampler->next2D();
        BSDFQueryRecord bsdfQueryRecord(its.toLocal(-ray.d));

        Color3f bsdfSample;
        Vector3f oRay;

        if (its.mesh->isEmitter())
            return its.mesh->getEmitter()->getRadiance();

        bsdfSample = bsdf->sample(bsdfQueryRecord, sample);
        oRay = its.toWorld(bsdfQueryRecord.wo);

        //throughput *= bsdfSample * std::abs(Frame::cosTheta(-ray.d)) / bsdf->pdf(bsdfQueryRecord);

        //continueProb = std::min(throughput.maxComponent() * bsdfQueryRecord.eta * bsdfQueryRecord.eta, 0.99f);

        if (sampler->next1D() < continueProb) {
            result = (1.0 / continueProb) * bsdfSample * Li(scene, sampler, Ray3f(p + 1e-4f * oRay, oRay));
        }
        else {
            result = Color3f(0.0);
        }

        return result;
    }
    */

    // Loop
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f L(0.0f);
        Color3f throughput = Color3f(1.f);
        float eta_prod = 1.0f;

        Ray3f currentRay = ray;
        int bounce = 0;

        while (true) {
            Intersection its;

            bool hit = scene->rayIntersect(currentRay, its);

            if (!hit || (its.mesh->getBSDF()->isDiffuse() && its.shFrame.n.dot(-currentRay.d) <= 0))
                break;

            if (its.mesh->isEmitter()) {
                L += throughput * its.mesh->getEmitter()->getRadiance();
            }

            const BSDF* bsdf = its.mesh->getBSDF();

            if (!bsdf)
                break;

            BSDFQueryRecord bRec(its.toLocal(-currentRay.d));
            Point2f u = sampler->next2D();

            float pdf = 0.0f;
            Color3f f = bsdf->sample(bRec, u);
            if (f.isZero())
                break;

            pdf = bsdf->pdf(bRec);

            if (pdf < 0.0f || f.isZero())
                break;

            throughput *= f;

            eta_prod *= bRec.eta;

            if (bounce >= 3) {
                float p = std::min(throughput.maxCoeff() * eta_prod * eta_prod, 0.99f);

                if (sampler->next1D() > p)
                    break;

                throughput /= p;
            }

            Vector3f wo_world = its.toWorld(bRec.wo);
            currentRay = Ray3f(its.p, wo_world);

            bounce++;
        }

        return L;
    }

    std::string toString() const {
        return "PATHMATSIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(PATHMATSIntegrator, "path_mats");
NORI_NAMESPACE_END