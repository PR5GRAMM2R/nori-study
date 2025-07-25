
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

    std::string toString() const {
        return "PATHMATSIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(PATHMATSIntegrator, "path_mats");
NORI_NAMESPACE_END