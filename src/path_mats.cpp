
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

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f result(0.0, 0.0, 0.0);
        Color3f throughput(1.0, 1.0, 1.0);
        float productedEta = 1.0f;
        Ray3f currentRay = ray;

        int rayBounceCount = 0;
        while (true) {
            Intersection its;
            bool hit = scene->rayIntersect(currentRay, its);

            if (!hit) {
                break;
            }

            if (its.mesh->isEmitter()) {
                result += throughput * its.mesh->getEmitter()->getRadiance();
            }

            const BSDF* bsdf = its.mesh->getBSDF();
            if (!bsdf) break;

            BSDFQueryRecord bRec(its.toLocal(-currentRay.d));
            Point2f sample = sampler->next2D();
            Color3f f = bsdf->sample(bRec, sample);
            if (f.isZero())
                break;

            float cosTheta = std::abs(Frame::cosTheta(bRec.wo));
            float pdf = bsdf->pdf(bRec);
            if (pdf == 0.f || f.isZero())
                break;
            throughput *= f * cosTheta / pdf;

            productedEta *= bRec.eta;

            if (rayBounceCount >= 3) {
                float pCont = std::min(throughput.maxComponent() * productedEta * productedEta, 0.99f);
                if (sampler->next1D() > pCont)
                    break;
                throughput /= pCont;
            }

            Vector3f wo_world = its.toWorld(bRec.wo);
            currentRay = Ray3f(its.p, wo_world);

            rayBounceCount++;
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