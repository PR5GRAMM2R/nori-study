
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <pcg32.h>
#include <chrono>

NORI_NAMESPACE_BEGIN

#define MAX_SAMPLES 10

class PATHEMSIntegrator : public Integrator {
public:
    PATHEMSIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f L(0.0f);
        Color3f throughput = Color3f(1.f);
        float eta_prod = 1.0f;

        Ray3f currentRay = ray;
        int bounce = 0;

        std::vector<Mesh*> meshes = scene->getMeshes();
        Mesh** lightMeshes = new Mesh * [meshes.size()];

        int lightCount = 0;
        for (Mesh* mesh : meshes) {
            if (mesh->isEmitter()) {
                lightMeshes[lightCount++] = mesh;
            }
        }

        while (true) {
            Intersection its;

            bool hit = scene->rayIntersect(currentRay, its);

            if (!hit || (its.mesh->getBSDF()->isDiffuse() && its.shFrame.n.dot(-currentRay.d) <= 0))
                break;

            Point3f p = its.p;
            Normal3f n = its.shFrame.n;

            if (its.mesh->isEmitter()) {
                L += throughput * its.mesh->getEmitter()->getRadiance();

                break;
            }

            if (its.mesh->getBSDF()->isDiffuse()) {
                Mesh* selectedLightMesh = lightMeshes[(int)(sampler->next1D() * lightCount)];

                // 선택된 Light Mesh 에서 point, normal, pd 를 랜덤하게 뽑아오기
                Normal3f lightNormal;
                Point3f lightPoint;
                float lightPD = selectedLightMesh->samplePosition(sampler, lightPoint, lightNormal);

                Normal3f dirFromLight = p - lightPoint;
                float dirToLightDist = dirFromLight.norm();
                dirFromLight = dirFromLight.normalized();

                Intersection lightIts;
                Ray3f rayFromIntersection(lightPoint, dirFromLight, 1e-4f, dirToLightDist - 1e-4f);
                float visible = (scene->getAccel()->rayIntersect(rayFromIntersection, lightIts, true)) ? 0.0 : 1.0;

                float geometry = visible * (std::fmax(0, n.dot(-dirFromLight)) * std::fmax(0, lightNormal.dot(dirFromLight)) / std::pow(dirToLightDist, 2));

                BSDFQueryRecord bsdfQueryRecord(its.toLocal(-dirFromLight), its.toLocal(-ray.d), ESolidAngle);
                Color3f bsdf = its.mesh->getBSDF()->eval(bsdfQueryRecord);

                if (lightPD != 0) {
                    L += throughput * bsdf * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
                }
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
        return "PATHEMSIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(PATHEMSIntegrator, "path_ems");
NORI_NAMESPACE_END