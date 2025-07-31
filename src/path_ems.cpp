
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

            if (its.mesh->isEmitter() && bounce == 0) {
                L += throughput * its.mesh->getEmitter()->getRadiance();

                break;
            }
            else if (its.mesh->isEmitter() && bounce > 0) {
                L += throughput * 0.5 * its.mesh->getEmitter()->getRadiance();

                break;
            }

            const BSDF* bsdf = its.mesh->getBSDF();

            if (!bsdf)
                break;

            /////////////////////

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

            BSDFQueryRecord bRecNEE(its.toLocal(-dirFromLight), its.toLocal(-currentRay.d), ESolidAngle);
            Point2f u1 = sampler->next2D();

            Color3f fNEE = its.mesh->getBSDF()->eval(bRecNEE); //sample(bRecNEE, u1);

            /////////////////////

            BSDFQueryRecord bRecBRDF(its.toLocal(-currentRay.d));
            Point2f u2 = sampler->next2D();

            Color3f fBRDF = bsdf->sample(bRecBRDF, u2);
            if (fBRDF.isZero())
                break;

            Vector3f wo_world = its.toWorld(bRecBRDF.wo);
            currentRay = Ray3f(its.p, wo_world);

            Intersection tempIts;

            if (lightPD != 0) {
                L += throughput * 0.5 * fNEE * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
            }

            throughput *= fBRDF;

            eta_prod *= bRecBRDF.eta;

            if (bounce >= 3) {
                float p = std::min(throughput.maxCoeff() * eta_prod * eta_prod, 0.99f);

                if (sampler->next1D() > p)
                    break;

                throughput /= p;
            }

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