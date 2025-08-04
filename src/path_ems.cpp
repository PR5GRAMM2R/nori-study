
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

    /*
    // Recursive
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        float continueProb = 0.95;

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f normal = its.shFrame.n;
        Point3f point = its.p;
        Color3f result(0.0f);

        static bool firstBounce = true;

        const BSDF* bsdf = its.mesh->getBSDF();

        BSDFQueryRecord bRecBRDF(its.toLocal(-ray.d));
        bsdf->sample(bRecBRDF, sampler->next2D());
        Ray3f newRay(point + 1e-4f * its.toWorld(bRecBRDF.wo), its.toWorld(bRecBRDF.wo));

        ///////////////////////////////////////////////////////////////////////////////////////////////

        if (!bsdf->isDiffuse()) {                       // Mirror or Dielectric
            if (sampler->next1D() < continueProb) {
                return (1 / continueProb) * Li(scene, sampler, newRay);
            }
            else
                return Color3f(0.0f);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////

        if (its.mesh->isEmitter() && its.shFrame.n.dot(-ray.d) > 0) {
            //if (firstBounce) {
                return its.mesh->getEmitter()->getRadiance();
            //    firstBounce = false;
            //}
            //else
            //    return Color3f(0.0f);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////

        Color3f resultNEE(0.0);

        std::vector<Mesh*> meshes = scene->getMeshes();
        Mesh** lightMeshes = new Mesh * [meshes.size()];

        int lightCount = 0;
        for (Mesh* mesh : meshes) {
            if (mesh->isEmitter()) {
                lightMeshes[lightCount++] = mesh;
            }
        }

        Mesh* selectedLightMesh = lightMeshes[int(sampler->next1D() * lightCount)];

        Normal3f lightNormal;
        Point3f lightPoint;
        float lightPD = (1.0 / (float)lightCount) * selectedLightMesh->samplePosition(sampler, lightPoint, lightNormal);

        Vector3f dirToLight = lightPoint - point;
        float dirToLightDist = dirToLight.norm();
        dirToLight = dirToLight.normalized();

        bool visible = false;
        Intersection lightIts;
        Ray3f rayToLight(point + 1e-4f * dirToLight, dirToLight);
        if (scene->rayIntersect(rayToLight, lightIts) && lightIts.mesh == selectedLightMesh) {
            visible = true;
        }

        float geometry = (std::fmax(0, normal.dot(dirToLight)) * std::fmax(0, lightNormal.dot(-dirToLight)) / std::pow(dirToLightDist, 2));

        BSDFQueryRecord bRecNEE(its.toLocal(-ray.d), its.toLocal(dirToLight), ESolidAngle);
        Color3f evalNEE = bsdf->eval(bRecNEE);
        Color3f albedoNEE = bsdf->sample(bRecNEE, sampler->next2D());

        if (visible) {
            resultNEE = geometry * selectedLightMesh->getEmitter()->getRadiance() * evalNEE / lightPD;
        }
        else {
            resultNEE = Color3f(0.0);
        }

        if (sampler->next1D() < continueProb) {
            if(visible)
                return (1 / continueProb) * (resultNEE + albedoNEE * Li(scene, sampler, newRay)) * 0.5;
            else
                return (1 / continueProb) * albedoNEE * Li(scene, sampler, newRay);
        }
        else 
            return Color3f(0.f);
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

            if (!hit)
                break;

            Normal3f normal = its.shFrame.n;
            Point3f point = its.p;

            const BSDF* bsdf = its.mesh->getBSDF();

            if (!bsdf)
                break;

            BSDFQueryRecord bRecBRDF(its.toLocal(-currentRay.d));
            Point2f u = sampler->next2D();

            float pdf = bsdf->pdf(bRecBRDF);
            Color3f f;

            if (bsdf->isDiffuse() && pdf > 0) {
                f = bsdf->sample(bRecBRDF, u) / pdf;
                if (f.isZero())
                    break;
            }
            else {
                f = bsdf->sample(bRecBRDF, u);
                if (f.isZero())
                    break;
            }

            //////////////////////

            if (its.mesh->isEmitter()) {
                L += throughput * its.mesh->getEmitter()->getRadiance();

                break;
            }

            currentRay = Ray3f(point, its.toWorld(bRecBRDF.wo));

            if (bsdf->isDiffuse()) {
                std::vector<Mesh*> meshes = scene->getMeshes();
                Mesh** lightMeshes = new Mesh * [meshes.size()];

                int lightCount = 0;
                for (Mesh* mesh : meshes) {
                    if (mesh->isEmitter()) {
                        lightMeshes[lightCount++] = mesh;
                    }
                }

                Mesh* selectedLightMesh = lightMeshes[int(sampler->next1D() * lightCount)];

                // 선택된 Light Mesh 에서 point, normal, pd 를 랜덤하게 뽑아오기
                Normal3f lightNormal;
                Point3f lightPoint;
                float lightPD = selectedLightMesh->samplePosition(sampler, lightPoint, lightNormal);

                Normal3f dirToLight = lightPoint - point;
                float dirToLightDist = dirToLight.norm();
                dirToLight = dirToLight.normalized();

                Intersection lightIts;
                Ray3f rayFromIntersection(point, dirToLight, 1e-4f, dirToLightDist + 1e-4f);
                //float visible = (selectedLightMesh->rayIntersect(rayFromIntersection)) ? 1.0 : 0.0;
                float visible = (scene->getAccel()->rayIntersect(rayFromIntersection, lightIts, false) && lightIts.mesh == selectedLightMesh) ? 1.0 : 0.0;

                float geometry = visible * (std::fmax(0, normal.dot(dirToLight)) * std::fmax(0, lightNormal.dot(-dirToLight)) / std::pow(dirToLightDist, 2));

                BSDFQueryRecord bRecNEE(its.toLocal(-currentRay.d), its.toLocal(dirToLight), ESolidAngle);
                Color3f fNEE = bsdf->eval(bRecNEE);

                //throughput *= 0.5;

                Intersection tempIts;
                bool isHitLight = scene->getAccel()->rayIntersect(currentRay, tempIts, false) && tempIts.mesh == selectedLightMesh;

                //printf("%f \n", visible);
                //printf("%f \n", bsdf->pdf(bRecNEE));
                printf("%f %f %f \n", fNEE.x(), fNEE.y(), fNEE.z());

                if (lightPD > 0 && !isHitLight) {
                    Color3f Ltemp = throughput * fNEE * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
                    L += Ltemp;
                }
            }
            
            //////////////////////

            throughput *= f;

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