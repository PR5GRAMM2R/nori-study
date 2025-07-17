
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <pcg32.h>
#include <chrono>

NORI_NAMESPACE_BEGIN

#define MAX_SAMPLES 1

class WHITTEDIntegrator : public Integrator {
public:
    WHITTEDIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        Point3f p = its.p;

        if (its.mesh->getBSDF()->isDiffuse()) {
            std::vector<Mesh*> meshes = scene->getMeshes();
            Mesh** lightMeshes = new Mesh * [meshes.size()];

            int lightCount = 0;
            for (Mesh* mesh : meshes) {
                if (mesh->isEmitter()) {
                    lightMeshes[lightCount++] = mesh;
                }
            }

            if (lightCount == 0)
                return Color3f(0.0f);

            Color3f result(0.0, 0.0, 0.0);

            int count = 0;
            for (int i = 0; i < MAX_SAMPLES; i++)
            {
                //Mesh* selectedLightMesh = lightMeshes[int(sampler->next1D() * (float)lightCount)];
                for (int l = 0; l < lightCount; l++) {
                    Mesh* selectedLightMesh = lightMeshes[l];

                    // 선택된 Light Mesh 에서 point, normal, pd 를 랜덤하게 뽑아오기
                    Normal3f lightNormal;
                    Point3f lightPoint;// = selectedLightMesh->sample(sampler, lightNormal, lightPD);
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
                        result += bsdf * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
                        count++;
                    }
                }
            }

            count /= lightCount;

            if (count > 0)
                result /= float(count);
            else
                result = Color3f(0.0f);

            if (its.mesh->isEmitter())
                return its.mesh->getEmitter()->getRadiance() + result;
            else
                return result;
        }
        else {
            Vector2f sample = sampler->next2D();

            BSDFQueryRecord bsdfQueryRecord(its.toLocal(-ray.d));
            Color3f c = its.mesh->getBSDF()->sample(bsdfQueryRecord, sample);

            Vector3f oRay = its.toWorld(bsdfQueryRecord.wo);
            Color3f li;

            if (sampler->next1D() < 0.95) {
                li = (1.0 / 0.95) * c * Li(scene, sampler, Ray3f(p + 1e-4f * oRay, oRay));
            }
            else{
                li = Color3f(0.0);
            }

            return li;
        }
    }

    std::string toString() const {
        return "WHITTEDIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(WHITTEDIntegrator, "whitted");
NORI_NAMESPACE_END



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <pcg32.h>
#include <chrono>

NORI_NAMESPACE_BEGIN

#define MAX_SAMPLES 10

class WHITTEDIntegrator : public Integrator {
public:
    WHITTEDIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        Point3f p = its.p;

        std::vector<Mesh*> meshes = scene->getMeshes();
        Mesh** lightMeshes = new Mesh * [meshes.size()];

        int lightCount = 0;
        for (Mesh* mesh : meshes) {
            if (mesh->isEmitter()) {
                lightMeshes[lightCount++] = mesh;
            }
        }

        if (lightCount == 0)
            return Color3f(0.0f);

        Color3f result(0.0, 0.0, 0.0);

        int count = 0;
        for (int i = 0; i < MAX_SAMPLES; i++)
        {
            //Mesh* selectedLightMesh = lightMeshes[int(sampler->next1D() * (float)lightCount)];
            for (int l = 0; l < lightCount; l++) {
                Mesh* selectedLightMesh = lightMeshes[l];

                // 선택된 Light Mesh 에서 point, normal, pd 를 랜덤하게 뽑아오기
                Normal3f lightNormal;
                Point3f lightPoint;// = selectedLightMesh->sample(sampler, lightNormal, lightPD);
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
                    result += bsdf * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
                    count++;
                }
            }
        }

        count /= lightCount;

        if (count > 0)
            result /= float(count);
        else
            result = Color3f(0.0f);

        if (its.mesh->isEmitter())
            return its.mesh->getEmitter()->getRadiance() + result;
        else
            return result;
    }

    std::string toString() const {
        return "WHITTEDIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(WHITTEDIntegrator, "whitted");
NORI_NAMESPACE_END
*/