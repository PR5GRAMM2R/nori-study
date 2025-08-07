
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

    // Recursive
    /*
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        return Li(scene, sampler, ray, 0, Color3f(1.0f), 1.0f, false);
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray, int depth, Color3f& throughput, float eta_prod, bool lastWasSpecular) const {
        Color3f L(0.0f);

        Intersection its;

        bool hit = scene->rayIntersect(ray, its);

        if (!hit)
            return Color3f(0.0f);

        Normal3f normal = its.shFrame.n;
        Point3f point = its.p;

        const BSDF* bsdf = its.mesh->getBSDF();

        if (!bsdf)
            return Color3f(0.0f);

        if (its.mesh->isEmitter()) {
            if (depth == 0 || lastWasSpecular) {
                L += throughput * its.mesh->getEmitter()->getRadiance();
            }
        }

        ////////// BSDF //////////

        //Color3f resultBSDF(0.0f);
        Color3f sampleBSDF(0.0f);
        float pdfBSDF = 0;

        BSDFQueryRecord bRecBSDF(its.toLocal(-ray.d));
        sampleBSDF = bsdf->sample(bRecBSDF, sampler->next2D());

        pdfBSDF = bsdf->pdf(bRecBSDF);

        if (!bsdf->isDiffuse()) {
            lastWasSpecular = true;
        }
        else {
            lastWasSpecular = false;
        }

        //////////////////////////

        ////////// NEE //////////

        Color3f resultNEE(0.0f);
        float pdfNEE = 0;

        if (!lastWasSpecular) {
            std::vector<Mesh*> meshes = scene->getMeshes();
            Mesh** lightMeshes = new Mesh * [meshes.size()];
            Mesh* selectedLightMesh = nullptr;

            int lightCount = 0;
            for (Mesh* mesh : meshes) {
                if (mesh->isEmitter()) {
                    lightMeshes[lightCount++] = mesh;
                }
            }

            selectedLightMesh = lightMeshes[int(sampler->next1D() * lightCount)];

            Normal3f lightNormal;
            Point3f lightPoint;
            float lightPD = (1.0 / (float)lightCount) * selectedLightMesh->samplePosition(sampler, lightPoint, lightNormal);

            Vector3f dirToLight = lightPoint - point;
            float dirToLightDist = dirToLight.norm();
            dirToLight = dirToLight.normalized();

            Intersection lightIts;
            Ray3f rayToLight(point, dirToLight, 1e-4f, dirToLightDist + 1e-4f);
            float visible = (scene->rayIntersect(rayToLight, lightIts) && lightIts.mesh->getName() == selectedLightMesh->getName()) ? 1.0 : 0.0;

            float geometry = visible * (std::fmax(0, normal.dot(dirToLight)) * std::fmax(0, lightNormal.dot(-dirToLight)) / std::pow(dirToLightDist, 2));

            BSDFQueryRecord bRecNEE(its.toLocal(-ray.d), its.toLocal(dirToLight), ESolidAngle);
            Color3f evalNEE = bsdf->eval(bRecNEE);
            pdfNEE = bsdf->pdf(bRecNEE);

            if (lightPD != 0) {
                resultNEE = throughput * evalNEE * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
            }
        }

        /////////////////////////

        if (pdfBSDF < 0.0f || sampleBSDF.isZero())
            return Color3f(0.0f);

        throughput *= sampleBSDF;
        eta_prod *= bRecBSDF.eta;

        if (depth >= 3) {
            float p = std::min(throughput.maxCoeff() * eta_prod * eta_prod, 0.99f);

            if (sampler->next1D() > p)
                return Color3f(0.0f);

            throughput /= p;
        }

        Ray3f newRay(point, its.toWorld(bRecBSDF.wo));

        if (lastWasSpecular) {
            L += Li(scene, sampler, newRay, depth + 1, throughput, eta_prod, lastWasSpecular);
        }
        else {
            L += resultNEE + Li(scene, sampler, newRay, depth + 1, throughput, eta_prod, lastWasSpecular);
        }

        return L;
    }
    */
    
    // Loop
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Color3f L(0.0f);
        Color3f throughput = Color3f(1.f);
        float eta_prod = 1.0f;

        Ray3f currentRay = ray;
        int depth = 0;

        bool lastWasSpecular = false;

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

            if (its.mesh->isEmitter()) {
                if (depth == 0 || lastWasSpecular) {
                    L += throughput * its.mesh->getEmitter()->getRadiance();
                }
            }

            lastWasSpecular = false;

            ////////// BSDF //////////

            Color3f resultBSDF(0.0f);
            Color3f sampleBSDF(0.0f);
            float pdfBSDF = 0;

            BSDFQueryRecord bRecBSDF(its.toLocal(-currentRay.d));
            sampleBSDF = bsdf->sample(bRecBSDF, sampler->next2D());

            pdfBSDF = bsdf->pdf(bRecBSDF);

            if (!bsdf->isDiffuse()) {
                lastWasSpecular = true;
            }
            else {
                lastWasSpecular = false;
            }

            //////////////////////////

            ////////// NEE //////////

            Color3f resultNEE(0.0f);
            float pdfNEE = 0;

            if (!lastWasSpecular) {
                std::vector<Mesh*> meshes = scene->getMeshes();
                Mesh** lightMeshes = new Mesh * [meshes.size()];
                Mesh* selectedLightMesh = nullptr;

                int lightCount = 0;
                for (Mesh* mesh : meshes) {
                    if (mesh->isEmitter()) {
                        lightMeshes[lightCount++] = mesh;
                    }
                }

                selectedLightMesh = lightMeshes[int(sampler->next1D() * lightCount)];

                Normal3f lightNormal;
                Point3f lightPoint;
                float lightPD = (1.0 / (float)lightCount) * selectedLightMesh->samplePosition(sampler, lightPoint, lightNormal);

                Vector3f dirToLight = lightPoint - point;
                float dirToLightDist = dirToLight.norm();
                dirToLight = dirToLight.normalized();

                Intersection lightIts;
                Ray3f rayToLight(point, dirToLight, 1e-4f, dirToLightDist + 1e-4f);
                float visible = (scene->rayIntersect(rayToLight, lightIts) && lightIts.mesh->getName() == selectedLightMesh->getName()) ? 1.0 : 0.0;

                float geometry = visible * (std::fmax(0, normal.dot(dirToLight)) * std::fmax(0, lightNormal.dot(-dirToLight)) / std::pow(dirToLightDist, 2));

                BSDFQueryRecord bRecNEE(its.toLocal(-currentRay.d), its.toLocal(dirToLight), ESolidAngle);
                Color3f evalNEE = bsdf->eval(bRecNEE);
                pdfNEE = bsdf->pdf(bRecNEE);

                if (lightPD != 0) {
                    resultNEE = evalNEE * geometry * selectedLightMesh->getEmitter()->getRadiance() / lightPD;
                }

                L += throughput * resultNEE;
            }

            /////////////////////////

            if (pdfBSDF < 0.0f || sampleBSDF.isZero())
                break;

            throughput *= sampleBSDF;
            eta_prod *= bRecBSDF.eta;

            if (depth >= 3) {
                float p = std::min(throughput.maxCoeff() * eta_prod * eta_prod, 0.99f);

                if (sampler->next1D() > p)
                    break;

                throughput /= p;
            }

            currentRay = Ray3f(point, its.toWorld(bRecBSDF.wo));

            depth++;
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