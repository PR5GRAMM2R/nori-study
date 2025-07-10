
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <pcg32.h>
#include <chrono>

NORI_NAMESPACE_BEGIN

#define MAX_SAMPLES 10

class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList& props) {
        /*sample = new Point2f[MAX_SAMPLES];
        rayFromIntersectionDir = new Vector3f[MAX_SAMPLES];
        rayFromIntersectionDirProb = new float[MAX_SAMPLES];

        for (int i = 0; i < MAX_SAMPLES; i++) {
            sample[i] = Point2f(rng.nextFloat(), rng.nextFloat());
            rayFromIntersectionDir[i] = Warp::squareToCosineHemisphere(sample[i]);
            rayFromIntersectionDirProb[i] = Warp::squareToCosineHemispherePdf(rayFromIntersectionDir[i]);
        }*/
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        //auto start = std::chrono::high_resolution_clock::now();                                    ////////////////////

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        Point3f p = its.p;

        float visible;
        Vector3f rayFromIntersectionDir;
        int count = 0;

        //pcg32 rng;
        Vector3f result(0.0, 0.0, 0.0);

        for (int i = 0; i < MAX_SAMPLES; i++) {
            rayFromIntersectionDir = Warp::squareToUniformSphere(sampler->next2D());
            float rayFromIntersectionDirProb = Warp::squareToUniformHemispherePdf(rayFromIntersectionDir);
            //Vector3f rayFromIntersectionDirWorld = its.shFrame.toWorld(rayFromIntersectionDir).normalized();
            //Ray3f rayFromIntersection(p + 1e-4f * rayFromIntersectionDirWorld, rayFromIntersectionDirWorld);
            rayFromIntersectionDir = (n.dot(Normal3f(0, 0, 1)) >= 0) ? rayFromIntersectionDir : -rayFromIntersectionDir;
            Ray3f rayFromIntersection(p + 1e-4f * rayFromIntersectionDir, rayFromIntersectionDir);
            visible = (scene->getAccel()->rayIntersect(rayFromIntersection, its, true)) ? 0 : 1;
            //float cosTheta = its.shFrame.cosTheta(rayFromIntersectionDir);
            if (rayFromIntersectionDirProb != 0.0) {
                result += Vector3f(visible) * std::fmax(0, n.dot(rayFromIntersectionDir)) / M_PI / rayFromIntersectionDirProb;
                //result += Vector3f(visible) * cosTheta / M_PI / rayFromIntersectionDirProb;
                count++;
            }
        }

        if(count > 0)
            result /= float(count);

        //auto end = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        //printf("11111 Li %f => %lld us.\n", ray.o.x(), duration.count());

        return Color3f(result.x(), result.y(), result.z());
    }

    std::string toString() const {
        return "AOIntegrator[]";
    }

protected:
    /*pcg32 rng;
    Point2f* sample;
    Vector3f* rayFromIntersectionDir;
    float* rayFromIntersectionDirProb;*/
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END
