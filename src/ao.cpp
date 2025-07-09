
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

#define MAX_SAMPLES 10

class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        Point3f p = its.p;

        float hitValue;
        Vector3f rayFromIntersectionDir;
        Point2f sample;

        pcg32 rng;
        Color3f result(0.0, 0.0, 0.0);

        for (int i = 0; i < MAX_SAMPLES; i++) {
            sample = Point2f(rng.nextFloat(), rng.nextFloat());
            rayFromIntersectionDir = Warp::squareToUniformSphere(sample);
            Ray3f rayFromIntersection(p, rayFromIntersectionDir, 1e-4f, INFINITY);
            hitValue = (scene->getAccel()->rayIntersect(rayFromIntersection, its, true)) ? 1.0 : 0.0;
            result += hitValue * (std::fmax(0, n.dot(rayFromIntersectionDir)) / M_PI);
        }

        return result;
    }

    std::string toString() const {
        return "AOIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END
