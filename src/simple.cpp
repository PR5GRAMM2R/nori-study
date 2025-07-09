
#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

struct Light {
    Point3f position;
    Color3f energy;
};

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList& props) {
        light.position = props.getPoint("position");
        light.energy = props.getColor("energy");
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        Point3f p = its.p;
        Vector3f dirToLight = light.position - p;
        float dirToLightDistance = dirToLight.norm();
        dirToLight = dirToLight.normalized();

        float visible;
        Ray3f rayFromIntersection(p, dirToLight, 1e-4f, dirToLightDistance);
        visible = (scene->getAccel()->rayIntersect(rayFromIntersection, its, true)) ? 0.0 : 1.0;

        return Color3f((light.energy / (4.0 * std::pow(M_PI, 2))) *
            (std::fmax(0, n.dot(dirToLight) / std::pow(dirToLightDistance, 2))) *
            visible);
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }

protected:
    Light light;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END
