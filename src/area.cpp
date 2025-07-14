#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList& props) {
        m_radiance = props.getColor("radiance");
    }

    std::string toString() const {
        return tfm::format(
            "Area Light[\n"
            "  radiance = %s\n"
            "]", m_radiance.toString());
    }
private:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END