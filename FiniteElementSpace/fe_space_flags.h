#ifndef FESPACETEST_FE_SPACE_FLAGS_H
#define FESPACETEST_FE_SPACE_FLAGS_H

namespace chi_math::finite_element
{
  enum class FESpaceFlags : unsigned int
  {
    NO_FLAGS_SET                          = 0,
    PREBUILD_CELL_MAPPINGS                = (1 << 0),
    PRECOMPUTE_IntV_gradShapeI_gradShapeJ = (1 << 1),
    PRECOMPUTE_IntV_shapeI_gradShapeJ     = (1 << 2),
    PRECOMPUTE_IntV_shapeI_shapeJ         = (1 << 3),
    PRECOMPUTE_IntV_shapeI                = (1 << 4),
    PRECOMPUTE_IntV_gradShapeI            = (1 << 5),

    PRECOMPUTE_IntS_shapeI_shapeJ         = (1 << 6),
    PRECOMPUTE_IntS_shapeI                = (1 << 7),
    PRECOMPUTE_IntS_shapeI_gradshapeJ     = (1 << 8),

    PRECOMPUTE_FACE_NODE_MAPPINGS         = (1 << 9)
  };

  inline FESpaceFlags
  operator|(const FESpaceFlags f1, const FESpaceFlags f2)
  {
    return static_cast<FESpaceFlags>(static_cast<unsigned int>(f1) |
                                     static_cast<unsigned int>(f2));
  }
}//namespace chi_math

#endif //FESPACETEST_FE_SPACE_FLAGS_H
