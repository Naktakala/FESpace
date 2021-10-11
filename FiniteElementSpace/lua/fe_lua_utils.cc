#include "fe_lua_utils.h"

void chi_math::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lua_register(L, "chiFESpaceTest",
               chi_math::lua_utils::chiFESpaceTest);
}