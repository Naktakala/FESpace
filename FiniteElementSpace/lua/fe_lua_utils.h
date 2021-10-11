#ifndef FINITE_ELEMENT_LUA_UTILS_H
#define FINITE_ELEMENT_LUA_UTILS_H

#include "ChiLua/chi_lua.h"

namespace chi_math::lua_utils
{
  int chiFESpaceTest(lua_State* L);
  void RegisterLuaEntities(lua_State* L);
}

#endif //FINITE_ELEMENT_LUA_UTILS_H