#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

#include "FiniteElementSpace/lua/fe_lua_utils.h"

int main(int argc, char* argv[])
{
  ChiLog&     log     = ChiLog::GetInstance();

  log.Log(LOG_0) << "FESpaceTest - Test run started";

  ChiTech::Initialize(argc,argv);

  auto& lua_console = ChiConsole::GetInstance();

  chi_math::lua_utils::RegisterLuaEntities(lua_console.consoleState);

  ChiTech::RunBatch(argc,argv);

  ChiTech::Finalize();

  log.Log(LOG_0) << "FESpaceTest - Execution finished";

  return 0;
}