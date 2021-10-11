CHI_TECH_DIR = os.getenv("CHI_TECH_DIR")

--############################################### Setup mesh
chiMeshHandlerCreate()

chiUnpartitionedMeshFromWavefrontOBJ(CHI_TECH_DIR.."/"..
        "ChiResources/TestObjects/TriangleMesh2x2Cuts.obj")

region1 = chiRegionCreate()

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED);

chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)



chiFESpaceTest()
