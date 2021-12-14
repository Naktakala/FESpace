CHI_TECH_DIR = os.getenv("CHI_TECH_DIR")

--############################################### Setup mesh
chiMeshHandlerCreate()

if (N==nil) then N=10; end
dx = 1.0/N
nodesx={}
for k=1,(N+1) do
    nodesx[k] = 0.0 + (k-1)*dx
end
umesh = chiMeshCreateUnpartitioned2DOrthoMesh(nodesx,nodesx)

region1 = chiRegionCreate()

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED);

chiSurfaceMesherExecute();
chiVolumeMesherExecute();

chiDestroyUnpartitionedMesh(umesh)

chiRegionExportMeshToVTK(region1, "ZMesh")

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)



chiFESpaceTest()
