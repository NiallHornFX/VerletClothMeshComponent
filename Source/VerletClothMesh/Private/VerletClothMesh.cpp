// VerletClothMesh.cpp - VerletClothMeshComponent Plugin - Niall Horn 2020. 
#include "Modules/ModuleInterface.h"
#include "Modules/ModuleManager.h"

class FVerletClothMeshModule : public IModuleInterface
{
public:

	virtual void StartupModule() override;
	virtual void ShutdownModule() override;
};

IMPLEMENT_MODULE(FVerletClothMeshModule, VerletClothMesh)

void FVerletClothMeshModule::StartupModule()
{
	
}

void FVerletClothMeshModule::ShutdownModule()
{

}
	
