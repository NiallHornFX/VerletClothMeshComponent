#pragma once

#include "CoreMinimal.h"
#include "GenericPlatform/GenericPlatform.h"
#include "UObject/ObjectMacros.h"
#include "Engine/EngineTypes.h"

class UVerletClothMeshComponent;
struct FVerletClothParticle;

// Private Classes, for internal plugin use only, not exposed to edtior. 

// Hash Grid Class to support Spatial Hashing of cloth particles, of some finite grid_size into some set cell count. To Acclerate Self Collisions.
class hash_grid  
{
public:
	hash_grid(UVerletClothMeshComponent *cloth, UWorld *world, int32 c_dim, float g_size, bool draw);

	~hash_grid();

	void particle_hash(); 

private:

	uint32 C_dim, C_count;
	float C_size, G_size;
	bool draw_grid;

	UVerletClothMeshComponent *Cloth;
	UWorld *World;

	TArray<FVerletClothParticle*> **hashgrid;
	TArray<bool> hashgrid_mask;

	friend class UVerletClothComponent;
};