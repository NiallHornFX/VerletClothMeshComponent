
// Implements
#include "HashGrid.h"

#include "VerletClothMeshComponent.h" // FVerletClothParticle

#include "Math/RandomStream.h"
#include "DrawDebugHelpers.h"

//#include "Stats/Stats.h"
//DECLARE_STATS_GROUP(TEXT("VerletCloth"), STATGROUP_VerletClothComponent, STATCAT_Advanced);
//DECLARE_CYCLE_STAT(TEXT("VerletCloth Grid Hash"), STAT_VerletCloth_GridHashTime, STATGROUP_VerletClothComponent);

hash_grid::hash_grid(UVerletClothMeshComponent *cloth, UWorld *world, int32 c_dim, float g_size, bool draw) 
	: Cloth(cloth) 
	, World(world)
	, C_dim(c_dim)
	, G_size(g_size)
	, draw_grid(draw)
{
	C_count = C_dim * C_dim * C_dim;
	C_size = G_size / (float)C_dim;

	// Allocate Outer HashTable/Grid bounds. 
	hashgrid = new TArray<FVerletClothParticle*>*[C_count];
	for (uint32 i = 0; i < C_count; ++i) { hashgrid[i] = nullptr; }
	hashgrid_mask.AddDefaulted(C_count);
}

hash_grid::~hash_grid()
{
	for (uint32 i = 0; i < C_count; ++i)
	{
		if (auto &ta = hashgrid[i])
		{
			delete ta; ta = nullptr;
		}
	}
	delete[] hashgrid;
}

void hash_grid::particle_hash()
{
	//SCOPE_CYCLE_COUNTER(STAT_VerletCloth_GridHashTime);

	float C_size_rec = 1.0f / C_size;
	// Particle Position Key, hash --> 1D grid cell index. Output == 0-(bck_cnt-1) Indices. 
	auto hashPt = [&](const FVector &PPos) -> uint32
	{
		int32 xx = static_cast<int32>(PPos.X * C_size_rec), yy = static_cast<int32>(PPos.Y * C_size_rec), zz = static_cast<int32>(PPos.Z * C_size_rec); // int3 key.
		int32 n = 73856093 * xx ^ 19349663 * yy ^ 83492791 * zz;
		return static_cast<uint32>(n % C_count); // 1D cell idx. 
	};

	for (int32 i = 0; i < Cloth->Particles.Num(); ++i)
	{
		FVerletClothParticle &curPt = Cloth->Particles[i];
		uint32 h_idx = hashPt(curPt.Position);

		// Alloc TArray Particle List at resulting cell in hashgrid (if not alloc'd), and append curPt ptr to that list
		if (hashgrid[h_idx] == nullptr) hashgrid[h_idx] = new TArray<FVerletClothParticle*>;
		hashgrid[h_idx]->Push(&curPt);

		hashgrid_mask[h_idx] = true;
		curPt.C_idx = h_idx; // Store cell on pt. 

		// Viz - ! Colour per hash cell.
		if (draw_grid)
		{
			FRandomStream rng_r(h_idx), rng_g(h_idx + 123), rng_b(h_idx + 4567);
			int32 r = static_cast<int32>(rng_r.FRandRange(0, 255)), g = static_cast<int32>(rng_g.FRandRange(0, 255)), b = static_cast<int32>(rng_b.FRandRange(0, 255));
			curPt.Col = FColor(r, g, b, 255);
			#ifdef DBG_DRAW_GRIDPTS
			DrawDebugSphere(World, curPt.Position, Cloth->ParticleRadius, 3, FColor(r, g, b, 255), false, Cloth->SubstepTime);
			#endif
		}
	}
}