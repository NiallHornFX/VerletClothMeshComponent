// Implements
#include "VerletClothMeshComponent.h"

// UE4 Headers
#include "Engine/World.h"
#include "DrawDebugHelpers.h"
#include "Math/RandomStream.h"
#include "WorldCollision.h"
#include "Engine/CollisionProfile.h"
#include "Stats/Stats.h"

#define DEBUG_PRINT_LOG
#define DEBUG_DRAW_CONSTRAINTS

DECLARE_STATS_GROUP(TEXT("VerletCloth"),                STATGROUP_VerletCloth,               STATCAT_Advanced);
DECLARE_CYCLE_STAT(TEXT("VerletCloth Sim"),             STAT_VerletCloth_SimTime,            STATGROUP_VerletCloth);
DECLARE_CYCLE_STAT(TEXT("VerletCloth Constraints"),     STAT_VerletCloth_ConTime,            STATGROUP_VerletCloth);
DECLARE_CYCLE_STAT(TEXT("VerletCloth World Collision"), STAT_VerletCloth_WorldCollisionTime, STATGROUP_VerletCloth);
DECLARE_CYCLE_STAT(TEXT("VerletCloth Self Collision"),  STAT_VerletCloth_SelfCollisionTime,  STATGROUP_VerletCloth);
DECLARE_CYCLE_STAT(TEXT("VerletCloth Integrate"),       STAT_VerletCloth_IntegrateTime,      STATGROUP_VerletCloth);

UVerletClothMeshComponent::UVerletClothMeshComponent(const FObjectInitializer& ObjectInitializer)
	: Super(ObjectInitializer)
{
	PrimaryComponentTick.bCanEverTick = true; bTickInEditor = true;

	// SM Init
	sm = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("ClothStaticMesh"));

	// UPropery Members Init
	bShowStaticMesh = true;
	bSimulate = false;
	ConstraintIterations = 4;
	bUseBendConstraints = true;
	bShow_Constraints = true;
	bWorldCollision = true;
	bSelfCollision = false;
	SubstepTime = 0.02f;
	StiffnessCoefficent = 0.5f; 
	CollisionFriction = 0.2f; 
	ClothForce = FVector(0.0f);
	ClothGravityScale = 0.1f; 
	bUse_Sleeping = false;
	bShow_Sleeping = false;
	Sleep_DeltaThreshold = 0.025f; 

	ParticleMass = 1.0f; ParticleRadius = 2.5f;

	// smData init
	smData.vb = nullptr, smData.cvb = nullptr, smData.smvb = nullptr, smData.ib = nullptr;
}

void UVerletClothMeshComponent::OnRegister()
{
	Super::OnRegister();

	// Recreate Editor Cloth State, at PIE Time.
	UWorld *world = GetWorld();
	if (world->IsPlayInEditor())
	{
		BuildClothState(); 
	}

	// Prop Updates
	sm->SetVisibility(bShowStaticMesh);
}

void UVerletClothMeshComponent::TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	Dt = DeltaTime;
	St = FMath::Max(SubstepTime, 0.01f); // Clamp to min substeptime. 

	if (bSimulate)
	{
		At += Dt; 
		while (At > St)
		{
			SCOPE_CYCLE_COUNTER(STAT_VerletCloth_SimTime);
			// Solve
			Integrate(St);
			EvalClothConstraints();
			ClothCollisionWorld(); 
			
			At -= St;
		}

		TickUpdateCloth();
	}
	else if (!bSimulate)
	{
		// Sleep
	}

}

// Reset Cloth to Inital State defined by Cloth Static Mesh.
void UVerletClothMeshComponent::ResetToInitalState()
{
	// Check this was not called before BuildClothState();
	if (!clothStateExists) { UE_LOG(LogType, Error, TEXT("VerletClothMeshComponent: No Inital Cloth State Exists; BuildClothState First.")); return; }

	for (int32 i = 0; i < particleCount; ++i)
	{
		FVector P = smData.Pos[i] + GetComponentLocation(); // Pts Pos, with Component Translation Offset.
		Particles[i].Position = P, Particles[i].PrevPosition = P;
	}
	TickUpdateCloth();
}

// Get Static Mesh Component, Static Mesh Data and Build Procedual Mesh Component Version of Inital State. Also Pass to and initalize Particles.
void UVerletClothMeshComponent::BuildClothState()
{
	UStaticMesh *usm = sm->GetStaticMesh();
	if (sm == nullptr) { UE_LOG(LogTemp, Error, TEXT("ERR::VerletCloth::No Static Mesh Set"));  return; }
	
	// Store SMLOD0 Buffer Pointers
	FStaticMeshLODResources *lod0 = *(usm->RenderData->LODResources.GetData()); 
	smData.vb   = &(lod0->VertexBuffers.PositionVertexBuffer); // Pos
	smData.smvb = &(lod0->VertexBuffers.StaticMeshVertexBuffer); // Static Mesh Buffer
	smData.cvb  = &(lod0->VertexBuffers.ColorVertexBuffer); // Colour
	smData.ib   = &(lod0->IndexBuffer); // Tri Inds

	smData.vert_count = smData.vb->GetNumVertices(); 
	smData.ind_count = smData.ib->GetNumIndices(); 
	smData.tri_count = smData.ind_count / 3;
	particleCount = smData.vert_count;
		
	#ifdef DEBUG_PRINT_LOG
	UE_LOG(LogTemp, Error, TEXT("STATIC MESH HAS ADJACNEY INFO = %d"), static_cast<int32>(lod0->bHasAdjacencyInfo));
	UE_LOG(LogTemp, Warning, TEXT("DBG::Static Mesh Vertex Count == %d | Index Count = %d"), smData.vert_count, smData.ind_count);
	#endif

	smData.Pos.AddDefaulted(smData.vert_count); smData.Col.AddDefaulted(smData.vert_count); smData.Normal.AddDefaulted(smData.vert_count); smData.Tang.AddDefaulted(smData.vert_count); smData.UV.AddDefaulted(smData.vert_count); 
	smData.Ind.AddDefaulted(smData.ind_count); smData.Tris.AddDefaulted(smData.tri_count);
	Particles.AddDefaulted(particleCount);

	// Need to add checks to delete previous procedual mesh data if exists.
	ClearAllMeshSections();
	
	// SmData Buffer --> Array Deserialization.
	for (int32 i = 0; i < smData.vert_count; ++i)
	{
		// SMesh-ProcMesh Init
		smData.Pos[i] = smData.vb->VertexPosition(i); // Pass Verts Without Component Location Offset initally.
		smData.Normal[i] = smData.smvb->VertexTangentZ(i);
		smData.Tang[i] = FProcMeshTangent(FVector(smData.smvb->VertexTangentX(i).X, smData.smvb->VertexTangentX(i).Y, smData.smvb->VertexTangentX(i).Z), false);
		lod0->bHasColorVertexData == true ?  smData.Col[i] = smData.cvb->VertexColor(i) : smData.Col[i] = FColor(255, 255, 255);
		smData.smvb->GetNumTexCoords() != 0 ? smData.UV[i] = smData.smvb->GetVertexUV(i, 0) : smData.UV[i] = FVector2D(0.0f); // Only support 1 UV Channel fnow.
		
		// Particle Init
		FVector vertPtPos = GetComponentLocation() + smData.vb->VertexPosition(i); // Pts With Component Location Offset
		Particles[i].Position = vertPtPos, Particles[i].PrevPosition = vertPtPos; 
		Particles[i].ID = i;
		lod0->bHasColorVertexData == true ? Particles[i].Col = smData.cvb->VertexColor(i) : Particles[i].Col = FColor(255, 255, 255);
	}
	// Indices 
	for (int32 i = 0; i < smData.ind_count; ++i) smData.Ind[i] = static_cast<int32>(smData.ib->GetIndex(i));

	// Build Cloth Mesh Section
	CreateMeshSection(0, smData.Pos, smData.Ind, smData.Normal, smData.UV, smData.Col, smData.Tang, false);
	bShowStaticMesh = false; sm->SetVisibility(bShowStaticMesh);
	clothStateExists = true;

	// Build Tri Array and Per Vert Shared Tri Array
	BuildTriArrays();

	// Build Constraints
	BuildClothConstraints();
}

void UVerletClothMeshComponent::BuildTriArrays()
{
	// Store Tris as 3 Indices within IntVector. 
	for (int32 t = 0, i = 0; t < smData.tri_count; ++t)
	{
		smData.Tris[t] = FIntVector(smData.Ind[i], smData.Ind[i + 1], smData.Ind[i + 2]);
		i += 3;
	}

	// Per Vert, store array of Tris (implictlly as indices) that i'm a part of. ( O(n^2) - ok for now as only done once. ) 
	if (smData.vtris) { delete smData.vtris; smData.vtris = nullptr; }
	else { smData.vtris = new TArray<FIntVector*>[smData.vert_count]; }
	for (int32 v = 0; v < smData.vert_count; ++v)
	{
		for (int32 t = 0; t < smData.tri_count; ++t)
		{
			FIntVector *t_tri = &(smData.Tris[t]);
			if (v == t_tri->X || v == t_tri->Y || v == t_tri->Z) smData.vtris[v].Add(t_tri);
		}
		#ifdef DEBUG_PRINT_LOG
		UE_LOG(LogTemp, Warning, TEXT("Vertex %d Tri Count = %d"), v, smData.vtris[v].Num());
		#endif
	}
}

// Builds Cloth Constraints and Stores them, oppose to evaulating per frame. This should be done only once, and not per tick unless topology is changing/tearing (tbd).
void UVerletClothMeshComponent::BuildClothConstraints()
{
	UWorld *world = GetWorld();
	// Particle : Vertex Index Mapping is 1:1 

	// Clear Previous Constraints
	Constraints.Empty();

	// Structual Constraints Only For now ...

	// For Each Particle(Vert) get triangles i'm a part of,
	for (int32 p = 0; p < particleCount; ++p)
	{
		FVerletClothParticle &curPt = Particles[p];
		TArray<FIntVector*> &tris = smData.vtris[p]; // Tris as Vert/Pt Indices within FIntVector.
		curPt.conCount = 0;

		// Each Tri im part of,
		for (int32 t = 0; t < tris.Num(); ++t)
		{
			FIntVector *TriInd = tris[t];
			// Each Vert/Particle Index of that Tri
			for (int32 i = 0; i < 3; ++i)
			{
				bool is_copy = false;
				int32 tvi = (*TriInd)[i];
				FVerletClothParticle &triPt = Particles[tvi];
				// First Check if this Constraint Pair (ID) Already Exists to avoid double constraints. 
				int32 Cur_Tri_conID = curPt.ID * triPt.ID;
				for (int32 c = 0; c < Constraints.Num(); ++c) if (Constraints[c].conID == Cur_Tri_conID) is_copy = true;
				// Or If Current Tri Pt, is Self CurPt, avoid self constraints.
				if (triPt.ID == curPt.ID || is_copy == true) continue;
				else
				{
					// Append New Constraint Of Particle Pair (CurPt and TriPt)
					Constraints.Emplace(curPt, triPt, this);
					curPt.conCount++;
				}
			}
		}
	}
}

// Rebuild Mesh Section 0, within tick with updated particle positions.  ToDo Update Normals and Tangents ...
void UVerletClothMeshComponent::TickUpdateCloth()
{
	check(particleCount == smData.vert_count);

	TArray<FVector> UpdtPos; UpdtPos.AddDefaulted(particleCount);
	for (int32 i = 0; i < particleCount; ++i)
	{
		UpdtPos[i] = Particles[i].Position - GetComponentLocation(); // Subtract Comp Translation Off as is added to ProcMesh Verts internally. 
	}
	UpdateMeshSection(0, UpdtPos, smData.Normal, smData.UV, smData.Col, smData.Tang);
}


/* Use Verlet Integration to Integrate Postion x(n+1) from Acceleration and implicit velocity using x(n) - x(n-1) */
void UVerletClothMeshComponent::Integrate(float i_St)
{
	SCOPE_CYCLE_COUNTER(STAT_VerletCloth_IntegrateTime);

	const float SubstepTimeSqr = i_St * i_St;
	const FVector Gravity = FVector(0, 0, GetWorld()->GetGravityZ()) * ClothGravityScale;

	for (int32 pt = 0; pt < particleCount; pt++)
	{
		FVerletClothParticle& Particle = Particles[pt];

		// Cloth Accel x''(t) = f/m (+ g) 
		const FVector Accel = Gravity + (ClothForce / ParticleMass);

		// x(n+1) = 2x(n) - x(n-1) + a(x) * dt^2
		// Integrate x''(n) to x(n+1) = x(n) + (x(n) - x(n-1)) + (a(x) * dt^2)
		const FVector NewPosition = Particle.Position + (Particle.Position - Particle.PrevPosition) + (Accel * SubstepTimeSqr);
		Particle.PrevPosition = Particle.Position; Particle.Position = NewPosition;
	}
}

// Cloth Particle and World Collision 
void UVerletClothMeshComponent::ClothCollisionWorld()
{
	SCOPE_CYCLE_COUNTER(STAT_VerletCloth_WorldCollisionTime);

	world_collided = false;
	UWorld *World = GetWorld();

	if (World && GetCollisionEnabled() != ECollisionEnabled::NoCollision)
	{
		// Get collision settings from component
		FCollisionQueryParams Params(SCENE_QUERY_STAT(VerletClothCollision));
		ECollisionChannel TraceChannel = GetCollisionObjectType();
		FCollisionResponseParams ResponseParams(GetCollisionResponseToChannels());

		for (int32 pt = 0; pt < particleCount; ++pt)
		{
			FVerletClothParticle &Particle = Particles[pt];
			//if (bPinTop && Particle.state == 0) continue; // Pinned Pt. 
			FHitResult Result;
			bool bHit = World->SweepSingleByChannel(Result, Particle.PrevPosition, Particle.Position, FQuat::Identity, TraceChannel, FCollisionShape::MakeSphere(ParticleRadius), Params, ResponseParams);
			if (bHit)
			{
				if (Result.bStartPenetrating)
				{
					Particle.Position += (Result.Normal * Result.PenetrationDepth);
					world_collided = true;
				}
				else
				{
					Particle.Position = Result.Location;
				}

				// Zero out any positive restitution velocity.
				FVector Delta = Particle.Position - Particle.PrevPosition;
				float NormalDelta = Delta | Result.Normal;
				FVector PlaneDelta = Delta - (NormalDelta * Result.Normal);
				Particle.PrevPosition += (NormalDelta * Result.Normal);

				if (CollisionFriction > 1e-04)
				{
					FVector ScaledPlaneDelta = PlaneDelta * CollisionFriction;
					Particle.PrevPosition += ScaledPlaneDelta;
				}
			}
		}
	}
}


// Solve a single distance constraint between a pair of particles by modifying particle postions to minimize CurDist-RestLength delta. 
void UVerletClothMeshComponent::SolveDistanceConstraint(FVerletClothParticle& ParticleA, FVerletClothParticle& ParticleB, float RestLength, float StiffnessCoeff)
{
	FVector PPos_A = ParticleA.Position, PPos_B = ParticleB.Position;

	// Find current vector between particles
	FVector Delta = ParticleB.Position - ParticleA.Position;
	float CurrentDistance = Delta.Size();
	float ErrorFactor = (CurrentDistance - RestLength) / CurrentDistance;

	// Project Pairs to Minimize Delta (If active state)
	if (ParticleA.state == 1) ParticleA.Position += ErrorFactor * StiffnessCoeff * Delta;
	if (ParticleB.state == 1) ParticleB.Position -= ErrorFactor * StiffnessCoeff * Delta;
}

// Solve Pre-Build Constraints on the Cloth Particles for k number of iterations.
void UVerletClothMeshComponent::EvalClothConstraints()
{
	UWorld *world = GetWorld();

	for (int32 k = 0; k < ConstraintIterations; ++k)
	{
		SCOPE_CYCLE_COUNTER(STAT_VerletCloth_ConTime);
		for (FVerletClothConstraint &con : Constraints)
		{
			SolveDistanceConstraint(con.Pt0, con.Pt1, con.restLength, StiffnessCoefficent);
			if (bShow_Constraints) {
				if (k == ConstraintIterations - 1) DrawDebugLine(world, con.Pt0.Position, con.Pt1.Position, FColor(255, 0, 0), false, Dt);
			}
		}
	}
}

// Debug Function to show Particle Radii in the Editor.
void UVerletClothMeshComponent::DBG_ShowParticles()
{
	UWorld *world = GetWorld();
	for (const FVerletClothParticle &pt : Particles)
	{
		DrawDebugSphere(world, pt.Position, ParticleRadius, 3, FColor(255, 0, 0, 1), false, 3.0f);
	}
}

void UVerletClothMeshComponent::DBG_ShowTris()
{
	// For Now just show Indv Spheres per Tri, to validate storing Tris implicilty by there vert indices only within a FIntVector ...
	UWorld *world = GetWorld();
	for (int32 i = 0; i < smData.tri_count; ++i)
	{
		FRandomStream rng_r(i), rng_g(i + 123), rng_b(i + 4567);
		int32 r = static_cast<int32>(rng_r.FRandRange(0, 255)), g = static_cast<int32>(rng_g.FRandRange(0, 255)), b = static_cast<int32>(rng_b.FRandRange(0, 255));
		FVector cPos = (smData.Pos[smData.Tris[i].X] + smData.Pos[smData.Tris[i].Y] + smData.Pos[smData.Tris[i].Z]) * 0.333f; 
		cPos += GetComponentLocation();
		DrawDebugSphere(world, cPos, ParticleRadius * 1.25f, 3, FColor(r, g, b, 255), false, 5.0f);
	}
}