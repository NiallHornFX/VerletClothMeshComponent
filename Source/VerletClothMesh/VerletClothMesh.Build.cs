// Verlet Cloth Mesh Plugin: VerletClothMesh.Build.cs 

using UnrealBuildTool;

public class VerletClothMesh : ModuleRules
{
	public VerletClothMesh(ReadOnlyTargetRules Target) : base(Target)
	{
		PublicDependencyModuleNames.AddRange(
		new string[]
		{
            "Core",
            "CoreUObject",
            "Engine",
			"ProceduralMeshComponent",
            "RHI"
		}
		);
		
		PrivateIncludePaths.Add("VerletClothMesh/Private/Headers"); // Private Headers.
	}
}
