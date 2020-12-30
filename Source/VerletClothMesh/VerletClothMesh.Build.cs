// VerletClothMesh.Build.cs - VerletClothMeshComponent Plugin - Niall Horn 2020. 

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
