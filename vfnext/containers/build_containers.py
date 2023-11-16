import subprocess
import os

containers = [
    "pangolin_latest.sif",
    "singularity_snpeff.sif"]

container_commands = [
    "singularity build -F --fakeroot --sandbox pangolin_latest.sif Singularity_pangolin",
    "singularity build -F --fakeroot --sandbox singularity_snpeff.sif Singularity_snpEff"
    ]

failed_containers = []
already_built = []

def container_exists(container):
    return os.path.isdir(container)

def build_container(container, command):
    print(f"@ Building {container}...")
    
    if container_exists(container):
        print(f"  > Already exists <")
        already_built.append(container)
        return True
    
    try:
        subprocess.check_call(command, shell=True)
        print(f"  > Done <")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  > Failed <")
        print(f"Error: {e}")
        failed_containers.append((container, command))
        return False

print("Building containers:")

for container, command in zip(containers, container_commands):
    if not build_container(container, command):
        pass

print("\nSummary:")

if failed_containers:
    print("\nSome containers failed to build. Try to build the containers on 'ViralFlow/vfnext/containers/' directory. Here are the details:")
    for container, command in failed_containers:
        print(f"\nContainer {container} failed to build.")
        print(f"Suggested command to build {container}:")
        print(f"  {command}")
else:
    print("\nAll containers were successfully built.")

print("\nExecuting additional steps:")

print("  > Loading sars-cov2 nextclade dataset...")
nextclade_command = "singularity exec -B nextclade_dataset/sars-cov-2:/tmp nextclade:2.4.sif nextclade dataset get --name 'sars-cov-2' --output-dir '/tmp'"
try:
    subprocess.check_call(nextclade_command, shell=True)
    print("    > Done <")
except subprocess.CalledProcessError as e:
    print("    > Failed <")
    print(f"Error: {e}")

print("  > Downloading snpeff database catalog...")
snpeff_command = "singularity exec singularity_snpeff.sif snpEff databases > snpEff_DB.catalog"
try:
    subprocess.check_call(snpeff_command, shell=True)
    print("    > Done <")
except subprocess.CalledProcessError as e:
    print("    > Failed <")
    print(f"Error: {e}")
