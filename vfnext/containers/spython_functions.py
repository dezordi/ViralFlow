from spython.main import Client
import os


def get_repository(repository_list):
    containers_name_list = []
    with open(repository_list, "r") as containers:
        file = containers.readlines()
        for container in file:
            repository_project_name = container.split("/")[1]
            container_version = container.split("/")[2].replace("\n", "")
            containers_name_list.append((repository_project_name, container_version))
    return containers_name_list


def container_pull(containers_dir, containers_name_list):
    for container in containers_name_list:
        print(f"Downloading container {container[1]}. This could be take a while. Please Wait ...")
        os.system(f"singularity pull -F {container[1]}.sif library://wallaulabs/viralflow/{container[1]}")  



def check_containers(containers_name_list, downloaded_list):
    download_list = os.scandir(downloaded_list)
    container_download_list = [container.name for container in download_list]
    containers_name = [(cont_name[0], f"{cont_name[1]}.sif") for cont_name in containers_name_list]
    missing_containers = []
    for i in containers_name:
        if i[1] not in container_download_list:
            missing_containers.append((i[0], i[1].split(".sif")[0]))
    return missing_containers


def containers_routine_pull(missing_containers_list, containers_dir, containers_names_list):
    lost_containers = missing_containers_list
    attempts = len(lost_containers) * 3
    while len(lost_containers) > 0 or attempts > 0:
        container_pull(containers_dir, missing_containers_list)
        lost_containers = check_containers(containers_names_list, containers_dir)
        attempts -= 1
        if len(lost_containers) == 0:
            break

