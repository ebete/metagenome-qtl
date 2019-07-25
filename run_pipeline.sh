#!/usr/bin/env bash

# input arguments
SNAKEMAKE_PARAMS=${@:1}

SLURM_CFG="${HOME}/.config/snakemake/slurm/"

# set work dir to this script location
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd "${dir}"

check_updates() {
	UPSTREAM=${1:-'@{u}'}
	LOCAL=$(git rev-parse @)
	REMOTE=$(git rev-parse "${UPSTREAM}")
	BASE=$(git merge-base @ "${UPSTREAM}")

	REPO_STATUS="\033[97;41;1m diverged"
	if [[ -z "${LOCAL}" ]]; then
		REPO_STATUS="\033[97;41;1m no git workspace found"
	elif [[ "${LOCAL}" == "${REMOTE}" ]]; then
		REPO_STATUS="\033[32;1m up-to-date"
	elif [[ "${LOCAL}" == "${BASE}" ]]; then
		REPO_STATUS="\033[33;1m need to pull"
	elif [[ "${REMOTE}" == "${BASE}" ]]; then
		REPO_STATUS="\033[35;1m need to push"
	fi

	printf "\n\tGit workspace status: %b \033[0m\n\n" "${REPO_STATUS}"
}

# print the status of the git workspace
check_updates 2>/dev/null

# when running on the SLURM master node
if [[ "$(hostname -s)" == "nioo0004" ]]; then
	if [[ ! -d "${SLURM_CFG}" ]]; then
		echo "SLURM configuration not found at '${SLURM_CFG}'. Exiting ..."
		exit 1
	fi

	SNAKEMAKE_PARAMS="--cluster-config \"slurm.yml\" \
		--latency-wait 10 \
		--local-cores 1 \
		--profile \"${SLURM_CFG}\" \
		--cluster \"${SLURM_CFG}/slurm-submit.py --account {cluster.account} --job-name {cluster.job-name} --nodes {cluster.nodes}\" \
		${SNAKEMAKE_PARAMS}"
fi

# print and execute command
CMD="snakemake --use-conda --printshellcmds --restart-times 0 ${SNAKEMAKE_PARAMS}"
echo "${CMD}"
eval "${CMD}"
