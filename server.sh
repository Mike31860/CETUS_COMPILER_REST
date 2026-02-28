#!/bin/bash
#SBATCH --partition=large-mem
#SBATCH --cpus-per-task=16
#SBATCH --job-name=SERVER_REST_COMPILER
#SBATCH --output=SERVER_REST_COMPILER.out
#SBATCH --error=SERVER_REST_COMPILER.err

source /opt/shared/anaconda/2024.02/etc/profile.d/conda.sh
conda activate CETUS_REST

echo "STARTING_SERVER"

./mvnw clean compile -U

echo "SETTING ENVIRONMENT VARIABLES..."

export C_INCLUDE_PATH="/lustre/parot/MiguelRosas/Section_Level_Tuning/TuningSystemJasonFinal-Section/SNU_NPB-1.0.3/NPB3.3-SER-C/CG"

echo "STARTING_SERVER_BUILD"
# Clean and compile first to ensure no stale classes
./mvnw clean compile

echo "RUNNING_SPRING_BOOT"
# Run the Spring Boot application. 
# We use -Dspring-boot.run.arguments to pass custom ports if needed.
./mvnw spring-boot:run