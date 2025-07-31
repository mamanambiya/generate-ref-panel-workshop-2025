#!/bin/bash

# Multi-Platform Cromwell Execution Script
# GA4GH Hackathon 2025 - African Genomics Team
#
# This script demonstrates different approaches to handle Docker platform
# arguments at the Cromwell engine level instead of in WDL files

set -e

# Configuration
WDL_FILE="federated_imputation_pipeline.wdl"
INPUT_FILE="test_input.json"
CROMWELL_JAR="cromwell.jar"

echo "🚀 Multi-Platform Cromwell Execution Options"
echo "============================================="

# Detect platform
PLATFORM_DETECTED=""
if [[ "$OSTYPE" == "darwin"* ]] && [[ "$(uname -m)" == "arm64" ]]; then
    PLATFORM_DETECTED="apple-silicon"
    echo "🍎 Detected: Apple Silicon (ARM64)"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM_DETECTED="linux"
    echo "🐧 Detected: Linux"
else
    PLATFORM_DETECTED="other"
    echo "💻 Detected: Other platform"
fi

echo ""
echo "📋 Available Execution Methods:"
echo ""

# Method 1: Using Cromwell Configuration File
echo "1️⃣  Method 1: Cromwell Configuration File (Recommended)"
echo "   Uses cromwell.conf for automatic platform detection"
echo ""
echo "   Command:"
echo "   java -Dconfig.file=cromwell.conf -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE}"
echo ""

# Method 2: Environment Variables
echo "2️⃣  Method 2: Environment Variables"
echo "   Set Docker platform via environment variable"
echo ""
echo "   Commands:"
if [[ "$PLATFORM_DETECTED" == "apple-silicon" ]]; then
    echo "   export CROMWELL_DOCKER_PLATFORM='--platform linux/amd64'"
else
    echo "   export CROMWELL_DOCKER_PLATFORM=''"
fi
echo "   java -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE}"
echo ""

# Method 3: Backend-specific execution
echo "3️⃣  Method 3: Backend-Specific Execution"
echo "   Use different backends for different platforms"
echo ""
echo "   Local (with platform detection):"
echo "   java -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE} --backend Local"
echo ""
echo "   Azure (platform args removed automatically):"
echo "   java -jar ${CROMWELL_JAR} run ${WDL_FILE} -i test_input_azure.json --backend Azure"
echo ""

# Method 4: Docker environment preset
echo "4️⃣  Method 4: Docker Environment Preset"
echo "   Pre-configure Docker for platform compatibility"
echo ""
echo "   Commands:"
if [[ "$PLATFORM_DETECTED" == "apple-silicon" ]]; then
    echo "   docker context create multiarch --docker \"host=unix:///var/run/docker.sock,platform=linux/amd64\""
    echo "   docker context use multiarch"
fi
echo "   java -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE}"
echo ""

# Interactive selection
echo "🎯 Choose execution method:"
echo "1) Cromwell config file (recommended)"
echo "2) Environment variables"
echo "3) Backend-specific"
echo "4) Docker context"
echo "5) Exit"
echo ""

read -p "Enter your choice (1-5): " choice

case $choice in
    1)
        echo ""
        echo "🔧 Executing with Cromwell configuration file..."
        if [[ ! -f "cromwell.conf" ]]; then
            echo "❌ cromwell.conf not found! Please ensure the configuration file exists."
            exit 1
        fi
        java -Dconfig.file=cromwell.conf -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE}
        ;;
    2)
        echo ""
        echo "🔧 Executing with environment variables..."
        if [[ "$PLATFORM_DETECTED" == "apple-silicon" ]]; then
            export CROMWELL_DOCKER_PLATFORM='--platform linux/amd64'
            echo "Set CROMWELL_DOCKER_PLATFORM='--platform linux/amd64'"
        else
            export CROMWELL_DOCKER_PLATFORM=''
            echo "Set CROMWELL_DOCKER_PLATFORM=''"
        fi
        java -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE}
        ;;
    3)
        echo ""
        echo "🔧 Choose backend:"
        echo "1) Local (with platform detection)"
        echo "2) Azure (platform args removed)"
        read -p "Backend choice (1-2): " backend_choice
        
        if [[ "$backend_choice" == "1" ]]; then
            echo "Executing with Local backend..."
            java -Dconfig.file=cromwell.conf -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE} --backend Local
        elif [[ "$backend_choice" == "2" ]]; then
            echo "Executing with Azure backend..."
            java -Dconfig.file=cromwell.conf -jar ${CROMWELL_JAR} run ${WDL_FILE} -i test_input_azure.json --backend Azure
        else
            echo "❌ Invalid backend choice"
            exit 1
        fi
        ;;
    4)
        echo ""
        echo "🔧 Setting up Docker context..."
        if [[ "$PLATFORM_DETECTED" == "apple-silicon" ]]; then
            echo "Creating multiarch Docker context..."
            docker context create multiarch --docker "host=unix:///var/run/docker.sock,platform=linux/amd64" || true
            docker context use multiarch
            echo "Docker context set to multiarch"
        else
            echo "Using default Docker context (no changes needed)"
        fi
        java -jar ${CROMWELL_JAR} run ${WDL_FILE} -i ${INPUT_FILE}
        ;;
    5)
        echo "👋 Exiting..."
        exit 0
        ;;
    *)
        echo "❌ Invalid choice"
        exit 1
        ;;
esac

echo ""
echo "✅ Execution completed!"
echo ""
echo "💡 Pro Tips:"
echo "   - Use cromwell.conf for automatic platform detection"
echo "   - Environment variables work for simple overrides"
echo "   - Backend-specific configs handle cloud vs local differences"
echo "   - Docker contexts provide global platform settings"
echo ""
echo "🎉 Ready for GA4GH Hackathon 2025!" 