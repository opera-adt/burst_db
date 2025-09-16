#!/bin/bash

# Script to automate the creation of a new release
# Usage: ./create-release.sh v0.13.0

set -e  # Exit on error

# Color codes for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check if version argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <version>"
    echo "Example: $0 v0.13.0"
    exit 1
fi

VERSION_TAG=$1
# Strip 'v' prefix if present to get just the number
VERSION_NUM=${VERSION_TAG#v}

# Ensure version tag starts with 'v'
if [[ ! $VERSION_TAG == v* ]]; then
    VERSION_TAG="v${VERSION_TAG}"
fi

echo -e "${GREEN}Creating release for version: ${VERSION_TAG} (${VERSION_NUM})${NC}"
echo "========================================="

# Step 1: Create git tag
if git rev-parse "$VERSION_TAG" >/dev/null 2>&1; then
    echo -e "${YELLOW}⚠ Git tag ${VERSION_TAG} already exists (skipping)${NC}"
else
    echo "Creating git tag ${VERSION_TAG}..."
    git tag ${VERSION_TAG}
    echo -e "${GREEN}✓ Git tag created${NC}"
fi

# Step 2: Install package in editable mode
echo "Installing package in editable mode..."
pip install -e . > /dev/null 2>&1
echo -e "${GREEN}✓ Package installed${NC}"

# Step 3: Create test directory
# Extract version numbers for directory naming (e.g., 0.13.0 -> 013)
VERSION_DIR=$(echo ${VERSION_NUM} | sed 's/\.//' | sed 's/\..*//')
TEST_DIR="test_${VERSION_DIR}"

if [ -d "${TEST_DIR}" ]; then
    echo -e "${YELLOW}⚠ Test directory ${TEST_DIR} already exists${NC}"
    read -p "Overwrite? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 1
    fi
    rm -rf ${TEST_DIR}
fi

echo "Creating test directory: ${TEST_DIR}..."
mkdir -p ${TEST_DIR}
cd ${TEST_DIR}
echo -e "${GREEN}✓ Test directory created${NC}"

# Step 4: Look for CMR survey from previous release
# Try to find the most recent test directory
PREV_TEST_DIR=$(ls -dt ../test_* 2>/dev/null | grep -v "${TEST_DIR}" | head -n1)

if [ -n "${PREV_TEST_DIR}" ] && [ -d "${PREV_TEST_DIR}" ]; then
    echo "Looking for CMR survey in: ${PREV_TEST_DIR}..."
    CMR_SURVEY=$(ls -t ${PREV_TEST_DIR}/cmr_survey*.csv.tar.gz 2>/dev/null | head -n1)
    if [ -n "$CMR_SURVEY" ]; then
        echo "Copying CMR survey from previous release..."
        cp "$CMR_SURVEY" .
        echo -e "${GREEN}✓ CMR survey copied${NC}"
    else
        echo -e "${YELLOW}⚠ No CMR survey found in ${PREV_TEST_DIR}${NC}"
        echo "You may need to obtain a new survey from SDS"
    fi
else
    echo -e "${YELLOW}⚠ No previous test directory found${NC}"
    echo "You'll need to obtain the CMR survey from SDS"
fi

# Step 5: Run make with VERSION parameter
echo "Running make to generate release files..."
echo "========================================="
make -f ../Makefile VERSION=${VERSION_NUM}

echo ""
echo "========================================="
echo -e "${GREEN}✓ Release ${VERSION_TAG} created successfully!${NC}"
echo "========================================="
echo ""
echo "Test directory: ${PWD}"
echo ""
echo "Next steps:"
echo "  1. Verify the generated files"
echo "  2. Push the tag: git push origin ${VERSION_TAG}"
echo "  3. Create GitHub release if everything looks good"
echo ""
echo "Generated files:"
ls -1 | head -20
