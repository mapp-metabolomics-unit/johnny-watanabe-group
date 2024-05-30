#!/bin/bash

echo "Starting project initialization script..."

# we print the PWD to see if we are in /tmp/ or not
echo "Current working directory: $PWD"


# Initialize git repository only if not already initialized
git init
if [ $? -ne 0 ]; then
    echo "Failed to initialize git repository"
    exit 1
fi
echo "Git repository initialized."

# Add remote origin only if not already added
git remote add origin "https://github.com/mapp-metabolomics-unit/johnny-watanabe-group.git"
if [ $? -ne 0 ]; then
    echo "Failed to add remote origin"
    exit 1
fi
echo "Remote origin added."

# Add all files to git
git add .
if [ $? -ne 0 ]; then
    echo "Failed to add files to git"
    exit 1
fi
echo "Files added to git."

# Commit the files
git commit -m "Initial commit"
if [ $? -ne 0 ]; then
    echo "Failed to commit files"
    exit 1
fi
echo "Files committed."

# Push to the remote repository
git push -u origin main
if [ $? -ne 0 ]; then
    echo "Failed to push to remote repository"
    exit 1
fi
echo "Pushed to remote repository."

echo "Project initialization script completed."
