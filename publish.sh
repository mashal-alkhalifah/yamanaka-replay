#!/usr/bin/env bash
# Publish this repository to GitHub under your own account.
#
#   ./publish.sh <github-username>
#
# Before running, create an EMPTY repository named yamanaka-replay at
# https://github.com/new -- no README, no .gitignore, no licence. Those files
# already exist here, and letting GitHub add its own makes the first push
# conflict.
#
# The script replaces the OWNER placeholder in README.md and CITATION.cff with
# your account name, amends the last commit so the published history carries the
# right links, pushes, and sets the repository topics.
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "usage: ./publish.sh <github-username>" >&2
  exit 1
fi

OWNER="$1"
cd "$(dirname "$0")"

sed -i.bak "s|github.com/OWNER/yamanaka-replay|github.com/${OWNER}/yamanaka-replay|g" \
  README.md CITATION.cff paper/main.tex && rm -f README.md.bak CITATION.cff.bak paper/main.tex.bak

git add README.md CITATION.cff paper/main.tex
git commit -q --amend --no-edit

git remote remove origin 2>/dev/null || true
git remote add origin "https://github.com/${OWNER}/yamanaka-replay.git"
git branch -M main
git push -u origin main

# Topics make the repository findable; failure here is not fatal.
gh api -X PUT "repos/${OWNER}/yamanaka-replay/topics" \
  -f "names[]=group-testing" \
  -f "names[]=experimental-design" \
  -f "names[]=active-learning" \
  -f "names[]=cellular-reprogramming" \
  -f "names[]=perturbation-biology" \
  -f "names[]=reproducible-research" >/dev/null 2>&1 || true

echo
echo "Published: https://github.com/${OWNER}/yamanaka-replay"
