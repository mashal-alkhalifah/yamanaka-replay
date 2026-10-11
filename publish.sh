#!/usr/bin/env bash
# One-shot publish helper.
#
#   ./publish.sh <github-username>
#
# Replaces the OWNER placeholder in README.md and CITATION.cff with the real
# account name, amends the initial commit, and pushes to
# https://github.com/<github-username>/yamanaka-replay
#
# The remote repository must already exist and be empty (no README, no license
# added by GitHub at creation time).
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "usage: ./publish.sh <github-username>" >&2
  exit 1
fi

OWNER="$1"
cd "$(dirname "$0")"

sed -i "s|github.com/OWNER/yamanaka-replay|github.com/${OWNER}/yamanaka-replay|g" \
  README.md CITATION.cff

git add README.md CITATION.cff
git commit --amend --no-edit

git remote remove origin 2>/dev/null || true
git remote add origin "https://github.com/${OWNER}/yamanaka-replay.git"
git branch -M main
git push -u origin main

echo
echo "Published: https://github.com/${OWNER}/yamanaka-replay"
