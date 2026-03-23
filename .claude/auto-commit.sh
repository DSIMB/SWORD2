#!/usr/bin/env bash
# Auto-commit hook: stages and commits changes when Claude finishes a response
# Triggered by the Stop hook event

set -euo pipefail

# Check if we're in a git repo
git rev-parse --git-dir >/dev/null 2>&1 || exit 0

# Check if there are any changes (staged, unstaged, or untracked)
if git diff --quiet HEAD 2>/dev/null && [ -z "$(git status --porcelain 2>/dev/null)" ]; then
    exit 0
fi

# Build commit message from changed files
CHANGED=$(git diff --name-only HEAD 2>/dev/null; git ls-files --others --exclude-standard 2>/dev/null)
NUM_FILES=$(echo "$CHANGED" | grep -c . || true)

if [ "$NUM_FILES" -eq 0 ]; then
    exit 0
fi

if [ "$NUM_FILES" -eq 1 ]; then
    FILE_BASE=$(basename "$(echo "$CHANGED" | head -1)")
    MSG="auto: update $FILE_BASE"
elif [ "$NUM_FILES" -le 3 ]; then
    FILES=$(echo "$CHANGED" | xargs -I{} basename {} | sort -u | paste -sd ', ')
    MSG="auto: update $FILES"
else
    # Summarize by common directory or extension
    DIRS=$(echo "$CHANGED" | xargs -I{} dirname {} | sort -u | head -2 | paste -sd ', ')
    MSG="auto: update ${NUM_FILES} files in $DIRS"
fi

# Truncate message to 72 chars
MSG="${MSG:0:72}"

# Stage all changes and commit
git add -A
git commit -m "$MSG" --no-verify >/dev/null 2>&1 || true