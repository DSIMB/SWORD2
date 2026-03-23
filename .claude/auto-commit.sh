#!/usr/bin/env bash
# Auto-commit hook: stages and commits changes after tool use
# Receives hook JSON on stdin

set -euo pipefail

# Read stdin JSON
INPUT=$(cat)

# Check if we're in a git repo
git rev-parse --git-dir >/dev/null 2>&1 || exit 0

# Check if there are any changes (staged or unstaged or untracked)
if git diff --quiet HEAD 2>/dev/null && [ -z "$(git status --porcelain 2>/dev/null)" ]; then
    exit 0
fi

# Extract tool info for commit message
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // "unknown"')

case "$TOOL_NAME" in
    Edit|Write)
        FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // .tool_response.filePath // "files"')
        FILE_BASE=$(basename "$FILE_PATH" 2>/dev/null || echo "$FILE_PATH")
        MSG="auto: update $FILE_BASE"
        ;;
    Bash)
        CMD=$(echo "$INPUT" | jq -r '.tool_input.command // "command"' | head -c 80)
        MSG="auto: after bash - ${CMD}"
        ;;
    *)
        MSG="auto: changes after $TOOL_NAME"
        ;;
esac

# Truncate message to 72 chars
MSG="${MSG:0:72}"

# Stage all changes and commit
git add -A
git commit -m "$MSG" --no-verify >/dev/null 2>&1 || true
