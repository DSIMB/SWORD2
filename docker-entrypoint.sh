#!/bin/bash
set -e

# Default UID and GID
USER_ID=${USER_ID:-1000}
GROUP_ID=${GROUP_ID:-1000}

# Create a group with the specified GID if it doesn't exist
if ! getent group "swordgroup" >/dev/null; then
    groupadd -g "$GROUP_ID" swordgroup
fi

# Create a user with the specified UID and GID if it doesn't exist
if ! id "sworduser" >/dev/null 2>&1; then
    useradd -m -u "$USER_ID" -g "swordgroup" sworduser
fi

# Change ownership of the app directory to the new user
chown -R sworduser:swordgroup /app

# Execute the main application as the new user
exec gosu sworduser python SWORD2.py "$@"
