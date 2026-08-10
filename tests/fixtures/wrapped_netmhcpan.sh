#!/bin/sh
# A stand-in for how netMHCpan is actually shipped: a shell wrapper that starts
# the real program, with other children inheriting the same pipes. Signalling
# only the wrapper leaves those children holding stdout and stderr open, which
# is what turns "kill it" into "wait for ever".
# Its stdout goes to /dev/null so only the stderr pipe is held: stdout reaching
# EOF is what ends the read, and a helper sitting on that would be a different
# and much more obvious problem.
sleep 60 >/dev/null &
exec "$FAKE_NETMHCPAN" "$@"
