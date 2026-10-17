# Verification satchel for the length-eight note

Hello again from Yang Shifu. I have been following sequence coding, rank-metric
questions, and finite-field geometry for a long time, so I tend to carry a
small satchel of scripts whenever an equation starts looking suspicious.

The ZIP files in this folder are the finite verification packages used by the
length-eight proof note. They cover the cubic and sextic isolated-subline
closures, the determinant-rank four/five/six branches, Hermitian tangent and
small-field blocks, projection capacity, and the non-Hermitian graph checks.
Each package contains its own README, source files, field conventions,
coverage description, and replay command.

The file general_q_finite_replay_audit_20261009.json records the local replay
audit. Binary arrays match their recorded SHA-256 values; text certificates
are compared after the CRLF/LF conversion used by Windows and Unix tools.

The public proof PDF carries the mathematical argument. These packages carry
the computational luggage needed to retrace the finite branches without
guessing which counterexample a mysterious success flag had in mind.
