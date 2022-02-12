cmna
====

_A proof-of-concept C MNA implementation_

The core of this library is a LAPACK-dependent, _extremely_ bare LTI circuit
solver. With some effort, it could likely be included in an abstraction which
used multi-step linear estimates of nonlinear components.

Documentation is forthcoming; for now, look at the `app/cmna/main.c` source
file for usage.

Documentation of the method (modified nodal analysis) is best deferred to the
truly-excellent [QUCS Technical Papers][qucstp].

[qucstp]: http://qucs.sourceforge.net/tech/technical.html
