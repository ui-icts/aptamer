`HAB_STUDIO_SUP=false hab studio enter`

So that we don't start the supervisor and
accidentally pull in a newer copy of a dependency

The habitat/plan.sh file links to specific versions
of dependencies in order to avoid having to get 
ViennaRNA to recompile with these newer versions
of GCC

```
build
source results/last_build.env
hab pkg upload results/$pkg_artifact
```

You can use the habitat GUI to promote to stable


