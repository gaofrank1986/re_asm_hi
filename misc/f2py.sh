f2py -m spy -c ElemIntgl0.f90 mesh.f90 extend_mesh.f90 mod_func.f90 tripole_mod.f90 ./hi_intgl/hi_integral.f90 ./hi_intgl/hi_funcs.f90 
f2py -m test -c mesh.f90 tripole_mod.f90 extend_mesh.f90 mod_func.f90 ./hi_intgl/hi_integral.f90 ./hi_intgl/hi_funcs.f90 ElemIntgl0.f90
