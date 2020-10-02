FC=pgf90
FCFLAGS=-O3 
HOMEPATH=$(PWD)
MACHINE="CRC"



msd_rot_calc: src/fortran/msd_rot_calc.f90 src/fortran/visc_calc.f90 src/fortran/flux_side.f90
	@echo "f=open('machine.name','w')" > src/python/path.py
	@echo "f.write('$(MACHINE)\n')" >> src/python/path.py
	@echo "f.write('$(HOMEPATH)\n')" >> src/python/path.py
	@echo "f.close()" >> src/python/path.py
	test -f module/path.include || echo "module use $(HOMEPATH)/module" >> ~/.bash_profile
	test -f module/path.include || source ~/.bash_profile
	@echo "prepend_path('PATH', '$(HOMEPATH)/src/exec/')" > module/path.include
	@echo "prepend_path('PATH', '$(HOMEPATH)')" >> module/path.include
	@echo "prepend_path('PATH', '$(HOMEPATH)/bin/')" >> module/path.include
	cat module/direct_calc_header module/path.include > module/Dir_Calc_Fluct.lua
	$(FC) $(FCFLAGS) -o src/exec/msd_rot_calc src/fortran/msd_rot_calc.f90
	$(FC) $(FCFLAGS) -o src/exec/matom_msd_rot_calc src/fortran/multi-atom_msd_rot_calc.f90
	$(FC) $(FCFLAGS) -o src/exec/visc_calc    src/fortran/visc_calc.f90
	$(FC) $(FCFLAGS) -o src/exec/flux_side    src/fortran/flux_side.f90
	mkdir -p bin/
	touch bin/test
	rm bin/*
	ln -s $(HOMEPATH)/src/python/vel_reselect.py bin/
	ln -s $(HOMEPATH)/src/python/grab_press.py bin/
	ln -s $(HOMEPATH)/src/python/init_segments.py bin/
	ln -s $(HOMEPATH)/src/python/combine_segments.py bin/
	ln -s $(HOMEPATH)/src/python/combine_segments2.py bin/
	ln -s $(HOMEPATH)/src/python/msd_fit.py bin/
	ln -s $(HOMEPATH)/src/python/reor_fit.py bin/
	ln -s $(HOMEPATH)/src/python/set_msd_calcs.py bin/
	ln -s $(HOMEPATH)/src/python/grab_flucts.py bin/
	ln -s $(HOMEPATH)/src/python/jump_rot.py bin/
	ln -s $(HOMEPATH)/src/python/parse_fit_results.py bin/
	ln -s $(HOMEPATH)/src/python/combine_weighted.py bin/
	ln -s $(HOMEPATH)/src/python/jump_fit.py bin/
	ln -s $(HOMEPATH)/src/python/msd.py bin/
	ln -s $(HOMEPATH)/src/python/cn_rot_calc_special.py bin/
	ln -s $(HOMEPATH)/src/python/parse_and_combine.py bin/
	ln -s $(HOMEPATH)/src/python/dat_to_pickle.py bin/
	ln -s $(HOMEPATH)/src/python/spectral_diffusion.py bin/
	ln -s $(HOMEPATH)/src/python/calc_eta.py bin/
	chmod 777 backbone.py
	chmod 777 bin/*
	chmod 777 src/exec



