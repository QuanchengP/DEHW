# DEHW
Design and analysis of double-enveloping hourglass worm (DEHW) drive with planar generatrix. Some adopted approaches are interval Newton method, virtual element method, dual mortar method, et al.

In case of any problem in code compilation or execution, please feel free to contact us by email or other available communications.

1.Platform: Cygwin on Windows or Linux (only one of them is needed)
> + Cygwin link: <https://www.cygwin.com/>
> + WSL link: <https://learn.microsoft.com/zh-cn/windows/wsl/basic-commands>

2.Prerequisite package: make, g++

3.Prerequisite library: PETSc, Eigen, filib++ (PETSc is only needed for Linux and not needed for Cygwin)
> + PETSc link: <https://petsc.org/release/>
> + Eigen link: <https://gitlab.com/libeigen/eigen/>
>   + or: <https://eigen.tuxfamily.org/index.php?title=Main_Page>
> + filib++ link: <http://www2.math.uni-wuppertal.de/wrswt/software/filib.html>
>   + The source code of filib++ is slightly modified ("throw(interval_io_exception)" to "noexcept(false)").

4.Compiling command: $ make
> + For Cygwin: the variables “INCL, LIB” in “makefile” need to be modified to specify the directory of prerequisite library.
> + For Linux: the variables “INCL1, INCL2 INCL3, LIB1, LIB2” in “makefile” need to be modified to specify the directory of prerequisite library.

5.Execution: the minimum required RAM is 64GB
> + The Cygwin supports execution of all numerical examples, whereas the efficiency is slow at several numerical examples because of the adoption of IncompleteLUT preconditioner of BiCGSTAB in Eigen. The efficiency is improved by KSPBICG in PETSc.
> + The 23/40-th step of I/R/F2/T1/C1, the 24/40-th step of III/R/F2/T1/C1 encounter nonconvergence of active set strategy.
