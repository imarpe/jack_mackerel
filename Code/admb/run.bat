jjm -nox -ind %1.ctl -iprint 100
copy jjm.par arc\%1.par
copy jjm.rep arc\%1.rep
copy jjm.std arc\%1.std
:: copy jjm.cor arc\%1.cor
copy fprof.yld arc\%1.yld
copy for_r.rep arc\%1_R.rep
call cleanad.bat
