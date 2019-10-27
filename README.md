# OpenFOAM installation
OpenFOAM and utilities installation is very complete and has lots of possibilities but can be hard to master

## Paraview full support
for complete Paraview support, run
`./makeParaview -python`
remark: keep all options in headers set as false  (they will change later) OTHERWISE -mpi and -mesa support does not work. Moreovr -mpi does not work for me

remark: read comment in script headers, and use variables with
`whereis`
and do not forget to install time
`sudo apt install time`

My personnal header looks like:
```

```

./makeParaView  -cmake "/usr/bin/cmake" -qmake "/usr/bin/qmake" -mesa -mesa-include "/usr/include/GL" -mesa-lib  "/usr/lib/x86_64-linux-gnu/libOSMesa.so" -python -python-lib "/usr/lib/x86_64-linux-gnu/libpython2.7.so.1" -python-include "/usr/include/python2.7"

-mpi