# OpenFOAM installation
OpenFOAM and utilities installation is very complete and has lots of possibilities but can be hard to master

## Paraview full support
Third Party provides you with paraview stuff but it does not work so here is how to do.

Server side (Linux):

    Setup (using `uname -a`): `Linux debian 4.19.0-6-amd64 #1 SMP Debian 4.19.67-2+deb10u1 (2019-09-20) x86_64 GNU/Linux`
    IP:        XX.XX.XX.XX
    Install official paraview package using `apt install paraview`. Currently, `paraview -V` gives ` paraview version 5.4.1`. Additional packages such as MPI may be installed.

Client side (Windows):

    Install official paraview binaries from `https://www.paraview.org/download/`; As I want better performances I choose the parallel version corresponding to the one I installed for Linux, namely `ParaView-5.4.1-Qt5-OpenGL2-MPI-Windows-64bit.zip` (Aug 21 2017).
    Install official Microsoft MPI package. Paraview versions have to be equal so they both need MPI support. Just go to `https://www.microsoft.com/en-us/download/details.aspx?id=57467` and click `Download`.

Connection (Windows client to Linux server):

    On the linux machine, start the paraview server with `mpirun -np 4 pvserver` and you should see the `nameofmachine:portnumber` (port is 11111 by default which I suggest you should not change) and a message saying it is waiting for connection.
    Fire up paraview (a blank cmd.exe window may appear, just leave it) and check version is correct in Help>About.
    Click on File>Connect then `Add Server` then give a name, set `Server Type` as `Client / Server`, set Host as Ip of server XX.XX.XX.XX and set Port as the port displayed in the console of the server (again, 11111 by default). Click 1Configure` and leave the settings to `Manual` and click Ok. Then double click on the server you just created and you should be good to go.

Possible issues: As said online, main issues come from versioning so make sure versions correspond to each other. Also make sure good ports are open and check firewall options if necessary.

What to do next: Note that the Pipeline Browser now features `remote(cs:XX.XX.XX.XX:yourport)`. If the server is purely graphical (that's the case for me) you might have a warning box `Server DISPLAY not accessible .... Remote rendering will be disabled`. I seek performance on server side so it looks bad but I suspect it might vanish as soon as I configure a proper environment on the server. Nevertheless you should be able to browse the remote server when clicking on `Open`.

Concerning OpenFOAM: the original `paraFoam` script (for local OpenFOAM post-processing) recognizes 2 kinds of data structure:

    test_case.OpenFOAM: OpenFOAM reader module created by OpenFOAM for paraview. It comes with `makeParaview` (installation of local version of paraview, included in Third-Party-x package). I had problem installing this on the server so I don't know nothing about that. But if we assume that this local paraview has been built with same options that our two equivalent paraview (which I assume require most of the flags of the script to be used, e.g. -python -python-lib etc..), there should exist a local version of pvserver as well that we can use the same way we use pvserver of official paraview package (in that scenario, you would change directory to local paraview of OpenFOAM and run `mpirun -np 4 ./pvserver`).
    test_case.foam: The official OpenFOAM reader from paraview (KitWare) itself. This one will work with our remote client-server thingy in any case.

I don't know which module has been created first but I know that the first one works well if the local paraview (the one from OF, in Third-Party-x) compiled  and the second one works in every case, as long as the paraview version compiled without major flaw.  However, the second options leads to a recurrent error stating `vtkMultiBlockDataSet (0x556b0eee10a0): Structure does not match. You must use CopyStructure before calling this method.`  This is really annoying and this is how to circumvent it:

    (On Server side) Load OpenFOAM environment and change directory to your test case. Run the command `foamToVTK` (can be run in parallel) which will create a VTK structure out of OpenFOAM data, in a folder `./VTK/`.
    (On Client side) Connect to server as explained above and instead of test_case.foam, open `test_case_..vtk`. Click `Apply` in Pipeline Browser and Paraview should not complain anymore.

