# xForms-based Scope

## How to run and troubleshooting

To compile OAI, build with scope support:
```
./build_oai --build-lib nrscope ...
```

To use the scope, run the xNB or the UE with option `-d`. If you receive the
error
```
In fl_initialize() [flresource.c:995]: 5G-gNB-scope: Can't open display :0
In fl_bgn_form() [forms.c:347]: Missing or failed call of fl_initialize()
```
you can allow root to open X windows:
```
xhost +si:localuser:root
```

### Usage in gdb with Scope enabled

In gdb, when you break, you can refresh immediately the scope view by calling the display function.
The first parameter is the graph context, nevertheless we keep the last value for a dirty call in gdb (so you can use '0')

Example with no variable known
```
phy_scope_nrUE(0, PHY_vars_UE_g[0][0], 0, 0, 0)
```
or
```
phy_scope_gNB(0, phy_vars_gnb, phy_vars_ru, UE_id)
```

# Qt-based Scope

## Building Instructions
For the new qt-based scope designed for NR, please consider the following:

1. run the gNB or the UE with the option '--dqt'.
2. make sure to install the Qt5 packages before running the scope. Otherwise, the scope will NOT be displayed! Note that Qt6 does NOT work.
3. To build the new scope, add 'nrqtscope' after the '--build-lib' option. So, the complete command would be

   ```
   ./build_oai --gNB -w USRP --nrUE --build-lib nrqtscope
   ```

## New Features

1. New KPIs for both gNB and UE, e.g., BLER, MCS, throughout, and number of scheduled RBs.
2. For each of the gNB and UE, a main widget is created with a 3x2 grid of sub-widgets, each to display one KPI.
3. Each of the sub-widgets has a drop-down list to choose the KPI to show in that sub-widget.
4. Both of the gNB and UE scopes can be resized using the mouse movement.

## Troubleshoot

Similar as with the Xforms-based scope, you should allow root to open X
windows if using Xorg:
```
xhost +si:localuser:root
```

If you run on Wayland (e.g., Gnome on Ubuntu 22), you might get this error:
```
Warning: Ignoring XDG_SESSION_TYPE=wayland on Gnome. Use QT_QPA_PLATFORM=wayland to run on Wayland anyway.
Authorization required, but no authorization protocol specified
qt.qpa.xcb: could not connect to display :0
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
```

We have to force the use of Wayland. Before that, make sure that the Qt-Wayland
plugin is installed. Also, there seems to be a problem with the qt5-gtk
platform theme, so remove it:
```
sudo apt-get install qtwayland5
sudo apt-get remove qt5-gtk-platformtheme
```

Now, force the Wayland plugin; also make sure to pass `-E` to sudo:
```
sudo -E QT_QPA_PLATFORM=wayland ./nr-softmodem --sa -O <conf> --dqt
```

Furthermore, if you get the following warning, _followed by a segfault_ (i.e.,
if it runs, you can ignore this)
```
QStandardPaths: wrong ownership on runtime directory /run/user/1000, 1000 instead of 0
```
You can set the XDG runtime directory like one of the following options (try in
order, keep `QT_QPA_PLATFORM` if necessary)
```
sudo -E XDG_RUNTIME_DIR=/run/user/0 ./nr-softmodem ...
sudo -E XDG_RUNTIME_DIR=/tmp/runtime-root ./nr-softmodem ...
sudo -E XDG_RUNTIME_DIR= ./nr-softmodem ...
```
