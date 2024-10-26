# Viscoelastic3D
This is the 3D extension of Basilisk C viscoelastic solver


## running the codes

There are two ways to run the codes:

1. Using the vanilla basilisk method:

```shell
qcc -O2 -Wall -disable-dimensions {NameOfFile}.c -o {NameOfFile} -lm 
./{NameOfFile}
```

2. Using the makefile (can be interactively run using bview browser):

```shell
CFLAGS=-DDISPLAY=-1 make {NameOfFile}.tst
```

Check the localhost on {NameOfFile}/display.html. something like: [http://basilisk.fr/three.js/editor/index.html?ws://localhost:7100](http://basilisk.fr/three.js/editor/index.html?ws://localhost:7100) and run interactively.
