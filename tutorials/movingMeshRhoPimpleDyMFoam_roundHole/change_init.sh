#!/bin/sh
for f in 0/* ; do sed -i -e "s/calculated/fixedValue/"  $f ; done;
