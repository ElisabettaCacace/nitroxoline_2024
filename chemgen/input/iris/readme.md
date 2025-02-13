# TLDR

java -Xmx1g -jar iris.jar "Colony growth" "folder/where/jpg_files/are"

# Iris GUI

To run Iris as a window (GUI) application, open a command line window, cd into
the folder the Iris JAR file is in, and run the following command:

java -Xmx2g -jar iris.jar

Note: the Iris jar file can also be run on all systems by double clicking. Iris
will attempt to set the heap space automatically. If you run into “out of
memory” error messages, please make sure you have permanently set the maximum
Java heap size to 2g (2 Gigabyte). Google on how to do that in your particular
system/setup.

# Iris console (e.g. for scripting)

To run Iris as a console application (no GUI window), which can be rather
useful for scripting, you need to provide Iris with a profile name and a folder
that pictures are found in:

java -Xmx1g -jar iris.jar "profile name" "folder location"

NOTE: please keep the quotes when inputting the profile name

Valid profile names are the following:

“Colony growth”,
“Colony growth inverted”,
“B.subtilis sporulation”,
“CPRG profile”,
“Biofilm formation”,
“Morphology profile”,
“Morphology&Color profile”,
“Xgal assay”
