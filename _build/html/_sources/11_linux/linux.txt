==============
Linux Commands
==============

Here is a list of some commands in Linux

.. contents::
   :local:

How to remove a non-empty directory?
====================================

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Command
     - Meaning
   * - ::

           $ rm -rf directory-name

     - Removes directory-name forcefully (`f`) and recursively (`r`)

How to unzip tar.gz files?
==========================

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Command
     - Meaning
   * - ::

           $ tar -zxvf file.tar.gz

     - Unzips a tar.gz file

How to move files from Downloads to a new folder?
=================================================

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Command
     - Meaning
   * - ::

           $ mv $HOME/filename.txt $HOME/new-folder/filenname.txt

     - Moves one file from $HOME to $HOME/new-folder 
   * - ::

           $ mv $HOME/*.* $HOME/new-folder

     - Moves all files with names and extensions from $HOME to $HOME/new-folder



