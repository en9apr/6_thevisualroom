==========================
How to Use Version Control
==========================

.. contents::
   :local:

Commands
========

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ git clone https://user@bitbucket.org/username/file.git

     - Clone the repository to the required directory. This will create a directory called "file" with the repository inside it
   * - ::

           $ git fetch origin

     - Fetch any changes from origin
   * - ::

           $ git merge origin/master

     - Merge the origin with the master branch
   * - ::

           $ git pull origin master

     - Fetch any changed from origin **AND** merge the origin with the master branch
   * - ::

           $ git pull

     - Fetch any changed from origin **AND** merge the origin with the master branch (origin and master are defaults)

How to configure git for a new virtual machine?
===============================================

* This is done once per machine

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ git config --global user.name "Andrew Roberts"
           $ git config --global user.email en9apr@hotmail.co.uk

     - This identifies the user of git


How to download a repository?
=============================

* cd into the directory you want the repository downloaded to
* It will create a subfolder with the same name as the repository on github

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ git clone https://github.com/en9apr/sphinx.git

     - This will clone the .git repository and create a new subdirectory called sphinx

How to put a directory under version control?
=============================================

* Firstly create a repository in GitHub (it must be public - unless you want to pay)
* cd into the directory you want under version control

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ gedit README.md

     - Create a readme file (add some text to this)
   * - ::

           $ git init

     - Initilises git
   * - ::

           $ git add README.md

     - Create a dummy file. Add all changes
   * - ::

           $ git commit -m "First commit"

     - Commit the changes
   * - ::

           $ git remote add origin https://github.com/en9apr/VM_Backup.git

     - Adds the origin to a github respository
   * - ::

           $ git push -u origin master

     - Pushes the changes from the master to the origin

You can now add new files and repeat the add, commit and push stages.

How to change the origin?
=========================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ git remote rm origin

     - Removes the current origin
   * - ::

           $ git remote set-url origin https://github.com/en9apr/sphinx.git

     - Changes the origin to github respository
   * - ::

           $ git add -A

     - Add all changes
   * - ::

           $ git push -u origin master

     - Pushes the changes to the origin (username and password required)


How to commit changes to git and push to respository?
=====================================================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ git add -A

     - Add all changes
   * - ::

           $ git commit -m "The new changes are identified by this sentence"

     - Commit the changes
   * - ::

           $ git push

     - Push the changes to the origin (will have to enter username and password for the origin)
