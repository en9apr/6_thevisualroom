===
Git
===

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

How to configure git for a machine?
===================================

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
* Alternatively GitLab can be private and you don't have to pay https://gitlab.com
* cd into the directory you want under version control

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ cd existing_folder

     - Change to existing folder
     
   * - ::

           $ gedit README.md

     - Create a readme file (add some text to this)
   * - ::

           $ git init

     - Initilises git
   * - ::

           $ git remote add origin https://gitlab.com/apr207/twoPhaseEulerSrivastava.git

     - Adds the origin to a github respository  
     
   * - ::

           $ git add .

     - Create a dummy file. Add all changes
   * - ::

           $ git commit -m "Solver before any changes are made"

     - Commit the changes

   * - ::

           $ git push -u origin master

     - Pushes the changes from the master to the origin

You can now add new files and repeat the add, commit and push stages.

Example of how to create a project
==================================

Exact steps after creating a project twoPhaseEulerSrivastavaFrictionFoam:

::

    cd twoPhaseEulerSrivastavaFrictionFoam
    git init
    git remote add origin https://gitlab.com/apr207/twoPhaseEulerSrivastavaFrictionFoam.git
    git add .
    git commit -m "16/05/2018 solver working with explicit friction, but massflow is too high"
    git push -u origin master



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

           $ git remote set-url origin https://gitlab.com/apr207/twoPhaseEulerSrivastava.git

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

           $ git add .

     - Add all changes
   * - ::

           $ git commit -m "The new changes are identified by this sentence"

     - Commit the changes
   * - ::

           $ git push

     - Push the changes to the origin (will have to enter username and password for the origin)
     
     
How to commit three folders to git and push to respository?
===========================================================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           $ git add 0 constant system

     - Add all changes
   * - ::

           $ git commit -m "The new changes are identified by this sentence"

     - Commit the changes
   * - ::

           $ git push -u origin master

     - Push the changes to the origin (will have to enter username and password for the origin)
          
     
How to change a project from private to public in gitlab?
=========================================================
     
Projects > Your Projects > (click on project) 

Settings > Permissions (expand) > Project visibility = Private
     
How to create a new repository
==============================

Sign into github as en9apr:

https://github.com/

Click New button

Public repository

Create repository

cd existing_folder

echo "# darcyFoam" >> README.md

git init

git add .

git commit -m "first commit"

git remote add origin https://github.com/en9apr/darcyFoam.git

git push -u origin master

username
password

If I change something?
======================

git add .

git commit -m "second commit"

In gitk: File > Reload

If I change something else?
===========================

git add .

git commit -m "third commit"

In gitk: File > Reload

Use git k to revert commits (you can see the tree)
==================================================

Right click > revert commit

(can only take you one step backwards)

how to remove reverts?
======================

hard reset. 

In git k: Right click > Reset master branch to here (select Hard reset)

In git k: File > Reload (should remove all commits)

How to create a branch?
=======================

In git k: right click > new branch > name_of_branch

How to delete a branch?
=======================

git branch -d name_of_branch

How to view branches?
=====================

git branch

How to branch off to previous state?
====================================

Create a branch where the state is - In git k: right click > new branch > name_of_branch

git checkout name_of_branch

can still make commits to this branch:

git add .

git commit -m "third commit"

In gitk: File > Reload

How to go back to master branch?
================================

git checkout master     
     
     
     
     
     
     
