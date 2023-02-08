###################################
*RNAlysis* quick-start video guide
###################################


Open the program
=================
If you installed *RNAlysis* from PyPi, you can open the graphical user interface by either executing the command `rnalysis-gui` from your terminal, or by typing the following code into a Python console::

    >>> from rnalysis import gui
    >>> gui.run_gui()

Load a table
============
Choose a file from your computer, and click 'start' to load it into *RNAlysis*.

.. image:: ../../rnalysis/gui/videos/load_table.webp

Examine your table
==================
You will now be able to see an overview of your data, including the table's name, type, and dimensions. Click the 'View full table' button to see your table in its entirety.

.. image:: ../../rnalysis/gui/videos/view_table.webp

Filter your table
=================
Choose a filtering function, set your desired parameters, and click 'Apply' to filter your table. The changes you make will not affect your original file until you save them.

.. image:: ../../rnalysis/gui/videos/filter_table.webp

Undo the operations you applied to your table
=============================================
At any moment, you can use the 'Command history' window to undo or redo an operation you applied to your table.

.. image:: ../../rnalysis/gui/videos/undo_actions.webp

Apply your operations 'in-place' or apply on a new table
========================================================
Instead of applying operations 'in-place', you can choose to apply the operation to a copy of your table in a new tab. The table in the original tab won't be modified.

.. image:: ../../rnalysis/gui/videos/apply_inplace.webp

Save the changes you made to your table
=======================================
To save result of your filtering operations, click the 'Save table' button and choose where to save the modified table.

.. image:: ../../rnalysis/gui/videos/save_table.webp

Work on multiple tables at the same time
========================================
You can work on multiple tables at the same time by opening a new tab and loading another table.

.. image:: ../../rnalysis/gui/videos/new_tab.webp

Different types of tables offer different ways to filter and analyze your data
==============================================================================
When loading a table, you can specify its type. Different types of tables support different types of functions: for example, count matrices support clustering analysis.

.. image:: ../../rnalysis/gui/videos/table_types.webp

Create and save graphs
======================
Some functions can generate graphs of your data. You can resize those graphs, and save them to your computer in multiple file formats.

.. image:: ../../rnalysis/gui/videos/generate_graphs.webp

Quickly look-up genes in your database of choice
==================================================
Easily access information about your genes with a simple right-click.
Select from a range of biological databases such as GeneCards, NCBI Genes, UniProtKB and others,
which can be configured in the settings menu to fit your specific needs.
.. image:: ../../rnalysis/gui/videos/quick_search.webp

Sort your tabs and change their icons
=====================================
To help organize your workspace, you can sort tabs by right-clicking a tab and choosing a sorting method. You can also change specific tabs' colors, to help you differentiate them.

.. image:: ../../rnalysis/gui/videos/sort_tabs.webp

Restore tabs you closed
=======================
If you accidentally closed one of your tabs - don't worry! You can restore closed tabs through the 'Edit' menu.

.. image:: ../../rnalysis/gui/videos/restore_tabs.webp

Import lists of genes as Gene Sets
==================================
In addition to tables, *RNAlysis* can also import lists of genes as Gene Sets. We will soon review what we can do with those gene sets.

.. image:: ../../rnalysis/gui/videos/import_gene_sets.webp

Visualize the intersections between your tables and gene sets
=============================================================
In the 'Visualize Gene Sets' window you can create Venn diagrams and UpSet plots that will display the various intersections between your tables and gene sets.

.. image:: ../../rnalysis/gui/videos/visualize_gene_sets.webp

Apply set operations to your tables and gene sets
=================================================
In the 'Set Operations' window you can extract specific subsets from your data. Either use predefined set operations, or click on specific subsets in the preview pane to select them.

.. image:: ../../rnalysis/gui/videos/set_operations.webp

Perform enrichment analysis on your tables and gene sets
========================================================
In the 'Enrichment Analysis' window, you can perform various types of enrichment analysis on the tables and gene sets you filtered.

.. image:: ../../rnalysis/gui/videos/enrichment_analysis.webp

Create Pipelines to streamline your data analysis
=================================================
You can group multiple operations in a specific order and with specific parameters into a Pipeline. Just add those functions to the Pipeline in the order you choose.

.. image:: ../../rnalysis/gui/videos/create_pipeline.webp

Apply Pipelines to one or more of your tables
=============================================
You can apply a Pipeline to a group of tables through the 'Pipelines' menu. Using Pipelines to analyze multiple datasets can make your workflow faster and less error-prone.

.. image:: ../../rnalysis/gui/videos/apply_pipeline.webp

Export and share Pipelines to make your analysis more reproducible
==================================================================
Pipelines you export can be imported from any computer, and can be shared with others to help make your analysis easier to understand and more reproducible.

.. image:: ../../rnalysis/gui/videos/export_pipeline.webp
