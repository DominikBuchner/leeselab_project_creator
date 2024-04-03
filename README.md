# leeselab_project_creator
 Python script to generate plate layout lists as well as worklist for laboratory projects in the Aquatic Ecosystem Research Group.
 Automatically suggests a pooling strategy as well as optimal primer order for library generation.

## How to install and run
Installation via pip:

`pip install leeselab_project_creator`

How to update:

`pip install --upgrade leeselab_project_creator`

How to run:

`leeselab_project_creator` or `python -m leeselab_project_creator`

## How to use
After starting the GUI will look like this:  
![alt text](image-1.png)

The project creator needs an Excel list with sample names as input. 
You can set the Path to the sample list, the sheet in which it can be found as well as the cell of the first sample.
Select between 1 and 500 extraction replicates.
Select between 1 and 500 PCR replicates.
Select what primers to use, seperated by comma (e.g. fwh2, rbcl, 16S).
Select which 1st step primer tags are available for the primer sets, seperated by comma (e.g. 1,2,3,4,5,6,7,8).
Enter a project name and select a folder where to save the plate layout.
Hit "Generate sample list".

## Example
As an example sample list a list with 12 samples was selected:  
![alt text](image.png)  

The example setup looks like this when it is completly filled:  
![alt text](image-2.png)  

The project creator will generate a plate layout according to the settings.
This will include the technical replicates to control for cross-contamination as well as the selected amount of extraction replicates (indicated by _1, _2 ...)

![alt text](image-3.png)

Two more sheet will be generated, the extraction worklist and the PCR worklist:
![alt text](image-4.png)

In the PCR worklist a pooling strategy for the librarys will be calculated as well as an optimal use of primers to maximize library diversity.
The selected PCR replicates will be generated as seperate plates in the PCR worklist.
![alt text](image-5.png)

All sheets can be printed and taken into the lab to document all required data to sufficiently describe the analysis in the end.

