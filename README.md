MAST-U Shot Design Tool
=======================

To install and run the tool, you will need the following:
    - Python (Version 3.9 or Higher)


Quickstart Guide
----------------

On Windows 10 or 11:

    - Make sure you have Python version 3.9 to 3.12. Press the Windows key and type Python - it it's installed it will suggest Apps named by version number.

    - On this webpage (https://github.com/cboswel/musd/) click the green Code button near the top right of the page. Select Download as Zip. Save it somewhere sensible then extract it.

    - Use file manager to navigate to the extracted folder. Try to find install.bat and double click it. This should install lots of Python libraries then run the tool. The run.bat file can be used if you use install.bat first, exit out, and then want to open it again without installing all the libraries from scratch.


Installation Notes
------------------

Just for fun, we could do all of the above from the command line. Linux and Mac users will need to do it this way too because the scripts were only made to work for Windows. Using the command line opens up access to lots of powerful tools that will become very useful in your computing career. So give this a try first, but it should be possible to just double click the right file as mentioned later.

Let's start by bringing up a terminal. Command prompt or powershell should work fine for Windows. You can access these by searching in the start menu or by right clicking in a folder you want to work in and selecting "Open in Terminal". On Windows 10, you need to hold the shift key as you right click and navigate that menu.

Check which version of Python you have installed by typing `python --version` into the command prompt and pressing enter. It should be between 3.9 and 3.12. If it isn't, download and install a version in that range.

Now we need to download all our code. Click the green "code" button on the top right of this git repository and choose Download Zip. Extract it. Alternatively, paste the following command into the terminal to download the tool into the current directory (You might need to install git for this method):

git clone https://github.com/CharlieBoswell/musd.git

In your terminal, navigate to the extracted folder. You can do this by writing `cd <name of directory you want to go into>` into the command prompt and pressing enter. After writing `cd `, you can press tab to see all the possible things you could write. If you need to go back you can use `cd ../`. To check all the files in the current directory, use `ls`. Alternatively, you can get there from File Explorer and right clicking as mentioned earlier.

If you're using Windows, you need to run the install.bat script like this: `.\install.bat`. On Linux, `./install.sh`. This will create something called a virtual environment in Python so we can keep everything in the project contained in one corner. Then all of the libraries needed for the tool are downloaded and installed.

Running
-------

We have to run it from inside our virtual environment. If you just set up the environment you probably don't need to do anything, but you can get back into an environment with `.\venv\Scripts\activate` on Windows or `source venv/bin/activate` on Linux. You could even just run the install script again.

To boot up the tool, execute the run script with .\run.bat on Windows or ./run.sh on Linux. Have fun exploring!
