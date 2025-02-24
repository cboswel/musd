MAST-U Shot Design Tool
=======================

To install and run the tool, you will need the following:
    - Python (Version 3.9 or Higher)

Installation
------------

First of all, bring up a terminal. Command prompt or powershell should work fine for Windows. You can access these by searching in the start menu or by right clicking in a folder you want to work in and selecting "Open in Terminal".

Check which version of Python you have installed by typing `python --version` into the command prompt and pressing enter. It should be between 3.9 and 3.12. If it isn't, download and install a version in that range.

Now we need to download all our code. Click the green "code" button on the top right of this git repository and choose Download Zip. Extract it. Alternatively, paste the following command into the terminal to download the tool into the current directory (You might need to install git for this method):

git clone https://github.com/CharlieBoswell/musd.git

In your terminal, navigate to the extracted folder. You can do this with `cd <name of directory you want to go into>`. Press tab to see all the possibilities. If you need to go back you can use `cd ../`. To check all the files in the current directory, use `ls`. Alternatively, you can get there from File Explorer and right clicking as mentioned earlier.

If you're using Windows, you need to run the install.bat script like this: `.\install.bat`. On Linux, `./install.sh`. This will create something called a virtual environment in Python so we can keep everything in the project contained in one corner. Then all of the libraries needed for the tool are downloaded and installed.

Running
-------

We have to run it from inside our virtual environment. If you just set up the environment you probably don't need to do anything, but you can get back into an environment with `.\venv\Scripts\activate` on Windows or `source venv/bin/activate` on Linux. You could even just run the install script again.

To boot up the tool, execute the run script with .\run.bat on Windows or ./run.sh on Linux. Have fun exploring!
