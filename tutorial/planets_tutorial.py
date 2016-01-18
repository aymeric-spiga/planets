
# coding: utf-8

# # PLANETS tutorial
# 
# *Author: [Aymeric SPIGA](http://www.lmd.jussieu.fr/~aslmd)*
# 
# This tutorial explains how to use the `planets` utilities written in `python`. This package would be useful to anyone interested in an easy way to get and manipulate planetary constants in `python` scripts and interactive use with `ipython`. Note that the `planets` package is inspired from one of the package in R. Pierrehumbert `python` suite in support of the book *Principles of Planetary Climates*.

# <small>NB: If you have `ipython notebook` installed on your computer, download and use this tutorial interactively with the command `ipython notebook planets_tutorial.ipynb`</small>

# First and foremost, import the `planets` package by the command

# In[1]:

import planets


# ## Side use: physical constants

# It is always useful to have a handful of physical constants around: this is not the goal of `planets` - still you will appreciate this possibility in further calculations.

# In[2]:

print planets.h # Planck's constant
print planets.c # Speed of light


# In[3]:

print planets.k # Boltzman thermodynamic constant
print planets.sigma # Boltzman thermodynamic constant
print planets.G # Gravitational constant
print planets.N_avogadro # Avogadro's number
print planets.Rstarkilo # Universal gas constant


# ## Main use: planetary constants and calculations 

# It is pretty straightforward to import planetary constants for any planet. Let us choose Venus as an example.

# In[4]:

myplanet = planets.Venus


# <small> Note that the command above + the `import` line can be equivalently obtained with the pretty intuitive alternative
# 
# ~~~python
# from planets import Venus as myplanet
# ~~~
# 
# Nevertheless this method is less flexible than the one indicated above, in which you have access to the whole `planets` variables and functions.</small>
# 
# OK, now you have loaded a `Planet` python object named `myplanet` and its attributes are set to the preset Venus values. The list of parameters in each `Planet` object can be obtained with the `show()` method.

# In[5]:

myplanet.show()


# Now each of those parameter is easy to get from the attribute variable indicated in the left column. For instance, get the mean radius of the planet.

# In[6]:

print myplanet.a


# If you are unsure about the units, you can still get the description by typing

# In[7]:

print planets.desc["a"]


# So now it is fairly easy to work out a calculation ! One example. Compute the equivalent temperature
# $$T_{eq} = \left( \frac{(1-A)F}{4\sigma} \right)^{\frac{1}{4}}$$

# In[8]:

teq = ( ( (1-myplanet.albedo)*myplanet.L ) / (4*planets.sigma) )**0.25
print teq


# Actually this is already coded in `planets` because it is a quantity planetary scientists often use. Thus to compute equivalent temperature, you can simply use the method `eqtemp()`

# In[9]:

teq = myplanet.eqtemp()
print teq


# The additional cool point of setting up a method is that it is easy to display values for a bunch of preset planets

# In[10]:

print planets.Earth.eqtemp()
print planets.Mars.eqtemp()
print planets.Venus.eqtemp()
print planets.Saturn.eqtemp()


# New methods are often added to `planets`. Check for the source for this. Other examples of methods available at the time of writing are

# In[11]:

print myplanet.H() # atmospheric scale height
print myplanet.N2() # Brunt-Väisälä frequency


# Those methods have additional arguments. For instance you can change the temperature (in Kelvin) for the atmospheric scale height

# In[14]:

print myplanet.H(T0=273.)


# This works also with an array of temperature values `temp`. Hence `myplanet.H(T0=temp)` will be an array of scale height values with the same dimensions as `temp`.

# ## Customize planetary constants

# The values of preset planetary constants can be modified (because, well, those are constants to a certain point, or measured with a given accuracy, or subject to change with new missions) in the text files located in the `planet` folder.
# 
# The use of text files to set planetary constants allows for a very flexible use. You might even set your own `Planet` object with your user-defined planetary constants, which a pretty useful feature in an universe in which we get to discover more and more exoplanets. Here is how to do that. Let us assume you want to create a set of planetary constant for HD189733b
# 
#     1. Open a file HD189733b.txt
#     2. Fill it with constants (look the existing files in the `planet` folder for templates)
#     3. Put this file in the planet folder (or have it where you work out your computations)
#     
# The new planet is now available in `planets`! Simply set a `Planet` object (name it the way you want) and load the values for the chosen planet with the `ini()` method. Then you can use all the available variables and methods associated to the object.
# 
# ~~~python
# a_hot_jupiter = Planet()
# a_hot_jupiter.ini("HD189733b")
# print a_hot_jupiter.L
# print a_hot_jupiter.eqtemp()
# ~~~
# 
# Note that if one particular planetary constant is missing in the text file, it will be set to `None` by default.

# ## Concluding remark

# `planets` is still a work in progress. Should you have bug fixes - or useful txt files for missing planetary environments, please feel free to drop me an email.
