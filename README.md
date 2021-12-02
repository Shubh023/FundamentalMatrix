# VISION - TP Fundamental Initial

### **Student :** *PATEL Shubhamkumar at M2 IMA*


This program is meant to take as input two images and calculated a Fundamental matrix (the best we can obtain using RANSAC algorithm).

## **Build the project** :

#### *Generate an executable binary Fundamental for your distribution* :

- ***``` mkdir build ```***

- ***``` cd build ```***

- ***``` cmake ../CMakeLists.txt -B . ```***

- ***``` make ```***

- ***``` ./Fundamental ```***

**PS** : *For my local setup i had to add ```set(OpenGL_GL_PREFERENCE GLVND)``` in the **CMakeLists.txt**. So in case there is problem related to this setting please remove this line and try to build with the same process as above.*

## **Calculating the fundamental Matrix** :
#### Once the matches are found using the SIFT, you need to click on any of the two images to calculate the fundamental matrix. In the termainal/Console the F matrix will be printed. Once that is done click once again on any of the images. If you get less inlier just re-execute the program until you get enough inliers.

#### When you will click on the screen points will appear either in BLUE or in RED according to the image you interacted with. And a line will be drawn to display the result of the matching point to image you clicked on.
#### When you click on the right image points will appear in RED on the right image and corresponding RED line will be drawn on the left image. And similarliy the points and lines will appear BLUE if the interaction was on the left screen.

#### I've added multiple results in the ***Results*** directory with files **3-Result.png** and **4-Result.png** which show all points and the corresponding lines that i was able to draw using the calculated Fundamental Matrix. 

