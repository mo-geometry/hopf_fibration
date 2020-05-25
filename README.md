# hopf fibration

Python scripts to animate the hopf fibration using matplotlib.

![plot3d](https://user-images.githubusercontent.com/62537514/79318089-7b152e00-7efe-11ea-83f9-4e4e50f67e39.png)

Unfortunately the 3d animation isn't the best as mutually overlapping objects are not handled by matplotlib (v2.1.1).
Lines & points are plotted in layers according to programming order - display graphics is a series of 2d layers which replicate a 3d renderer. For example

matlab - python comparison:

![exampleError](https://user-images.githubusercontent.com/62537514/79315662-2b813300-7efb-11ea-829e-c57a19ef4efe.png)
