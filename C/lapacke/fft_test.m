x = (0:3);
y = x;
z = x;

[y_3D,x_3D,z_3D] = meshgrid(x,y,z);
o = fftn(x_3D);