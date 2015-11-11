function mat = angles2rotMat(angles)
angle_x = angles(1);
angle_y = angles(2);
angle_z = angles(3);

    A       = cos(angle_x);
    B       = sin(angle_x);
    C       = cos(angle_y);
    D       = sin(angle_y);
    E       = cos(angle_z);
    F       = sin(angle_z);
    AD      =   A * D;
    BD      =   B * D;
    mat(1,1)  =   C * E;
    mat(1,2)  =  -C * F;
    mat(1,3)  =   D;

    mat(2,1)  =  BD * E + A * F;
    mat(2,2)  = -BD * F + A * E;
    mat(2,3)  =  -B * C;
    
    mat(3,1)  = -AD * E + B * F;
    mat(3,2) =  AD * F + B * E;
    mat(3,3) =   A * C;
    
    