import numpy as np

def rotate(vec, theta=0.0, phi=0.0, alpha=0.0):
    theta = np.radians(theta)
    phi = np.radians(phi)
    alpha = np.radians(alpha)
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(alpha), -np.sin(alpha)],
                   [0, np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                   [0, 1, 0],
                   [-np.sin(theta), 0, np.cos(theta)]])
    Rz = np.array([[np.cos(phi), -np.sin(phi), 0],
                   [np.sin(phi), np.cos(phi), 0],
                   [0, 0, 1]])
    rotated_vec = np.dot(Rz, np.dot(Ry, np.dot(Rx, vec)))
    return rotated_vec

if __name__ == '__main__':
    Ex = np.array([1.0, 1.0, 0.0])
    Ex /= np.linalg.norm(Ex)
    Ey = np.array([0.0, 1.0, 1.0])
    Ey /= np.linalg.norm(Ey)

    # angle between Ex and Ey (should be 60 deg)
    dotproduct = np.dot(Ex, Ey)
    print(dotproduct)
    angle_rad = np.arccos(dotproduct)
    angle_deg = np.rad2deg(angle_rad)
    print(angle_rad, angle_deg)

    # cross product
    cross = np.cross(Ex, Ey)
    cross /= np.linalg.norm(cross)
    print(cross)

    #angle betweem cros and z
    z_vector = np.array([0, 0, 1])
    dot2 = np.dot(cross, z_vector)
    angle_rad2 = np.arccos(dot2)
    angle_deg2 = np.rad2deg(angle_rad2)
    print(f'theta = {angle_rad2} rad,  {angle_deg2} deg, {angle_rad2 / np.pi}*pi rad')


    ks = rotate(z_vector, phi=-45, theta=angle_deg2)
    print(ks, np.linalg.norm(ks))

    kp = np.array([0.5, 0, 0.5])
    kp /= np.linalg.norm(kp)
    print(kp)

    #angle between kp and ks
    dot3 = np.dot(kp,ks)
    angle_rad3 = np.arccos(dot3)
    angle_deg3 = np.rad2deg(angle_rad3)
    print(angle_rad3, angle_deg3)

    #By ģeometrija
    print('By geometry')

    Ex = np.array([1.0, 0.0, 1.0])
    Ex /= np.linalg.norm(Ex)
    Ey = np.array([1.0, 1.0, 0.0])
    Ey /= np.linalg.norm(Ey)

    # angle between Ex and Ey (should be 60 deg)
    dotproduct = np.dot(Ex, Ey)
    print(dotproduct)
    angle_rad = np.arccos(dotproduct)
    angle_deg = np.rad2deg(angle_rad)
    print(angle_rad, angle_deg)

    # cross product
    cross = np.cross(Ex, Ey)
    cross /= np.linalg.norm(cross)
    print(cross)

    # angle betweem cros and z
    z_vector = np.array([0, 0, 1])
    dot2 = np.dot(cross, z_vector)
    angle_rad2 = np.arccos(dot2)
    angle_deg2 = np.rad2deg(angle_rad2)
    print(f'theta = {angle_rad2} rad,  {angle_deg2} deg, {angle_rad2 / np.pi}*pi rad')

    ks = rotate(z_vector, phi=90+45, theta=angle_deg2)
    print('ks = ',ks, np.linalg.norm(ks))

    kp = np.array([0., 1., 1.])
    kp /= np.linalg.norm(kp)
    print(kp)

    # angle between kp and ks
    dot3 = np.dot(kp, ks)
    angle_rad3 = np.arccos(dot3)
    angle_deg3 = np.rad2deg(angle_rad3)
    print(angle_rad3, angle_deg3)



