import numpy as np

def rotate(vec, theta=0.0, phi=0.0, alpha=0.0):
    '''

    :param vec:
    :param theta: angle in deg
    :param phi: angle in deg
    :param alpha: angle in deg
    :return:
    '''
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

def rotate_about_axis(vector, axis, angle):
    '''

    :param vector:
    :param axis:
    :param angle: angle in deg
    :return:
    '''
    angle = np.radians(angle)
    N_mat = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
    rot_mat = np.eye(3) + np.sin(angle) * N_mat + (1.0 - np.cos(angle)) * np.dot(N_mat, N_mat)
    rotated_vector = np.dot(rot_mat, vector)
    return rotated_vector

if __name__ == '__main__':
    Ex = np.array([1.0, 1.0, 0.0])
    Ex /= np.linalg.norm(Ex)
    Ey = np.array([0.0, 1.0, 1.0])
    Ey /= np.linalg.norm(Ey)

    # angle between Ex and Ey (should be 60 deg)
    dotproduct = np.dot(Ex, Ey)
    print('angle between Ex and Ey (should be 60 deg)')
    print(f'dotproduct={dotproduct}')
    angle_rad = np.arccos(dotproduct)
    angle_deg = np.rad2deg(angle_rad)
    print(f'angle = {angle_rad} rad, angle = {angle_deg} deg')

    # cross product
    print('cross product np.cross(Ex, Ey)')
    cross = np.cross(Ex, Ey)
    cross /= np.linalg.norm(cross)
    print(f'ks = {cross}')

    #angle between cros and z
    print('angle between ks and z')
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
    print('angle between kp and ks')
    dot3 = np.dot(kp, ks)
    angle_rad3 = np.arccos(dot3)
    angle_deg3 = np.rad2deg(angle_rad3)
    print(f'{angle_rad3} rad, {angle_deg3} deg')

    #perpendicular polarization
    print('Ex perpendicular')
    # cross2 = np.cross(Ex, ks)
    cross2 = np.cross(ks, Ex)
    cross2 /= np.linalg.norm(cross2)
    print(f'Ex perp = {cross2}')

    print('angle between Ex and Ex perp')
    dot2 = np.dot(cross2, Ex)
    angle_rad2 = np.arccos(dot2)
    angle_deg2 = np.rad2deg(angle_rad2)
    print(f'theta = {angle_rad2} rad,  {angle_deg2} deg, {angle_rad2 / np.pi}*pi rad')

    print('angle between Ex perp and z')
    dot2 = np.dot(cross2, z_vector)
    angle_rad2 = np.arccos(dot2)
    angle_deg2 = np.rad2deg(angle_rad2)
    print(f'theta = {angle_rad2} rad,  {angle_deg2} deg, {angle_rad2 / np.pi}*pi rad')
    Ex2 = rotate(z_vector, phi=-45, theta=angle_deg2)
    print(Ex2, np.linalg.norm(Ex2))

    # rotate about axis
    print('rotate Ex about ks')
    x_vector = np.array([1, 0, 0])
    print(f'x_vector = {x_vector[:2]}')
    for rot_angle in [0, 45, 60, 90]:
        print()
        print(f'Ex = {Ex}')
        print(f'ks = {ks}')
        Ex_rot = rotate_about_axis(Ex, ks, rot_angle)
        print(f'angle = {rot_angle} deg = {np.radians(rot_angle)} rad')
        print(f'Ex_rot = {Ex_rot}, {np.linalg.norm(Ex_rot)}')
        #calculate theta
        dot2 = np.dot(Ex_rot, z_vector)
        theta_rad2 = np.arccos(dot2)
        theta_deg2 = np.rad2deg(theta_rad2)
        print(f'theta = {theta_rad2} rad,  {theta_deg2} deg, {theta_rad2 / np.pi}*pi rad')
        #calculate phi
        dot_phi = np.dot(Ex_rot[:2]/np.linalg.norm(Ex_rot[:2]), x_vector[:2]/np.linalg.norm(x_vector[:2]))
        phi_rad2 = np.arccos(dot_phi)
        phi_deg2 = np.rad2deg(phi_rad2)
        print(f'phi = {phi_rad2} rad,  {phi_deg2} deg, {phi_rad2 / np.pi}*pi rad')

        Ex2 = rotate(z_vector, phi=phi_deg2, theta=theta_deg2)
        print(f'Ex2 = {Ex2}, {np.linalg.norm(Ex2)}')
        print(np.round(Ex_rot,5) == np.round(Ex2,5))

    # #By ģeometrija
    # print('By geometry')
    #
    # Ex = np.array([1.0, 0.0, 1.0])
    # Ex /= np.linalg.norm(Ex)
    # Ey = np.array([1.0, 1.0, 0.0])
    # Ey /= np.linalg.norm(Ey)
    #
    # # angle between Ex and Ey (should be 60 deg)
    # dotproduct = np.dot(Ex, Ey)
    # print(dotproduct)
    # angle_rad = np.arccos(dotproduct)
    # angle_deg = np.rad2deg(angle_rad)
    # print(angle_rad, angle_deg)
    #
    # # cross product
    # cross = np.cross(Ex, Ey)
    # cross /= np.linalg.norm(cross)
    # print(cross)
    #
    # # angle betweem cros and z
    # z_vector = np.array([0, 0, 1])
    # dot2 = np.dot(cross, z_vector)
    # angle_rad2 = np.arccos(dot2)
    # angle_deg2 = np.rad2deg(angle_rad2)
    # print(f'theta = {angle_rad2} rad,  {angle_deg2} deg, {angle_rad2 / np.pi}*pi rad')
    #
    # ks = rotate(z_vector, phi=90+45, theta=angle_deg2)
    # print('ks = ',ks, np.linalg.norm(ks))
    #
    # kp = np.array([0., 1., 1.])
    # kp /= np.linalg.norm(kp)
    # print(kp)
    #
    # # angle between kp and ks
    # dot3 = np.dot(kp, ks)
    # angle_rad3 = np.arccos(dot3)
    # angle_deg3 = np.rad2deg(angle_rad3)
    # print(angle_rad3, angle_deg3)
    #


