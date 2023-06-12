import numpy as np
import pandas as pd
import cv2
from matplotlib import pyplot as plt
import csv
import math
import os
import glob
import sys
import re

dir = '/Users/konrai/Dropbox/Quercus/Quercus_画像/調査-鳥取大演習林/original'
dir_csv = '/Users/konrai/Dropbox/Quercus/Quercus_全体共有/調査-鳥取大演習林/EFD/csv'
#os.chdir('/Users/konrai/Dropbox/Konrai_exp/konara_mizunara/Threshold')
os.chdir(dir)
### Select image ###
thresh = 170 # Threshold value (arbitrary value)
files = glob.glob('*.jpg') # All scan image ara loaded
#files_compare = glob.glob('../csv/*.csv')
files_compare = glob.glob(dir_csv+'/*.csv')
files_compare = [s.replace('csv', 'jpg') for s in files_compare]
files_compare = [s.split('/')[-1] for s in files_compare]

print('Scan list:')
for x in files:
    if x in files_compare:
        print(x, 'Done')
    else:
        print(x, 'Undone')


file_src = input('Which image?:')

### Get Contour Information ###
#file_src = '/Users/konrai/Dropbox/Konrai_exp/konara_mizunara/Threshold/449-konara-img023.jpg'
file_name = file_src.split('/')[-1] # Get file name
file_name = file_name.split('.')[0] # remove Filename extention (拡張子)

img_src = cv2.imread(file_src, 0) # load as gray scale
img_src = img_src[300:img_src.shape[0]-300, 300:img_src.shape[1]-300]
height = img_src.shape[0]
width = img_src.shape[1]
img_blank = np.zeros((height, width, 3), np.uint8)
#plt.imshow(img_src, cmap='gray', vmin=0, vmax=255)
#plt.show()

# Rabelling
threshold_area = 100 # Threshold area (removing noizes)
l_info = [['label', 'area', 'perimeter', 'roundness', 'angle']]
nlabel, img_lab = cv2.connectedComponents(image = img_src)

for n in range(1, nlabel):
    img_tmp = cv2.compare(img_lab, n, cv2.CMP_EQ) # Extract image which is labeled n
    # Get contours coordinate
    contours_tmp = cv2.findContours(img_tmp, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[1]
    area = cv2.contourArea(contours_tmp[0]) # area
    perimeter = cv2.arcLength(np.array(contours_tmp[0]), True)
    if area != 0 and perimeter != 0:
        roundness = 4 * np.pi * area / perimeter / perimeter
    else:
        roundness = 0

    if threshold_area < area:
        id = len(l_info)
        # Center of gravity, COG
        m = cv2.moments(img_tmp)
        x_g = m['m10'] / m['m00']
        y_g = m['m01'] / m['m00']
        ang = 0.5 * math.atan2(2.0 * m['mu11'], m['mu20'] - m['mu02'])
        img_col = cv2.cvtColor(img_tmp, cv2.COLOR_GRAY2RGB)
        img_blank = cv2.add(img_blank, img_col)

        r = 100
        cv2.circle(img=img_blank, center=(int(x_g), int(y_g)), radius=10, color=(156,200,115), thickness=-1) # plot COG
        cv2.line(img=img_blank, pt1=(int(x_g+r*math.cos(ang)), int(y_g+r*math.sin(ang))), pt2=(int(x_g-r*math.cos(ang)), int(y_g-r*math.sin(ang))), color=(156,200,115), thickness=5) # Major axis
        cv2.line(img=img_blank, pt1=(int(x_g+r*math.cos(math.pi/2+ang)), int(y_g+r*math.sin(math.pi/2+ang))), pt2=(int(x_g-r*math.cos(math.pi/2+ang)), int(y_g-r*math.sin(math.pi/2+ang))), color=(156,200,115), thickness=5) # Minor axis

        ### Rotate as to normalize angle (Get img_rot)
        # Bounding rencangle (外接長方形)
        x, y, w, h = cv2.boundingRect(img_tmp)
        cv2.rectangle(img_blank, (x,y), (x+w,y+h), color=(85,91,241), thickness=5)

        img_crp = img_tmp[y:y+h, x:x+w] # Extract ROI
        cv2.putText(img_blank, text=str(id), org=(int(x), int(y-15)), fontFace=cv2.FONT_HERSHEY_SIMPLEX, fontScale=5, color=(85,91,241), thickness=10, lineType=cv2.LINE_AA) # Add ID label

        # Calculate image size after rotaiton
        w_rot = int(np.round(h*np.absolute(np.sin(ang+math.radians(10)))+w*np.absolute(np.cos(ang+math.radians(10)))))
        h_rot = int(np.round(h*np.absolute(np.cos(ang+math.radians(10)))+w*np.absolute(np.sin(ang+math.radians(10)))))
        size_rot = (w_rot, h_rot)
        margin = 500
        img_blk_rot = np.zeros((h_rot+margin, w_rot+margin), np.uint8) # image will be pasted on this blank image. Because contour must be closed, margins are neede.

        # Rotate around center of the original image
        center = (w/2, h/2)
        scale = 1.0
        rotation_matrix = cv2.getRotationMatrix2D(center, math.degrees(ang)+10, scale)

        # Add parallel shift (rotation + translation)
        affine_matrix = rotation_matrix.copy()
        affine_matrix[0][2] = affine_matrix[0][2] -w/2 + w_rot/2
        affine_matrix[1][2] = affine_matrix[1][2] -h/2 + h_rot/2
        img_rot = cv2.warpAffine(img_crp, affine_matrix, size_rot, flags=cv2.INTER_CUBIC)

        img_blk_rot[int(margin/2):h_rot+int(margin/2), int(margin/2):w_rot+int(margin/2)] = img_rot

        plt.imshow(img_blk_rot, cmap='gray', vmin=0, vmax=255)
        #plt.colorbar()
        plt.show()

        dir = input('Is direction OK? Correct: 1, Opposite: 2 :')
        dir = int(dir)

        if dir ==2:
            center = (w/2, h/2)
            rotation_matrix = cv2.getRotationMatrix2D(center, math.degrees(ang)+190, scale)

            affine_matrix = rotation_matrix.copy()
            affine_matrix[0][2] = affine_matrix[0][2] -w/2 + w_rot/2
            affine_matrix[1][2] = affine_matrix[1][2] -h/2 + h_rot/2
            img_rot = cv2.warpAffine(img_crp, affine_matrix, size_rot, flags=cv2.INTER_CUBIC)
            img_blk_rot[int(margin/2):h_rot+int(margin/2), int(margin/2):w_rot+int(margin/2)] = img_rot

        else:
            print('Direction is correct!')


        nlabel_rot, img_lab_rot = cv2.connectedComponents(image = img_blk_rot)
        for i in range(1, nlabel_rot):
            img_tmp_rot = cv2.compare(img_lab_rot, i, cv2.CMP_EQ) # Extract image which is labeled n
            # Get contours coordinate
            contours_tmp_rot = cv2.findContours(img_tmp_rot, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[1]
            area = cv2.contourArea(contours_tmp_rot[0]) # area
            perimeter = cv2.arcLength(np.array(contours_tmp[0]), True)
            if area != 0 and perimeter != 0:
                roundness = 4 * np.pi * area / perimeter / perimeter
            else:
                roundness = 0

            if threshold_area < area:
                # Get information of contours
                contours = cv2.findContours(img_tmp_rot, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[1]
                area = cv2.contourArea(contours[0]) # area
                perimeter = cv2.arcLength(np.array(contours[0]), True)

                print('ID:', id)
                print('Label number:', n)
                print('Area:', area)
                print('Perimeter:', perimeter)
                print('Roundness:', roundness)
                print('angle:', math.degrees(ang), '(deg)')

                l_tmp = [id, area, perimeter, roundness, math.degrees(ang)]
                l_info.append(l_tmp)

                # Save images
                #cv2.imwrite(f"/Users/konrai/Dropbox/Konrai_exp/konara_mizunara/label/{file_name}_{id}_original.jpg", img_tmp)
                cv2.imwrite(f"/Users/konrai/Dropbox/Quercus/Quercus_画像/調査-鳥取大演習林/label/{file_name}_{id}_original.jpg", img_tmp)
                cv2.imwrite(f"/Users/konrai/Dropbox/Quercus/Quercus_画像/調査-鳥取大演習林/label/{file_name}_{id}_crop.jpg", img_crp)
                cv2.imwrite(f"/Users/konrai/Dropbox/Quercus/Quercus_画像/調査-鳥取大演習林/label/{file_name}_{id}_normalized.jpg", img_tmp_rot)

                l_contours = list()
                for i in range(len(contours[0][:])):
                    #print(contours[0][i][0])
                    l_contours.append(contours[0][i][0])
                    with open(f"/Users/konrai/Dropbox/Quercus/Quercus_画像/調査-鳥取大演習林/contour/contour_{file_name}_{id}.csv", 'w') as f:
                        writer = csv.writer(f) # Save as csv file
                        writer.writerows(l_contours)

df = pd.DataFrame(l_info)
df.to_csv(f"/Users/konrai/Dropbox/Quercus/Quercus_全体共有/調査-鳥取大演習林/EFD/csv/{file_name}.csv")
print('Dataframe is made as below;\n', df)
cv2.imwrite(f"/Users/konrai/Dropbox/Quercus/Quercus_画像/調査-鳥取大演習林/label/{file_name}.jpg", img_blank)
