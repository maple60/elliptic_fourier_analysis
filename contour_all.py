#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 07:52:32 2022

@author: konrai
今までは1枚ずつやっていたが、まとめて葉の輪郭座標を取得したい!
一気に全画像の輪郭を抽出する
"""

#%% import libraries
import numpy as np
import pandas as pd
import cv2
#from matplotlib import pyplot as plt
import csv
import math
import os
import glob

#%% define function
# pixelをmmに変換する関数
dpi = 300
def cvtmm(x, dpi=300): # xはピクセル数
    mm = x * 25.4 / dpi
    return mm

# 300ピクセルは1インチだから25.4mmになる
#cvtmm(300)

#%% directory setting
# dir = '/Users/konrai/Dropbox/Quercus' # original scan images
dir = '/Users/konrai/Library/CloudStorage/Dropbox/hawaii_leaf_personal' # original scan images

#site = '調査-宮城'
os.chdir(dir) # set working directory
# dir_csv = './Quercus_全体共有/' + site +'EFD/csv' # 使われていない変数
#files = glob.glob('./Quercus_画像/' + site + '/original/*.jpg') # All scan image are loaded
dir_img = '2014_leaf_scan' # フォルダ単位での処理
files = glob.glob('./' + dir_img + '/*.jpg') # All scan image are loaded

#%% main processing
thresh = 170 # Threshold value (arbitrary value)
margin = 10 # トリミング画像の四方の余白
r = 100 # 何を表すか忘れた
thresh_area = 300 # 面積の閾値 (単位は多分pixel?) ハワイフトモモだとこれの値が小さくなることに注意
thresh_area_up = 1000000 # もし上限を設定したかったらこのパラメータをいじる。
y_min = 150 # スケールを省く：上から指定することに注意

for file_src in files:
    # file_name = file_src.split(sep='/')[4].split(sep='.')[0] # ID(アルファベット+4ケタ数字)のみを取り出し
    file_name = file_src.split(sep='/')[2].split(sep='.')[0] # 拡張子を取り除いたファイル名を取得
    print(file_name)

    # グレースケールで画像を読み込む
    img_src = cv2.imread(file_src, 0)
    height = img_src.shape[0]
    width = img_src.shape[1]

    # スケールを省く
    img_src = img_src[y_min:height, 0:width]

    ret, im_bin = cv2.threshold(img_src, thresh, 255, cv2.THRESH_BINARY_INV) # binary image

    contours = cv2.findContours(im_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[1] # get contours
    area_list = list()
    for ctr in range(len(contours)):
        area = cv2.contourArea(contours[ctr])
        area_list.append(area)

    area_leaf = list(filter(lambda x: thresh_area_up > x > thresh_area, area_list))
    leaf_ID = list()

    for LA in area_leaf:
        leaf_ID.append(area_list.index(LA))
        print('ID:', area_list.index(LA), ', Area: ', cvtmm(LA), ' mm2')

    #plt.rcParams['figure.figsize'] = [20, 10]
    img_blank = np.zeros((height, width, 3), np.uint8)
    img_blank_1ch = np.zeros((height, width, 1), np.uint8)
    img_blank_3ch = np.zeros((height, width, 3), np.uint8)
    img_id_col_base = img_blank_3ch.copy()
    img_id_all = img_blank_3ch.copy()
    d_csv = [['ID1', 'ID2', 'area', 'perimeter', 'roundness', 'width', 'length']]
    img_lbl = cv2.imread(file_src, cv2.IMREAD_COLOR)
    img_lbl = img_lbl[y_min:height, 0:width]

    ID = 1
    for LID in leaf_ID:
        img_cnt = cv2.imread(file_src, cv2.IMREAD_COLOR)
        img_cnt = img_cnt[y_min:height, 0:width]
        cv2.drawContours(img_cnt, contours, LID, color=(255,0,0), thickness=10) # 元画像の上に赤色で輪郭を描画
        img_id = img_blank.copy()
        cv2.drawContours(img_id, contours, LID, color=(255,255,255), thickness=-1) # 黒のブランク画像の上に輪郭を白塗りつぶしで描画

        x, y, w, h = cv2.boundingRect(contours[LID])
        img_br = cv2.rectangle(img_lbl,(x,y),(x+w,y+h),(255,0,0),thickness=10)
        cv2.putText(img_br, text=str(leaf_ID.index(LID)+1), org=(int(x), int(y-15)), fontFace=cv2.FONT_HERSHEY_SIMPLEX, fontScale=5, color=(85,91,241), thickness=15, lineType=cv2.LINE_AA) # Add ID label
        save_name = './trim/' + dir_img + '/' + file_name # トリミングした画像を保存
        if LID == leaf_ID[len(leaf_ID)-1]:
            save_name_label = save_name + '_label.jpg'
            cv2.imwrite(save_name_label, img_br)

        # それぞれのラベル画像を1枚ずつトリミングして保存
        img_trm = cv2.imread(file_src, cv2.IMREAD_COLOR)
        img_trm = img_trm[y_min:height, 0:width]
        img_trm = img_trm[y-margin : y+h+margin, x-margin : x+w+margin] # 少し余白を取って保存
        #img_trm = img_trm[y-margin : y+h+margin, x: x+w+margin] # marginとって右側見切れているときはこっちを使う。
        save_name_trimming = save_name + '_' + str(ID) + '.jpg'
        cv2.imwrite(save_name_trimming, img_trm) # ここはplt.imsaveでも良い?

        ### 面積などの情報を取得する
        area = cvtmm(cv2.contourArea(contours[LID])) # 面積(mm2)
        perimeter = cvtmm(cv2.arcLength(np.array(contours[LID]), True)) # 周長(mm)
        roundness = 4 * np.pi * area / perimeter / perimeter # 真円度, 1に近いほど真円に近い

        img_id = img_blank_1ch.copy() # black background image
        cv2.drawContours(img_id, contours, LID, color=(255,255,255), thickness=-1) # filled contour image

        # Center of gravity, COG
        m = cv2.moments(img_id)
        x_g = m['m10'] / m['m00']
        y_g = m['m01'] / m['m00']
        ang = 0.5 * math.atan2(2.0 * m['mu11'], m['mu20'] - m['mu02'])

        # plot COG
        img_id_col = cv2.cvtColor(img_id, cv2.COLOR_GRAY2RGB) # 重心とか外接長方形はカラーで描きたいので、3チャンネルに変換
        img_id_col = cv2.add(img_id_col_base, img_id_col)
        cv2.circle(img=img_id_col, center=(int(x_g), int(y_g)), radius=10, color=(156,200,115), thickness=-1) # plot COG
        cv2.line(img=img_id_col, pt1=(int(x_g+r*math.cos(ang)), int(y_g+r*math.sin(ang))), pt2=(int(x_g-r*math.cos(ang)), int(y_g-r*math.sin(ang))), color=(156,200,115), thickness=5) # Major axis
        cv2.line(img=img_id_col, pt1=(int(x_g+r*math.cos(math.pi/2+ang)), int(y_g+r*math.sin(math.pi/2+ang))), pt2=(int(x_g-r*math.cos(math.pi/2+ang)), int(y_g-r*math.sin(math.pi/2+ang))), color=(156,200,115), thickness=5) # Minor axis

        #  Bounding rencangle (外接長方形)
        x, y, w, h = cv2.boundingRect(img_id)
        cv2.rectangle(img_id_col, (x,y), (x+w,y+h), color=(85,91,241), thickness=5)

        # Add label
        cv2.putText(img_id_col, text=str(ID), org=(int(x), int(y-15)), fontFace=cv2.FONT_HERSHEY_SIMPLEX, fontScale=5, color=(85,91,241), thickness=10, lineType=cv2.LINE_AA) # Add ID label

        img_id_all = cv2.add(img_id_all, img_id_col)

        d_tmp = [file_name, ID, area, perimeter, roundness, cvtmm(w), cvtmm(h)]
        d_csv.append(d_tmp)

        ### 輪郭の座標をcsvで保存する
        d_contour = list()
        #save_name_contour = './Quercus_画像/' + site + '/'+ 'contour/contour_' + file_name + '_' + str(ID) + '.csv'
        save_name_contour = './contour/' + dir_img + '/contour_' + file_name + '_' + str(ID) + '.csv'
        for j in range(len(contours[LID][:])):
            d_contour.append(contours[LID][j][0])
            with open(save_name_contour, 'w') as f:
                writer = csv.writer(f) # Save as csv file
                writer.writerows(d_contour)
        print('ID', ID, 'is saved')
        ID = ID + 1

    ### 重心やラベルを描画した画像を保存
    print('Labeled as below image')
    img_id_all_plt = img_id_all.copy()
    img_id_all_plt = cv2.cvtColor(img_id_all, cv2.COLOR_BGR2RGB)
    save_name_label_all = './label/' + dir_img + '/' + file_name +'_label.jpg'
    cv2.imwrite(save_name_label_all, img_id_all)

    ### 面積などの情報をcsvで保存
    d_csv = pd.DataFrame(d_csv)
    save_name_csv = '/leaf_info/' + dir_img + '/' + file_name + '_leaf_info.csv'
    #d_csv.to_csv(f"/Users/konrai/Dropbox/hawaii_leaf_personal/Quercus_全体共有/調査-宮城/EFD/csv/{file_name}.csv")
    d_csv.to_csv(f"/Users/konrai/Dropbox/hawaii_leaf_personal/leaf_info/{dir_img}/{file_name}_leaf_info.csv")

    print('A dataframe have been made as below;\n', d_csv)
    print(file_name, ' finished')

print('All processing finished')