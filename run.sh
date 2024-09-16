cd data
python ../prep_ims_2D.py 0-5-6-7-12.xyz
python ../get_ims_data.py 0-5-6-7-12.xyz
cd ../
python ./main_plot.py data 0-5-6-7-12

