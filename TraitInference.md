<a name="Machine Learning"/>

# Machine Learning for MAG Trait Inference

Need to install genephene repo:
```
cd ~/repos/
git clone https://github.com/chrisquince/genephene.git
cd genephene
cd datafiles
gunzip taxid_ko_matrix_all_full.csv.gz 
cd ..
```

```
pip3 install jupyter --user
sudo pip3 install pandas
sudo pip3 install scikit-learn
sudo pip3 install matplotlib

export PATH=/home/ubuntu/.local/bin/:$PATH
```

Log on to the virtual machine. Launch a Jupyter Notebook instance by typing:

```
jupyter notebook --no-browser --port=8889 --ip=127.0.0.1
```

Then go to a local terminal and type:

```
 ssh -N -f -L localhost:8887:localhost:8889 ubuntu@xxx.yyy.zzz.kkk
 ```
 
 This allows you to access the notebook server from your local machine

 Finally, open a web browser window and type ‘localhost:8887’ into the URL bar

 Open the 'train_classifier.ipynb' file and then clear output and restart kernel.
 
  