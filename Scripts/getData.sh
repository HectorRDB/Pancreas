loc="./"
curl https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/baron-human.rds > baron.rds
curl https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/muraro.rds > muraro.rds
curl https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/xin.rds > xin.rds
curl https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/segerstolpe.rds > segerstolpe.rds
mv baron-human.rds $loc/
mv muraro.rds $loc/
mv xin.rds $loc/
mv segerstolpe.rds $loc/
