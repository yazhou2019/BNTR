#Read_the_data
chname=c('ECoG_ch1.mat','ECoG_ch2.mat','ECoG_ch3.mat','ECoG_ch4.mat',
         'ECoG_ch5.mat','ECoG_ch6.mat','ECoG_ch7.mat','ECoG_ch8.mat',
         'ECoG_ch9.mat','ECoG_ch10.mat','ECoG_ch11.mat','ECoG_ch12.mat',
         'ECoG_ch13.mat','ECoG_ch14.mat','ECoG_ch15.mat','ECoG_ch16.mat',
         'ECoG_ch17.mat','ECoG_ch18.mat','ECoG_ch19.mat','ECoG_ch20.mat',
         'ECoG_ch21.mat','ECoG_ch22.mat','ECoG_ch23.mat','ECoG_ch24.mat',
         'ECoG_ch25.mat','ECoG_ch26.mat','ECoG_ch27.mat','ECoG_ch28.mat',
         'ECoG_ch29.mat','ECoG_ch30.mat','ECoG_ch31.mat','ECoG_ch32.mat',
         'ECoG_ch33.mat','ECoG_ch34.mat','ECoG_ch35.mat','ECoG_ch36.mat',
         'ECoG_ch37.mat','ECoG_ch38.mat','ECoG_ch39.mat','ECoG_ch40.mat',
         'ECoG_ch41.mat','ECoG_ch42.mat','ECoG_ch43.mat','ECoG_ch44.mat',
         'ECoG_ch45.mat','ECoG_ch46.mat','ECoG_ch47.mat','ECoG_ch48.mat',
         'ECoG_ch49.mat','ECoG_ch50.mat','ECoG_ch51.mat','ECoG_ch52.mat',
         'ECoG_ch53.mat','ECoG_ch54.mat','ECoG_ch55.mat','ECoG_ch56.mat',
         'ECoG_ch57.mat','ECoG_ch58.mat','ECoG_ch59.mat','ECoG_ch60.mat',
         'ECoG_ch61.mat','ECoG_ch62.mat','ECoG_ch63.mat','ECoG_ch64.mat'
         
)

chlist=list()
library("R.matlab")


path<-('20100802S1_Epidural-ECoG+Food-Tracking_B_Kentaro+Shimoda_mat_ECoG64-Motion6')



for(i in 1:64){
  pathname<-file.path(path,chname[i])
  chlist[[i]]<-readMat(pathname)
}


n=length(chlist[[1]][[1]])

chall=matrix(0,64,n)
for(i in 1:64){
  chall[i,]=chlist[[i]][[1]]
}


write.table(chall,'chall_big_monkey1.csv',row.names = F,col.names = F)
