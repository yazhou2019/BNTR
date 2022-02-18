% Motion.mat is the original motion file from the source
motion=load('Motion.mat');
timeindex=15001:25000;
usedindex=floor(1000*motion.MotionTime(timeindex));
