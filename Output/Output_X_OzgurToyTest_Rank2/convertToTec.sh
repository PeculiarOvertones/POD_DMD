args=("$@")
i=1
for i in `seq 0 1 ${args[0]}`;
do
	python post.py keep POD $i
done

for i in `seq 0 1 ${args[0]}`;
do
	python post.py keep DirSVD $i
done

for i in `seq 0 1 ${args[0]}`;
do
	python post.py keep DMD_Real $i
done

for i in `seq 0 1 ${args[0]}`;
do
	python post.py keep DMD_Imag $i
done
