echo 
echo "RAND DATASET"
for i in `seq 10 10 100`;
do
	./s_m istanze/${i}_rand_dataset
done

echo "CLUSTER DATASET"
for i in `seq 10 10 100`;

do ./s_m istanze_cl/${i}_cl_dataset
done

echo "LINE DATASET"
for i in `seq 10 10 100`;

do ./s_m istanze_line/${i}_line_dataset
done

echo "CIRCLE DATASET"
for i in `seq 10 10 100`;

do ./s_m istanze_circle/${i}_circle_dataset
done
