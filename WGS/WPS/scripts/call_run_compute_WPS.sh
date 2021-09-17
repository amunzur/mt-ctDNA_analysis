scriptpath="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/scripts/run_compute_WPS.sh";

window_size="15";

for bed_file in /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/data/beds/*.bed; do
sbatch ${scriptpath} ${bed_file} ${window_size};
done
