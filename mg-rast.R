# GLUSEEN MG-RAST download ####

library(devtools)
install_github("djeppschmidt/mgrastr")
library(mgrastr)
# script executed in bash/python ####
# found here: https://github.com/MG-RAST/MG-RAST-Tools/blob/master/scripts/mg-compare-functions.py

# download instructions: 
#git clone http://github.com/MG-RAST/MG-RAST-Tools
#cd MG-RAST-Tools
#python setup.py build
#sudo python setup.py install

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812553.3,mgm4812534.3,mgm4812503.3,mgm4812482.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Function/mgRast_ontologyLVL3_list66.txt

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812451.3,mgm4812551.3,mgm4812526.3,mgm4812531.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Function/mgRast_ontologyLVL3_list67.txt

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812503.3,mgm4812531.3,mgm4812526.3,mgm4812553.3,mgm4812451.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list01.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812482.3,mgm4812534.3,mgm4812551.3,mgm4681020.3,mgm4681131.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list02.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681059.3,mgm4679356.3,mgm4681106.3,mgm4679401.3,mgm4679358.3,mgm4681135.3,mgm4679309.3,mgm4735876.3,mgm4815619.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Function/mgRast_ontologyLVL3_list03.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679349.3,mgm4681079.3,mgm4681103.3,mgm4679341.3,mgm4681097.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list04.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679370.3,mgm4679316.3,mgm4679300.3,mgm4681090.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list05.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679304.3,mgm4681101.3,mgm4679416.3,mgm4681136.3,mgm4679328.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list06.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681128.3,mgm4679312.3,mgm4679383.3,mgm4681054.3,mgm4679371.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list07.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681052.3,mgm4681089.3,mgm4681034.3,mgm4679343.3,mgm4679376.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list08.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681104.3,mgm4681025.3,mgm4679338.3,mgm4679308.3,mgm4679394.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list09.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681039.3,mgm4679314.3,mgm4679325.3,mgm4679296.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list10.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679298.3,mgm4681124.3,mgm4681061.3,mgm4681069.3,mgm4679355.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list11.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679388.3,mgm4681030.3,mgm4679297.3,mgm4679320.3,mgm4679417.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list12.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679346.3,mgm4679403.3,mgm4679378.3,mgm4681088.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list13.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679380.3,mgm4679339.3,mgm4679390.3,mgm4679301.3,mgm4679382.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list14.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679345.3,mgm4679311.3,mgm4735862.3,mgm4735939.3,mgm4679367.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list15.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679336.3,mgm4681023.3,mgm4679318.3,mgm4679387.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list16.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679414.3,mgm4679332.3,mgm4679360.3,mgm4679374.3,mgm4681064.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list17.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679377.3,mgm4679303.3,mgm4679381.3,mgm4679407.3,mgm4679368.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list18.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679413.3,mgm4684216.3,mgm4679397.3,mgm4681112.3,mgm4681041.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list19.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679323.3,mgm4679395.3,mgm4679342.3,mgm4679361.3,mgm4679334.3,mgm4679384.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list20.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679363.3,mgm4681035.3,mgm4679348.3,mgm4679335.3,mgm4679359.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list21.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681066.3,mgm4679295.3,mgm4679354.3,mgm4679357.3,mgm4679404.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list22.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681140.3,mgm4679409.3,mgm4681126.3,mgm4681019.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list23.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679419.3,mgm4679362.3,mgm4681095.3,mgm4679366.3,mgm4679307.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list24.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679347.3,mgm4681127.3,mgm4679322.3,mgm4679299.3,mgm4681077.3,mgm4681110.3,mgm4679294.3,mgm4681130.3,mgm4681040.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list25.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679305.3,mgm4681080.3,mgm4679411.3,mgm4679351.3,mgm4681072.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list26.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679329.3,mgm4679364.3,mgm4681037.3,mgm4681121.3,mgm4679330.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list27.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681129.3,mgm4679340.3,mgm4679344.3,mgm4679302.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list28.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679352.3,mgm4679324.3,mgm4684214.3,mgm4679333.3,mgm4681092.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list29.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681082.3,mgm4679319.3,mgm4679398.3,mgm4681102.3,mgm4681043.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list30.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679406.3,mgm4681022.3,mgm4681087.3,mgm4679420.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list31.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679321.3,mgm4681081.3,mgm4679393.3,mgm4681060.3,mgm4679353.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list32.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681018.3,mgm4679410.3,mgm4684215.3,mgm4681119.3,mgm4679402.3,mgm4679415.3,mgm4679396.3,mgm4679385.3,mgm4681042.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list33.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681139.3,mgm4679379.3,mgm4679331.3,mgm4681027.3,mgm4679372.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list34.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679412.3,mgm4681100.3,mgm4679317.3,mgm4679327.3,mgm4679391.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list35.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679408.3,mgm4679375.3,mgm4679313.3,mgm4679389.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list36.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679392.3,mgm4679400.3,mgm4681141.3,mgm4679421.3,mgm4681114.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list37.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679306.3,mgm4681138.3,mgm4679315.3,mgm4679369.3,mgm4679399.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list38.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4681137.3,mgm4679405.3,mgm4679386.3,mgm4679310.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list39.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4679337.3,mgm4679350.3,mgm4812491.3,mgm4812559.3,mgm4812480.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list40.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812570.3,mgm4812477.3,mgm4812550.3,mgm4812486.3,mgm4812487.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list41.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812452.3,mgm4812528.3,mgm4812498.3,mgm4812515.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list42.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812461.3,mgm4812524.3,mgm4812508.3,mgm4812518.3,mgm4812517.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list43.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812459.3,mgm4812469.3,mgm4812563.3,mgm4812495.3,mgm4812573.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list44.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812544.3,mgm4812554.3,mgm4812471.3,mgm4812571.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list45.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812542.3,mgm4812535.3,mgm4812557.3,mgm4812460.3,mgm4812473.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list46.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812455.3,mgm4812541.3,mgm4812496.3,mgm4812453.3,mgm4812575.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list47.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812506.3,mgm4812525.3,mgm4812505.3,mgm4812468.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list48.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812465.3,mgm4812476.3,mgm4812462.3,mgm4812529.3,mgm4812485.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list49.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812565.3,mgm4812537.3,mgm4812467.3,mgm4812456.3,mgm4812549.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list50.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812478.3,mgm4812543.3,mgm4812561.3,mgm4812479.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list51.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812545.3,mgm4812556.3,mgm4812538.3,mgm4812502.3,mgm4812509.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list52.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812450.3,mgm4812555.3,mgm4812499.3,mgm4812527.3,mgm4812500.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list53.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812512.3,mgm4812558.3,mgm4812519.3,mgm4812457.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list54.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812472.3,mgm4812466.3,mgm4812454.3,mgm4812489.3,mgm4812501.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list55.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812514.3,mgm4812504.3,mgm4812516.3,mgm4812458.3,mgm4812484.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list56.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812530.3,mgm4812490.3,mgm4812494.3,mgm4812560.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list57.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812536.3,mgm4812562.3,mgm4812572.3,mgm4812546.3,mgm4812567.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list58.txt #
 
#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812539.3,mgm4812493.3,mgm4812522.3,mgm4812475.3,mgm4812521.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list59.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812523.3,mgm4812483.3,mgm4812464.3,mgm4812492.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list60.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812513.3,mgm4812507.3,mgm4812449.3,mgm4812463.3,mgm4812548.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list61.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812540.3,mgm4812564.3,mgm4812569.3,mgm4812488.3,mgm4812568.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list62.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812510.3,mgm4812566.3,mgm4812511.3,mgm4812520.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list63.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812470.3,mgm4812532.3,mgm4812574.3,mgm4812481.3,mgm4812552.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list64.txt #

#mg-compare-functions.py --token "gbdwiyCZqACamvn8fG59aTs3Z" --ids "mgm4812547.3,mgm4812497.3,mgm4812533.3,mgm4812474.3" --level level3 --source Subsystems --format text --evalue 10 > ~/Desktop/PhD/Metagenome/Download/mgRast_ontologyLVL3_list65.txt #

# get organisms ####
library(httr)
library(jsonlite)
base<-"https://api-ui.mg-rast.org/metagenome/"
sample<-"703dc7e20f6d676d343733353834342e33"
end<-"?verbosity=stats&detail=taxonomy&auth=gbdwiyCZqACamvn8fG59aTs3Z"

call1<-paste(base,sample,end, sep="")
sample_003tax<-GET(call1)
sample_003tax_text<-content(sample_003tax, "text")
sample_003tax_json<-fromJSON(sample_003tax_text, flatten = TRUE)

downloadTaxa<-function(sample, auth){
  require(httr)
  require(jsonlite)
  base<-"https://api-ui.mg-rast.org/metagenome/"
  end<-"?verbosity=stats&detail=taxonomy&auth="
  auth="gbdwiyCZqACamvn8fG59aTs3Z"
  s.tax<-fromJSON(content(GET(paste(base,sample,end,auth,sep="")), "text"), flatten=T)
  s.tax
}

# test taxa download ####

s003Tax<-downloadTaxa(sample) #it works!!!

ag.Phylum(l.taxa$`003_R1_022417`)

# make list of samples ####
O.key<-read.csv("~/Desktop/PhD/Metagenome/SampleKey.csv")
O.key<-data.frame(O.key)
O.key$Seq_ID<-as.character(O.key$Seq_ID)
O.key$Private_ID<-as.character(O.key$MG_ID)

l.sID<-as.list(O.key$MG_ID)
names(l.sID)<-O.key$Seq_ID

l.taxa<-lapply(l.sID, downloadTaxa, auth="gbdwiyCZqACamvn8fG59aTs3Z")
#l.taxa2<-download.O(l.sID[1:10], auth="gbdwiyCZqACamvn8fG59aTs3Z")
#saveRDS(l.taxa, "~/Desktop/PhD/Metagenome/Organism/Organism_Download.RDS")
#l.taxa<-readRDS("~/Desktop/PhD/Metagenome/Organism/Organism_Download.RDS")


# proof of concept table parsing####
library(reshape2)
phylum<-ldply(l.taxa, ag.Phylum)
names(phylum)<-c("ID", "taxa", "values")
t.phylum<-dcast(phylum, ID~taxa)
t.phylum[is.na(t.phylum)]<-0


a<-c(1,2,3,4,5)
b<-c(1,1,1,1,1)
c<-c(3,3,3,3,3)
df<-data.frame(a, "c"=b+c)
df

# table parsing ####
t.domain<-ag.tax(l.taxa, ag.Domain)
t.phylum<-ag.tax(l.taxa, ag.Phylum)
t.class<-ag.tax(l.taxa, ag.Class)
t.order<-ag.tax(l.taxa, ag.Order)
t.family<-ag.tax(l.taxa, ag.family)
t.genus<-ag.tax(l.taxa, ag.Genus)
t.species<-ag.tax(l.taxa, ag.Species)

sample_key<-read.csv("~/Desktop/PhD/Metagenome/SampleKey.csv")
sample_key<-data.frame(sample_key)
sample_key$Seq_ID<-as.character(sample_key$Seq_ID)
rownames(sample_key)<-sample_key$Seq_ID

sam<-as.data.frame(read.csv("~/Desktop/PhD/Metagenome/GLUSEEN_HM2_metadata.csv"))
sam$Sample_ID2<-as.character(sam$Sample_ID2)
rownames(sam)<-sam$Sample_ID2


ct.domain<-condense(t.domain, samkey=sample_key, sam=sam)
ct.phylum<-condense(t.phylum, samkey=sample_key, sam=sam)
ct.class<-condense(t.class, samkey=sample_key, sam=sam)
ct.order<-condense(t.order, samkey=sample_key, sam=sam)
ct.family<-condense(t.family, samkey=sample_key, sam=sam)
ct.genus<-condense(t.genus, samkey=sample_key, sam=sam)
ct.species<-condense(t.species, samkey=sample_key, sam=sam)

saveRDS(ct.domain, "~/Desktop/PhD/Metagenome/Organism/ct_domain.RDS")
saveRDS(ct.phylum, "~/Desktop/PhD/Metagenome/Organism/ct_phylum.RDS")
saveRDS(ct.class, "~/Desktop/PhD/Metagenome/Organism/ct_class.RDS")
saveRDS(ct.order, "~/Desktop/PhD/Metagenome/Organism/ct_order.RDS")
saveRDS(ct.family, "~/Desktop/PhD/Metagenome/Organism/ct_family.RDS")
saveRDS(ct.genus, "~/Desktop/PhD/Metagenome/Organism/ct_genus.RDS")
saveRDS(ct.species, "~/Desktop/PhD/Metagenome/Organism/ct_species.RDS")

# I don't actually need all these for the final analysis, just the sample sums for the bacterial domain; but I'll need these later ;)

#<-sample_sums(subset_taxa(ct.domain, Domain=="Bacteria"))



# deprecated now that I have an R version working ####
path1<-"~/Desktop/PhD/Metagenome/Function/"

library(plyr)
library(phyloseq)

F.files<-list.files(path1, pattern=".txt")
F.filesLong<-paste(path1, F.files, sep="")
F.list.df<-lapply(F.filesLong, read_files)


ps.F<-merge_phyloseq(F.list.df[[1]], F.list.df[[2]],F.list.df[[3]],F.list.df[[4]],F.list.df[[5]],F.list.df[[6]],F.list.df[[7]],F.list.df[[8]],F.list.df[[9]],F.list.df[[10]],F.list.df[[11]], F.list.df[[12]],F.list.df[[13]],F.list.df[[14]],F.list.df[[15]],F.list.df[[16]],F.list.df[[17]],F.list.df[[18]],F.list.df[[19]],F.list.df[[20]],F.list.df[[21]], F.list.df[[22]],F.list.df[[23]],F.list.df[[24]],F.list.df[[25]],F.list.df[[26]],F.list.df[[27]],F.list.df[[28]],F.list.df[[29]],F.list.df[[30]],F.list.df[[31]], F.list.df[[32]],F.list.df[[33]],F.list.df[[34]],F.list.df[[35]],F.list.df[[36]],F.list.df[[37]],F.list.df[[38]],F.list.df[[39]],F.list.df[[40]],F.list.df[[41]], F.list.df[[42]],F.list.df[[43]],F.list.df[[44]],F.list.df[[45]],F.list.df[[46]],F.list.df[[47]],F.list.df[[48]],F.list.df[[49]],F.list.df[[50]],F.list.df[[51]], F.list.df[[52]],F.list.df[[53]],F.list.df[[54]],F.list.df[[55]],F.list.df[[56]],F.list.df[[57]],F.list.df[[58]],F.list.df[[59]],F.list.df[[60]],F.list.df[[61]], F.list.df[[62]],F.list.df[[63]],F.list.df[[64]],F.list.df[[65]],F.list.df[[66]],F.list.df[[67]])

# Find missing accessions

key<-read.csv("~/Desktop/PhD/Metagenome/SampleKey.csv")
key<-data.frame(key)
setdiff(key$MG_ID,colnames(otu_table(ps.F)))
rownames(key)<-key$MG_ID
# output: "mgm4681135.3" "mgm4679309.3" "mgm4735876.3" "mgm4815619.3"

condense.F<-function(df, samkey, sam){
  require(phyloseq)
  ps<-phyloseq(df, sample_data(samkey))
  tax<-matrix(data=NA, nrow=length(taxa_names(ps)), ncol=1)
  rownames(tax)<-taxa_names(ps)
  tax[,1]<-taxa_names(ps)
  colnames(tax)<-"Function"
  tax_table(ps)<-tax
  ps2<-merge_samples(ps, "Sample_ID", fun = sum)
  sample_data(ps2)<-sample_data(sam)
  ps2
}

ps.F<-condense.F(ps.F, key, sam)

#test.file1<-read_files("~/Desktop/PhD/Metagenome/Function/mgRast_ontologyLVL3_list01.txt")
#test.file2<-read_files("~/Desktop/PhD/Metagenome/Function/mgRast_ontologyLVL3_list07.txt")
#test.ps<-merge_phyloseq(test.file1, test.file2)
#otu_table(test.ps)[1,]
#taxa_names(test.file1)

ncol(otu_table(ps.F))

# function download with R ####
#https://api-ui.mg-rast.org/metagenome/f6d2be46306d676d343831323437372e33?verbosity=stats&detail=ontology

#https://api.mg-rast.org//matrix/function?group_level=level3&source=Subsystems&evalue=10&identity=60&length=15&version=1&result_type=abundance&asynchronous=1&id=mgm4812540.3&id=mgm4812564.3&id=mgm4812569.3&id=mgm4812488.3&id=mgm4812568.3
O.key<-read.csv("~/Desktop/PhD/Metagenome/SampleKey.csv")
O.key<-data.frame(O.key)
O.key$Seq_ID<-as.character(O.key$Seq_ID) # unique sample IDs
O.key$Private_ID<-as.character(O.key$MG_ID) # unique accession codes in mgm... .3 format

l.sID<-as.character(O.key$MG_ID) # make list of accession numbers
names(l.sID)<-O.key$Seq_ID # name list according to unique sample IDs


url<-loadFunc(l.sID,  auth="gbdwiyCZqACamvn8fG59aTs3Z", ont="Subsystems", level="3", E="20", length="15", id="60", parallel=FALSE)

dt<-downloadFunc(url, ont="Subsystems", level=3)
dt<-condense(dt, sample_key, sam)

url2<-loadFunc(l.sID,  auth="gbdwiyCZqACamvn8fG59aTs3Z", ont="Subsystems", level="1", E="20", length="15", id="60", parallel=FALSE)
dt.l1<-downloadFunc(url2, ont="Subsystems", level=1)
dt.l1<-condense(dt.l1, sample_key, sam)
saveRDS(dt.l1, "~/Desktop/PhD/Metagenome/Function/dt20L1.RDS")


url3<-loadFunc(l.sID,  auth="gbdwiyCZqACamvn8fG59aTs3Z", ont="Subsystems", level="4", E="20", length="15", id="60", parallel=FALSE)
dt.l4<-downloadFunc(url3, ont="Subsystems", level=4)
dt.l4<-condense(dt.l4, sample_key, sam)
saveRDS(dt.l1, "~/Desktop/PhD/Metagenome/Function/dt20L1.RDS")


url.15.4<-loadFunc(l.sID, auth="gbdwiyCZqACamvn8fG59aTs3Z", ont="Subsystems", level="4", E="15", length="15", id="60", parallel=FALSE)#why does this keep screwing up?

url.15.4<-lapply(l.sID, mg.load, auth="gbdwiyCZqACamvn8fG59aTs3Z", ont="Subsystems", level="4", E="15", length="15", id="60")

dt15.l1<-downloadFunc(url.15.4, ont="Subsystems", level=1)
dt15.l1<-condense(dt15.l1, sample_key, sam)
dt15.l2<-downloadFunc(url.15.4, ont="Subsystems", level=2)
dt15.l2<-condense(dt15.l2, sample_key, sam)
dt15.l3<-downloadFunc(url.15.4, ont="Subsystems", level=3)
dt15.l3<-condense(dt15.l3, sample_key, sam)
dt15.l4<-downloadFunc(url.15.4, ont="Subsystems", level=4)
dt15.l4<-condense(dt15.l4, sample_key, sam)
saveRDS(dt15.l1, "~/Desktop/PhD/Metagenome/Function/dt15l1.RDS")
saveRDS(dt15.l2, "~/Desktop/PhD/Metagenome/Function/dt15l2.RDS")
saveRDS(dt15.l3, "~/Desktop/PhD/Metagenome/Function/dt15l3.RDS")
saveRDS(dt15.l4, "~/Desktop/PhD/Metagenome/Function/dt15l4.RDS")
#explore

dt20<-downloadFunc(url, ont="Subsystems", level=3)
dt20<-condense(dt20, sample_key, sam)
hist(sample_sums(dt20), main="dt20")


dt10<-downloadFunc(url.10, ont="Subsystems", level=3)
dt10<-condense(dt10, sample_key, sam)
hist(sample_sums(dt10), main="dt10")

dt5<-downloadFunc(url.5, ont="Subsystems", level=3)
dt5<-condense(dt5, sample_key, sam)
hist(sample_sums(dt5), main="dt5")

# going to keep dt20

saveRDS(dt20, "~/Desktop/PhD/Metagenome/Function/dt20.RDS")
     