import os
output = open("data/stats/summary.txt", "w") 

for file_name in os.listdir("data/stats"): 
    if file_name.endswith(".stats.txt"): 
        stats_file = open("data/stats/" + file_name)
        for line in stats_file:
            if "+ 0 mapped" in line: 
                output.write (file_name.replace(".stats.txt", "") + 
                       "," + 
                       line.rsplit(" ")[4].replace("(","").replace("%","") +
                             "\n")
                
output.close()
