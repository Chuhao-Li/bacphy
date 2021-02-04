rp_list = ["30S ribosomal protein S" + str(i) for i in range(1, 22)] + \
          ["50S ribosomal protein L" + str(i) for i in range(1, 37) if i not in [7, 8, 12, 26]] + \
          ['50S ribosomal protein L7/L12']

for i in rp_list:
    print(i)
