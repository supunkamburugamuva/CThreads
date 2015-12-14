__author__ = 'supun'

import subprocess

def exec_main():
    iter = [100]
    rows = [87, 174, 348, 696, 1392, 2088]
    cols = [2000, 1000, 4000]
    threads = [1, 2, 4, 8, 16, 24]

    for i in iter:
        for r in rows:
            for c in cols:
                for t in threads:
                    p = subprocess.Popen(['./a.out', str(t), str(i), str(r * c / 100000), str(c)], stdout=subprocess.PIPE);
                    out, err = p.communicate()
                    print out + "\n"

if __name__ == "__main__":
    exec_main()
