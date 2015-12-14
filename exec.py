__author__ = 'supun'

import subprocess

def exec_main():
    iter = [100]
    rows = [87, 174, 348, 696, 1392, 2088]
    cols = [200000, 100000, 400000]
    threads = [1, 2, 4, 8, 16, 24]

    for i in iter:
        for c in cols:
            r_in = 87 * c / 100000
            for t in threads:
                r = r_in * t
                p = subprocess.Popen(['./a.out', str(t), str(i), str(r), str(c)], stdout=subprocess.PIPE);
                out, err = p.communicate()
                print out

if __name__ == "__main__":
    exec_main()
