
import pandas as pd
import os
import sys

def screen_finder_func(screen_bed_directory, o_bed_path):
    # Create an empty dataframe to store the concatenated results
    all_screen_data = pd.DataFrame()

    # Iterate through each .bed file in the directory
    for screen_bed_file in os.listdir(screen_bed_directory):
        if screen_bed_file.endswith('.bed'):
            # Read the current .bed file
            screen_bed_path = os.path.join(screen_bed_directory, screen_bed_file)
            screen = pd.read_csv(screen_bed_path, sep='\t', header=None, index_col=False,
                                 names=['chr', 'start', 'end', 'wv1', 'wv2', 'wv3'])
            screen.sort_values(by=['chr', 'start'], inplace=True)

            # Save the sorted screen data to a CSV file
            screen_csv_name = f"{os.path.splitext(screen_bed_file)[0]}_screen_sorted.csv"
            screen.to_csv(screen_csv_name, sep='\t', index=False)

            # Append the current screen data to the concatenated dataframe
            all_screen_data = pd.concat([all_screen_data, screen])

    chr_array = all_screen_data['chr'].to_numpy()
    chrs = {}
    inter = [0]

    o_bed_files = [f for f in os.listdir(o_bed_path) if f.endswith('.bed')]
    l_strip = 'lifted_'
    bed_strip = '.bed'

    for i in range(len(chr_array) - 1):
        if chr_array[i] != chr_array[i + 1]:
            inter.append(i + 1)  # beginning of the next chromosome

    inter.append(len(chr_array))

    for i in range(len(inter) - 1):
        chrs[chr_array[inter[i]]] = (inter[i], inter[i + 1] - 1)  # the interval is left-closed right-open [ )

    print(chrs)

    counter = 0
    overlap_num = []
    threshold = 1

    for file in o_bed_files:
        counter += 1
        res = []
        line_num = 0
        with open(os.path.join(o_bed_path, file)) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                line_num += 1
                line = line.rstrip('\n').split('\t')
                chr_name = line[0]
                chr_start = chrs[chr_name][0]
                chr_end = chrs[chr_name][1]
                start = int(line[1])  # start of the target
                end = int(line[2])  # end of the target

                idx_start = all_screen_data.iloc[chr_start:chr_end + 1]['start'].searchsorted(start, side='left') + chr_start
                idx_end = all_screen_data.iloc[chr_start:chr_end + 1]['end'].searchsorted(end, side='right') + chr_start

                if idx_start > chr_start and start < all_screen_data.iloc[idx_start - 1]['end']:
                    idx_start -= 1
                if idx_end <= chr_end and end > all_screen_data.iloc[idx_end]['start']:
                    idx_end += 1

                for i in range(idx_start, idx_end):
                    template_start = all_screen_data.iloc[i]['start']
                    template_end = all_screen_data.iloc[i]['end']
                    collect = True
                    percent = 1

                    if template_start < start and end < template_end:
                        collect = True
                        percent = (end - start + 1) / (template_end - template_start)
                    elif template_start < start:
                        percent = (template_end - start + 1) / (template_end - template_start + 1)
                        if percent < threshold:
                            collect = False
                    elif end < template_end:
                        percent = (end - template_start + 1) / (template_end - template_start + 1)
                        if percent < threshold:
                            collect = False

                    if collect:
                        res.append([chr_name, start, end, template_start, template_end, percent])

        output = pd.DataFrame(res, columns=['chrNO', 'my_start', 'my_end', 'screen_start', 'screen_end', 'overlap %'])
        output.to_csv(file.replace(l_strip, '').replace(bed_strip, '') + '.csv', header=True, index=False)
        overlap_num.append([file.replace(l_strip, '').replace(bed_strip, ''), len(output)])
        print('file {0} has {1} overlapping gene {2}/{3}...'.format(file, len(output), counter, len(o_bed_files)))

    lxf_df = pd.DataFrame(overlap_num, columns=['Cell_type', 'overlap_num'])
    lxf_df.to_csv('overlapping.csv', header=True, index=False)

    print('Processed in total {0} files'.format(counter))
    print('Programme reached the end_hello')
    print


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <screen_bed_directory> <o_bed_path>")
    else:
        screen_bed_directory = sys.argv[1]
        o_bed_path = sys.argv[2]
        screen_finder_func(screen_bed_directory, o_bed_path)
