# BSD License 2.0

# Copyright (c) 2020, Janka Uryga <lolzatu2@gmail.com>
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the organization nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# example row:

# {
#     '#pattern name': 'MapolyY_B0017',
#     'Family': 'BBR-BPC',
#     'sequence name': 'cle1',
#     'strand': '-',
#     'start': 31,
#     'stop': 67,
#     'score': float,
#     'p-value': float,
#     'q-value': float,
#     'matched sequence': 'CTCTCTCTCTCTCTCTCTCTC',
# }


COLNAMES = [
    '#pattern name',
    'Family',
    'sequence name',
    'start',
    'stop',
    'strand',
    'score',
    'p-value',
    'q-value',
    'matched sequence',
]


def bounds(r):
    return (r['start'], r['stop'])



def range_union(a, b):
    (a_start, a_stop) = a
    (b_start, b_stop) = b
    return (min(a_start, b_start), max(a_stop, b_stop))



def group_by(xs, f):
    grouped = {}
    for x in xs:
        fx = f(x)
        if fx in grouped:
            grouped[fx].append(x)
        else:
            grouped[fx] = [x]
    return grouped



def group_overlapping_ranges(ranges: 'Sequence[dict]'):
    if not ranges:
        return []
    ranges = sorted(ranges, key=lambda r: r['start'])
    first = ranges[0]
    start, stop = bounds(first)
    grouped = [[first]]
    for r in ranges[1:]:
        (r_start, r_stop) = bounds(r)
        if (start <= r_start <= stop):
            grouped[-1].append(r)
            start, stop = range_union((start, stop), (r_start, r_stop))
        else:
            grouped.append([r])
            start, stop = (r_start, r_stop)
    return grouped


import re
is_repeated_2mer = lambda seq: re.match(r'^[A-Z]?([A-Z]{2})\1+[A-Z]?$', seq)



import csv

def csv_write_dicts(rows: 'Sequence[dict]', file):
    w = csv.DictWriter(file, COLNAMES)
    w.writerow({name: name for name in COLNAMES})
    w.writerows(rows)



def main():
    import sys

    IN_NAME, OUT_NAME = sys.argv[1:]
    assert OUT_NAME.endswith('.csv')
    REMOVED_NAME = OUT_NAME.replace('.csv', '-removed.csv')

    f = open(IN_NAME)
    # f = open('fimo-1.tsv')
    rows = list(csv.DictReader(f, delimiter='\t'))

    # parse
    orig_rows = rows = [
        {
            **row,
            'start': int(row['start']),
            'stop':  int(row['stop']),
            'score':   float(row['score']),
            'p-value': float(row['p-value']),
            'q-value': float(row['q-value']),
            'matched sequence': row['matched sequence'].upper(), # normalize to uppercase
        }
        for row in rows
    ]

    # remove 1-letter sequences and repeated 2-meres
    rows = [
        row for row in rows
        if not(
            len(set(row['matched sequence'])) == 1 or
            is_repeated_2mer(row['matched sequence'])
        )
    ]

    removed = [
        row for row in orig_rows
        if (
            len(set(row['matched sequence'])) == 1 or
            is_repeated_2mer(row['matched sequence'])
        )
    ]

    with open(REMOVED_NAME, 'w', newline='') as out:
        csv_write_dicts(removed, file=out)


    no_overlaps = [
        # max(name_row_group, key=lambda row: row['stop'] - row['start'])
        min(name_row_group, key=lambda row: row['p-value'])
        for name, name_rows in group_by(rows, lambda row: (row['#pattern name'], row['sequence name'])).items()
        for name_row_group in group_overlapping_ranges(name_rows)
    ]

    # write output
    with open(OUT_NAME, 'w', newline='') as out:
        csv_write_dicts(no_overlaps, file=out)


if __name__ == '__main__':
    main()
