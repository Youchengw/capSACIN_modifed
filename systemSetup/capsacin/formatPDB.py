import pandas as pd

def formatPDB(df):
    df['x'] = df['x'].astype('float').apply('{0:.3f}'.format)
    df['y'] = df['y'].astype('float').apply('{0:.3f}'.format)
    df['z'] = df['z'].astype('float').apply('{0:.3f}'.format)
    df['b'] = df['b'].astype('float').apply('{0:.2f}'.format)
    df['occ'] = df['occ'].astype('float').apply('{0:.2f}'.format)
    # Keep the record and coordinate columns fixed beyond 99,999 atoms.
    # This matches MDAnalysis PDBWriter's five-digit serial convention.
    df['idx'] = df['idx'].astype(int) % 100000
    df = df.astype('string')
    df.atom = df.atom.str.strip().str.ljust(6, ' ')
    df.idx = df.idx.str.rjust(5, ' ').str.ljust(7, ' ')
    df.name = df.name.str.ljust(4, ' ') # add 1-3 spaces to right
    df.resname = df.resname.str.ljust(4, ' ') # add 1-3 spaces to right
    df.chain = df.chain.str.ljust(2, ' ') # add 1-2 spaces to right
    df.resids = df.resids.str.ljust(6, ' ')
    df.x = df.x.str.rjust(8, ' ')
    df.y = df.y.str.rjust(8, ' ')
    df.z = df.z.str.rjust(8, ' ')
    df.occ = df.occ.str.rjust(6, ' ')
    df.b = df.b.str.rjust(6, ' ')
    df.type = df.type.str.rjust(12, ' ')

    return df
