import pacmap
import pandas as pd

def run_pacmap(in_csv, out_csv, metric):
    df = pd.read_csv(in_csv)
    mat = df.iloc[:,1:].to_numpy()
    embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0, distance=metric)
    transformed = embedding.fit_transform(mat)
    pd.DataFrame(transformed, index=df.ID, columns=['Component 1', 'Component 2']).to_csv(out_csv)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('in_csv')
    parser.add_argument('out_csv')
    parser.add_argument('metric')

    args = parser.parse_args()

    run_pacmap(args.in_csv, args.out_csv, args.metric)