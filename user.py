from pathlib import Path


def convert_df(df, name):
	filepath = Path(name)
	return df.to_csv(filepath)
