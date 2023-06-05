import streamlit as st


def _convert_df(df, name):
	print(name)
	return df.to_csv(name)


def download_df_as_csv(df, name):
	_convert_df(df, name)
	label = '{0}{1}{2}'.format("Download ", name, " as CSV")
	st.download_button(
		label=label,
		data=name,
		file_name=name)
