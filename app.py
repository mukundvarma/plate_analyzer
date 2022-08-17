import numpy as np
import os
import pandas as pd
import streamlit as st
from analysis import analyze, calculate_standard_curve, collate_results, VALID_WELLS
from plotting import plot_rfu_panel, plot_standard_curve, plot_platemap

st.set_page_config(layout='wide')

@st.cache
def read_input(file, index_col='Cycle'):
	filename = file.name
	ext = filename.split(".")[-1]
	assert ext in ['csv', 'xlsx', 'xls'], "Unsupported format. Please upload a csv or xlsx file"
	read_cmd = {"csv": pd.read_csv, "xlsx": pd.read_excel, "xls": pd.read_excel}[ext]
	try:
		df = read_cmd(file, index_col = index_col)
	except KeyError:
		if ext in ["xls", "xlsx"]:
			raise Exception("You might be using an export from an old and unsupported version of excel. Try uploading a csv export instead.")
	if df.shape[1] == 97:
		df = df.iloc[:, 1:]

	if df.shape[1] != 96:
		raise ValueError("Uploaded file has the wrong number of columns for a 96-well plate. Please check the input.")
	return df


@st.cache
def convert_df(df):
 # IMPORTANT: Cache the conversion to prevent computation on every rerun
 return df.to_csv().encode('utf-8')

def main():
	st.title("Plate analyzer")
	st.sidebar.header("Input file upload")

	uploaded_file = st.sidebar.file_uploader("Choose a file")
	index_column = st.sidebar.text_input("Index column name", "Cycle")

	if uploaded_file is not None:

		rfu_tab, sc_tab, res_tab, raw_tab = st.tabs(['RFU plots', 'Standard Curve', 'Results', 'Raw Data'])

		df = read_input(uploaded_file, index_column)
		plate_map = pd.DataFrame(df.columns.astype(str).values.reshape(8, 12), index=np.arange(1, 9), columns=np.arange(1, 13))
		df.columns = pd.Categorical(df.columns, categories=VALID_WELLS, ordered=True)		

		merged, slopes, intercepts = analyze(df)


		results_table = collate_results(slopes, intercepts, df)


		with rfu_tab:
			
			if sum(slopes < 0) > 0:
				st.warning("Some rates infered from the linear fit were found to be negative.")
			fig1 = plot_rfu_panel(merged, ret=True)
			st.plotly_chart(fig1)

		with sc_tab:
			fold_dilutions = st.number_input("Fold Dilutions", 2, )
			min_nz = st.number_input("Min non-zero DNA Amount", value=1.0, min_value=0.0, step=1e-6, format='%.6f')
			means = calculate_standard_curve(merged, min_nz, fold_dilutions)
			standard_slope, r, fig2 = plot_standard_curve(means)
			st.write(f"Standard Slope: {standard_slope:.3f}")
			st.write(f"\nR-squared: {r ** 2:.3f}")
			st.plotly_chart(fig2)
		with res_tab:
			col1, col2 = st.columns([1, 3])
			try:
				results_table['Rate (nmol NTP/cycle)'] = results_table['Rate (RFU/Cycle)'] / standard_slope 
				results_table.index = pd.Categorical(results_table.index, categories=VALID_WELLS, ordered=True)
			except:
				st.write("Error in calculating conversion factor.")

			with col1:
				st.write(results_table.style.format("{:.2f}"))
				csv = convert_df(results_table)

				st.download_button(
				     label="Download data as CSV",
				     data=csv,
				     file_name='results.csv',
				     mime='text/csv',
				 )
			with col2:
				for column in results_table.columns:
					with st.expander(column):
						hm = plot_platemap(results_table[column], plate_map)
						st.plotly_chart(hm)

		with raw_tab:
			st.dataframe(df.style.format("{:.2f}"))




if __name__== "__main__":
	main()

