#!/usr/bin/env python
import io
import sys
import polars as pl
import streamlit as st
import plotly.express as px
# import pandas as pd
# import numpy as np

TABLE_SCHEMA = {
    "Contig": pl.String,
    "Pos": pl.UInt32,
    "End": pl.UInt32,
    "GeneID": pl.String,
    "GeneName": pl.String,
    "TPM": pl.Float16,
    "Donor": pl.String,
    "TissueID": pl.String,
    "Tissue": pl.String,
    "TR_RNA_overlap": pl.UInt16,
    "TRID": pl.UInt32,
    "Motif": pl.String,
    "STRID": pl.String
}


@st.cache_data(show_spinner="Loading parquet file...")
def read_parquet(filename: str) -> pl.DataFrame:
    df = pl.read_parquet(filename, schema=TABLE_SCHEMA)
    return df


def apply_filters(
        df: pl.DataFrame,
        donors: list[str],
        tissues: list[str],
        tissues_id: list[str],
        contigs: list[str],
        pos_str: str,
        end_str: str,
        trid_str: str,
        gene_id: str,
        gene_name: str,
    ) -> tuple[pl.DataFrame, bool]:
    exprs = []
    pos: int
    end: int
    trid: int

    if donors:
        exprs.append(pl.col("Donor").is_in(donors))
    if tissues:
        exprs.append(pl.col("Tissue").is_in(tissues))
    if tissues_id:
        exprs.append(pl.col("TissueID").is_in(tissues_id))
    if contigs:
        exprs.append(pl.col("Contig").is_in(contigs))
    if pos_str:
        try:
            pos = int(pos_str)
        except ValueError:
            pos = 1
        exprs.append(pl.col("Pos") == pos)
    if end_str:
        try:
            end = int(end_str)
        except ValueError:
            end = 10000
        exprs.append(pl.col("End") == end)
    if trid_str:
        try:
            trid = int(trid_str)
        except ValueError:
            trid = 0
        exprs.append(pl.col("TRID") == trid)
    if gene_id.strip():
        exprs.append(
            pl.col("GeneID").str.to_lowercase().str.contains(gene_id.strip().lower(), literal=True)
        )
    if gene_name.strip():
        exprs.append(
            pl.col("GeneName").str.to_lowercase().str.contains(gene_name.strip().lower(), literal=True)
        )

    filtered = df.filter(exprs) if exprs else df
    return filtered, bool(exprs)


def plotter(displayed):
    st.title("Interactive Plotly Boxplot")

    # Create the boxplot
    genex = px.box(displayed.sort(by=("TissueID")), x="GeneName", y="TPM", color="Tissue", title="Expression by gene")
    # Display in Streamlit
    st.plotly_chart(genex, width='stretch')

    # Create the boxplot
    donox = px.box(displayed.sort(by=("Donor")), x="GeneName", y="TPM", color="Donor", title="Expression by donor")
    # Display in Streamlit
    st.plotly_chart(donox, width='stretch')


st.set_page_config(page_title="Gene Expression Explorer", layout="wide")
st.title("Gene Expression Explorer")

uploaded = ""
if len(sys.argv) > 1:
    _, uploaded = sys.argv
else:
    uploaded = "/run/media/lpaz/drive2/smaht/donors/rna/tr_vs_rna/p25_SMaHT_TR_v_RNA.pq"

df: pl.DataFrame = pl.DataFrame({})
print(f'loading: {uploaded} ...')
try:
    df = read_parquet(uploaded)
except Exception as exc:
    st.error(f"Could not read the parquet file: {exc}")
    st.stop()

with st.sidebar:
    st.header("Filters")
    donors = st.multiselect("Donor", sorted(df["Donor"].drop_nulls().unique().to_list()))
    tissues = st.multiselect("Tissue", sorted(df["Tissue"].drop_nulls().unique().to_list()))

    trid = st.text_input("TR ID", placeholder=0)
    gene_name = st.text_input("GeneName", placeholder="Substring match")
    
    contigs = st.multiselect("Contig", sorted(df["Contig"].drop_nulls().unique().to_list()))
    pos = st.text_input("Pos", placeholder=1)
    end = st.text_input("End", placeholder=1000)

    tissues_id = st.multiselect("TissueID", sorted(df["TissueID"].drop_nulls().unique().to_list()))
    gene_id = st.text_input("GeneID", placeholder="Substring match")

filtered, filters_active = apply_filters(
    df,
    donors,
    tissues,
    tissues_id,
    contigs,
    pos,
    end,
    trid,
    gene_id,
    gene_name,
)

if filters_active:
    displayed = filtered
    st.caption(f"Showing all {displayed.height:,} matching rows.")
else:
    displayed = df.head(20)
    st.caption(f"Showing the first 20 of {df.height:,} rows. Use the filters to fetch data.")


st.dataframe(displayed, width='stretch', hide_index=True, height='stretch')

csv_buffer = io.BytesIO()
displayed.write_csv(csv_buffer)
st.download_button(
    "Download shown data as CSV",
    data=csv_buffer.getvalue(),
    file_name="gene_expression.csv",
    mime="text/csv",
)

if st.button("Show plots"):
    plotter(displayed)
    # if fig is not None:
    #     st.plotly_chart(fig, wi=True)
