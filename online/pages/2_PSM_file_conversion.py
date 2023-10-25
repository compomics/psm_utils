"""psm_utils Streamlit-based web server."""

import os
from pathlib import Path
from tempfile import NamedTemporaryFile

import streamlit as st
from _base import StreamlitPage

from psm_utils.io import FILETYPES, READERS, WRITERS, _infer_filetype, convert


class StreamlitPageConvert(StreamlitPage):
    """Conversion page for psm_utils online Streamlit web server."""

    def _main_page(self):
        st.markdown(
            """
            ## PSM file conversion
            """
        )
        with st.form(key="main_form"):
            st.session_state["convert_input_file"] = st.file_uploader(
                label="Input PSM file",
            )
            st.session_state["convert_input_filetype"] = st.selectbox(
                label="Input PSM file type",
                options=["infer"] + list(READERS.keys()),
                index=0,
            )
            st.session_state["convert_output_filetype"] = st.selectbox(
                label="Output PSM file type",
                options=WRITERS.keys(),
                index=2,
            )

            submitted = st.form_submit_button("Convert")

        if submitted:
            if not self.state["convert_input_file"]:
                st.error("Input PSM file required.")
            else:
                self._convert()

    def _convert(self):
        with st.spinner("Converting PSM file..."):
            # Infer filetype if required
            if st.session_state["convert_input_filetype"] == "infer":
                st.session_state["convert_input_filetype"] = _infer_filetype(
                    st.session_state["convert_input_file"].name
                )

            try:
                # Write uploaded file to disk and convert
                with NamedTemporaryFile(mode="wb", delete=False) as tmp_file:
                    tmp_file.write(st.session_state["convert_input_file"].getvalue())
                    tmp_file.flush()
                    tmp_file.close()
                    convert(
                        input_filename=tmp_file.name,
                        output_filename="output_filename",
                        input_filetype=st.session_state["convert_input_filetype"],
                        output_filetype=st.session_state["convert_output_filetype"],
                    )
            except Exception as e:
                st.exception(e)
            else:
                st.success("PSM file successfully converted!")

                # Construct output filename with new extension
                output_filename = (
                    Path(st.session_state["convert_input_file"].name).stem
                    + FILETYPES[st.session_state["convert_output_filetype"]]["extension"]
                )

                # Open converted file in memory for download button
                with open("output_filename", "rb") as file:
                    is_downloaded = st.download_button(
                        "Download",
                        data=file,
                        file_name=output_filename,
                        mime="text/plain",
                    )
                os.remove("output_filename")


if __name__ == "__main__":
    StreamlitPageConvert()
