import anndata as ad
import numpy as np
import pandas as pd
import pytest
import zarr

ad.settings.allow_write_nullable_strings = False

N_OBS = 10
N_VAR = 10


class TestClass:
    @pytest.fixture(scope="class")
    def adata(self):

        obs = pd.DataFrame(
            {
                "id": pd.Series([f"cell{i}" for i in range(N_OBS)]).astype("object"),
                "integer": pd.Series([i for i in range(N_OBS)]).astype("int32"),
                "float": pd.Series([float(i) for i in range(N_OBS)]).astype("float32"),
                "boolean": pd.Series(np.random.rand(N_OBS) > 0.5).astype("bool"),
                "categorical": pd.Series([f"cat{i%3}" for i in range(N_OBS)])
                .astype("str")
                .astype("object"),
            },
        )
        # Avoid pandas string array
        obs.index = obs.index.astype("str").astype("object")
        obs["integer_with_nan"] = obs["integer"]
        obs.loc["0", "integer_with_nan"] = np.nan
        obs["float_with_nan"] = obs["float"]
        obs.loc["0", "float_with_nan"] = np.nan
        obs["categorical_with_nan"] = obs["categorical"]
        obs.loc["0", "categorical_with_nan"] = np.nan

        var = pd.DataFrame(
            {
                "id": pd.Series([f"gene{i}" for i in range(N_VAR)]).astype("object"),
            },
        )
        # Avoid pandas string array
        var.index = var.index.astype("str").astype("object")

        uns = {
            "categorical_colors": ["#ff0000", "#00ff00", "#0000ff"],
        }

        adata = ad.AnnData(np.ones((N_OBS, N_VAR)), obs=obs, var=var, uns=uns)

        return adata

    @pytest.fixture(scope="class")
    def anndata_zarr_v2(self, adata, tmp_path_factory):
        ad.settings.zarr_write_format = 2
        fn = tmp_path_factory.mktemp("data") / "anndata.zarr"
        adata.write_zarr(fn)
        return fn

    @pytest.fixture(scope="class")
    def anndata_zarr_v3(self, adata, tmp_path_factory):
        ad.settings.zarr_write_format = 3
        fn = tmp_path_factory.mktemp("data") / "anndata.zarr"
        adata.write_zarr(fn)
        return fn

    def test_get_obs_col_names(self, anndata_zarr_v2, anndata_zarr_v3):
        from cherita.dataset.metadata import get_obs_col_names

        for zarr_path in [anndata_zarr_v2, anndata_zarr_v3]:
            adata_group = zarr.open_group(zarr_path, mode="r")
            obs_col_names = get_obs_col_names(adata_group)
            assert set(obs_col_names) == {
                "id",
                "integer",
                "float",
                "boolean",
                "categorical",
                "integer_with_nan",
                "float_with_nan",
                "categorical_with_nan",
            }

    def test_get_obs_col_metadata(self, anndata_zarr_v2, anndata_zarr_v3):
        from cherita.dataset.metadata import get_obs_cols_metadata

        for zarr_path in [anndata_zarr_v2, anndata_zarr_v3]:
            adata_group = zarr.open_group(zarr_path, mode="r")
            obs_cols_metadata = get_obs_cols_metadata(adata_group)
            metadata = {m["name"]: m for m in obs_cols_metadata}
            assert len(obs_cols_metadata) == 8
            assert metadata["id"]["type"] == "discrete"
            assert metadata["integer"]["type"] == "continuous"
            assert metadata["float"]["type"] == "continuous"
            assert metadata["boolean"]["type"] == "boolean"
            assert metadata["categorical"]["type"] == "categorical"
            assert len(metadata["id"]["values"]) == metadata["id"]["n_values"] == N_OBS
            for col in ["integer_with_nan", "float_with_nan", "categorical_with_nan"]:
                assert "undefined" in metadata[col]["codes"]
                assert metadata[col]["codes"]["undefined"] == -1
                assert "-1" in metadata[col]["codesMap"]
                assert metadata[col]["codesMap"]["-1"] == "undefined"

    def test_obs_bin_data(self, anndata_zarr_v2, anndata_zarr_v3):
        from cherita.dataset.metadata import get_obs_col_metadata

        for zarr_path in [anndata_zarr_v2, anndata_zarr_v3]:
            adata_group = zarr.open_group(zarr_path, mode="r")
            metadata = get_obs_col_metadata(adata_group, col="integer")
            assert metadata["name"] == "integer"
            assert metadata["type"] == "continuous"
            assert metadata["bins"]["nBins"] == len(metadata["bins"]["binEdges"]) == 5
            assert metadata["mean"] == sum(range(N_OBS)) / N_OBS
            assert metadata["n_unique"] == N_OBS

    def test_obs_colors(self, anndata_zarr_v2, anndata_zarr_v3):
        from cherita.dataset.metadata import get_obs_col_metadata

        for zarr_path in [anndata_zarr_v2, anndata_zarr_v3]:
            adata_group = zarr.open_group(zarr_path, mode="r")
            metadata = get_obs_col_metadata(adata_group, "categorical")
            assert metadata["colors"] == ["#ff0000", "#00ff00", "#0000ff"]

    def test_get_var_col_names(self, anndata_zarr_v2, anndata_zarr_v3):
        from cherita.dataset.metadata import get_var_col_names

        for zarr_path in [anndata_zarr_v2, anndata_zarr_v3]:
            adata_group = zarr.open_group(zarr_path, mode="r")
            var_col_names = get_var_col_names(adata_group)
            assert set(var_col_names) == {"id"}

    def test_search_var_names(self, anndata_zarr_v2, anndata_zarr_v3):
        from cherita.dataset.search import search_var_names

        for zarr_path in [anndata_zarr_v2, anndata_zarr_v3]:
            adata_group = zarr.open_group(zarr_path, mode="r")
            results = search_var_names(adata_group, col="id", text="gene1")
            assert len(results) == 1
            assert results[0]["name"] == "gene1"
