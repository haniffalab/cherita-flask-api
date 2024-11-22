from __future__ import annotations
import numpy as np
import zarr
from dataclasses import dataclass
from typing import Union, Callable
from cherita.utils.adata_utils import get_group_index, get_index_in_array
from cherita.resources.errors import (
    InvalidKey,
    InvalidVar,
)


@dataclass
class Marker:
    """
    Represents a marker in the data.

    Attributes:
        index (Union[int, list[int]]): The index or indices of the marker.
        name (str): The name of the marker.
        matrix_index (Union[int, list[int]]): The matrix index or indices of the marker.
        adata_group (zarr.Group): The zarr group containing the data.
        isSet (bool): Indicates whether the marker is a set of markers or not.
        _aggregation_function (Callable[[np.ndarray]]): The aggregation function to
            apply to the data.
    """

    index: Union[int, list[int]]
    name: str
    matrix_index: Union[int, list[int]]
    adata_group: zarr.Group
    isSet: bool = False
    _aggregation_function: Callable[[np.ndarray]] = None

    @classmethod
    def from_any(
        cls,
        adata_group: zarr.Group,
        index_or_varset: Union[int, str, dict],
        var_names_col: str = None,
        aggregation_func: Callable[[np.ndarray], np.ndarray] = lambda x: np.mean(
            x, axis=0
        ),
    ):
        """
        Initialize a Marker instance by determining the type of input and calling the
        appropriate class method.

        Args:
            adata_group (zarr.Group): The zarr group containing the data.
            index_or_varset (Union[int, str, dict]): The index, name, or varset of
                the marker.
            var_names_col (str, optional): The name of the column containing the
                variable names. Defaults to None.
            aggregation_func (Callable[[np.ndarray], np.ndarray], optional): The
                aggregation function to apply to the data. Defaults to np.mean.
        """
        if isinstance(index_or_varset, (int, str)):
            return cls.from_index(adata_group, index_or_varset, var_names_col)
        elif isinstance(index_or_varset, dict):
            return cls.from_varset(
                adata_group, index_or_varset, var_names_col, aggregation_func
            )
        else:
            raise InvalidVar(f"Invalid input type: {type(index_or_varset)}")

    @classmethod
    def from_index(
        cls, adata_group: zarr.Group, index: Union[int, str], var_names_col: str = None
    ) -> Marker:
        """
        Create a Marker instance from an index.

        Args:
            adata_group (zarr.Group): The zarr group containing the data.
            index (Union[int, str]): The index or name of the marker.
            var_names_col (str, optional): The name of the column containing the
                variable names. Defaults to None.

        Returns:
            Marker: An instance of the Marker class.

        Raises:
            InvalidVar: If the index is invalid.
        """
        if isinstance(index, int):
            matrix_index = index
            try:
                index = get_group_index(adata_group.var)[matrix_index]
            except:
                raise InvalidVar(f"Invalid feature index {index}")
        elif isinstance(index, str):
            index = index
            try:
                matrix_index = get_index_in_array(
                    get_group_index(adata_group.var), index
                )
            except InvalidKey:
                raise InvalidVar(f"Invalid feature index {index}")
        else:
            raise InvalidVar(f"Invalid feature type {type(index)}")

        return cls(
            index=index,
            name=(
                adata_group.var[var_names_col][matrix_index] if var_names_col else index
            ),
            matrix_index=matrix_index,
            isSet=False,
            adata_group=adata_group,
        )

    @classmethod
    def from_varset(
        cls,
        adata_group: zarr.Group,
        varset: dict,
        var_names_col: str = None,
        aggregation_func: Callable[[np.ndarray]] = lambda x: np.mean(x, axis=0),
    ) -> Marker:
        """
        Create a Marker instance from a set of indices.

        Args:
            adata_group (zarr.Group): The zarr group containing the data.
            varset (dict): A dictionary containing a string `name` and an array
                `indices` containing the indices of multiple vars.
            var_names_col (str, optional): The name of the column containing the
                variable names. Defaults to None.
            aggregation_func (Callable[[np.ndarray]], optional): The aggregation
                function to apply to the data. Defaults to np.mean.

        Returns:
            Marker: An instance of the Marker class.
        """
        name = varset["name"]
        markers = [
            Marker.from_index(adata_group, i, var_names_col) for i in varset["indices"]
        ]

        marker_instance = cls(
            index=[m.index for m in markers],
            name=name,
            matrix_index=[m.matrix_index for m in markers],
            isSet=True,
            adata_group=adata_group,
        )
        marker_instance._aggregation_function = aggregation_func

        return marker_instance

    def get_X_at(self, indices: list[int] = None) -> np.ndarray:
        """
        Get the data associated with the marker at the specified indices.

        Args:
            indices (list[int], optional): The indices to get the data for.
                Defaults to None.

        Returns:
            np.ndarray: The data associated with the marker at the specified indices.
        """
        if indices is not None and not len(indices):
            return np.array([])
        if self.isSet:
            # return all data for each marker in the set instead of aggregated data
            return [
                self.adata_group.X[indices if indices is not None else slice(None), i]
                for i in self.matrix_index
            ]
        else:
            return self.adata_group.X[
                indices if indices is not None else slice(None), self.matrix_index
            ]

    @property
    def X(self) -> np.ndarray:
        """
        Get the data associated with the marker.

        Returns:
            np.ndarray: The data associated with the marker.
        """
        if self.isSet:
            if self._aggregation_function is None:
                raise ValueError(
                    "Aggregation function is not set for a set of markers."
                )
            return self._aggregation_function(
                [self.adata_group.X[:, i] for i in self.matrix_index]
            )
        else:
            return self.adata_group.X[:, self.matrix_index]
