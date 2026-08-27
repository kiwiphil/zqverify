"""Oriented finite cubical lattices with exact incidence data."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from itertools import combinations, product
from typing import Iterable, Iterator, Sequence

Vertex = tuple[int, ...]


@dataclass(frozen=True, order=True)
class Edge:
    """An edge oriented in the positive ``axis`` direction from ``base``."""

    base: Vertex
    axis: int

    @property
    def head(self) -> Vertex:
        point = list(self.base)
        point[self.axis] += 1
        return tuple(point)


@dataclass(frozen=True, order=True)
class Plaquette:
    """A positively oriented coordinate plaquette with axes ``i < j``."""

    base: Vertex
    i: int
    j: int

    def __post_init__(self) -> None:
        if self.i >= self.j:
            raise ValueError("Plaquette axes must satisfy i < j")


@dataclass(frozen=True)
class OrientedBox:
    """A free-boundary cubical box.

    ``lengths[a]`` is the number of unit cells along coordinate axis ``a``.
    Vertices therefore have coordinates from 0 through ``lengths[a]``.
    """

    lengths: tuple[int, ...]

    def __post_init__(self) -> None:
        if len(self.lengths) < 2:
            raise ValueError("The lattice dimension must be at least 2")
        if any(length < 1 for length in self.lengths):
            raise ValueError("Every box length must be positive")

    @property
    def dimension(self) -> int:
        return len(self.lengths)

    @cached_property
    def vertices(self) -> tuple[Vertex, ...]:
        ranges = [range(length + 1) for length in self.lengths]
        return tuple(product(*ranges))

    @cached_property
    def edges(self) -> tuple[Edge, ...]:
        result: list[Edge] = []
        for axis, length in enumerate(self.lengths):
            ranges = [range(value + 1) for value in self.lengths]
            ranges[axis] = range(length)
            for base in product(*ranges):
                result.append(Edge(tuple(base), axis))
        return tuple(sorted(result))

    @cached_property
    def plaquettes(self) -> tuple[Plaquette, ...]:
        result: list[Plaquette] = []
        for i, j in combinations(range(self.dimension), 2):
            ranges = [range(value + 1) for value in self.lengths]
            ranges[i] = range(self.lengths[i])
            ranges[j] = range(self.lengths[j])
            for base in product(*ranges):
                result.append(Plaquette(tuple(base), i, j))
        return tuple(sorted(result))

    @cached_property
    def vertex_index(self) -> dict[Vertex, int]:
        return {vertex: index for index, vertex in enumerate(self.vertices)}

    @cached_property
    def edge_index(self) -> dict[Edge, int]:
        return {edge: index for index, edge in enumerate(self.edges)}

    @cached_property
    def plaquette_index(self) -> dict[Plaquette, int]:
        return {plaquette: index for index, plaquette in enumerate(self.plaquettes)}

    def _shift(self, point: Vertex, axis: int, amount: int = 1) -> Vertex:
        shifted = list(point)
        shifted[axis] += amount
        return tuple(shifted)

    def plaquette_boundary(self, plaquette: Plaquette) -> tuple[tuple[Edge, int], ...]:
        """Return the four oriented edges in ``partial plaquette``.

        For ``i < j`` the convention is
        ``+i, +j, -i, -j`` around the plaquette.
        """

        x = plaquette.base
        i, j = plaquette.i, plaquette.j
        return (
            (Edge(x, i), +1),
            (Edge(self._shift(x, i), j), +1),
            (Edge(self._shift(x, j), i), -1),
            (Edge(x, j), -1),
        )

    def edge_boundary(self, edge: Edge) -> tuple[tuple[Vertex, int], ...]:
        return ((edge.base, -1), (edge.head, +1))

    @cached_property
    def boundary2_matrix(self) -> tuple[tuple[int, ...], ...]:
        """Integer matrix of ``partial_2: C_2 -> C_1``.

        Rows are edges and columns are plaquettes.
        """

        matrix = [[0 for _ in self.plaquettes] for _ in self.edges]
        for p_index, plaquette in enumerate(self.plaquettes):
            for edge, sign in self.plaquette_boundary(plaquette):
                matrix[self.edge_index[edge]][p_index] = sign
        return tuple(tuple(row) for row in matrix)

    @cached_property
    def boundary1_matrix(self) -> tuple[tuple[int, ...], ...]:
        """Integer matrix of ``partial_1: C_1 -> C_0``.

        Rows are vertices and columns are edges.
        """

        matrix = [[0 for _ in self.edges] for _ in self.vertices]
        for e_index, edge in enumerate(self.edges):
            for vertex, sign in self.edge_boundary(edge):
                matrix[self.vertex_index[vertex]][e_index] = sign
        return tuple(tuple(row) for row in matrix)

    def incident_plaquette_indices(self, edge_index: int) -> tuple[int, ...]:
        row = self.boundary2_matrix[edge_index]
        return tuple(index for index, value in enumerate(row) if value != 0)

    def planar_filling_indices(
        self,
        origin: Sequence[int],
        axes: tuple[int, int],
        sizes: tuple[int, int],
    ) -> tuple[int, ...]:
        """Indices of a rectangular planar filling."""

        origin_tuple = tuple(origin)
        if len(origin_tuple) != self.dimension:
            raise ValueError("Origin has the wrong dimension")
        i, j = sorted(axes)
        if i == j:
            raise ValueError("Rectangle axes must be distinct")
        size_i, size_j = sizes if axes == (i, j) else (sizes[1], sizes[0])
        if size_i < 1 or size_j < 1:
            raise ValueError("Rectangle sizes must be positive")
        if origin_tuple[i] + size_i > self.lengths[i]:
            raise ValueError("Rectangle exceeds the box along its first axis")
        if origin_tuple[j] + size_j > self.lengths[j]:
            raise ValueError("Rectangle exceeds the box along its second axis")

        indices: list[int] = []
        for a in range(size_i):
            for b in range(size_j):
                base = list(origin_tuple)
                base[i] += a
                base[j] += b
                indices.append(self.plaquette_index[Plaquette(tuple(base), i, j)])
        return tuple(indices)

    def rectangular_loop(
        self,
        origin: Sequence[int],
        axes: tuple[int, int],
        sizes: tuple[int, int],
    ) -> tuple[int, ...]:
        """Return the oriented 1-chain bounding a planar rectangle.

        The chain is computed as the boundary of the planar filling, ensuring
        consistency with the incidence convention.
        """

        coefficients = [0] * len(self.plaquettes)
        for index in self.planar_filling_indices(origin, axes, sizes):
            coefficients[index] = 1
        return self.boundary2(coefficients)

    def boundary2(self, plaquette_coefficients: Sequence[int]) -> tuple[int, ...]:
        if len(plaquette_coefficients) != len(self.plaquettes):
            raise ValueError("Wrong number of plaquette coefficients")
        return tuple(
            sum(row[p] * plaquette_coefficients[p] for p in range(len(self.plaquettes)))
            for row in self.boundary2_matrix
        )

    def boundary1(self, edge_coefficients: Sequence[int]) -> tuple[int, ...]:
        if len(edge_coefficients) != len(self.edges):
            raise ValueError("Wrong number of edge coefficients")
        return tuple(
            sum(row[e] * edge_coefficients[e] for e in range(len(self.edges)))
            for row in self.boundary1_matrix
        )

    def signed_incident_types(self, edge_index: int) -> tuple[tuple[int, int], ...]:
        """Deterministic order of incident signed plaquette types.

        Each item is ``(plaquette_index, sign)`` with sign ``+1`` listed before
        sign ``-1`` for each plaquette.
        """

        result: list[tuple[int, int]] = []
        for p_index in self.incident_plaquette_indices(edge_index):
            result.append((p_index, +1))
            result.append((p_index, -1))
        return tuple(result)
