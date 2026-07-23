import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from obstacles import Circle, Rectangle


def metric_scale_field(manifold, bounds, res=300):
    x0, x1, y0, y1 = bounds
    xs = np.linspace(x0, x1, res)
    ys = np.linspace(y0, y1, res)
    field = np.full((res, res), np.nan)
    
    for j in range(len(ys)):
        for i in range(len(xs)):
            p = np.array([xs[i], ys[j]])
            if manifold.in_domain(p):
                d = np.linalg.det(manifold.metric(p))
                if d > 0:
                    field[j, i] = np.sqrt(d)
    return field


def plot_background(ax, manifold, bounds, res=300, cmap="viridis", alpha=1.0, colorbar=False):
    field = metric_scale_field(manifold, bounds, res)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(alpha=0.0)
    img = ax.imshow(np.ma.masked_invalid(field), origin="lower", extent=[bounds[0], bounds[1], bounds[2], bounds[3]], cmap=cmap_obj, alpha=alpha, aspect="equal", zorder=0)
    
    if hasattr(manifold, "rmax"):
        ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec="black", lw=1.5, zorder=1))
        ax.add_patch(plt.Circle((0, 0), manifold.rmax, fill=False, ec="black", lw=1.0, ls="--", alpha=0.6, zorder=1))
    if colorbar:
        cb = ax.figure.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\sqrt{\det g}$")
    return img


def plot_obstacles(ax, obstacles, color="0.15", alpha=0.85, hatch=None):
    for obs in obstacles:
        if isinstance(obs, Circle):
            ax.add_patch(plt.Circle((obs.cx, obs.cy), obs.r, facecolor=color, edgecolor="white", lw=1.0, alpha=alpha, hatch=hatch, zorder=5))
        elif isinstance(obs, Rectangle):
            ax.add_patch(plt.Rectangle((obs.sx, obs.sy), obs.ex - obs.sx, obs.ey - obs.sy, facecolor=color, edgecolor="white", lw=1.0, alpha=alpha, hatch=hatch, zorder=5))


def plot_graph(ax, points, adj, node_size=4, node_color="white", edge_color="0.6", edge_alpha=0.25, edge_width=0.5, show_edges=True):
    points = np.asarray([np.asarray(p, dtype=float) for p in points])
    if show_edges:
        edges = [[points[v], points[u]] for v in range(len(adj)) for u in adj[v] if u > v]
        lines = LineCollection(edges, colors=edge_color, alpha=edge_alpha, linewidths=edge_width, zorder=2)
        ax.add_collection(lines)
    ax.scatter(points[:, 0], points[:, 1], s=node_size, c=node_color, zorder=3, edgecolors="none")


def plot_components(ax, points, comp, cmap="tab10", node_size=12):
    points = np.asarray([np.asarray(p, dtype=float) for p in points])
    scatter = ax.scatter(points[:, 0], points[:, 1], c=np.asarray(comp), cmap=cmap, s=node_size, zorder=3)
    return scatter


def plot_path(ax, points, path, color="red", lw=2.5, label=None, marker_ends=True, zorder=6):
    if not path:
        return
    points = np.asarray([np.asarray(points[i], dtype=float) for i in path])
    ax.plot(points[:, 0], points[:, 1], color=color, lw=lw, zorder=zorder, label=label, solid_capstyle="round", solid_joinstyle="round")
    
    if marker_ends:
        ax.scatter([points[0, 0]], [points[0, 1]], c="green", s=110, zorder=zorder + 1, edgecolors="black", linewidths=1.2)
        ax.scatter([points[-1, 0]], [points[-1, 1]], c="red", s=160, marker="*", zorder=zorder + 1, edgecolors="black", linewidths=1.0)


def plot_field(ax, xs, ys, field, manifold=None, cmap="magma", colorbar=False, label="diff field"):
    field = np.asarray(field, dtype=float).T
    if manifold is not None:
        keep = np.ones_like(field, dtype=bool)
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                if not manifold.in_domain(np.array([x, y])):
                    keep[j, i] = False
        field = np.ma.masked_where(~keep, field)
    
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(alpha=0.0)
    im = ax.imshow(field, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap=cm, aspect="equal", zorder=0)
    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(label)
    return im


def plot_separators(ax, index, max_depth=4, show_portals=True, cmap="turbo"):
    segs, depths = [], []
    for sid in range(index.num_seps):
        axis, sv, depth = index.sep_meta[sid]
        if depth > max_depth:
            continue
        x0, x1, y0, y1 = index.sep_bbox[sid]
        segs.append([(sv, y0), (sv, y1)] if axis == 0 else [(x0, sv), (x1, sv)])
        depths.append(depth)
        
    depths = np.array(depths)
    dmax = max(1, depths.max()) if len(depths) else 1
    colors = plt.get_cmap(cmap)(depths / dmax)
    lines = LineCollection(segs, colors=colors, linewidths=np.clip(3.0 - 0.5 * depths, 0.5, 3.0), zorder=4)
    ax.add_collection(lines)
    
    if show_portals:
        ids = [index.sep_portals[s] for s in range(index.num_seps) if index.sep_meta[s][2] <= max_depth]
        if ids:
            allp = np.unique(np.concatenate(ids))
            ax.scatter(index.points[allp, 0], index.points[allp, 1], s=16, c="red", zorder=5, edgecolors="black", linewidths=0.3)


def plot_polyline(ax, coords, color="black", lw=2.0, ls="--", label=None, zorder=6):
    coords = np.asarray(coords, dtype=float)
    ax.plot(coords[:, 0], coords[:, 1], color=color, lw=lw, ls=ls, label=label, zorder=zorder)


def visualize_query(manifold, bounds, points, adj=None, obstacles=(), path=None, title=None, background=True, show_graph=True, figsize=(7, 7), ax=None, colorbar=False, **bg_kw):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    if background:
        plot_background(ax, manifold, bounds, colorbar=colorbar, **bg_kw)
    if obstacles:
        plot_obstacles(ax, obstacles)
    if show_graph and adj is not None:
        plot_graph(ax, points, adj)
    if path:
        plot_path(ax, points, path)
        
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if title:
        ax.set_title(title)
    return fig, ax
