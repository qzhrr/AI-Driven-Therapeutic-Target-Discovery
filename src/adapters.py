"""
adapters.py — portable ``host`` shim for running the pipeline outside Claude Science.

The pipeline was originally executed inside Claude Science, where a ``host`` object
provided two capabilities that are not plain Python:

  * ``host.llm(request | [requests], max_concurrency=N)`` — batched LLM calls used
    for weight selection (step 9), target nomination (step 10) and the gene-masked
    ablation (step 11a).
  * ``host.mcp(server, method, **kwargs)`` — connector queries used during data
    acquisition (GTEx expression) and tractability (ChEMBL, Open Targets).

This module reproduces that surface so the reconstructed scripts run unchanged.

Design
------
* ``host.llm`` is implemented faithfully against the public Anthropic Messages API
  (``pip install anthropic``; set ``ANTHROPIC_API_KEY``). The request/response
  shapes match what the original code consumes: each request is a dict with
  ``prompt``, ``system``, ``tools``, ``tool_choice``, ``model``, ``max_tokens``;
  each result is ``{"text", "tool_use": {"name", "input"}, "content", "model",
  "usage", "stop_reason"}`` or ``{"error": ...}``.

* ``host.mcp`` is connector-mediated on the original platform. Those connectors
  are not public packages, so live connector mode is not reproduced here. Every
  connector-dependent step is written cache-first (see ``common.use_cache``): with
  the committed ``data/cache/`` responses the pipeline reproduces fully offline.
  Calling ``host.mcp`` in live mode raises a clear, actionable error.

Only ``host.llm`` and ``host.current_model`` are needed to *recompute* the
stochastic LLM steps live; everything else reproduces from cache.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from typing import Any


# ---------------------------------------------------------------------------
# Model selection
# ---------------------------------------------------------------------------
# The original runs used the session model (Claude Opus). Downstream code only
# needs *a* capable model id; make it configurable and default to a current one.
_DEFAULT_MODEL = os.environ.get("ANTHROPIC_MODEL", "claude-opus-4-20250514")


class _Host:
    """Minimal re-implementation of the Claude Science ``host`` object."""

    # -- model ------------------------------------------------------------
    def current_model(self) -> str:
        return _DEFAULT_MODEL

    def reasoning_model(self) -> str:
        return _DEFAULT_MODEL

    # -- LLM --------------------------------------------------------------
    def llm(self, request_or_list, max_concurrency: int = 8):
        """
        Run one request (dict) or many (list of dicts) against the Anthropic API.

        Mirrors ``host.llm``: a single request returns a single result dict; a
        list returns a positionally-matched list. Failures are returned in-band
        as ``{"error": "..."}`` so batch callers can count errors without the
        whole batch aborting.
        """
        single = isinstance(request_or_list, dict)
        requests = [request_or_list] if single else list(request_or_list)

        try:
            import anthropic  # deferred: only needed for live LLM calls
        except ImportError as e:  # pragma: no cover
            raise ImportError(
                "host.llm requires the 'anthropic' package. Install it "
                "(pip install anthropic) and set ANTHROPIC_API_KEY, or run the "
                "pipeline in cached mode (the default) to reproduce the "
                "committed LLM outputs without any API calls."
            ) from e

        client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from the env

        def _one(req: dict) -> dict:
            try:
                return self._call_one(client, req)
            except Exception as exc:  # noqa: BLE001 - surface as in-band error
                return {"error": f"{type(exc).__name__}: {exc}"}

        if len(requests) == 1:
            results = [_one(requests[0])]
        else:
            with ThreadPoolExecutor(max_workers=max(1, int(max_concurrency))) as ex:
                results = list(ex.map(_one, requests))

        return results[0] if single else results

    @staticmethod
    def _call_one(client, req: dict) -> dict:
        """Translate one host.llm request dict into an Anthropic Messages call."""
        # A host.llm request accepts either 'prompt' (str) or 'messages' (list).
        messages = req.get("messages")
        if messages is None:
            messages = [{"role": "user", "content": req["prompt"]}]

        kwargs: dict[str, Any] = {
            "model": req.get("model", _DEFAULT_MODEL),
            "max_tokens": req.get("max_tokens", 1024),
            "messages": messages,
        }
        if req.get("system"):
            kwargs["system"] = req["system"]
        if req.get("tools"):
            kwargs["tools"] = req["tools"]
        if req.get("tool_choice"):
            kwargs["tool_choice"] = req["tool_choice"]
        if req.get("temperature") is not None:
            kwargs["temperature"] = req["temperature"]

        msg = client.messages.create(**kwargs)

        # Normalise the response to the shape the original code consumes.
        text = ""
        tool_use = None
        content = []
        for block in msg.content:
            btype = getattr(block, "type", None)
            if btype == "text":
                text += block.text
                content.append({"type": "text", "text": block.text})
            elif btype == "tool_use":
                tool_use = {"name": block.name, "input": block.input}
                content.append(
                    {"type": "tool_use", "name": block.name, "input": block.input}
                )

        return {
            "text": text,
            "tool_use": tool_use,
            "content": content,
            "model": msg.model,
            "usage": {
                "input_tokens": msg.usage.input_tokens,
                "output_tokens": msg.usage.output_tokens,
            },
            "stop_reason": msg.stop_reason,
        }

    # -- MCP connectors ---------------------------------------------------
    def mcp(self, server: str, method: str, **kwargs):
        """
        Connector-mediated queries are platform-specific and not reproduced here.

        Every connector-dependent step is cache-first, so the committed
        ``data/cache/`` responses reproduce these steps offline. This method only
        fires in live mode (GBM_LIVE=1) when a cache is missing.
        """
        raise NotImplementedError(
            f"host.mcp({server!r}, {method!r}) is a Claude Science connector and "
            "is not reproduced in this standalone repo. The committed cache in "
            "data/cache/ reproduces every connector-dependent step offline (the "
            "default cached mode). To recompute this specific query live, replace "
            "this call with the connector's public REST equivalent "
            "(GTEx portal API / EBI ChEMBL API / Open Targets GraphQL)."
        )


# The singleton the pipeline scripts import: ``from adapters import host``.
host = _Host()
