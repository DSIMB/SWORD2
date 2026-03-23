"""
Web server for SWORD3 predictions.

Provides a REST API and simple HTML visualization for protein domain
partitioning from sequence.

Usage:
    python -m sword2_dl.server --checkpoint best.pt --port 8000

API endpoints:
    POST /predict  - Submit sequence, get domain predictions
    GET  /health   - Health check
    GET  /         - Interactive web interface
"""

import io
import json
import logging
import os
from pathlib import Path

import numpy as np
import torch

logger = logging.getLogger(__name__)


def create_app(checkpoint_path: str, config_path: str | None = None):
    """Create FastAPI application with model loaded."""
    try:
        from fastapi import FastAPI, HTTPException
        from fastapi.responses import HTMLResponse, JSONResponse
        from pydantic import BaseModel, Field
    except ImportError:
        logger.error("FastAPI not installed. Install with: pip install fastapi uvicorn")
        raise

    from .config import Config
    from .pipeline import load_model, predict_from_embeddings, compute_embeddings_esm2

    # Load model at startup
    config = Config.from_yaml(config_path) if config_path else Config()
    model, config, device = load_model(checkpoint_path, config)
    logger.info("Model loaded successfully")

    app = FastAPI(
        title="SWORD3",
        description="Protein domain partitioning from sequence",
        version="1.0.0",
    )

    class PredictionRequest(BaseModel):
        sequence: str = Field(..., min_length=10, max_length=2000,
                              description="Protein sequence (amino acid letters)")
        protein_id: str = Field(default="query", description="Protein identifier")
        min_confidence: float = Field(default=0.1, ge=0, le=1)
        min_domain_size: int = Field(default=20, ge=5)

    class DomainSegment(BaseModel):
        start: int
        end: int

    class Domain(BaseModel):
        segments: list[DomainSegment]

    class Partitioning(BaseModel):
        domains: list[Domain]
        num_domains: int
        confidence: float

    class PredictionResponse(BaseModel):
        protein_id: str
        sequence_length: int
        partitionings: list[Partitioning]
        contact_map: list[list[float]] | None = None

    @app.get("/health")
    async def health():
        return {"status": "ok", "model": "sword3"}

    @app.post("/predict", response_model=PredictionResponse)
    async def predict(request: PredictionRequest):
        """Predict domain partitionings from protein sequence."""
        # Validate sequence
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        seq = request.sequence.upper()
        if not all(c in valid_aa for c in seq):
            raise HTTPException(400, "Invalid amino acid characters in sequence")

        try:
            # Compute embeddings
            embeddings = compute_embeddings_esm2(
                {request.protein_id: seq}, device=device
            )
            emb = embeddings[request.protein_id]

            # Pad to expected dim if needed
            expected_dim = config.model.total_embed_dim
            if emb.shape[-1] < expected_dim:
                pad = torch.zeros(emb.shape[0], expected_dim - emb.shape[-1])
                emb = torch.cat([emb, pad], dim=-1)

            # Predict
            result = predict_from_embeddings(
                model, emb, device,
                min_confidence=request.min_confidence,
                min_domain_size=request.min_domain_size,
            )

            # Format response
            partitionings = []
            for part in result["partitionings"]:
                domains = []
                for dom in part["domains"]:
                    segments = [DomainSegment(start=s["start"], end=s["end"]) for s in dom]
                    domains.append(Domain(segments=segments))
                partitionings.append(Partitioning(
                    domains=domains,
                    num_domains=part["num_domains"],
                    confidence=part["confidence"],
                ))

            contact_map = None
            if "contact_map" in result:
                contact_map = result["contact_map"].tolist()

            return PredictionResponse(
                protein_id=request.protein_id,
                sequence_length=len(seq),
                partitionings=partitionings,
                contact_map=contact_map,
            )

        except Exception as e:
            logger.exception("Prediction failed")
            raise HTTPException(500, f"Prediction failed: {str(e)}")

    @app.get("/", response_class=HTMLResponse)
    async def index():
        """Simple web interface for predictions."""
        return _get_html_interface()

    return app


def _get_html_interface() -> str:
    """Generate the HTML interface for the web server."""
    return """<!DOCTYPE html>
<html>
<head>
    <title>SWORD3: Protein Domain Prediction</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }
        h1 { color: #333; margin-bottom: 10px; }
        .subtitle { color: #666; margin-bottom: 20px; }
        .input-section { background: white; padding: 20px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
        textarea { width: 100%; height: 100px; font-family: monospace; padding: 10px; border: 1px solid #ddd; border-radius: 4px; resize: vertical; }
        button { background: #2563eb; color: white; padding: 10px 24px; border: none; border-radius: 4px; cursor: pointer; font-size: 16px; margin-top: 10px; }
        button:hover { background: #1d4ed8; }
        button:disabled { background: #93c5fd; cursor: not-allowed; }
        .results { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); display: none; }
        .partitioning { margin: 15px 0; padding: 15px; background: #f9fafb; border-radius: 4px; border-left: 4px solid #2563eb; }
        .domain { margin: 5px 0; padding: 5px 10px; border-radius: 3px; display: inline-block; margin-right: 8px; font-size: 14px; }
        .confidence { color: #666; font-size: 13px; }
        .domain-bar { height: 30px; display: flex; margin: 10px 0; border-radius: 4px; overflow: hidden; }
        .domain-segment { display: flex; align-items: center; justify-content: center; color: white; font-size: 11px; font-weight: bold; }
        .error { color: #dc2626; padding: 10px; background: #fef2f2; border-radius: 4px; }
        .loading { text-align: center; padding: 20px; color: #666; }
        canvas { max-width: 100%; margin-top: 15px; }
        label { font-weight: 600; display: block; margin-bottom: 5px; }
    </style>
</head>
<body>
    <h1>SWORD3</h1>
    <p class="subtitle">Predict protein domain partitionings from sequence</p>

    <div class="input-section">
        <label for="sequence">Protein Sequence:</label>
        <textarea id="sequence" placeholder="Enter protein sequence (e.g., MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH...)"></textarea>
        <button id="submit" onclick="predict()">Predict Domains</button>
    </div>

    <div id="results" class="results"></div>

    <script>
    const COLORS = ['#2563eb','#dc2626','#059669','#d97706','#7c3aed','#db2777','#0891b2','#65a30d','#ea580c','#6366f1'];

    async function predict() {
        const seq = document.getElementById('sequence').value.trim().replace(/[^A-Za-z]/g, '');
        if (seq.length < 10) { alert('Sequence too short (min 10 residues)'); return; }

        const btn = document.getElementById('submit');
        const results = document.getElementById('results');
        btn.disabled = true;
        btn.textContent = 'Predicting...';
        results.style.display = 'block';
        results.innerHTML = '<div class="loading">Computing embeddings and running model...</div>';

        try {
            const resp = await fetch('/predict', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({sequence: seq})
            });
            if (!resp.ok) throw new Error((await resp.json()).detail || 'Prediction failed');
            const data = await resp.json();
            renderResults(data, seq.length);
        } catch(e) {
            results.innerHTML = `<div class="error">Error: ${e.message}</div>`;
        } finally {
            btn.disabled = false;
            btn.textContent = 'Predict Domains';
        }
    }

    function renderResults(data, seqLen) {
        const el = document.getElementById('results');
        let html = `<h2>Results (${data.sequence_length} residues, ${data.partitionings.length} partitioning(s))</h2>`;

        data.partitionings.forEach((part, i) => {
            html += `<div class="partitioning">`;
            html += `<strong>Partitioning ${i+1}</strong> &mdash; ${part.num_domains} domain(s) `;
            html += `<span class="confidence">(confidence: ${part.confidence.toFixed(3)})</span>`;

            // Domain bar
            html += `<div class="domain-bar">`;
            part.domains.forEach((dom, di) => {
                dom.segments.forEach(seg => {
                    const width = ((seg.end - seg.start + 1) / seqLen * 100).toFixed(1);
                    const left = (seg.start / seqLen * 100).toFixed(1);
                    html += `<div class="domain-segment" style="width:${width}%;margin-left:${di===0&&seg===dom.segments[0]?left:'0'}%;background:${COLORS[di%COLORS.length]}" title="Domain ${di+1}: ${seg.start+1}-${seg.end+1}">D${di+1}</div>`;
                });
            });
            html += `</div>`;

            // Domain details
            part.domains.forEach((dom, di) => {
                const segs = dom.segments.map(s => `${s.start+1}-${s.end+1}`).join(' + ');
                html += `<span class="domain" style="background:${COLORS[di%COLORS.length]}20;border:1px solid ${COLORS[di%COLORS.length]}40;color:${COLORS[di%COLORS.length]}">Domain ${di+1}: ${segs}</span>`;
            });
            html += `</div>`;
        });

        // Contact map canvas
        if (data.contact_map) {
            html += `<h3 style="margin-top:20px">Predicted Contact Map</h3>`;
            html += `<canvas id="contactCanvas" width="400" height="400"></canvas>`;
        }

        el.innerHTML = html;

        // Draw contact map
        if (data.contact_map) {
            const canvas = document.getElementById('contactCanvas');
            const ctx = canvas.getContext('2d');
            const L = data.contact_map.length;
            const scale = Math.min(400 / L, 4);
            canvas.width = L * scale;
            canvas.height = L * scale;
            for (let i = 0; i < L; i++) {
                for (let j = 0; j < L; j++) {
                    const v = data.contact_map[i][j];
                    ctx.fillStyle = `rgb(${255-v*255},${255-v*255},${255})`;
                    ctx.fillRect(j*scale, i*scale, scale, scale);
                }
            }
        }
    }
    </script>
</body>
</html>"""


def main():
    """CLI entry point for web server."""
    import argparse

    parser = argparse.ArgumentParser(description="SWORD3 Web Server")
    parser.add_argument("--checkpoint", required=True, help="Model checkpoint")
    parser.add_argument("--config", default=None, help="Config YAML")
    parser.add_argument("--host", default="0.0.0.0", help="Host to bind")
    parser.add_argument("--port", type=int, default=8000, help="Port")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    app = create_app(args.checkpoint, args.config)

    import uvicorn
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
