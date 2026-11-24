// NanoFusion page controller: parameter form, jobs, file open/save, frame player.
import { ScatterViewer } from './viewer.js';

const $ = (selector) => document.querySelector(selector);
const form = $('#params');
const viewer = new ScatterViewer($('#plot'));
const ui = {
  run: $('#run-btn'), errors: $('#form-errors'), count: $('#particle-count'),
  job: $('#job'), jobLabel: $('#job-label'), jobProgress: $('#job-progress'), cancel: $('#cancel-btn'),
  fusion: $('#fusion'), fusionDetail: $('#fusion-detail'), gotoFusion: $('#goto-fusion'),
  open: $('#open-btn'), openInput: $('#open-input'), save: $('#save-link'),
  frameSlider: $('#frame-slider'), frameLabel: $('#frame-label'),
  restart: $('#restart'), prev: $('#prev'), next: $('#next'), last: $('#last'), animate: $('#animate'),
  fps: $('#fps'), fpsOut: $('#fps-out'), zslice: $('#zslice'), zOut: $('#z-out'),
  simTime: $('#sim-time'), compTime: $('#comp-time'), empty: $('#empty-hint'), status: $('#status'),
  about: $('#about'), aboutBody: $('#about-body'),
  report: $('#report'), reportBody: $('#report-body'), autofit: $('#autofit'),
  impIons: $('#imp-ions'), impSizeLabel: $('#imp-size-label'),
};

let meta = null;        // metadata of the result being shown
let summary = null;     // latest {id, name, canSave, saved} from /api/state
let current = 0;        // frame index
let pollTimer = null;
let lastJobReported = null;
const anim = { running: false, last: 0 };

// --- API --------------------------------------------------------------------

async function api(path, options = {}) {
  const response = await fetch(path, options);
  if (!response.ok) {
    let error = new Error(`${response.status} ${response.statusText}`);
    if ((response.headers.get('Content-Type') || '').includes('json')) {
      const body = await response.json();
      error = new Error(body.error);
      error.errors = body.errors;
    }
    throw error;
  }
  return response;
}

const getJSON = async (path) => (await api(path)).json();
const post = (path, body, headers = {}) =>
  api(path, { method: 'POST', body, headers: { 'X-NanoFusion': '1', ...headers } });

// --- status line -------------------------------------------------------------

function setStatus(text, kind = '') {
  ui.status.textContent = text;
  ui.status.className = `status ${kind}`;
}

// --- form ----------------------------------------------------------------------

function setForm(values) {
  for (const [name, value] of Object.entries(values)) {
    const field = form.elements.namedItem(name);
    if (!field) continue;
    if (field.type === 'checkbox') field.checked = Boolean(value);
    else field.value = value;
  }
  updateFormState();
}

function readForm() {
  const values = {};
  for (const field of form.elements) {
    if (!field.name) continue;
    if (field.type === 'radio') {
      if (field.checked) values[field.name] = field.value;
    } else {
      values[field.name] = field.type === 'checkbox' ? field.checked : field.value;
    }
  }
  return values;
}

const currentModel = () => form.elements.model.value;

function updateFormState() {
  const model = currentModel();
  for (const el of form.querySelectorAll('[data-model]')) el.hidden = el.dataset.model !== model;
  form.elements.impCore.disabled = form.elements.impProfile.value !== 'expanded';
  form.elements.impExpansion.disabled = form.elements.impProfile.value !== 'expanded';
  form.elements.impGap.disabled = form.elements.impGeometry.value !== 'pair';
  const sphere = form.elements.impShape.value === 'sphere';
  ui.impSizeLabel.textContent = sphere ? 'Radius (nm)' : 'Ions per edge';
  const size = Number.parseFloat(form.elements.impSize.value);
  if (size > 0) {
    // ~ ions in a sphere of LiD: volume / (2.04 Å)^3
    const pair = form.elements.impGeometry.value === 'pair' ? 2 : 1;
    const ions = pair * (sphere ? Math.round((4 / 3) * Math.PI * (size * 10 / 2.04) ** 3) : Math.round(size) ** 3);
    ui.impIons.value = `${sphere ? '≈ ' : ''}${ions.toLocaleString()}${ions > 20000 ? ' (slow: minutes)' : ''}`;
  } else {
    ui.impIons.value = '–';
  }
  const logTime = form.elements.logTime.checked;
  form.elements.transition.disabled = !logTime;
  form.elements.damping.disabled = !logTime;
  const n = ['sizeX', 'sizeY', 'sizeZ'].map((k) => Number.parseInt(form.elements[k].value, 10));
  ui.count.value = n.every((v) => v > 0) ? (n[0] * n[1] * n[2]).toLocaleString() : '–';
}

function showErrors(errors) {
  ui.errors.replaceChildren();
  ui.errors.hidden = errors.length === 0;
  if (!errors.length) return;
  const list = document.createElement('ul');
  for (const message of errors) {
    const item = document.createElement('li');
    item.textContent = message;
    list.append(item);
  }
  ui.errors.append(list);
}

form.addEventListener('input', updateFormState);
form.addEventListener('submit', async (event) => {
  event.preventDefault();
  showErrors([]);
  try {
    await post('/api/run', JSON.stringify(readForm()), { 'Content-Type': 'application/json' });
    poll();
  } catch (error) {
    showErrors(error.errors || [error.message]);
  }
});

// --- jobs and polling ---------------------------------------------------------------

function renderJob(job) {
  const running = job?.state === 'running';
  ui.job.hidden = !running;
  ui.run.disabled = running;
  ui.open.disabled = running;
  if (running) {
    ui.jobLabel.textContent = `${job.label}…`;
    if (job.total > 0) ui.jobProgress.value = Math.round((1000 * job.done) / job.total);
    else ui.jobProgress.removeAttribute('value');
    ui.cancel.disabled = false;
  } else if (job && job.message) {
    const key = `${job.label}|${job.state}|${job.elapsed}`;
    if (key !== lastJobReported) {
      lastJobReported = key;
      const kind = job.state === 'failed' ? 'error' : '';
      setStatus(job.message, kind);
    }
  }
}

async function poll() {
  clearTimeout(pollTimer);
  let state;
  try {
    state = await getJSON('/api/state');
  } catch {
    setStatus('Lost contact with the NanoFusion server. Is it still running?', 'error');
    pollTimer = setTimeout(poll, 2000);
    return;
  }
  renderJob(state.job);
  summary = state.result;
  updateSave();
  if (summary && summary.id !== meta?.id) await loadResult();
  if (state.job?.state === 'running') pollTimer = setTimeout(poll, 250);
}

ui.cancel.addEventListener('click', async () => {
  ui.cancel.disabled = true;
  await post('/api/cancel', '');
  poll();
});

// --- results ------------------------------------------------------------------------

async function loadResult() {
  const info = await getJSON('/api/result');
  const response = await api('/api/frames');
  if (response.headers.get('X-Result-Id') !== String(info.id)) return loadResult(); // changed meanwhile
  const firstLoad = meta === null;
  const positions = new Float32Array(await response.arrayBuffer());
  stopAnimation();
  meta = info;
  viewer.setData({
    positions,
    particles: info.particles,
    isLi: Uint8Array.from(info.isLi),
    layer: Int32Array.from(info.layer),
    lo: info.bounds.lo,
    hi: info.bounds.hi,
    name: info.name,
  });
  ui.empty.hidden = true;
  ui.frameSlider.max = String(info.frames);
  ui.frameSlider.disabled = false;
  ui.zslice.max = String(Math.max(0, info.layers - 1));
  ui.zslice.value = '0';
  viewer.setLayer(0);
  ui.animate.disabled = info.frames < 2;
  // Show the parameters behind the result, unless the user just ran it from this form.
  if (info.origin === 'file' || firstLoad) setForm(info.form);
  ui.autofit.checked = info.model === 'improved';
  viewer.setAutoFit(ui.autofit.checked);
  renderFusion();
  renderReport();
  setFrame(0);
  if (info.warning) setStatus(info.warning, 'warning');
}

const sci = (x, digits = 2) => {
  if (x === null || x === undefined || Number.isNaN(x)) return '–';
  if (x === 0) return '0';
  const a = Math.abs(x);
  if (a >= 0.01 && a < 1e5) return x.toPrecision(digits + 1).replace(/\.?0+$/, '');
  const [m, e] = x.toExponential(digits).split('e');
  return `${m}×10${superscript(Number(e))}`;
};
const superscript = (n) => String(n).replace(/[-0-9]/g, (c) => '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'['-0123456789'.indexOf(c)]);

function renderFusion() {
  const fusion = meta?.fusion;
  ui.gotoFusion.hidden = true;
  if (meta?.report) {
    const f = meta.report.fusion;
    const s = meta.report.scaled;
    ui.fusion.textContent = `≈ ${sci(f.total)}`;
    ui.fusion.className = f.total >= 0.01 ? 'yes' : 'neutral';
    ui.fusionDetail.textContent = `Expected fusion reactions in this crystal per pulse ` +
      `(inside the crystal + ions stopping in surrounding LiD). Scaled to a ${s.target_radius_nm * 2} nm ` +
      `crystal: ≈ ${sci(s.beam_target_D_Li + s.beam_target_D_D)}.`;
    return;
  }
  if (!meta) {
    ui.fusion.textContent = '–';
    ui.fusion.className = 'unknown';
    ui.fusionDetail.textContent = '';
  } else if (!fusion) {
    ui.fusion.textContent = 'unknown';
    ui.fusion.className = 'unknown';
    ui.fusionDetail.textContent = 'The lattice layout in this file is not recognised.';
  } else if (fusion.count > 0) {
    const first = fusion.events[0];
    ui.fusion.textContent = `YES (${fusion.count.toLocaleString()})`;
    ui.fusion.className = 'yes';
    ui.fusionDetail.textContent = `${fusion.count.toLocaleString()} close approaches below ${meta.threshold} m. ` +
      `First: frame ${first.frame + 1}, particles ${first.p1} and ${first.p2}, ${first.distance} m apart.`;
    ui.gotoFusion.hidden = false;
  } else {
    ui.fusion.textContent = 'NO';
    ui.fusion.className = 'no';
    ui.fusionDetail.textContent = `No neighbouring pair came within ${meta.threshold} m.`;
  }
}

ui.gotoFusion.addEventListener('click', () => {
  const first = meta?.fusion?.events[0];
  if (first) setFrame(first.frame);
});

// --- improved-model report --------------------------------------------------------

function renderReport() {
  const r = meta?.report;
  ui.report.hidden = !r;
  if (!r) return;
  const L = r.laser, V = r.validity, D = r.dynamics, F = r.fusion, E = r.energy, S = r.scaled;
  const yesNo = (ok) => `<span class="${ok ? 'good' : 'bad'}">${ok ? 'yes' : 'no'}</span>`;
  const table = (title, rows) => `<div><h4>${title}</h4><table>${rows.map(([k, v]) =>
    `<tr><td>${k}</td><td>${v}</td></tr>`).join('')}</table></div>`;
  ui.reportBody.innerHTML = `<div class="report-grid">${[
    table('Laser and validity', [
      ['Peak field', `${sci(L.peak_field_V_per_m)} V/m`],
      ['Period / cycles in FWHM', `${L.period_fs.toFixed(1)} fs / ${L.cycles_in_fwhm.toFixed(2)}`],
      ['Electron a₀ (≳1: relativistic)', sci(L.electron_a0)],
      ['Field needed to strip this crystal', `${sci(V.field_needed_to_strip_V_per_m)} V/m ${yesNo(V.full_stripping_plausible)}`],
      [`… a ${S.target_radius_nm * 2} nm crystal`, `${sci(V.intensity_needed_to_strip_target_W_per_cm2)} W/cm² ${yesNo(V.full_stripping_plausible_at_target_size)}`],
      ['Laser skin depth in stripped LiD', `${V.skin_depth_nm.toFixed(0)} nm`],
    ]),
    table(`Dynamics (${r.particles.toLocaleString()} ions, ${r.simulated_fs.toFixed(0)} fs)`, [
      ['Max / mean D energy', `${sci(D.max_D_keV)} / ${sci(D.mean_D_keV)} keV`],
      ['Max / mean Li energy', `${sci(D.max_Li_keV)} / ${sci(D.mean_Li_keV)} keV`],
      ['Coulomb-explosion theory, max D', `${sci(D.coulomb_explosion_theory_max_D_keV)} keV`],
      ['Net laser drift: D / Li', `${sci(D.D_drift_keV)} / ${sci(D.Li_drift_keV)} keV`],
      ['Closest D–Li approach', `${sci(D.closest_D_Li_approach_A)} Å at ${D.closest_D_Li_approach_time_fs.toFixed(1)} fs`],
      ['Highest D–Li pair energy', `${sci(D.max_D_Li_pair_energy_keV)} keV (2.2 MeV barrier)`],
    ]),
    table('Expected fusion reactions, this crystal', [
      ['Inside the crystal: D+Li / D+D', `${sci(F.in_crystal_D_Li)} / ${sci(F.in_crystal_D_D)}`],
      ['Ions stopping in cold LiD: D+Li / D+D / Li+D', `${sci(F.beam_target_D_Li)} / ${sci(F.beam_target_D_D)} / ${sci(F.beam_target_Li_D)}`],
      ['… in a tritide (LiT) catcher instead', sci(F.catchers.LiT.total)],
      ['Fusion energy', `${sci(F.fusion_energy_MeV)} MeV`],
      ['Laser energy through the crystal', `${sci(E.laser_energy_through_crystal_MeV)} MeV`],
      ['Energy gain', E.gain === null ? '–' : sci(E.gain)],
    ]),
    table(`Extrapolated to a ${S.target_radius_nm * 2} nm crystal`, [
      ['Deuterons', S.deuterons.toLocaleString()],
      ['Max D energy', `${sci(S.max_D_keV)} keV`],
      ['Reactions: D+Li / D+D / Li+D', `${sci(S.beam_target_D_Li)} / ${sci(S.beam_target_D_D)} / ${sci(S.beam_target_Li_D)}`],
      ['Energy gain', S.gain === null ? '–' : sci(S.gain)],
    ]),
  ].join('')}</div>
  <h4>By crystal size (one crystal, one pulse; gain at the lowest intensity that strips it)</h4>
  <table><tr><td>Radius</td><td>Intensity to strip</td><td>Max D</td><td>LiD catcher (gain)</td><td>LiT catcher (gain)</td></tr>
  ${(r.size_table || []).map((row) => `<tr><td>${row.radius_nm} nm</td><td>${sci(row.intensity_to_strip_W_per_cm2)} W/cm²</td>
    <td>${sci(row.max_D_keV)} keV</td><td>${sci(row.LiD.total)} (${sci(row.LiD.gain_at_strip_threshold)})</td>
    <td>${sci(row.LiT.total)} (${sci(row.LiT.gain_at_strip_threshold)})</td></tr>`).join('')}</table>`;
}

// --- frame player --------------------------------------------------------------------

function setFrame(index) {
  if (!meta) return;
  current = Math.max(0, Math.min(meta.frames - 1, index));
  viewer.setFrame(current);
  ui.frameSlider.value = String(current + 1);
  ui.frameLabel.textContent = `Frame ${(current + 1).toLocaleString()} of ${meta.frames.toLocaleString()}`;
  ui.simTime.textContent = `Simulation time: ${meta.times[current]} s`;
  ui.compTime.textContent = `Computing time: ${meta.executionMs[current].toLocaleString()} ms`;
  updateNav();
}

function updateNav() {
  const atStart = !meta || current === 0;
  const atEnd = !meta || current >= meta.frames - 1;
  ui.restart.disabled = ui.prev.disabled = atStart || anim.running;
  ui.next.disabled = ui.last.disabled = atEnd || anim.running;
}

ui.frameSlider.addEventListener('input', () => setFrame(Number(ui.frameSlider.value) - 1));
ui.restart.addEventListener('click', () => setFrame(0));
ui.prev.addEventListener('click', () => setFrame(current - 1));
ui.next.addEventListener('click', () => setFrame(current + 1));
ui.last.addEventListener('click', () => setFrame(Infinity));

function startAnimation() {
  if (!meta || meta.frames < 2) return;
  if (current >= meta.frames - 1) setFrame(0);
  anim.running = true;
  anim.last = performance.now();
  ui.animate.textContent = 'Stop';
  updateNav();
  requestAnimationFrame(tick);
}

function stopAnimation() {
  anim.running = false;
  ui.animate.textContent = 'Animate';
  updateNav();
}

function tick(now) {
  if (!anim.running) return;
  const interval = 1000 / Number(ui.fps.value);
  if (now - anim.last >= interval) {
    anim.last = now;
    if (current >= meta.frames - 1) {
      stopAnimation();
      return;
    }
    setFrame(current + 1);
  }
  requestAnimationFrame(tick);
}

ui.animate.addEventListener('click', () => (anim.running ? stopAnimation() : startAnimation()));
ui.fps.addEventListener('input', () => { ui.fpsOut.value = ui.fps.value; });

for (const radio of document.querySelectorAll('input[name="view"]')) {
  radio.addEventListener('change', () => {
    const mode = document.querySelector('input[name="view"]:checked').value;
    viewer.setMode(mode);
    ui.zslice.disabled = mode !== '2d' || !meta;
  });
}
ui.autofit.addEventListener('change', () => viewer.setAutoFit(ui.autofit.checked));
ui.zslice.addEventListener('input', () => {
  ui.zOut.value = ui.zslice.value;
  viewer.setLayer(Number(ui.zslice.value));
});

document.addEventListener('keydown', (event) => {
  const target = event.target;
  if (!meta || event.altKey || event.ctrlKey || event.metaKey) return;
  if (target.closest('input, textarea, select, button, dialog, a')) return;
  if (event.key === 'ArrowRight') setFrame(current + 1);
  else if (event.key === 'ArrowLeft') setFrame(current - 1);
  else if (event.key === 'Home') setFrame(0);
  else if (event.key === 'End') setFrame(Infinity);
  else if (event.key === ' ') anim.running ? stopAnimation() : startAnimation();
  else return;
  event.preventDefault();
});

// --- open / save ---------------------------------------------------------------------

const unsaved = () => Boolean(summary?.canSave && !summary.saved);

function updateSave() {
  const enabled = Boolean(summary?.canSave);
  ui.save.setAttribute('aria-disabled', String(!enabled));
  const what = meta?.report ? 'the report' : 'the results';
  ui.save.title = enabled ? (summary.saved ? `Download ${what} again` : `Download ${what} as JSON`) :
    'Run a simulation first';
  document.title = `${unsaved() ? '• ' : ''}${summary ? `${summary.name} – ` : ''}NanoFusion Simulation`;
}

ui.save.addEventListener('click', (event) => {
  if (ui.save.getAttribute('aria-disabled') === 'true') {
    event.preventDefault();
    return;
  }
  setTimeout(poll, 1500); // the server marks the result saved once it has sent it
});

ui.open.addEventListener('click', () => {
  if (unsaved() && !confirm('The current results are not saved. Discard them?')) return;
  ui.openInput.click();
});

ui.openInput.addEventListener('change', async () => {
  const file = ui.openInput.files[0];
  ui.openInput.value = '';
  if (!file) return;
  setStatus(`Uploading ${file.name}…`);
  try {
    await post('/api/open', file, {
      'Content-Type': 'application/octet-stream',
      'X-Filename': encodeURIComponent(file.name),
    });
    poll();
  } catch (error) {
    setStatus(error.message, 'error');
  }
});

window.addEventListener('beforeunload', (event) => {
  if (unsaved()) event.preventDefault();
});

// --- about --------------------------------------------------------------------------

$('#about-btn').addEventListener('click', async () => {
  if (!ui.aboutBody.childElementCount) {
    const html = await (await api('/static/about.html')).text();
    ui.aboutBody.innerHTML = html; // our own static file
  }
  ui.about.showModal();
});

// --- start-up -----------------------------------------------------------------------

(async () => {
  setForm(await getJSON('/api/defaults'));
  await poll();
  if (!meta && !ui.status.textContent) setStatus('Set the parameters and press Run.');
})();
