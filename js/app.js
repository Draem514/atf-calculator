/* ============================================================
   ATF Cladding Calculator — 核心计算引擎
   ============================================================ */

function $(id) { return document.getElementById(id); }
function val(id) { return parseFloat($(id).value) || 0; }

// 数值积分（Simpson）
function simpson(f, a_, b_, n) {
  n = n || 200;
  if (n % 2) n++;
  const h = (b_ - a_) / n;
  let s = f(a_) + f(b_);
  for (let i = 1; i < n; i++) s += (i % 2 ? 2 : 4) * f(a_ + i * h);
  return s * h / 3;
}

// 主计算函数
function calcSingle(p) {
  // 单位转换：mm→m, GPa→Pa, MPa→Pa, 厚度 μm→m
  const a = p.a * 1e-3, b = p.b * 1e-3, c = p.c * 1e-3;
  const pi = p.pi * 1e6, po = p.po * 1e6;
  const E1 = p.E1 * 1e9, nu1 = p.nu1, alpha1 = p.alpha1, k1 = p.k1, Sy1 = p.Sy1 * 1e6;
  const E2 = p.E2 * 1e9, nu2 = p.nu2, alpha2 = p.alpha2, k2 = p.k2, Sy2 = p.Sy2 * 1e6;
  const tau_int = p.tau_int * 1e6;
  const Ti = p.Ti, To = p.To;
  const eps_irr_r1 = p.eps_irr_r1, eps_irr_r2 = p.eps_irr_r2;
  const eps_irr_z1 = p.eps_irr_z1, eps_irr_z2 = p.eps_irr_z2;

  // === 温度场（双层稳态导热）===
  const c2 = (Ti - To) / ((k2 / k1) * Math.log(b / a) + Math.log(c / b));
  const c1 = (k2 / k1) * c2;
  const D1 = Ti + c1 * Math.log(a), D2 = To + c2 * Math.log(c);
  function T1(r) { return D1 - c1 * Math.log(r); }
  function T2(r) { return D2 - c2 * Math.log(r); }
  const Ti_interface = T1(b);

  // === 热-辐照特解积分 ===
  function I1(r) {
    // ∫_a^r [alpha1*(T1(ρ)-To) + eps_irr_r1] * ρ dρ
    const f = rho => (alpha1 * (T1(rho) - To) + eps_irr_r1) * rho;
    return simpson(f, a, r, 100);
  }
  function I2(r) {
    // ∫_b^r [alpha2*(T2(ρ)-To) + eps_irr_r2] * ρ dρ
    const f = rho => (alpha2 * (T2(rho) - To) + eps_irr_r2) * rho;
    return simpson(f, b, r, 100);
  }
  function epsInit1(r) { return alpha1 * (T1(r) - To) + eps_irr_r1; }
  function epsInit2(r) { return alpha2 * (T2(r) - To) + eps_irr_r2; }

  // === 求解 4x4 线性方程组 [A1,B1,A2,B2] ===
  const M = Array.from({length:4}, () => new Array(4).fill(0));
  const rhs = new Array(4).fill(0);
  // Eq1: r=a, sigma_r = -pi
  M[0][0] = 1; M[0][1] = 1 / (a * a);
  const sr_sp_a = -E1 / (1 - nu1) * I1(a) / (a * a);
  rhs[0] = -pi - sr_sp_a;
  // Eq2: r=c, sigma_r = -po
  M[1][2] = 1; M[1][3] = 1 / (c * c);
  const sr_sp_c = -E2 / (1 - nu2) * I2(c) / (c * c);
  rhs[1] = -po - sr_sp_c;
  // Eq3: r=b, sigma_r 连续
  M[2][0] = 1; M[2][1] = 1 / (b * b); M[2][2] = -1; M[2][3] = -1 / (b * b);
  const sr_b1 = -E1 / (1 - nu1) * I1(b) / (b * b);
  const sr_b2 = -E2 / (1 - nu2) * I2(b) / (b * b);
  rhs[2] = sr_b2 - sr_b1;
  // Eq4: r=b, u 连续
  M[3][0] = (1 + nu1) / E1 * (1 - 2 * nu1) * b;
  M[3][1] = -(1 + nu1) / E1 * 1 / b;
  M[3][2] = -(1 + nu2) / E2 * (1 - 2 * nu2) * b;
  M[3][3] = (1 + nu2) / E2 * 1 / b;
  const u_sp_b1 = (1 + nu1) / (1 - nu1) * I1(b) / b;
  const u_sp_b2 = (1 + nu2) / (1 - nu2) * I2(b) / b;
  rhs[3] = u_sp_b2 - u_sp_b1;

  // 高斯消元求解
  const A = M.map((row, i) => [...row, rhs[i]]);
  for (let col = 0; col < 4; col++) {
    let maxR = col;
    for (let r = col + 1; r < 4; r++) if (Math.abs(A[r][col]) > Math.abs(A[maxR][col])) maxR = r;
    [A[col], A[maxR]] = [A[maxR], A[col]];
    const piv = A[col][col];
    for (let j = col; j < 5; j++) A[col][j] /= piv;
    for (let r = 0; r < 4; r++) {
      if (r === col) continue;
      const f = A[r][col];
      for (let j = col; j < 5; j++) A[r][j] -= f * A[col][j];
    }
  }
  const [A1, B1, A2, B2] = [A[0][4], A[1][4], A[2][4], A[3][4]];

  // === 应力场 ===
  function sigmaR_Zr(r) {
    return A1 + B1 / (r * r) - E1 / (1 - nu1) * I1(r) / (r * r);
  }
  function sigmaTheta_Zr(r) {
    return A1 - B1 / (r * r) + E1 / (1 - nu1) * I1(r) / (r * r) - E1 * epsInit1(r) / (1 - nu1);
  }
  function sigmaR_Cr(r) {
    return A2 + B2 / (r * r) - E2 / (1 - nu2) * I2(r) / (r * r);
  }
  function sigmaTheta_Cr(r) {
    return A2 - B2 / (r * r) + E2 / (1 - nu2) * I2(r) / (r * r) - E2 * epsInit2(r) / (1 - nu2);
  }
  function sigmaZ_Zr(r) {
    const sr = sigmaR_Zr(r), st = sigmaTheta_Zr(r);
    return nu1 * (sr + st) - E1 * alpha1 * (T1(r) - To) / (1 - nu1) - E1 * eps_irr_z1 / (1 - nu1);
  }

  // 关键位置应力 (Pa)
  const st_inner = sigmaTheta_Zr(a);
  const st_b1 = sigmaTheta_Zr(b);
  const st_b2 = sigmaTheta_Cr(b);
  const st_outer = sigmaTheta_Cr(c);
  const sr_inner = sigmaR_Zr(a);
  const sz_inner = sigmaZ_Zr(a);
  // von Mises at inner wall (Zr side)
  const von_inner = Math.sqrt(0.5 * ((sr_inner - st_inner)**2 + (st_inner - sz_inner)**2 + (sz_inner - sr_inner)**2));

  // === 失效评估 ===
  const t_coat = c - b; // m
  // 涂层屈曲
  const sigma_cr = E2 / Math.sqrt(3 * (1 - nu2 ** 2)) * (t_coat / b);
  const safety_buckling = Math.abs(sigma_cr) / Math.abs(st_outer);
  // 界面剪切
  const Delta_st = Math.abs(st_b2 - st_b1);
  const tau_eq = Delta_st * t_coat / (2 * b);
  const safety_interface = tau_int / tau_eq;
  // 基体屈服
  const safety_yield = Sy1 / von_inner;

  // 裂纹断裂力学
  let K_I = null, K_ratio = null;
  if (p.a_crack > 0) {
    let r_tip, st_tip, KIC_mat;
    if (p.crackPos === 'inner') { r_tip = a + p.a_crack; KIC_mat = p.KIC_Zr; }
    else if (p.crackPos === 'interface') { r_tip = b; KIC_mat = p.KIC_Zr; }
    else { r_tip = c - p.a_crack; KIC_mat = p.KIC_Cr; }
    if (r_tip <= b) st_tip = sigmaTheta_Zr(r_tip);
    else st_tip = sigmaTheta_Cr(r_tip);
    if (st_tip > 0) {
      K_I = st_tip * Math.sqrt(Math.PI * p.a_crack) * 1.12; // Pa·sqrt(m)
      K_ratio = K_I / KIC_mat;
    }
  }

  return {
    T_interface: Ti_interface,
    sigmaTheta_inner: st_inner,
    sigmaTheta_interface_Zr: st_b1,
    sigmaTheta_interface_Cr: st_b2,
    sigmaTheta_outer: st_outer,
    vonMises_inner: von_inner,
    safety_buckling,
    safety_interface,
    safety_yield,
    K_I, K_ratio,
    t_coat_um: (p.c - p.b) * 1000 // mm直接转μm，不经过m
  };
}

// ===== 材料预设库 =====
const PRESET_MATERIALS = {
  zr: {
    zt4: { E: 75.7, nu: 0.30, alpha: 6.5e-6, k: 16.4, Sy: 200, label: 'Zircaloy-4' },
    zrnb: { E: 95.0, nu: 0.34, alpha: 7.0e-6, k: 22.0, Sy: 350, label: 'Zr-Nb 合金' }
  },
  cr: {
    purecr: { E: 280, nu: 0.22, alpha: 9.1e-6, k: 76, Sy: 1200, tau: 80, label: '纯 Cr' }
  }
};

// localStorage 键名
const LS_KEY_ZR = 'atf_custom_mat_zr';
const LS_KEY_CR = 'atf_custom_mat_cr';

// 从 localStorage 加载自定义材料
function loadCustomMaterials(side) {
  const key = side === 'zr' ? LS_KEY_ZR : LS_KEY_CR;
  try { return JSON.parse(localStorage.getItem(key)) || []; }
  catch { return []; }
}

// 保存自定义材料到 localStorage
function saveCustomMaterials(side, list) {
  const key = side === 'zr' ? LS_KEY_ZR : LS_KEY_CR;
  localStorage.setItem(key, JSON.stringify(list));
}

// 渲染自定义材料列表（可点击加载、可删除）
function renderCustomList(side) {
  const listEl = $(side === 'zr' ? 'custom_list_zr' : 'custom_list_cr');
  const items = loadCustomMaterials(side);
  if (items.length === 0) {
    listEl.classList.add('hidden');
    listEl.innerHTML = '';
    return;
  }
  listEl.classList.remove('hidden');
  listEl.innerHTML = '<div style="font-size:0.78rem;color:#5a7a9a;margin-bottom:0.3rem;">已保存的预设：</div>'
    + items.map((m, i) => `
      <div class="custom-item" data-idx="${i}">
        <span class="custom-item-name">${escHtml(m.name)} — E=${m.E} GPa, Sy=${m.Sy} MPa</span>
        <button class="custom-item-del" data-idx="${i}" title="删除">×</button>
      </div>`).join('');

  // 点击加载
  listEl.querySelectorAll('.custom-item').forEach(el => {
    el.addEventListener('click', e => {
      if (e.target.classList.contains('custom-item-del')) return;
      loadSavedMaterial(side, parseInt(el.dataset.idx));
    });
  });
  // 点击删除
  listEl.querySelectorAll('.custom-item-del').forEach(el => {
    el.addEventListener('click', e => {
      e.stopPropagation();
      deleteCustomMaterial(side, parseInt(el.dataset.idx));
    });
  });
}

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

// 加载已保存的自定义材料
function loadSavedMaterial(side, idx) {
  const items = loadCustomMaterials(side);
  const m = items[idx];
  if (!m) return;
  // 填充输入框
  const prefix = side === 'zr' ? 'in_' : 'in_';
  const numIdx = side === 'zr' ? '1' : '2';
  $(prefix + 'E' + numIdx).value = m.E;
  $(prefix + 'nu' + numIdx).value = m.nu;
  $(prefix + 'alpha' + numIdx).value = m.alpha;
  $(prefix + 'k' + numIdx).value = m.k;
  $(prefix + 'Sy' + numIdx).value = m.Sy;
  if (side === 'cr') $(prefix + 'tau').value = m.tau || 80;
  // 更新标签
  const nameEl = $(side === 'zr' ? 'matname_zr' : 'matname_cr');
  const lblIrr = $(side === 'zr' ? 'lbl_irr1a' : 'lbl_irr2a');
  nameEl.textContent = m.name + (side === 'zr' ? ' 基体' : ' 涂层');
  if (lblIrr) {
    const span = lblIrr.querySelector('.mat-name');
    if (span) span.textContent = m.name;
  }
}

// 删除已保存的自定义材料
function deleteCustomMaterial(side, idx) {
  const items = loadCustomMaterials(side);
  items.splice(idx, 1);
  saveCustomMaterials(side, items);
  renderCustomList(side);
}

// 填充预设材料参数到输入框并更新标签
function fillMaterial(side) {
  const selId = side === 'zr' ? 'sel_mat_zr' : 'sel_mat_cr';
  const btnId = side === 'zr' ? 'btn_fill_zr' : 'btn_fill_cr';
  const prefix = side === 'zr' ? 'in_' : 'in_';
  const numIdx = side === 'zr' ? '1' : '2';
  const sel = $(selId);
  const matKey = sel.value;
  const library = PRESET_MATERIALS[side];
  if (!library || !library[matKey]) return;
  const m = library[matKey];
  $(prefix + 'E' + numIdx).value = m.E;
  $(prefix + 'nu' + numIdx).value = m.nu;
  $(prefix + 'alpha' + numIdx).value = m.alpha;
  $(prefix + 'k' + numIdx).value = m.k;
  $(prefix + 'Sy' + numIdx).value = m.Sy;
  if (side === 'cr') $(prefix + 'tau').value = m.tau || 80;
  // 更新标签
  const nameEl = $(side === 'zr' ? 'matname_zr' : 'matname_cr');
  const lblIrr = $(side === 'zr' ? 'lbl_irr1a' : 'lbl_irr2a');
  nameEl.textContent = m.label + (side === 'zr' ? ' 基体' : ' 涂层');
  if (lblIrr) {
    const span = lblIrr.querySelector('.mat-name');
    if (span) span.textContent = m.label;
  }
  // 按钮反馈
  const btn = $(btnId);
  btn.textContent = '✓ 已加载';
  setTimeout(() => btn.textContent = '加载', 1000);
}

// 保存当前输入框参数为自定义材料
function saveCustomMaterial(side) {
  const prefix = side === 'zr' ? 'in_' : 'in_';
  const numIdx = side === 'zr' ? '1' : '2';
  const nameInput = $(side === 'zr' ? 'in_custom_name_zr' : 'in_custom_name_cr');
  const name = nameInput.value.trim();
  if (!name) { alert('请输入材料名称'); nameInput.focus(); return; }
  const item = {
    name,
    E: val(prefix + 'E' + numIdx),
    nu: val(prefix + 'nu' + numIdx),
    alpha: val(prefix + 'alpha' + numIdx),
    k: val(prefix + 'k' + numIdx),
    Sy: val(prefix + 'Sy' + numIdx)
  };
  if (side === 'cr') item.tau = val(prefix + 'tau');
  // 检查是否已存在同名材料
  const items = loadCustomMaterials(side);
  const existIdx = items.findIndex(it => it.name === name);
  if (existIdx >= 0) {
    if (!confirm(`"${name}" 已存在，是否覆盖？`)) return;
    items[existIdx] = item;
  } else {
    items.push(item);
  }
  saveCustomMaterials(side, items);
  renderCustomList(side);
  // 切换下拉框到该自定义项
  const sel = $(side === 'zr' ? 'sel_mat_zr' : 'sel_mat_cr');
  sel.value = 'custom';
  onCustomChange(side);
  // 按钮反馈
  const btn = $(side === 'zr' ? 'btn_save_zr' : 'btn_save_cr');
  btn.textContent = '✓ 已保存';
  setTimeout(() => btn.textContent = '保存', 1000);
}

// 当切换到"自定义"时的处理
function onCustomChange(side) {
  const sel = $(side === 'zr' ? 'sel_mat_zr' : 'sel_mat_cr');
  const fillBtn = $(side === 'zr' ? 'btn_fill_zr' : 'btn_fill_cr');
  const saveBtn = $(side === 'zr' ? 'btn_save_zr' : 'btn_save_cr');
  const nameRow = $(side === 'zr' ? 'custom_name_row_zr' : 'custom_name_row_cr');
  const nameEl = $(side === 'zr' ? 'matname_zr' : 'matname_cr');
  const lblIrr = $(side === 'zr' ? 'lbl_irr1a' : 'lbl_irr2a');
  const matLabel = side === 'zr' ? '基体' : '涂层';

  if (sel.value === 'custom') {
    fillBtn.classList.add('hidden');
    saveBtn.classList.remove('hidden');
    nameRow.classList.remove('hidden');
    nameEl.textContent = '自定义 ' + matLabel + '（请先填写参数再保存）';
    if (lblIrr) {
      const span = lblIrr.querySelector('.mat-name');
      if (span) span.textContent = '自定义';
    }
    renderCustomList(side);
  } else {
    fillBtn.classList.remove('hidden');
    saveBtn.classList.add('hidden');
    nameRow.classList.add('hidden');
  }
}

// ===== 读取页面参数 =====
function readParams() {
  return {
    a: val('in_a'), b: val('in_b'), c: val('in_c'),
    pi: val('in_pi'), po: val('in_po'),
    E1: val('in_E1'), nu1: val('in_nu1'), alpha1: val('in_alpha1'),
    k1: val('in_k1'), Sy1: val('in_Sy1'),
    E2: val('in_E2'), nu2: val('in_nu2'), alpha2: val('in_alpha2'),
    k2: val('in_k2'), Sy2: val('in_Sy2'), tau_int: val('in_tau'),
    Ti: val('in_ti'), To: val('in_to'),
    eps_irr_r1: val('in_epsirr_r1'), eps_irr_r2: val('in_epsirr_r2'),
    eps_irr_z1: val('in_epsirr_z1'), eps_irr_z2: val('in_epsirr_z2'),
    crackPos: $('in_crack_pos').value, a_crack: val('in_a_crack') * 1e-3, // um->mm for consistency
    KIC_Cr: val('in_KIC_Cr') * 1e6 * Math.sqrt(1e-3),
    KIC_Zr: val('in_KIC_Zr') * 1e6 * Math.sqrt(1e-3)
  };
}

// ===== 渲染单工况结果 =====
function fmt(v, d=4) { return Number(v.toFixed(d)).toString(); }
function cls(v) { return v >= 1 ? 'pass' : v < 0.8 ? 'fail' : 'warn'; }

function renderResult(res) {
  const rows = [
    ['Zr/Cr 界面温度', fmt(res.T_interface, 1) + ' K'],
    ['内壁环向应力 σ_θ (r=a)', fmt(res.sigmaTheta_inner / 1e6, 1) + ' MPa'],
    ['界面环向应力 σ_θ (r=b, Zr侧)', fmt(res.sigmaTheta_interface_Zr / 1e6, 1) + ' MPa'],
    ['界面环向应力 σ_θ (r=b, Cr侧)', fmt(res.sigmaTheta_interface_Cr / 1e6, 1) + ' MPa'],
    ['外壁环向应力 σ_θ (r=c)', fmt(res.sigmaTheta_outer / 1e6, 1) + ' MPa'],
    ['内壁 von Mises 应力', fmt(res.vonMises_inner / 1e6, 1) + ' MPa'],
  ];
  const sfRows = [
    ['涂层屈曲安全系数', res.safety_buckling, fmt(res.safety_buckling, 2)],
    ['界面剪切安全系数', res.safety_interface, fmt(res.safety_interface, 2)],
    ['基体屈服安全系数', res.safety_yield, fmt(res.safety_yield, 2)],
  ];
  let html = '<table class="result-table"><tr><th colspan="2">温度场与应力</th></tr>';
  rows.forEach(r => { html += `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`; });
  html += '<tr><th colspan="2">失效评估</th></tr>';
  sfRows.forEach(r => { html += `<tr><td>${r[0]}</td><td class="${cls(r[1])}">${r[2]}</td></tr>`; });
  if (res.K_I != null) {
    html += '<tr><th colspan="2">裂纹断裂力学</th></tr>';
    html += `<tr><td>K_I (Mode I)</td><td>${fmt(res.K_I / 1e6, 3)} MPa√m</td></tr>`;
    html += `<tr><td>K_I / K_IC</td><td class="${cls(1 / res.K_ratio)}">${fmt(res.K_ratio, 3)}</td></tr>`;
  }
  html += `<tr><td>涂层厚度</td><td>${fmt(res.t_coat_um, 1)} μm</td></tr></table>`;
  html += '<div class="sf-legend">'
    + '<span class="lg pass">■ ≥ 1.0 安全</span>'
    + '<span class="lg warn">■ 0.8–1.0 临界</span>'
    + '<span class="lg fail">■ &lt; 0.8 危险</span>'
    + '</div>';
  $('result-content').innerHTML = html;
  $('result-placeholder').classList.add('hidden');
  $('result-content').classList.remove('hidden');
}

// ===== 扫描 =====
function updateScanNote() {
  const param = $('scan_param').value;
  const start = val('scan_start'), end = val('scan_end'), steps = val('scan_steps');
  const unit = param === 't_coat' ? 'μm' : param === 'p_i' ? 'MPa' : 'K';
  const label = param === 't_coat' ? '涂层厚度' : param === 'p_i' ? '内压 p_i' : '内壁温度 T_inner';
  $('scan_note').textContent = `${label}: ${start} ~ ${end} ${unit}，共 ${steps} 步`;
}

function drawScanChart(labels, buckling, interf, yield_) {
  const canvas = $('scan-chart');
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const pad = { t: 40, r: 30, b: 50, l: 60 };
  const cW = W - pad.l - pad.r, cH = H - pad.t - pad.b;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#0d1b2a'; ctx.fillRect(0, 0, W, H);

  const allVals = [...buckling, ...interf, ...yield_].filter(v => isFinite(v));
  if (allVals.length === 0) return;
  const yMin = Math.min(...allVals) * 0.8, yMax = Math.max(...allVals) * 1.2;
  const yRange = yMax - yMin || 1;
  const n = labels.length;
  const xStep = n > 1 ? cW / (n - 1) : cW / 2;

  function toX(i) { return pad.l + i * xStep; }
  function toY(v) { return pad.t + cH - (v - yMin) / yRange * cH; }

  // 网格
  ctx.strokeStyle = '#1e3a5f'; ctx.lineWidth = 1;
  for (let i = 0; i <= 5; i++) {
    const y = pad.t + cH * i / 5;
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + cW, y); ctx.stroke();
    ctx.fillStyle = '#5a7a9a'; ctx.font = '11px Segoe UI'; ctx.textAlign = 'right';
    ctx.fillText((yMax - yRange * i / 5).toFixed(1), pad.l - 8, y + 4);
  }
  // x 标签
  ctx.textAlign = 'center';
  labels.forEach((l, i) => {
    if (i % Math.ceil(n / 10) === 0) {
      ctx.fillStyle = '#8ab4f8'; ctx.fillText(l, toX(i), H - pad.b + 20);
    }
  });
  // 失效临界线 y=1
  if (yMin < 1 && yMax > 1) {
    ctx.strokeStyle = '#f87171'; ctx.setLineDash([6, 4]); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(pad.l, toY(1)); ctx.lineTo(pad.l + cW, toY(1)); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#f87171'; ctx.font = 'bold 11px Segoe UI'; ctx.textAlign = 'left';
    ctx.fillText('失效临界', pad.l + cW - 60, toY(1) - 6);
  }
  // 三条曲线
  const lines = [{ data: buckling, color: '#f87171', label: '屈曲' },
                 { data: interf, color: '#60a5fa', label: '界面' },
                 { data: yield_, color: '#4ade80', label: '屈服' }];
  lines.forEach(({ data, color, label }) => {
    ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.beginPath();
    let started = false;
    data.forEach((v, i) => {
      if (!isFinite(v)) return;
      if (!started) { ctx.moveTo(toX(i), toY(v)); started = true; }
      else ctx.lineTo(toX(i), toY(v));
    });
    ctx.stroke();
    data.forEach((v, i) => {
      if (!isFinite(v)) return;
      ctx.fillStyle = color; ctx.beginPath(); ctx.arc(toX(i), toY(v), 3, 0, Math.PI * 2); ctx.fill();
    });
  });
  // 图例
  let lx = pad.l + 10;
  lines.forEach(({ color, label }) => {
    ctx.fillStyle = color; ctx.fillRect(lx, 12, 14, 10);
    ctx.fillStyle = '#e0e6f0'; ctx.font = '12px Segoe UI'; ctx.textAlign = 'left';
    ctx.fillText(label, lx + 18, 22); lx += 70;
  });
  // 轴标题
  ctx.fillStyle = '#c4a35a'; ctx.font = 'bold 12px Segoe UI'; ctx.textAlign = 'center';
  ctx.fillText('扫描参数值', pad.l + cW / 2, H - 5);
  ctx.save(); ctx.translate(14, pad.t + cH / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText('安全系数', 0, 0); ctx.restore();
}

function renderScanResult(paramsArr, resultsArr) {
  const labels = paramsArr.map(v => fmt(v, 1));
  const buckling = resultsArr.map(r => r.safety_buckling);
  const interf = resultsArr.map(r => r.safety_interface);
  const yield_ = resultsArr.map(r => r.safety_yield);

  const el = $('scan-result');
  let html = '<table class="result-table"><tr><th>扫描参数</th><th>屈曲SF</th><th>界面SF</th><th>屈服SF</th><th>界面温度(K)</th></tr>';
  resultsArr.forEach((r, i) => {
    html += `<tr><td>${labels[i]}</td>
      <td class="${cls(r.safety_buckling)}">${fmt(r.safety_buckling, 2)}</td>
      <td class="${cls(r.safety_interface)}">${fmt(r.safety_interface, 2)}</td>
      <td class="${cls(r.safety_yield)}">${fmt(r.safety_yield, 2)}</td>
      <td>${fmt(r.T_interface, 1)}</td></tr>`;
  });
  html += '</table>';
  html += '<div class="sf-legend">'
    + '<span class="lg pass">■ ≥ 1.0 安全</span>'
    + '<span class="lg warn">■ 0.8–1.0 临界</span>'
    + '<span class="lg fail">■ &lt; 0.8 危险</span>'
    + '</div>';
  el.innerHTML = html;
  el.classList.remove('hidden');
  $('scan-chart-area').classList.remove('hidden');
  drawScanChart(labels, buckling, interf, yield_);
}

// ===== 初始化 =====
document.addEventListener('DOMContentLoaded', () => {
  // 导航切换
  document.querySelectorAll('nav a').forEach(a => {
    a.addEventListener('click', e => {
      e.preventDefault();
      document.querySelectorAll('nav a').forEach(x => x.classList.remove('active'));
      a.classList.add('active');
      document.querySelectorAll('main').forEach(p => p.classList.remove('active'));
      document.getElementById('section-' + a.dataset.section).classList.add('active');
    });
  });

  // 涂层厚度实时预览
  ['in_b', 'in_c'].forEach(id => $(id).addEventListener('input', () => {
    $('disp_tcoat').textContent = ((val('in_c') - val('in_b')) * 1000).toFixed(1);
  }));

  // 扫描参数联动
  ['scan_param', 'scan_start', 'scan_end', 'scan_steps'].forEach(id =>
    $(id).addEventListener('input', updateScanNote));

  // 材料填充按钮
  $('btn_fill_zr').addEventListener('click', () => fillMaterial('zr'));
  $('btn_fill_cr').addEventListener('click', () => fillMaterial('cr'));
  // 保存按钮
  $('btn_save_zr').addEventListener('click', () => saveCustomMaterial('zr'));
  $('btn_save_cr').addEventListener('click', () => saveCustomMaterial('cr'));
  // 下拉框切换（预设 vs 自定义）
  $('sel_mat_zr').addEventListener('change', () => onCustomChange('zr'));
  $('sel_mat_cr').addEventListener('change', () => onCustomChange('cr'));
  // 初始化时渲染已有自定义材料列表
  renderCustomList('zr');
  renderCustomList('cr');

  // 计算按钮
  $('btn_calc').addEventListener('click', () => {
    try { renderResult(calcSingle(readParams())); }
    catch (e) { alert('计算出错: ' + e.message); console.error(e); }
  });

  // 扫描按钮
  $('btn_scan').addEventListener('click', () => {
    try {
      const param = $('scan_param').value;
      const start = val('scan_start'), end = val('scan_end');
      const steps = Math.max(2, Math.round(val('scan_steps')));
      const step = (end - start) / (steps - 1);
      const base = readParams();
      const paramsArr = [], resultsArr = [];
      for (let i = 0; i < steps; i++) {
        const v = start + i * step;
        const p = { ...base };
        if (param === 't_coat') { p.c = p.b + v * 1e-3; }
        else if (param === 'p_i') { p.pi = v; }
        else if (param === 'T_inner') { p.Ti = v; }
        paramsArr.push(v);
        resultsArr.push(calcSingle(p));
      }
      renderScanResult(paramsArr, resultsArr);
    } catch (e) { alert('扫描出错: ' + e.message); console.error(e); }
  });

  updateScanNote();
});
