import streamlit as st
import streamlit.components.v1 as components
import prompts
import chat_with_gpt
from string import Template
from rdkit import Chem
from rdkit.Chem import Draw
import requests

size = (250, 250)


# Set Page Configuration
st.set_page_config(page_title="QAI Model", layout="wide", page_icon="🧪")


# ---------------------------------------------------------------------------
# Design system
# ---------------------------------------------------------------------------
st.markdown("""
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
  :root{
    --ink:#0a0f14;
    --panel:#111a22;
    --panel-2:#16212b;
    --line:rgba(148,178,196,.16);
    --line-strong:rgba(148,178,196,.30);
    --text:#e7eef4;
    --muted:#8ea3b3;
    --accent:#2ed3c4;
    --accent-deep:#12867f;
    --warn:#f0b429;
    --danger:#ef6461;
    --radius:14px;
    --radius-sm:10px;
  }

  /* ---------- Shell ---------- */
  .stApp{
    background:
      radial-gradient(1100px 600px at 12% -10%, rgba(46,211,196,.10), transparent 60%),
      radial-gradient(900px 500px at 100% 0%, rgba(59,130,246,.08), transparent 55%),
      var(--ink);
    color:var(--text);
    font-family:'IBM Plex Sans', system-ui, -apple-system, 'Segoe UI', sans-serif;
  }
  .block-container{
    padding-top:1.6rem !important;
    padding-bottom:3.5rem !important;
    max-width:1180px;
  }
  header[data-testid="stHeader"]{background:transparent;}
  #MainMenu, footer{visibility:hidden;}

  h1,h2,h3,h4{font-family:'IBM Plex Sans', sans-serif; color:var(--text); letter-spacing:-.01em;}

  /* ---------- Masthead ---------- */
  .qai-mast{
    display:flex; align-items:center; gap:16px;
    padding:20px 24px;
    border:1px solid var(--line);
    border-radius:var(--radius);
    background:linear-gradient(180deg, rgba(46,211,196,.07), rgba(255,255,255,0)) , var(--panel);
    margin-bottom:22px;
  }
  .qai-mast .mark{
    flex:0 0 auto; width:46px; height:46px; border-radius:12px;
    display:grid; place-items:center; font-size:22px;
    background:linear-gradient(145deg, var(--accent-deep), rgba(46,211,196,.25));
    border:1px solid rgba(46,211,196,.45);
  }
  .qai-mast h1{font-size:1.5rem; font-weight:600; margin:0 0 2px 0;}
  .qai-mast p{margin:0; color:var(--muted); font-size:.92rem; line-height:1.45;}
  .qai-mast .stamp{
    margin-left:auto; text-align:right; color:var(--muted);
    font-family:'IBM Plex Mono', monospace; font-size:.74rem; line-height:1.5;
    border-left:1px solid var(--line); padding-left:16px;
  }
  .qai-mast .stamp b{color:var(--accent); font-weight:500;}

  /* ---------- Panels ---------- */
  .qai-panel{
    border:1px solid var(--line);
    border-radius:var(--radius);
    background:var(--panel);
    padding:22px 24px 8px 24px;
    margin-bottom:20px;
  }
  .qai-panel-head{
    display:flex; align-items:baseline; gap:10px;
    padding-bottom:14px; margin-bottom:18px;
    border-bottom:1px solid var(--line);
  }
  .qai-panel-head .t{font-size:1rem; font-weight:600;}
  .qai-panel-head .s{font-size:.82rem; color:var(--muted);}

  /* ---------- Widget labels ---------- */
  div[data-testid="stWidgetLabel"] label p,
  div[data-testid="stWidgetLabel"] label{
    color:var(--text) !important;
    font-size:.88rem !important;
    font-weight:500 !important;
  }

  /* ---------- Inputs ---------- */
  div[data-testid="stTextInput"] input,
  div[data-testid="stNumberInput"] input,
  div[data-testid="stTextArea"] textarea,
  div[data-testid="stSelectbox"] div[data-baseweb="select"] > div{
    background-color:var(--panel-2) !important;
    color:var(--text) !important;
    border:1px solid var(--line-strong) !important;
    border-radius:var(--radius-sm) !important;
    box-shadow:none !important;
    font-size:.92rem !important;
  }
  div[data-testid="stTextInput"] input,
  div[data-testid="stNumberInput"] input{padding:11px 13px !important;}
  div[data-testid="stTextArea"] textarea{
    padding:12px 13px !important;
    font-family:'IBM Plex Mono', monospace !important;
    font-size:.86rem !important;
    line-height:1.6 !important;
  }
  input::placeholder, textarea::placeholder{color:#6b8090 !important;}

  div[data-testid="stTextInput"] input:hover,
  div[data-testid="stNumberInput"] input:hover,
  div[data-testid="stTextArea"] textarea:hover,
  div[data-testid="stSelectbox"] div[data-baseweb="select"] > div:hover{
    border-color:rgba(46,211,196,.55) !important;
  }
  div[data-testid="stTextInput"] input:focus,
  div[data-testid="stNumberInput"] input:focus,
  div[data-testid="stTextArea"] textarea:focus,
  div[data-testid="stSelectbox"] div[data-baseweb="select"] > div:focus-within{
    border-color:var(--accent) !important;
    box-shadow:0 0 0 3px rgba(46,211,196,.18) !important;
    outline:none !important;
  }
  /* dropdown menu */
  div[data-baseweb="popover"] ul{background:var(--panel-2) !important; border:1px solid var(--line-strong) !important;}
  div[data-baseweb="popover"] li{color:var(--text) !important; font-size:.9rem !important;}
  div[data-baseweb="popover"] li:hover{background:rgba(46,211,196,.14) !important;}

  /* ---------- Checkbox ---------- */
  div[data-testid="stCheckbox"] label{color:var(--text) !important; font-size:.9rem !important;}
  div[data-testid="stCheckbox"] label span[data-baseweb="checkbox"] div:first-child{
    background:var(--panel-2) !important; border-color:var(--line-strong) !important;
  }

  /* ---------- Buttons ---------- */
  .stButton > button{
    width:100%;
    border-radius:var(--radius-sm);
    border:1px solid var(--line-strong);
    background:var(--panel-2);
    color:var(--text);
    font-family:'IBM Plex Sans', sans-serif;
    font-weight:500;
    font-size:.92rem;
    padding:.62rem 1.1rem;
    transition:border-color .16s ease, background .16s ease, color .16s ease;
  }
  .stButton > button:hover{
    border-color:var(--accent);
    background:rgba(46,211,196,.10);
    color:#ffffff;
  }
  .stButton > button:focus-visible{
    outline:2px solid var(--accent);
    outline-offset:2px;
  }
  /* primary action */
  .stButton > button[kind="primary"]{
    background:linear-gradient(100deg, var(--accent-deep), var(--accent));
    border:1px solid rgba(46,211,196,.7);
    color:#04201e;
    font-weight:600;
  }
  .stButton > button[kind="primary"]:hover{
    filter:brightness(1.08);
    color:#04201e;
  }

  /* ---------- Structure viewer ---------- */
  .qai-viewer{
    border:1px dashed var(--line-strong);
    border-radius:var(--radius-sm);
    background:var(--panel-2);
    padding:18px;
    text-align:center;
    color:var(--muted);
    font-size:.85rem;
    line-height:1.55;
  }
  div[data-testid="stImage"] img{
    background:#ffffff;
    border-radius:var(--radius-sm);
    border:1px solid var(--line-strong);
    padding:8px;
  }
  div[data-testid="stImage"] div[data-testid="caption"]{
    color:var(--muted) !important;
    font-family:'IBM Plex Mono', monospace;
    font-size:.78rem !important;
  }

  /* ---------- Spec list (result page) ---------- */
  .qai-spec{display:grid; gap:1px; background:var(--line); border:1px solid var(--line); border-radius:var(--radius); overflow:hidden;}
  .qai-spec .row{display:flex; justify-content:space-between; gap:18px; padding:14px 18px; background:var(--panel);}
  .qai-spec .k{color:var(--muted); font-size:.86rem;}
  .qai-spec .v{color:var(--text); font-family:'IBM Plex Mono', monospace; font-size:.88rem; text-align:right; word-break:break-word;}

  /* ---------- Report body ---------- */
  .qai-report{
    border:1px solid var(--line);
    border-radius:var(--radius);
    background:var(--panel);
    padding:8px 28px 18px 28px;
  }
  .qai-report p, .qai-report li{color:#d5e2ec; font-size:.95rem; line-height:1.72;}
  .qai-report h1{font-size:1.32rem;} .qai-report h2{font-size:1.14rem;} .qai-report h3{font-size:1rem;}
  .qai-report h1,.qai-report h2,.qai-report h3{margin-top:1.5rem;}
  .qai-report table{width:100%; border-collapse:collapse; margin:1rem 0; font-size:.88rem;}
  .qai-report th{background:var(--panel-2); color:var(--text); text-align:left; padding:10px 12px; border:1px solid var(--line);}
  .qai-report td{padding:10px 12px; border:1px solid var(--line); color:#d5e2ec;}
  .qai-report code{background:rgba(46,211,196,.10); color:var(--accent); padding:.12em .4em; border-radius:5px; font-family:'IBM Plex Mono', monospace;}

  /* ---------- Alerts / spinner ---------- */
  div[data-testid="stAlert"]{border-radius:var(--radius-sm); border:1px solid var(--line-strong);}
  div[data-testid="stSpinner"] p{color:var(--muted) !important; font-size:.88rem !important;}
  hr{border-color:var(--line) !important;}

  /* ---------- Responsive ---------- */
  @media (max-width:900px){
    .block-container{padding-left:1rem !important; padding-right:1rem !important;}
    .qai-mast{flex-wrap:wrap; gap:12px; padding:18px;}
    .qai-mast .stamp{margin-left:0; border-left:none; padding-left:0; text-align:left; width:100%;}
    .qai-panel{padding:18px 16px 6px 16px;}
    .qai-report{padding:6px 16px 14px 16px;}
    .qai-spec .row{flex-direction:column; gap:4px;}
    .qai-spec .v{text-align:left;}
  }
  @media (prefers-reduced-motion:reduce){
    *{transition:none !important; animation:none !important;}
  }
</style>
""", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
if "page" not in st.session_state:
    st.session_state.page = "form"
if "api_response" not in st.session_state:
    st.session_state.api_response = None
if "structure_fig" not in st.session_state:
    st.session_state.structure_fig = None
if "structure_label" not in st.session_state:
    st.session_state.structure_label = ""
if "structure_error" not in st.session_state:
    st.session_state.structure_error = ""

options = dict()


def rerun():
    """Rerun on both current and legacy Streamlit versions."""
    if hasattr(st, "rerun"):
        st.rerun()
    else:
        st.experimental_rerun()


# ---------------------------------------------------------------------------
# Data helpers (unchanged behaviour)
# ---------------------------------------------------------------------------
def get_cid_from_name(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        try:
            cids = response.json()["IdentifierList"]["CID"]
            return cids[0]  # Return the first matching CID
        except (KeyError, IndexError):
            return None
    else:
        return None


def get_pubchem_product_code(product_name):
    product_code_from_pubchem = ""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{product_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            smiles = response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            product_code_from_pubchem = smiles
        except (KeyError, IndexError):
            product_code_from_pubchem = "NO DRUG FOUND"
    else:
        product_code_from_pubchem = "NO DRUG FOUND"
    if product_code_from_pubchem == "NO DRUG FOUND":
        return ""
    else:
        return product_code_from_pubchem


def showStructure(product_name):
    product_code = ""
    product_code_from_pubchem = get_pubchem_product_code(product_name)
    if product_code_from_pubchem == "":
        product_code_prompt = prompts.STRUCTURE_PROMPT.substitute(product_name=product_name)
        print("Prompt is: " + product_code_prompt)
        product_code = chat_with_gpt.chatWithGpt(product_code_prompt)
        if product_code == "NO DRUG FOUND":
            return "", ""
    else:
        product_code = product_code_from_pubchem

    print("product code is: " + product_code)
    print("product code from pubchem: " + product_code_from_pubchem)

    product_code = (product_code or "").strip()
    m = Chem.MolFromSmiles(product_code)
    if m is None:
        return "", product_code
    fig = Draw.MolToImage(m, size=size)
    return fig, product_code


# ---------------------------------------------------------------------------
# Masthead
# ---------------------------------------------------------------------------
JURISDICTIONS = [
    "INDIAN PHARMACOPIEA",
    "BRITISH PHARMACOPIEA",
    "UNITED STATES PHARMACOPOEIA",
    "MARTINDALE-EXTRA PHARMACOPIEA",
    "COMPARE WITH ALL",
]

INFO_TYPES = [
    "METHOD OF PREPARATION",
    "CHARACTARIZATION/EVALUATION",
    "Both of above",
    "CHECK RESULTS",
]

st.markdown(
    """
    <div class="qai-mast">
      <div class="mark">🧪</div>
      <div>
        <h1>QAI Model</h1>
        <p>Generate pharmaceutical quality reports from monograph references and your own lab results.</p>
      </div>
      <div class="stamp">
        Reference sources<br>
        <b>PubChem</b> · <b>Pharmacopoeia</b>
      </div>
    </div>
    """,
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------------------
# FORM PAGE
# ---------------------------------------------------------------------------
if st.session_state.page == "form":

    left, right = st.columns([1.25, 1], gap="large")

    # ----- Product details -----
    with left:
        st.markdown(
            """
            <div class="qai-panel-head">
              <span class="t">Product details</span>
              <span class="s">All three fields are required</span>
            </div>
            """,
            unsafe_allow_html=True,
        )

        options["product_name"] = st.text_input(
            "💊 Product name",
            placeholder="e.g., Paracetamol",
            key="in_product_name",
        )

        c1, c2 = st.columns(2)
        with c1:
            options["quanOfMed"] = st.text_input(
                "📦 Quantity of medicine",
                placeholder="e.g., 1000 tablets",
                key="in_quan",
            )
        with c2:
            options["powerOfDrug"] = st.text_input(
                "⚡ Power of drug",
                placeholder="e.g., 500 mg",
                key="in_power",
            )

        options["jurisdiction"] = st.selectbox(
            "🌎 Jurisdiction",
            JURISDICTIONS,
            key="in_jurisdiction",
        )

        options["typeOfInfo"] = st.selectbox(
            "📊 Information required",
            INFO_TYPES,
            key="in_type_of_info",
        )

        if options["typeOfInfo"] == "CHECK RESULTS":
            options["resultsToCheck"] = st.text_area(
                "🔍 Your results",
                height=200,
                placeholder="Paste lab results here...",
                key="checkResults",
            )

        options["ftir_required"] = st.checkbox(
            "📡 Retrieve FTIR data with the report",
            key="in_ftir",
        )

        st.markdown("<div style='height:6px'></div>", unsafe_allow_html=True)

        submit_button = st.button(
            "🚀 Submit & Generate Report",
            type="primary",
            key="btn_submit",
        )

        st.markdown("<div style='height:10px'></div>", unsafe_allow_html=True)

    # ----- Structure preview -----
    with right:
        st.markdown(
            """
            <div class="qai-panel-head">
              <span class="t">Molecular structure</span>
              <span class="s">Optional</span>
            </div>
            """,
            unsafe_allow_html=True,
        )

        structure_button = st.button("🔬 Get structure", key="btn_structure")

        if structure_button:
            st.session_state.structure_fig = None
            st.session_state.structure_error = ""
            st.session_state.structure_label = ""
            if ("product_name" not in options) or ("product_name" in options and options["product_name"] == ""):
                st.session_state.structure_error = "Enter a product name first, then get the structure."
            else:
                with st.spinner("🛠️ Looking up the structure..."):
                    fig, code = showStructure(options["product_name"])
                if fig == "":
                    st.session_state.structure_error = (
                        "No structure found for that name. Check the spelling or try the generic name."
                    )
                else:
                    st.session_state.structure_fig = fig
                    st.session_state.structure_label = f"{options['product_name']} — {code}"

        if st.session_state.structure_error:
            st.error(f"⚠️ {st.session_state.structure_error}")
        elif st.session_state.structure_fig is not None:
            st.image(st.session_state.structure_fig, caption=st.session_state.structure_label)
        else:
            st.markdown(
                """
                <div class="qai-viewer">
                  The 2D structure appears here.<br>
                  Looked up on PubChem first, with the model as a fallback.
                </div>
                """,
                unsafe_allow_html=True,
            )

    # ----- Submit handling -----
    if submit_button:
        if not all([options["product_name"], options["quanOfMed"], options["powerOfDrug"]]):
            st.error("⚠️ Product name, quantity and power are all needed before generating a report.")
        else:
            prompt = prompts.getPromptForOptions(options)
            with st.spinner("🛠️ Processing... Please wait"):
                api_response = chat_with_gpt.chatWithGpt(prompt)
                st.session_state.api_response = api_response

            st.session_state.update(options)
            st.session_state.page = "result"
            rerun()


# ---------------------------------------------------------------------------
# RESULT PAGE
# ---------------------------------------------------------------------------
elif st.session_state.page == "result":

    nav_left, nav_right = st.columns([1, 3])
    with nav_left:
        if st.button("🔙 Go Back to Form", key="btn_back"):
            st.session_state.page = "form"
            rerun()

    st.markdown("<div style='height:18px'></div>", unsafe_allow_html=True)

    summary_col, report_col = st.columns([1, 2.1], gap="large")

    with summary_col:
        st.markdown(
            """
            <div class="qai-panel-head">
              <span class="t">📑 Submission summary</span>
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.markdown(
            f"""
            <div class="qai-spec">
              <div class="row"><span class="k">💊 Product name</span><span class="v">{st.session_state.product_name}</span></div>
              <div class="row"><span class="k">📦 Quantity of medicine</span><span class="v">{st.session_state.quanOfMed}</span></div>
              <div class="row"><span class="k">⚡ Power of drug</span><span class="v">{st.session_state.powerOfDrug}</span></div>
              <div class="row"><span class="k">🌎 Jurisdiction</span><span class="v">{st.session_state.get("jurisdiction", "—")}</span></div>
              <div class="row"><span class="k">📊 Information required</span><span class="v">{st.session_state.get("typeOfInfo", "—")}</span></div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        if st.session_state.structure_fig is not None:
            st.markdown("<div style='height:16px'></div>", unsafe_allow_html=True)
            st.image(st.session_state.structure_fig, caption=st.session_state.structure_label)

    with report_col:
        st.markdown(
            """
            <div class="qai-panel-head">
              <span class="t">📋 Generated report</span>
            </div>
            """,
            unsafe_allow_html=True,
        )

        if st.session_state.api_response:
            st.markdown('<div class="qai-report">', unsafe_allow_html=True)
            st.markdown(st.session_state.api_response)
            # components.html(st.session_state.api_response, height=1000, width=1000, scrolling=True)
            st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.warning("⚠️ No response received from API.")

        if st.session_state.get("ftir_required"):
            st.markdown("<div style='height:18px'></div>", unsafe_allow_html=True)
            st.markdown(
                """
                <div class="qai-panel-head">
                  <span class="t">🔬 FTIR data</span>
                </div>
                """,
                unsafe_allow_html=True,
            )
            with st.spinner("📡 Fetching FTIR Data..."):
                ftir_data = chat_with_gpt.get_ftir_from_gpt(st.session_state.product_name)
            st.markdown('<div class="qai-report">', unsafe_allow_html=True)
            st.write(ftir_data)
            st.markdown("</div>", unsafe_allow_html=True)
