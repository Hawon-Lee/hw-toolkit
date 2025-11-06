from playwright.sync_api import sync_playwright
import json

pdb_ids = ["5aqo", "5th4"]
base_url = "https://www.pdbbind-plus.org.cn/browse/{}"

def get_ligand_properties(page, pdb_id):
    url = base_url.format(pdb_id)
    page.goto(url, timeout=30000)  # 30초 타임아웃 설정

    # ligand-properties를 포함하는 엘리먼트가 나타날 때까지 기다림
    try:
        page.wait_for_selector("div.ligand-properties", timeout=15000)
    except:
        print(f"{pdb_id}: ligand properties 로딩 실패.")
        return {}

    properties = page.eval_on_selector_all(
        "div.ligand-properties tr",
        """
        rows => rows.map(row => {
            const cells = row.querySelectorAll('td');
            return cells.length >= 2 ? [cells[0].innerText.trim(), cells[1].innerText.trim()] : null;
        }).filter(item => item !== null)
        """
    )

    return dict(properties)

ligand_data = {}

with sync_playwright() as pw:
    browser = pw.chromium.launch(headless=True)
    context = browser.new_context()
    page = context.new_page()

    for pdb_id in pdb_ids:
        ligand_data[pdb_id] = get_ligand_properties(page, pdb_id)

    browser.close()

with open("ligand_properties.json", "w") as f:
    json.dump(ligand_data, f, indent=2)

print("Ligand properties 수집 완료.")
