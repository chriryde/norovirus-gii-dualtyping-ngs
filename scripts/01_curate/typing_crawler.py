from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError
from pathlib import Path

def run_norovirus_typing_tool(INPUT_FASTA, OUTPUT_DIR):
    URL = "https://mpf.rivm.nl/mpf/typingtool/norovirus/"

    INPUT_FASTA = Path(INPUT_FASTA).resolve()
    OUTPUT_DIR = Path(OUTPUT_DIR).resolve()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        page.goto(URL, wait_until="domcontentloaded")

        file_input_selector = 'input[type="file"][name="data"]'
        page.set_input_files(file_input_selector, str(INPUT_FASTA)) 

        page.wait_for_function(
            """() => {
                const el = document.querySelector('input[type="file"][name="data"]');
                return el && el.files && el.files.length > 0;
            }""",
            timeout=10000,
        )

        page.wait_for_function(
                """() => {
                    const el = document.querySelector('span.error-text.error');
                    if (!el) return false;
                    const style = window.getComputedStyle(el);
                    return style.display === 'none';
                }""",
                timeout=30000,
            )
        
        start_button_selector = 'button[id^="button-run_"]'
        page.locator(start_button_selector).click()

        page.wait_for_url("**/job/**", timeout=120000)

        job_url = page.url.rstrip("/")
        job_id = job_url.split("/")[-1]

        csv_link = page.locator('a[href*="results.csv"]')

        with page.expect_download(timeout=120000) as download_info:
            csv_link.click()

        download = download_info.value
        output_path = OUTPUT_DIR / f"{INPUT_FASTA.stem}_job_{job_id}_table.csv"
        download.save_as(output_path)

        print(f"Downloaded CSV to: {output_path}")
        browser.close()
        
        return output_path
    
    
