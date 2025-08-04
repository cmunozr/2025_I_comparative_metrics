from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException, TimeoutException
import time

# Helper Function for Automation
def automate_luke_download_simplified(year, theme_text, map_box_identifiers, email, accept_license=True):
    """
    Automates the download process on the Luke data service website.
    It navigates to the initial page to get the theme's numerical ID,
    then constructs and directly navigates to the download URL,
    and finally fills out the email and accepts the license on the download page.

    Args:
        year (int): The year to select for the data (e.g., 2019).
        theme_text (str): The exact visible text of the theme to select in the dropdown.
        map_box_identifiers (list): A list of map box identifiers (e.g., ["S4", "S5"])
                                     to include in the download request.
        email (str): The email address to enter for the download link.
        accept_license (bool, optional): Whether to check the license acceptance box. Defaults to True.
    """
    driver = webdriver.Chrome()
    driver.maximize_window()

    try:
        print(f"\n--- Initiating download for Theme: '{theme_text}' (Year: {year}) with {len(map_box_identifiers)} boxes ---")

        # Step 1: Navigate to the initial data selection page
        print("1. Navigating to the initial data selection page...")
        driver.get("https://kartta.luke.fi/opendata/valinta-en.html")
        time.sleep(3) # Give browser a moment to load initial page

        # Step 2: Wait for dropdowns to be present and clickable
        print("2. Waiting for 'Year' and 'Theme' dropdowns to be ready...")
        WebDriverWait(driver, 20).until(
            EC.element_to_be_clickable((By.ID, "vuosivalikko"))
        )
        WebDriverWait(driver, 20).until(
            EC.element_to_be_clickable((By.ID, "taso"))
        )
        print("   Dropdowns are ready.")
        time.sleep(1) # Small pause after dropdowns are ready

        # Step 3: Select the desired Year
        print(f"3. Selecting Year: {year}...")
        year_dropdown_element = Select(driver.find_element(By.ID, "vuosivalikko"))
        try:
            year_dropdown_element.select_by_value(str(year))
            print(f"   Selected year by value: {year}")
        except NoSuchElementException:
            year_dropdown_element.select_by_visible_text(str(year))
            print(f"   Selected year by visible text (fallback): {year}")
        time.sleep(2) # Allow page to update theme options after year selection

        # Step 4: Select the desired Theme and extract its numerical prefix (e.g., '362')
        print(f"4. Selecting Theme: '{theme_text}' and extracting its prefix...")
        theme_dropdown_element = Select(driver.find_element(By.ID, "taso"))
        try:
            theme_dropdown_element.select_by_visible_text(theme_text)
            print(f"   Selected theme: '{theme_text}'")
            # After selection, get the 'value' attribute of the currently selected option
            selected_option = theme_dropdown_element.first_selected_option
            theme_prefix = selected_option.get_attribute("value")
            print(f"   Extracted theme prefix: {theme_prefix}")

            if theme_prefix == "-1": # Check for the "Choose theme" default value
                print(f"Error: Theme '{theme_text}' was not correctly selected or found for year {year}. Aborting this request.")
                return # Exit the function for this specific request

        except NoSuchElementException:
            print(f"Error: Theme '{theme_text}' not found in dropdown for year {year}. Aborting this request.")
            return # Exit the function for this specific request

        # Step 5: Construct the direct download URL
        print("5. Constructing the direct download URL...")
        base_download_url = "https://kartta.luke.fi/lataus/palvelu.2021" # This part is fixed
        
        # Format each map box identifier with the extracted theme_prefix
        formatted_boxes = [f"{theme_prefix}+{box_id}" for box_id in map_box_identifiers]
        query_string = ";".join(formatted_boxes) # Join multiple boxes with semicolons
        
        final_download_url = f"{base_download_url}?{query_string}"
        print(f"   Constructed URL: {final_download_url}")

        # Step 6: Navigate directly to the constructed download URL
        print("6. Navigating directly to the constructed download URL (download page)...")
        driver.get(final_download_url)

        # Step 7: Wait for elements on the download page (email input, license checkbox) to be ready
        print("7. Waiting for 'Email' input and 'License' checkbox on the download page...")
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.XPATH, "//input[@type='email' or @type='text']"))
        )
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.XPATH, "//input[@type='checkbox']"))
        )
        print("   Download page elements are ready.")

        # Step 8: Fill in the email address
        if email:
            print(f"8. Attempting to enter email: '{email}'...")
            try:
                # Find the email input field (using type='email' or type='text')
                email_input = driver.find_element(By.XPATH, "//input[@type='email' or @type='text']")
                email_input.send_keys(email)
                print(f"   Entered email successfully.")
            except NoSuchElementException:
                print("   Email input field not found. Cannot enter email.")
            except Exception as e:
                print(f"   An error occurred while entering email: {e}")

        # Step 9: Accept the license
        if accept_license:
            print("9. Attempting to accept the license...")
            try:
                # Find the license checkbox
                license_checkbox = driver.find_element(By.XPATH, "//input[@type='checkbox']")
                if not license_checkbox.is_selected(): # Only click if not already selected
                    license_checkbox.click()
                    print("   Accepted license successfully.")
            except NoSuchElementException:
                print("   License checkbox not found. Cannot accept license.")
            except Exception as e:
                print(f"   An error occurred while accepting license: {e}")

        # Step 10: Click the "Send the link for download to my email" button (Finnish version)
        print("10. Clicking the 'Lähetä latausosoite sähköpostiin' button...")
        try:
            send_link_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.XPATH, "//input[@type='submit' and @value='Lähetä latausosoite sähköpostiin']"))
            )
            send_link_button.click()
            print("    Button clicked. Waiting for email processing...")
            time.sleep(20) # Give time for the server to process the request and send email
        except (NoSuchElementException, TimeoutException):
            print("    'Lähetä latausosoite sähköpostiin' button not found or not ready. Manual interaction might be needed for this request.")
        except Exception as e:
            print(f"    An error occurred while clicking 'Send link' button: {e}")

        time.sleep(5) # Brief pause before closing browser

    except Exception as e:
        print(f"\n!!! An unexpected error occurred during automation for this request: {e}")

    finally:
        print("Closing browser for this request.")
        driver.quit()

# --- Configuration ---

# All known map box identifiers from the entire region
all_map_boxes = [
    "X4", "X5", "W5", "W4", "W3", "V3", "V4", "V5", "U4", "U5",
    "T5", "T4", "S4", "S5", "R5", "R4", "Q3", "Q4", "Q5", "P3",
    "P4", "P5", "P6", "N3", "N4", "N5", "N6", "M3", "M4", "M5",
    "L2", "L3", "L4", "L5", "K2", "K3", "K4"
]

# The specific themes, with year-specific templates.
# Each key (common name) maps to a dictionary of year-specific templates or a 'default'.
desired_themes_template = {
    #"Stand mean diameter": {
    #    "2009": "Mean diameter of stand {} (cm)",
    #    "default": "Stand mean diameter of {} (cm)"
    #},
    #"Stand mean height": {
    #    "2009": "Mean height of stand {} (dm)", 
    #    "default": "Stand mean height {} (dm)"
    #},
    #"Canopy cover broad-leaved trees": {
    #    "2009": "Canopy cover of broad leaved trees {} (%)",
    #    "default": "Canopy cover of broad-leaved trees {} (%)"
    #},
    #"Canopy cover": {
    #    "2009": "Canopy cover {} (%)",
    #    "default": "Canopy cover {} (%)"
    #},
    #"Stand age": {
    #    "2009": "Stand age of growing stock {} (year)",
    #    "default": "Stand age {} (year)"
    #},
    #"Stand basal area": {
    #    "2009": "Stand basal area {} (m2/ha)", 
    #    "default": "Stand basal area {} (m2/ha)"
    #},
    "Volume, the growing stock": {
        "2009": "Volume of the growing stock {} (m3/ha)",
        "default": "Volume, the growing stock {} (m3/ha)"
    },
    "Volume, spruce": {
        "2009": "Spruce volume {} (m3/ha)",
        "default": "Volume, spruce {} (m3/ha)"
    }
}


# The years for which data is available and want to download (excluding 2006).
available_years = [2009, 2011, 2013, 2015, 2017, 2019, 2021]

# Maximum number of map boxes to request per single download URL.
batch_size = 19

# --- Helper Function for Batching ---

def chunks(lst, n):
    """Yield successive n-sized chunks from a list."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# --- Main Automation Execution Loop ---

if __name__ == "__main__":
    # Create batches of all map box identifiers
    map_box_batches = list(chunks(all_map_boxes, batch_size))

    # Outer loop: Iterate through each desired theme category
    for theme_key, year_templates in desired_themes_template.items():
        print(f"\n==================================================")
        print(f"*** Starting Processing for Theme Category: {theme_key} ***")
        print(f"==================================================")

        # Middle loop: Iterate through each available year
        for year in available_years:
            # Determine the correct theme text template based on the current year
            # .get(key, default_value) is used here.
            # It tries to find a template for `str(year)` (e.g., "2009").
            # If not found, it falls back to the "default" template.
            theme_format_string = year_templates.get(str(year), year_templates["default"])

            # Construct the exact theme text for the current year
            # For 2009, it will use the "Mean diameter {} (cm)" template.
            # For other years, it will use "Stand mean diameter of {} (cm)"
            current_theme_text = theme_format_string.format(year)
            
            print(f"\n--- Starting downloads for Theme: '{current_theme_text}' (Year: {year}) ---")

            # Inner loop: Iterate through each batch of map boxes for the current theme and year
            for i, batch in enumerate(map_box_batches):
                print(f"\n----- Processing Batch {i+1} of {len(map_box_batches)} for '{current_theme_text}' -----")
                
                # Construct a unique email address for this specific request for better organization
                email_suffix = f"{theme_key.replace(' ', '_').replace('(', '').replace(')', '').replace(',', '')}_{year}_batch{i+1}"
                actual_email = f"cmunozbiol+{email_suffix}@gmail.com"
                
                # Execute the automation helper function
                automate_luke_download_simplified(
                    year=year,
                    theme_text=current_theme_text,
                    map_box_identifiers=batch,
                    email=actual_email,
                    accept_license=True
                )
                
                # Pause between batches to be polite to the server
                if i < len(map_box_batches) - 1:
                    print(f"\n--- Waiting 10 seconds before next map box batch for '{current_theme_text}' ---")
                    time.sleep(10)
            
            print(f"\n--- Finished all map box batches for '{current_theme_text}' ---")
            # Pause between different years for the same theme
            print(f"\n--- Waiting 10 seconds before starting next year for {theme_key} ---")
            time.sleep(10)

        print(f"\n==================================================")
        print(f"*** Finished All Years for Theme Category: {theme_key} ***")
        print(f"==================================================")
        # Pause between different theme categories
        print(f"\n--- Waiting 10 seconds before starting the next theme category ---")
        time.sleep(10)

    print("\n==================================================")
    print("ALL SPECIFIED THEMES ACROSS ALL AVAILABLE YEARS HAVE BEEN PROCESSED.")
    print("Check email for download links.")
    print("==================================================")