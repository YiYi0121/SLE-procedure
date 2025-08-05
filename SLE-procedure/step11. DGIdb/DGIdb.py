import time
from selenium.webdriver.common.by import By
from selenium.webdriver import Chrome
from selenium.webdriver.common.keys import Keys
from selenium.webdriver import ActionChains
from selenium.webdriver.chrome.options import Options
import requests
import os
import csv


def login():
    url = 'https://www.dgidb.org/search_interactions'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches', ['enable-automation'])
    # opt.add_argument('--headless')
    web = Chrome(options=opt)
    web.get(url)
    return web


datalists = []

data = []
temp_list = []

def get_data(web):
    Name = web.find_element(By.XPATH, '//*[@id="STAT1-unique"]/div[1]/div/h3/a').text
    print(Name)
    ALL = web.find_element(By.XPATH,'//*[@id="STAT1-unique"]/div[2]/div/div/button[1]')
    web.execute_script("arguments[0].click();", ALL)
    time.sleep(5)
    table1 = web.find_elements(By.XPATH, '//*[@id="STAT1-STAT1"]/tbody[3]/tr')
    for tr in table1:
        Drug = tr.find_element(By.XPATH, './td[1]').text
        Interaction = tr.find_element(By.XPATH, './td[2]').text
        Score = tr.find_element(By.XPATH, './td[6]').text
        data.append(Name)
        data.append(Drug)
        data.append(Interaction)
        data.append(Score)
        temp_list = data[:]
        datalists.append(temp_list)
        del data[0:4]
        print(datalists)
    Name = web.find_element(By.XPATH, '//*[@id="JAK2-unique"]/div[1]/div/h3/a').text
    print(Name)
    ALL = web.find_element(By.XPATH, '//*[@id="JAK2-unique"]/div[2]/div/div/button[1]')
    web.execute_script("arguments[0].click();", ALL)
    time.sleep(5)
    table2 = web.find_elements(By.XPATH, '//*[@id="JAK2-JAK2"]/tbody[3]/tr')
    for tr in table2:
        Drug = tr.find_element(By.XPATH, './td[1]').text
        Interaction = tr.find_element(By.XPATH, './td[2]').text
        Score = tr.find_element(By.XPATH, './td[6]').text
        data.append(Name)
        data.append(Drug)
        data.append(Interaction)
        data.append(Score)
        temp_list = data[:]
        datalists.append(temp_list)
        del data[0:4]
        print(datalists)
    return datalists


def is_exit(web):
    try:
        web.find_element(By.XPATH,'//*[@id="disclaimer-modal"]/div')
        return True
    except:
        return False


data_list = []

def page_spider():
    f = open('./interactions.csv', mode='w', encoding='utf-8', newline='')
    writer = csv.writer(f)
    writer.writerow(['Name', 'Drug', 'Interaction Type & Directionality', 'Interaction Score'])
    web = login()
    time.sleep(10)
    if is_exit(web):
        close = web.find_element(By.XPATH,'//*[@id="disclaimer-modal"]/div/div/div[1]/button')
        web.execute_script("arguments[0].click();", close)
        time.sleep(2)
        search = 'STAT1, JAK2'
        web.find_element(By.XPATH, '//*[@id="identifiers"]').send_keys(search)
        time.sleep(5)
        enter = web.find_element(By.XPATH, '//*[@id="search-btn"]')
        web.execute_script("arguments[0].click();", enter)
        time.sleep(20)
        data_list = get_data(web)
        writer_list = []
        for detail in data_list:
            writer_list.append(detail[0])
            writer_list.append(detail[1])
            writer_list.append(detail[2])
            writer_list.append(detail[3])
            writer.writerow(writer_list)
            del writer_list[0:4]
    f.close()
    web.close()


if __name__ == '__main__':
    page_spider()
