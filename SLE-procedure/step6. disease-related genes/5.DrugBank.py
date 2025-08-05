import time
from selenium.webdriver.common.by import By
from selenium.webdriver import Chrome
from selenium.webdriver.chrome.options import Options
import csv


def login():
    url = 'https://www.drugbank.ca/'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches', ['enable-automation'])
    # opt.add_argument('--headless')
    opt.add_argument(
        'user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36')
    web = Chrome(options=opt)
    web.get(url)
    return web


def is_exit(web):
    try:
        web.find_element(By.XPATH, '//*[@id="DataTables_Table_0_next"]/a')
        return True
    except:
        return False


def get_page(web):
    number = web.find_element(By.XPATH, '/html/body/main/div/div[3]/div/div[2]/div/div[3]/div[1]/div').text
    number_str = number.split(' ', )[5]
    if (int(number_str) % 10) == 0:
        pageNumber = int(number_str) / 10
    else:
        pageNumber = int(int(number_str) / 10) + 1
    return pageNumber


def is_exit1(web):
    try:
        web.find_element(By.XPATH, '/html/body/main/div/div[2]/dl[1]/dd[4]/table/tbody/tr/td[3]/a')
        return True
    except:
        return False


Gene_list = []
templist = []
Genelist = []


def get_data(web):
    total_page = get_page(web)
    print(total_page)
    for i in range(1, total_page+1):
        tr_list = web.find_elements(By.XPATH, '/html/body/main/div/div[3]/div/div[2]/div/div[2]/div/table/tbody/tr')
        for tr in tr_list:
            target = tr.find_element(By.XPATH, './td[3]/a').get_attribute("href")
            js = "window.open('" + target + "');"
            web.execute_script(js)
            time.sleep(2)
            web.switch_to.window(web.window_handles[3])
            time.sleep(3)
            if is_exit1(web):
                web.find_element(By.XPATH, '/html/body/main/div/div[2]/dl[1]/dd[4]/table/tbody/tr/td[3]/a').click()
                time.sleep(4)
                GeneName = web.find_element(By.XPATH, '/html/body/main/div/div[2]/dl[1]/dd[3]').text
                time.sleep(3)
                web.close()
                time.sleep(2)
                web.switch_to.window(web.window_handles[2])
                time.sleep(3)
                Genelist.append(GeneName)
                templist = Genelist[:]
                Gene_list.append(templist)
                del Genelist[0]
            else:
                web.close()
                web.switch_to.window(web.window_handles[2])
        if is_exit(web):
            web.find_element(By.XPATH, '//*[@id="DataTables_Table_0_next"]/a').click()
            time.sleep(3)
        else:
            break
    web.close()
    return Gene_list




Genelists = []


def page_spider():
    f = open('./DrugBank.csv', mode='w', encoding='utf-8', newline='')
    writer = csv.writer(f)
    writer.writerow(['DrugBank'])
    web = login()
    time.sleep(4)
    print("请先在弹出的网页中勾选验证按钮，再执行下面内容！")
    search = input('请输入搜索内容(这里是：systemic lupus erythematosus)：')
    web.find_element(By.XPATH, '//*[@id="query"]').send_keys(search)
    time.sleep(2)
    web.find_element(By.XPATH, '/html/body/main/div/div[1]/div[1]/div[2]/form/div[2]/label[4]/span').click()
    time.sleep(3)
    web.find_element(By.XPATH, '/html/body/main/div/div[1]/div[1]/div[2]/form/div[1]/div/div/button').click()
    time.sleep(4)
    systemic = web.find_element(By.XPATH,
                                '/html/body/main/div/div[2]/div[2]/div/div/div[1]/div/h2/a').get_attribute("href")
    js = "window.open('" + systemic + "');"
    web.execute_script(js)
    time.sleep(2)
    web.switch_to.window(web.window_handles[1])
    time.sleep(2)
    drug = web.find_element(By.XPATH, '/html/body/main/div/div[3]/ul/li[2]/a').get_attribute('href')
    js = "window.open('" + drug + "');"
    web.execute_script(js)
    time.sleep(2)
    web.switch_to.window(web.window_handles[2])
    time.sleep(4)
    Genelists = get_data(web)
    writerlist = []
    for detail in Genelists:
        writerlist.append(detail[0])
        print(writerlist)
        writer.writerow(writerlist)
        del writerlist[0]
    f.close()


if __name__ == '__main__':
    page_spider()
