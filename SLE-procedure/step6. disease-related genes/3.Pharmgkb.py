import time
from selenium.webdriver.common.by import By
from selenium.webdriver import Chrome
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
import csv 


def login():
    url = 'https://www.pharmgkb.org/'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches',['enable-automation'])
    opt.add_argument('--headless')
    web = Chrome(options=opt)
    web.get(url)
    return web

def is_exit(web):
    try:
        web.find_element(By.XPATH,'//*[@id="app"]/div/div[1]/div/div/div/div[2]/nav[1]/a[2]')
        return True
    except:
        return False


def Searchpage(web):
    number = web.find_element(By.XPATH,'//*[@id="app"]/div/div[1]/div/div/div/div[2]/nav[1]/span[2]').text
    number_str = number.split(' ',)[3]
    print(int(int(number_str)))
    if(int(number_str) % 20) == 0:
        pageNumber = int (number_str) / 20
    else:
        pageNumber = int(int(number_str)/ 20) +1
    return pageNumber

Gene_list = []
templist=[]
Genelist = []
def get_gene():
    global Gene
    f = open('./pharmgkb.csv', mode='w', encoding='utf-8', newline='')
    writer = csv.writer(f)
    writer.writerow(['pharmgkb'])
    web =login()
    time.sleep(5)
    search = "systemic lupus erythematosus"
    web.find_element(By.XPATH, '//*[@id="app"]/div/div[1]/div/div[1]/div/div/div/input').send_keys(search, Keys.ENTER)
    time.sleep(2)
    web.find_element(By.XPATH,'//*[@id="geneChkControl"]').click()
    time.sleep(2)
    totalpage = Searchpage(web)
    print(totalpage)
    for i in range(1, totalpage+1):
        td_list = web.find_elements(By.XPATH, '//*[@id="app"]/div/div[1]/div/div/div/div[2]/div[2]/div')
        for div in td_list:
            Gene = div.find_element(By.XPATH, './div[1]/div[2]/div[2]/a').text
            Genelist.append(Gene)
            templist = Genelist[:]
            Gene_list.append(templist)
            del Genelist[0]
        if i != 1:
            if is_exit(web):
                web.find_element(By.XPATH, '//*[@id="app"]/div/div[1]/div/div/div/div[2]/nav[1]/a[2]').click()
                time.sleep(2)
            else:
                break
        else:
            web.find_element(By.XPATH, '//*[@id="app"]/div/div[1]/div/div/div/div[2]/nav[1]/a').click()
            time.sleep(2)
    writerlist = []
    for detail in Gene_list:
        writerlist.append(detail[0])
        print(writerlist)
        writer.writerow(writerlist)
        del writerlist[0]
    f.close()
    web.close()


if __name__ == '__main__':
    get_gene()
