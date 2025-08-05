import time
from selenium.webdriver.common.by import By
from selenium.webdriver import Chrome
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
import csv 


def login():
    url = 'https://db.idrblab.net/ttd/'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches',['enable-automation'])
    #opt.add_argument('--headless')
    web = Chrome(options=opt)
    web.get(url)
    return web

def is_exit(web):
    try:
        web.find_element(By.XPATH,'//*[@id="fixed-width-page"]/div/main/div[2]/nav[1]/ul/li[5]/a')
        return True
    except:
        return False


Gene_list = []
templist=[]
Genelist = []
def get_gene():
    global Gene
    f = open('./TTD.csv', mode='w', encoding='utf-8', newline='')
    writer = csv.writer(f)
    writer.writerow(['TTD'])
    web =login()
    time.sleep(5)
    search = "systemic lupus erythematosus"
    web.find_element(By.XPATH, '//*[@id="fixed-width-page"]/div/main/div[2]/div[1]/div/div/div[1]/div/form[1]/table/tbody/tr[2]/td[1]/input').send_keys(search, Keys.ENTER)
    time.sleep(2)
    page = 1
    while page < 5:
        tr_list = web.find_elements(By.XPATH, '//*[@id="fixed-width-page"]/div/main/div[2]/table')
        for tr in tr_list:
            tr.find_element(By.XPATH,'./tbody/tr[1]/th[1]/span/a').click()
            time.sleep(2)
            GeneName = web.find_element(By.XPATH,'//*[@id="table-target-general"]/tbody/tr[5]/td/div').text
            web.back()
            web.switch_to.window(web.window_handles[0])
            Genelist.append(GeneName)
            templist = Genelist[:]
            Gene_list.append(templist)
            del Genelist[0]
        if is_exit(web):
            web.find_element(By.XPATH, '//*[@id="fixed-width-page"]/div/main/div[2]/nav[1]/ul/li[5]/a').click()
            time.sleep(2)
        else:
            break
        page += 1
    writerlist = []
    for detail in Gene_list:
        writerlist.append(detail[0])
        writer.writerow(writerlist)
        del writerlist[0]
    f.close()
    web.close()




if __name__ == '__main__':
    get_gene()
