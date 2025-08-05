import time
from selenium.webdriver.common.by import By
from selenium.webdriver import Chrome
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
import requests
import os
import csv


def login():
    url = 'https://tcmsp-e.com/tcmspsearch.php?qs=herb_all_name&token=332a3150d47fc86ba5d7ee98903de4ff'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches', ['enable-automation'])
    opt.add_argument('--headless')
    web = Chrome(options=opt)
    web.get(url)
    return web


def login1():
    url1 = 'https://pubchem.ncbi.nlm.nih.gov'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches', ['enable-automation'])
    opt.add_argument('--headless')
    web1 = Chrome(options=opt)
    web1.get(url1)
    return web1


def login2():
    url2 = 'http://www.swissadme.ch/index.php'
    opt = Options()
    opt.page_load_strategy = 'eager'
    opt.add_experimental_option('excludeSwitches', ['enable-automation'])
    opt.add_argument('--headless')
    web2 = Chrome(options=opt)
    web2.get(url2)
    return web2


def is_exist_pubchemCID(web):
    try:
        web.find_element(By.XPATH, '//*[@id="container_lsp"]/div[2]/table/thead/tr[6]/td/a')
        return True
    except:
        return False


def is_exist(web1):
    try:
        web1.find_element(By.XPATH,
                          '//*[@id="collection-results-container"]/div/div/div[2]/ul/li/div/div/div[1]/div[2]/div[1]/a')
        return True
    except:
        return False


def is_exist1(web1):
    try:
        web1.find_element(By.XPATH,
                          '//*[@id="collection-results-container"]/div/div/div[2]/ul/li/div/div/div[1]/div[2]/div[1]/div/span/a')
        return True
    except:
        return False


def Searchpage(web):
    number = web.find_element(By.XPATH, '//*[@id="grid"]/div[3]/span ').text
    number_str = number.split(' ', )[4]
    if (int(number_str) % 15) == 0:
        pageNumber = int(number_str) / 15
    else:
        pageNumber = int(int(number_str) / 15) + 1
    return pageNumber


datalist = []
data = []


def get_data(web):
    frame_count = 0
    name1 = web.find_element(By.XPATH, '/html/body/div[2]/div[2]/div[2]/div[1]').text
    name = name1.split('(')[0]
    totalpage = Searchpage(web)
    for i in range(1, totalpage + 1):
        tr_list = web.find_elements(By.XPATH, '//*[@id="grid"]/div[2]/table/tbody/tr')
        for tr in tr_list:
            MOL_ID = tr.find_element(By.XPATH, './td[1]').text
            MOlecule = tr.find_element(By.XPATH, './td[2]/a').get_attribute('innerText')
            MOlecule_url = tr.find_element(By.XPATH, './td[2]/a').get_attribute('href')
            js = "window.open('" + MOlecule_url + "');"
            web.execute_script(js)
            time.sleep(2)
            web.switch_to.window(web.window_handles[1])
            time.sleep(2)
            InChIKey = web.find_element(By.XPATH, '//*[@id="container_lsp"]/div[2]/table/thead/tr[5]/td').text
            if is_exist_pubchemCID(web):
                web.find_element(By.XPATH, '//*[@id="container_lsp"]/div[2]/table/thead/tr[6]/td/a').click()
                time.sleep(6)
                SMILES = web.find_element(By.XPATH, '//*[@id="Canonical-SMILES"]/div[2]/div[1]').text
                web2 = login2()
                time.sleep(3)
                web2.find_element(By.XPATH, '//*[@id="smiles"]').send_keys(SMILES)
                time.sleep(2)
                web2.find_element(By.XPATH, '//*[@id="submitButton"]').click()
                time.sleep(6)
                GIabsorption = web2.find_element(By.XPATH,
                                                 '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[12]/td[2]').text
                if GIabsorption == "High":
                    Lipinski = web2.find_element(By.XPATH,
                                                 '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[22]/td[2]').text
                    Ghose = web2.find_element(By.XPATH,
                                              '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[23]/td[2]').text
                    Veber = web2.find_element(By.XPATH,
                                              '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[24]/td[2]').text
                    Egan = web2.find_element(By.XPATH,
                                             '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[25]/td[2]').text
                    Muegge = web2.find_element(By.XPATH,
                                               '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[26]/td[2]').text
                    Lipinski_str = Lipinski.split(';')[0]
                    Ghose_str = Ghose.split(';')[0]
                    Veber_str = Veber.split(';')[0]
                    Egan_str = Egan.split(';')[0]
                    Muegge_str = Muegge.split(';')[0]
                    text = [Lipinski_str, Ghose_str, Veber_str, Egan_str, Muegge_str]
                    sub = 'Yes'
                    if text.count(sub) >= 2:
                        swissADME = 1
                    else:
                        swissADME = 0
                else:
                    Lipinski = web2.find_element(By.XPATH,
                                                 '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[22]/td[2]').text
                    Ghose = web2.find_element(By.XPATH,
                                              '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[23]/td[2]').text
                    Veber = web2.find_element(By.XPATH,
                                              '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[24]/td[2]').text
                    Egan = web2.find_element(By.XPATH,
                                             '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[25]/td[2]').text
                    Muegge = web2.find_element(By.XPATH,
                                               '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[26]/td[2]').text
                    swissADME = 0
                    Lipinski_str = Lipinski.split(';')[0]
                    Ghose_str = Ghose.split(';')[0]
                    Veber_str = Veber.split(';')[0]
                    Egan_str = Egan.split(';')[0]
                    Muegge_str = Muegge.split(';')[0]
                web2.close()
                web.back()
                web.close()
                web.switch_to.window(web.window_handles[0])
                time.sleep(2)
            else:
                web1 = login1()
                time.sleep(4)
                web1.find_element(By.XPATH,
                                  '//*[@id="main-content"]/div[1]/div/div[2]/div/div[2]/form/div/div[1]/input').send_keys(
                    InChIKey, Keys.ENTER)
                time.sleep(15)
                if is_exist(web1):
                    web1.find_element(By.XPATH,
                                      '//*[@id="collection-results-container"]/div/div/div[2]/ul/li/div/div/div[1]/div[2]/div[1]/a').click()
                    time.sleep(5)
                    SMILES = web1.find_element(By.XPATH, '//*[@id="Canonical-SMILES"]/div[2]/div[1]').text
                    time.sleep(3)
                    web1.close()
                    web2 = login2()
                    time.sleep(2)
                    web2.find_element(By.XPATH, '//*[@id="smiles"]').send_keys(SMILES)
                    time.sleep(2)
                    web2.find_element(By.XPATH, '//*[@id="submitButton"]').click()
                    time.sleep(6)
                    GIabsorption = web2.find_element(By.XPATH,
                                                     '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[12]/td[2]').text
                    if GIabsorption == "High":
                        Lipinski = web2.find_element(By.XPATH,
                                                     '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[22]/td[2]').text
                        Ghose = web2.find_element(By.XPATH,
                                                  '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[23]/td[2]').text
                        Veber = web2.find_element(By.XPATH,
                                                  '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[24]/td[2]').text
                        Egan = web2.find_element(By.XPATH,
                                                 '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[25]/td[2]').text
                        Muegge = web2.find_element(By.XPATH,
                                                   '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[26]/td[2]').text

                        Lipinski_str = Lipinski.split(';')[0]
                        Ghose_str = Ghose.split(';')[0]
                        Veber_str = Veber.split(';')[0]
                        Egan_str = Egan.split(';')[0]
                        Muegge_str = Muegge.split(';')[0]
                        text = [Lipinski_str, Ghose_str, Veber_str, Egan_str, Muegge_str]
                        sub = 'Yes'
                        if text.count(sub) >= 2:
                            swissADME = 1
                        else:
                            swissADME = 0
                    else:
                        Lipinski = web2.find_element(By.XPATH,
                                                     '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[22]/td[2]').text
                        Ghose = web2.find_element(By.XPATH,
                                                  '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[23]/td[2]').text
                        Veber = web2.find_element(By.XPATH,
                                                  '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[24]/td[2]').text
                        Egan = web2.find_element(By.XPATH,
                                                 '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[25]/td[2]').text
                        Muegge = web2.find_element(By.XPATH,
                                                   '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[26]/td[2]').text
                        swissADME = 0
                        Lipinski_str = Lipinski.split(';')[0]
                        Ghose_str = Ghose.split(';')[0]
                        Veber_str = Veber.split(';')[0]
                        Egan_str = Egan.split(';')[0]
                        Muegge_str = Muegge.split(';')[0]
                    web2.close()
                    time.sleep(3)
                    web.back()
                    web.close()
                    web.switch_to.window(web.window_handles[0])
                    time.sleep(3)
                else:
                    if is_exist1(web1):
                        web1.find_element(By.XPATH,
                                          '//*[@id="collection-results-container"]/div/div/div[2]/ul/li/div/div/div[1]/div[2]/div[1]/div/span/a').click()
                        time.sleep(5)
                        SMILES = web1.find_element(By.XPATH, '//*[@id="Canonical-SMILES"]/div[2]/div[1]').text
                        web1.close()
                        web2 = login2()
                        time.sleep(2)
                        web2.find_element(By.XPATH, '//*[@id="smiles"]').send_keys(SMILES)
                        time.sleep(2)
                        web2.find_element(By.XPATH, '//*[@id="submitButton"]').click()
                        time.sleep(6)
                        GIabsorption = web2.find_element(By.XPATH,
                                                         '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[12]/td[2]').text
                        if GIabsorption == "High":
                            Lipinski = web2.find_element(By.XPATH,
                                                         '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[22]/td[2]').text
                            Ghose = web2.find_element(By.XPATH,
                                                      '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[23]/td[2]').text
                            Veber = web2.find_element(By.XPATH,
                                                      '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[24]/td[2]').text
                            Egan = web2.find_element(By.XPATH,
                                                     '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[25]/td[2]').text
                            Muegge = web2.find_element(By.XPATH,
                                                       '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[26]/td[2]').text
                            Lipinski_str = Lipinski.split(';')[0]
                            Ghose_str = Ghose.split(';')[0]
                            Veber_str = Veber.split(';')[0]
                            Egan_str = Egan.split(';')[0]
                            Muegge_str = Muegge.split(';')[0]
                            text = [Lipinski_str, Ghose_str, Veber_str, Egan_str, Muegge_str]
                            sub = 'Yes'
                            if text.count(sub) >= 2:
                                swissADME = 1
                            else:
                                swissADME = 0
                        else:
                            Lipinski = web2.find_element(By.XPATH,
                                                         '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[22]/td[2]').text
                            Ghose = web2.find_element(By.XPATH,
                                                      '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[23]/td[2]').text
                            Veber = web2.find_element(By.XPATH,
                                                      '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[24]/td[2]').text
                            Egan = web2.find_element(By.XPATH,
                                                     '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[25]/td[2]').text
                            Muegge = web2.find_element(By.XPATH,
                                                       '//*[@id="sib_body"]/div[11]/div[1]/div[4]/table/tbody/tr[26]/td[2]').text
                            swissADME = 0
                            Lipinski_str = Lipinski.split(';')[0]
                            Ghose_str = Ghose.split(';')[0]
                            Veber_str = Veber.split(';')[0]
                            Egan_str = Egan.split(';')[0]
                            Muegge_str = Muegge.split(';')[0]
                        time.sleep(2)
                        web2.close()
                        time.sleep(2)
                        web.back()
                        web.close()
                        web.switch_to.window(web.window_handles[0])
                        time.sleep(3)
                    else:
                        web1.close()
                        web.switch_to.window(web.window_handles[1])
                        time.sleep(1)
                        path = r"/桌面/结构图/"
                        if os.path.exists(path):
                            pass
                        else:
                            os.mkdir(path)
                        structure = web.find_element(By.XPATH,
                                                     '//*[@id="container_lsp"]/div[2]/table/thead/tr[3]/td/a/img')
                        frame_count += 1
                        img_path = structure.get_attribute('src')
                        img_name = path + "No." + name + str(frame_count) + "_" + MOlecule + "_" + img_path.split('/')[-1]
                        with open(img_name, 'wb') as f:
                            response = requests.get(img_path).content
                            f.write(response)
                            f.close()
                        SMILES = '-'
                        GIabsorption = '-'
                        Lipinski_str = '-'
                        Ghose_str = '-'
                        Veber_str = '-'
                        Egan_str = '-'
                        Muegge_str = '-'
                        swissADME = '-'
                        web.close()
                        web.switch_to.window(web.window_handles[0])
                        time.sleep(2)
            data.append(name)
            data.append(MOL_ID)
            data.append(MOlecule)
            data.append(MOlecule_url)
            data.append(InChIKey)
            data.append(SMILES)
            data.append(GIabsorption)
            data.append(Lipinski_str)
            data.append(Ghose_str)
            data.append(Veber_str)
            data.append(Egan_str)
            data.append(Muegge_str)
            data.append(swissADME)
            temp_list = data[:]
            datalist.append(temp_list)
            del data[0:13]
        web.find_element(By.XPATH,
                         '//*[@id="grid"]/div[@class="k-pager-wrap k-grid-pager k-widget"]/a[@title="Go to the next page"]').click()
        time.sleep(2)
    web.back()
    return datalist


def is_exist2(web):
    try:
        web.find_element(By.XPATH, '//*[@id="grid"]/div[2]/table/tbody/tr/td[3]/a')
        return True
    except:
        return False


data_list = []


def page_spider():
    f = open('./SwissADME.csv', mode='w', encoding='utf-8', newline='')
    writer = csv.writer(f)
    writer.writerow(
        ['drug', 'MOL ID', 'MOlecule', 'MOlecule_url', 'InChIKey', 'SMLIES', 'GIabsorption', 'Lipinski', 'Ghose',
         'Veber',
         'Egan', 'Muegge', 'SwissADME'])
    web = login()
    search = input('请输入搜索内容：')
    time.sleep(2)
    web.find_element(By.XPATH, '//*[@id="inputVarTcm"]').send_keys(search, Keys.ENTER)
    time.sleep(2)
    if is_exist2(web):
        web.find_element(By.XPATH, '//*[@id="grid"]/div[2]/table/tbody/tr/td[3]/a').click()
        time.sleep(2)
        data_list = get_data(web)
        print('page_data_list', data_list)
        writer_list = []
        for detail in data_list:
            writer_list.append(detail[0])
            writer_list.append(detail[1])
            writer_list.append(detail[2])
            writer_list.append(detail[3])
            writer_list.append(detail[4])
            writer_list.append(detail[5])
            writer_list.append(detail[6])
            writer_list.append(detail[7])
            writer_list.append(detail[8])
            writer_list.append(detail[9])
            writer_list.append(detail[10])
            writer_list.append(detail[11])
            writer_list.append(detail[12])
            writer.writerow(writer_list)
            del writer_list[0:13]
    else:
        print('数据库查不到此药物!')
    f.close()
    web.close()


if __name__ == '__main__':
    page_spider()
