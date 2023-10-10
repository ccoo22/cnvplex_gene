#!/home/genesky/software/python/3.9.4/bin/python3
import argparse
import glob
import os
import time

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(description="http://10.0.22.63:8888/soft/query.action?id= 下载探针设计结果", formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--input', '-i', type = str, required=True, help = "任务文件，一行一个任务，直接从提交的 http://10.0.22.63:8888/soft/probeDesign.action 页面复制文本即可")
    parser.add_argument('--output_dir', '-o', type = str, required=True, help = "结果输出目录, 例如 ./probe_file")
    args = parser.parse_args()

    return args

def read_input(file):
    tixi = {}
    with open(file, 'r') as fh:
        for line in fh:
            tixi_name, id_value = line.strip().split(':')
            tixi_name = tixi_name.replace('体系', "")
            tixi[tixi_name] = id_value
    return tixi

def download(tixi_id, output_file, sleep_time = 0):
    driver.get(f'http://10.0.22.63:8888/soft/query.action?id={tixi_id}')
    # 一定要等待
    time.sleep(sleep_time)
    # 弹窗accept
    driver.switch_to.alert.accept()
    # 获取文件下载地址
    address = driver.find_elements(By.TAG_NAME, 'a')
    address1 = address[0].get_property('href')
    # 本地保存路径
    
    # 下载
    print(f"    wget {address1} -O {output_file}")
    os.system(f"wget {address1} -O {output_file} -q")

args = set_and_parse_args()
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# 读入体系信息
tixi = read_input(args.input)

# 配置chrome浏览器的参数
options = webdriver.FirefoxOptions()
options.add_argument('--headless')  # 确保无头
options.add_argument('--disable-gpu')  # 无需要gpu加速
options.add_argument('--no-sandbox')  # 无沙箱
options.add_argument('window-size=1920x1080') # 4. 设置浏览器分辨率

print("初始化浏览器")
driver = webdriver.Firefox(options=options)

driver.implicitly_wait(30)

print("开始下载")

total_count = len(list(tixi.keys()))
count = 0
for tixi_name in tixi.keys():
    count += 1
    print(f'[process] {count} / {total_count}  {tixi_name}')

    
    retry = True
    try_round = 0
    sleep_time = 0.5 # 等待时间
    sleep_step = 0.5  # 重试时，增加等待时间
    while(retry):
        try:
            output_file = os.path.join(args.output_dir, tixi_name + ".xlsx")
            download(tixi[tixi_name], output_file, sleep_time)
            retry = False
        except Exception as err:
            try_round += 1
            sleep_time += sleep_step  # 增加等待时间
            print(f'异常，重试第 {try_round} 次')


driver.quit()
