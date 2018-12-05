# -*- coding: utf-8 -*-
import requests
import json
import sys
value = sys.argv

for num in [4, 3, 12]:
    print "num =" + str(num)
print "End"

requests.post('https://hooks.slack.com/services/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', data = json.dumps({
    'text': u'Test', # 投稿するテキスト
    'username': u'me', # 投稿のユーザー名
    'icon_emoji': u':ghost:', # 投稿のプロフィール画像に入れる絵文字
    'link_names': 1, # メンションを有効にする
}))
