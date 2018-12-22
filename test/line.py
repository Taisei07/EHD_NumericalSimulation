# -*- coding: utf-8 -*-

import json
import requests
sys.path.append('../../')
from module.LineAPI import line_notify_token


def LineNotify(A):
    line_notify_api = 'https://notify-api.line.me/api/notify'
    headers = {'Authorization': 'Bearer ' + line_notify_token}
    # メッセージ
    payload = {'message': str(A)}
    files = {"imageFile": open(A.png)}
    requests.post(line_notify_api, data=payload, headers=headers, files=files)

LineNotify()
