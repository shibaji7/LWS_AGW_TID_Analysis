import dropbox
from dropbox import DropboxOAuth2FlowNoRedirect

'''
This example walks through a basic oauth flow using the existing long-lived token type
Populate your app key and app secret in order to run this locally
'''
APP_KEY = "m5y6ql43dh80v0y"
APP_SECRET = "rs9eof2db5x53sx"
auth_code = "sl.BruRWiKd4vf9ITwE4UH7kD9YMyPuWh2gpBdGt1kfSg-6gU8oeF2L30Bki-SiwX8WRrUf9WowHQdLSNAUNUV2Je84VpZrH8ZXK-D7BSmaCOUO3Q7EnHq3ms-CtGF2TvUZBw4Re4lzrpYf"

auth_flow = DropboxOAuth2FlowNoRedirect(APP_KEY, APP_SECRET)

authorize_url = auth_flow.start()


try:
    oauth_result = auth_flow.finish(auth_code)
except Exception as e:
    print('Error: %s' % (e,))
    exit(1)

with dropbox.Dropbox(oauth2_access_token=oauth_result.access_token) as dbx:
    dbx.users_get_current_account()
    print("Successfully set up client!")