import os
import smtplib
import imghdr
from email.message import EmailMessage


def mail_function(subject,content):
    """
    this function will send a simpel mail to sacharens@gmail.com
    with a subject and short conotent
    :param subject: string
    :param content: string
    """
    EMAIL_ADDRESS = 'sinaipymail@gmail.com'
    EMAIL_PASSWORD = 'ftmzztowvqlwrwoq'

    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = 'sacharens@gmail.com'

    msg.set_content(content)

    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
        smtp.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
        smtp.send_message(msg)


if __name__ == "__main__":
    mail_function()

