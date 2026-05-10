"""Python translation of analysis/classes/Morpheus.m.

Morpheus sends status/error email notifications for jobs/job arrays.  This Python
version uses smtplib when SMTP settings are provided, or the local sendmail
binary if available.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from email.message import EmailMessage
from typing import Any, List, Optional
import smtplib
import subprocess
import shutil


@dataclass
class Morpheus:
    recipient: str
    verbose: bool = True
    name: str = ""
    maxNumErrors: float = float("inf")
    sender: str = "morpheus@localhost"
    smtpHost: Optional[str] = None
    smtpPort: int = 25
    smtpUser: Optional[str] = None
    smtpPassword: Optional[str] = None
    useTLS: bool = False
    numErrors: int = 0
    successListeners: List[Any] = field(default_factory=list)
    errorListeners: List[Any] = field(default_factory=list)

    def SendMessage(self, subject: str, body: str | list[str], sourceObj: Any = None) -> bool:
        lines = body if isinstance(body, list) else str(body).splitlines()
        if self.name:
            lines = [f"Message from morpheus: {self.name}"] + list(lines)
        if sourceObj is not None:
            lines.append(f"source: {sourceObj.__class__.__name__}")
        text = "\n".join(str(x) for x in lines)

        msg = EmailMessage()
        msg["Subject"] = str(subject)
        msg["From"] = self.sender
        msg["To"] = self.recipient
        msg.set_content(text)

        if self.smtpHost:
            with smtplib.SMTP(self.smtpHost, self.smtpPort, timeout=30) as smtp:
                if self.useTLS:
                    smtp.starttls()
                if self.smtpUser:
                    smtp.login(self.smtpUser, self.smtpPassword or "")
                smtp.send_message(msg)
            return True

        sendmail = shutil.which("sendmail")
        if sendmail:
            proc = subprocess.run([sendmail, "-t"], input=msg.as_string(), text=True, capture_output=True)
            return proc.returncode == 0

        if self.verbose:
            print(f"To: {self.recipient}\nSubject: {subject}\n{text}")
        return True

    send_message = SendMessage

    def AddSuccessListener(self, sourceObj: Any, eventName: str) -> None:
        self.successListeners.append((sourceObj, eventName))

    add_success_listener = AddSuccessListener

    def AddErrorListener(self, sourceObj: Any, eventName: str) -> None:
        self.errorListeners.append((sourceObj, eventName))

    add_error_listener = AddErrorListener

    def HandleSuccessMessageRequest(self, sourceObj: Any, event: Any = None) -> bool:
        if hasattr(sourceObj, "Status"):
            header, body = sourceObj.Status()
        else:
            header, body = "Success", str(sourceObj)
        return self.SendMessage(header, body, sourceObj)

    handle_success_message_request = HandleSuccessMessageRequest

    def HandleErrorMessageRequest(self, sourceObj: Any, event: Any = None) -> bool:
        self.numErrors += 1
        if hasattr(sourceObj, "Status"):
            header, body = sourceObj.Status()
        else:
            header, body = "Error", str(sourceObj)
        if self.numErrors <= self.maxNumErrors:
            return self.SendMessage(header, body, sourceObj)
        return False

    handle_error_message_request = HandleErrorMessageRequest
