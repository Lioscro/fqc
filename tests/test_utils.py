import os
import tempfile
import uuid
from unittest import TestCase

import fqc.utils as utils


class TestUtils(TestCase):

    def test_open_as_text_textfile(self):
        path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        with utils.open_as_text(path, 'w') as f:
            f.write('TESTING')
        self.assertTrue(os.path.exists(path))
        with utils.open_as_text(path, 'r') as f:
            self.assertEqual(f.read(), 'TESTING')

    def test_open_as_text_gzip(self):
        path = os.path.join(tempfile.gettempdir(), '{}.gz'.format(uuid.uuid4()))
        with utils.open_as_text(path, 'w') as f:
            f.write('TESTING')
        self.assertTrue(os.path.exists(path))
        with utils.open_as_text(path, 'r') as f:
            self.assertEqual(f.read(), 'TESTING')
