from django.db import models

class UploadedDataset(models.Model):
    name = models.CharField(max_length=255)
    file = models.FileField(upload_to='data/single_cell/')
    uploaded_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name
