from django.db import models
from django.contrib.auth.models import User

class Job(models.Model):
    STATUS_CHOICES = [
        ('PENDING', 'Pending'),
        ('RUNNING', 'Running'),
        ('SUCCESS', 'Success'),
        ('FAILURE', 'Failure'),
    ]
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    input_file = models.FileField(upload_to='uploads/')
    reference_genome = models.CharField(max_length=100)
    status = models.CharField(max_length=10, choices=STATUS_CHOICES, default='PENDING')
    created_at = models.DateTimeField(auto_now_add=True)
    log = models.TextField(blank=True, default='')

    def __str__(self):
        return f"Job {self.id} - {self.status}"
