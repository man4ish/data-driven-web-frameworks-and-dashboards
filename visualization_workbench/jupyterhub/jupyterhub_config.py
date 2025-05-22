c = get_config()

c.JupyterHub.bind_url = 'http://:8000'
c.Authenticator.admin_users = {'admin'}
c.Authenticator.allowed_users = {'admin', 'user1'}
c.LocalAuthenticator.create_system_users = True

c.JupyterHub.spawner_class = 'dockerspawner.DockerSpawner'
c.DockerSpawner.image = 'jupyter/scipy-notebook:latest'
c.DockerSpawner.network_name = 'visualization_workbench_default'
c.DockerSpawner.remove = True
c.DockerSpawner.volumes = {'jupyterhub-user-{username}': '/home/jovyan/work'}

