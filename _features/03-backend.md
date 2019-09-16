---
id: backend
name: Backend
heading: ReactomeGSA Backend<br>
subheading: A kubernetes application
image: "/assets/images/fig_backend.png"
---

The complete backend of the Reactome Gene Set Analysis service is developed
as one kubernetes application. It is designed so that it can easily be
deployed on any kubernetes cluster. Therefore, for cases where the public
analysis service cannot be used, it is straight forward to deploy one's own
analysis service. All details, as well as the source code are available at the
[project site](https://github.com/reactome/gsa-backend).
