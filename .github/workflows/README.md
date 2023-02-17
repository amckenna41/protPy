# Workflows used in protpy

* `build_test.yml` - build and test the protpy application, running all unit tests and uploading all test artifacts to workflow on GitHub.
* `deploy_test_pypi.yml` - after test workflow successful, deploy to test pypi server.
* `deploy_pypi.yml` - after deployment to test pypi server successful, deploy to pypi server.