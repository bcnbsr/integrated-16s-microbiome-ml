# Licensing and redistribution guidance

## Current status

This repository is not being presented as fully open-source yet. The two source repositories used by the workflow are public, but neither currently contains a license file:

- [Luna-Maan/CSthesis](https://github.com/Luna-Maan/CSthesis)
- [R-joumaa/diet-microbiome-depression-ml](https://github.com/R-joumaa/diet-microbiome-depression-ml)

Public visibility is not the same as a license. In the absence of an explicit license, copyright normally remains with the authors and redistribution rights should not be assumed.

## Recommended path

1. Ask both upstream authors for written permission to redistribute the copied source files in this combined repository.
2. Ask whether they want attribution, a particular license, or specific citation text.
3. Keep the copied code in clearly labelled upstream directories and retain the upstream repository URLs and commit references.
4. Add a `NOTICE` file with the agreed attribution.
5. License the original integration code and documentation separately if the ownership is clear.

## Recommended license for original work

I recommend **Apache License 2.0** for original integration code if the author wants an explicit patent grant and a clear contribution framework. **MIT** is a good alternative when a short, permissive license is preferred.

Do not add Apache-2.0 or MIT to the entire repository until the upstream redistribution permission is resolved. A single top-level license can incorrectly imply that it covers code the repository owner does not own.

## Public-release checklist

- [ ] Written redistribution permission received from both upstream authors.
- [ ] Any third-party dependencies and classifiers checked for their own licenses.
- [ ] Raw reads and participant-level files excluded.
- [ ] Notebook outputs, metadata, logs, and archives checked for identifiers and credentials.
- [ ] `NOTICE.md` updated with exact attribution and commit references.
- [ ] License scope stated explicitly for original and upstream components.
- [ ] A clean clone can be installed and run using the README.

This is project guidance rather than legal advice. For a thesis, publication, or public release with institutional data, confirm the final wording with the relevant supervisor or legal/technology-transfer office.
