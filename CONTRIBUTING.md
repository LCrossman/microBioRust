microBioRust is an open source project to create efficient bioinformatics tools in Rust for sustainable computing.  The current implementation offers a range of parsers for bioinformatics file formats, especially aimed at microbial genomics.  Microbial genomics is a rapidly evolving field, where we are often dealing with many thousands of short files.  In addition, we are offering total InterOp with Python (installable with pip) and soon-to-be R support.  

## First steps:
1. **Fork the repository:** 
   Please click the "Fork" button at the top right of the https://github.com/LCrossman/microBioRust GitHub page to create your own copy of microBioRust.
2. **Clone your fork locally:**
   ```bash
   git clone https://github.com/YOUR-USERNAME/microBioRust.git
   cd microBioRust
   ```
3. **Create a branch:**
   ```bash
   git checkout -b your-branch-name
   ```
   Set up your environment:

   Rust: Do ensure you have Rust installed and run cargo build.
   Python Bindings: Set up a virtual environment and pip install -r requirements.txt.

4. **Issues**: if you are planning big changes, please open an issue first to suggest your changes and mark it as ongoing.
   **Make your changes:** Write your code, add unit tests if applicable, and ensure everything runs smoothly (cargo test).

5. **Push and PR:** Push the branch to your fork, then navigate back to the main microBioRust repository to open your Pull Request.

## Ways to Contribute

We need lots of skills and expertise to improve microBioRust and keep current. You don't need to be a Rust expert or a bioinformatician! Here are the main areas we need help with:

* **Rust Core Engine:** Help to optimize functions, Don't-repeat-yourself (DRY), memory usage, fix bugs, add new parsers, or write tests.
* **Docs & Tutorials:** Docs are in a separate repo, [microBioRust-docs](https://github.com/microBioRust/microBioRust-docs). If you find a typo, think a function needs a better explanation, or want to write a tutorial on how you use the tool, please send a PR!
* **Biology/Bioinformatics:** We do need real-world edge cases. If you have an unusual genome file, such as an EMBL CONTIG JOIN, an opaque error we can clarify, or an idea for a new biological metric we should calculate, please open an issue
* **Downstream InterOp (Python, R):** The core is Rust, but users work in Python and R. Contributing to the InterOp scripts, data wrangling wrappers, or statistical workflows is definitely encouraged.

### Finding Something to Work On

If you are looking for a place to start, please check out the GitHub Issues. There are some specific labels to help you find the right task, please comment if you'd like to take something on so you can be assigned to it!

* **good first issue**: ideal starting point for new joiners. These are well-defined, isolated tasks that don't require deep knowledge of the full codebase.  Seqmetrics is a good starting place because it is relatively small and each metric task is separate and self-contained.
* **help wanted**: It would be great if a community member could take this on because the core maintainers lack the bandwidth or specific expertise.
* **documentation**: Tasks related to improving the README, wikis, or inline code comments.
* **Python / R / Rust**: Feel free to contribute in any of these languages.

If you can just comment on the issue saying if you would be assigned to that issue then we all know someone is working on it.
Otherwise, feel free to create a new issue or pull request!

Code of Conduct: This project follows the [Rust Code of Conduct](https://www.rust-lang.org/policies/code-of-conduct). By participating in the microBioRust community, we expect to uphold these standards and treat all contributors with respect

FAQ
Q: Should I contribute to https://github.com/LCrossman/microBioRust or https://github.com/microBioRust/microBioRust?
A: Please contribute to the original LCrossman/microBioRust repository. While the organisation repo (microBioRust/microBioRust) is the source for crates.io and both are kept in sync, it is technically a fork. Submitting PRs to the LCrossman repo keeps the Git history tidier.  This also applies to the docs folder, https://github.com/LCrossman/microBioRust-docs and https://github.com/microBioRust/docs.  If you're unsure, the simplest rule is: look at which repo has the most recent commit activity and submit there.
