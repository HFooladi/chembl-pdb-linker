# Contributing to ChEMBL-PDB Linker

Thank you for your interest in contributing to ChEMBL-PDB Linker! This document provides guidelines for contributing to the project.

## Development Setup

### Prerequisites

- Python 3.9 or higher
- [uv](https://github.com/astral-sh/uv) (recommended) or pip

### Installation

```bash
# Clone the repository
git clone https://github.com/HFooladi/chembl-pdb-linker.git
cd chembl-pdb-linker

# Create virtual environment and install dependencies
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

## Code Quality

This project uses several tools to maintain code quality:

### Linting and Formatting

We use [Ruff](https://github.com/astral-sh/ruff) for both linting and formatting:

```bash
# Check for linting issues
ruff check src/ tests/

# Auto-fix linting issues
ruff check --fix src/ tests/

# Check formatting
ruff format --check src/ tests/

# Apply formatting
ruff format src/ tests/
```

### Type Checking

We use [Pyright](https://github.com/microsoft/pyright) for static type checking:

```bash
pyright
```

### Running Tests

We use [pytest](https://pytest.org/) for testing:

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src/chembl_pdb_linker

# Run only unit tests
pytest tests/unit/

# Run only integration tests
pytest tests/integration/ -m integration

# Skip slow tests
pytest -m "not slow"
```

## Pull Request Process

1. **Fork the repository** and create a new branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the code style guidelines above.

3. **Run quality checks** before committing:
   ```bash
   ruff check src/ tests/
   ruff format src/ tests/
   pyright
   pytest
   ```

4. **Commit your changes** with a clear, descriptive message:
   ```bash
   git commit -m "feat: Add new feature description"
   ```

   We follow [Conventional Commits](https://www.conventionalcommits.org/) style:
   - `feat:` - New features
   - `fix:` - Bug fixes
   - `docs:` - Documentation changes
   - `test:` - Test additions or modifications
   - `refactor:` - Code refactoring
   - `chore:` - Maintenance tasks

5. **Push your branch** and open a Pull Request:
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Fill out the PR template** with:
   - Description of changes
   - Related issues (if any)
   - Checklist confirmation

## Reporting Issues

When reporting issues, please include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected vs actual behavior
- Python version and operating system
- Any relevant error messages or logs

## Code of Conduct

Please be respectful and constructive in all interactions. We welcome contributors of all experience levels.

## Questions?

If you have questions, feel free to:
- Open a [GitHub Issue](https://github.com/HFooladi/chembl-pdb-linker/issues)
- Check existing issues for similar questions
