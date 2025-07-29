# Contributing to Federated Imputation Pipeline

## GA4GH Hackathon 2025 - African Genomics Team

Thank you for your interest in contributing to the Federated Genotype Imputation Pipeline! This guide explains our development workflow, testing requirements, and CI/CD processes.

## Table of Contents

- [Development Workflow](#development-workflow)
- [CI/CD Pipeline](#cicd-pipeline)
- [Testing Requirements](#testing-requirements)
- [Pull Request Process](#pull-request-process)
- [Code Quality Standards](#code-quality-standards)
- [Security Guidelines](#security-guidelines)

## Development Workflow

### Branch Strategy
- **`main`**: Production-ready code with full CI/CD validation
- **`develop`**: Integration branch for feature development
- **`feature/*`**: Feature branches for new functionality
- **`hotfix/*`**: Critical bug fixes for production

### Getting Started
```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/prepare_minimac4_data.git
cd prepare_minimac4_data

# Create a feature branch
git checkout -b feature/your-feature-name

# Make your changes
# ... edit files ...

# Test locally (see Testing section below)
cd tests/
./run_unit_tests.sh

# Commit and push
git add .
git commit -m "feat: add your feature description"
git push origin feature/your-feature-name
```

## CI/CD Pipeline

Our repository uses GitHub Actions for automated testing and quality assurance. The pipeline consists of three main workflows:

### 1. Unit Tests Workflow
**File**: `.github/workflows/unit-tests.yml`
**Triggers**: Push/PR to `main`/`develop`

```yaml
Workflow: 🧬 Federated Imputation Pipeline - Unit Tests
├── Test Matrix (4 parallel jobs)
│   ├── ExtractRegion unit test
│   ├── QualityControl unit test  
│   ├── MinimacConversion unit test
│   └── Complete pipeline integration test
└── Test Summary & Reporting
```

**What it tests**:
- Individual WDL task functionality
- Input/output validation
- Error handling and edge cases
- Cross-platform Docker compatibility

### 2. Integration Tests Workflow
**File**: `.github/workflows/integration-tests.yml`
**Triggers**: Push to `main`, daily schedule, manual dispatch

```yaml
Workflow: 🚀 Federated Imputation Pipeline - Integration Tests
├── Matrix Strategy (2 configurations)
│   ├── Main pipeline execution
│   └── Enhanced test pipeline with validation
└── End-to-end validation
```

**What it tests**:
- Complete pipeline execution
- MSAV file generation and validation
- Performance benchmarks
- Output quality verification

### 3. Security Scan Workflow
**File**: `.github/workflows/security-scan.yml`
**Triggers**: Weekly schedule, manual dispatch

```yaml
Workflow: 🔒 Security Scan - Docker & Dependencies
├── Docker Image Security Scan
│   ├── Trivy vulnerability scanning
│   ├── Base image analysis
│   └── Dependency checks
└── Repository Security Analysis
    ├── Secret detection
    ├── WDL security patterns
    └── Security recommendations
```

## Testing Requirements

### Before Submitting a PR

#### 1. Local Unit Tests
```bash
# Run all unit tests
cd tests/
chmod +x run_unit_tests.sh
./run_unit_tests.sh

# Expected output:
# 🎯 Unit Test Suite Complete
# Results Summary:
#   Total Tests: 4
#   Passed: 4 ✅
#   Failed: 0 ❌
#   Success Rate: 100.0%
```

#### 2. Individual Component Testing
```bash
# Test specific components
java -jar cromwell.jar run tests/unit/test_extract_region.wdl -i tests/inputs/unit_test_config.json
java -jar cromwell.jar run tests/unit/test_quality_control.wdl -i tests/inputs/qc_unit_test_config.json
java -jar cromwell.jar run tests/unit/test_minimac_conversion.wdl -i tests/inputs/minimac_unit_test_config.json
```

#### 3. Integration Testing
```bash
# Test complete pipeline
java -jar cromwell.jar run workflows/federated_imputation_pipeline.wdl -i inputs/test_local.json
java -jar cromwell.jar run tests/test_pipeline.wdl -i inputs/test_local.json
```

### Test Coverage Requirements

All contributions must maintain or improve test coverage:

- ✅ **Unit Tests**: Every WDL task must have a corresponding unit test
- ✅ **Input Validation**: Test edge cases and invalid inputs
- ✅ **Output Verification**: Validate file formats, sizes, and content
- ✅ **Error Handling**: Test failure modes and error recovery
- ✅ **Documentation**: Update test documentation for new features

### Adding New Tests

When adding new WDL tasks or workflows:

1. **Create Unit Test**: `tests/unit/test_your_task.wdl`
2. **Add Input Config**: `tests/inputs/your_task_config.json`
3. **Update Test Runner**: Add to `tests/run_unit_tests.sh`
4. **Document**: Update `tests/README.md`

Example unit test structure:
```wdl
version 1.0
import "../../tasks/your_task.wdl" as YourTask

workflow TestYourTask {
    input {
        # Test inputs
    }
    
    call YourTask.YourTask { /* ... */ }
    call ValidateYourTask { /* ... */ }
    
    output {
        File validation_report = ValidateYourTask.validation_report
        Boolean test_passed = ValidateYourTask.test_passed
    }
}

task ValidateYourTask {
    # Validation logic
}
```

## Pull Request Process

### 1. Pre-PR Checklist
- [ ] All unit tests pass locally
- [ ] Integration tests pass locally
- [ ] Code follows WDL best practices
- [ ] Documentation updated (if applicable)
- [ ] No security vulnerabilities introduced

### 2. PR Requirements
- **Title**: Use conventional commits format (`feat:`, `fix:`, `docs:`, etc.)
- **Description**: Clear explanation of changes and motivation
- **Testing**: Include test results or explain testing approach
- **Breaking Changes**: Clearly document any breaking changes

### 3. Automated Checks
Your PR will automatically trigger:

- ✅ **Unit Tests**: All 4 unit tests must pass
- ✅ **Integration Tests**: Full pipeline validation
- ✅ **Security Scan**: Dependency and container security checks
- ✅ **Code Quality**: WDL syntax validation

### 4. Review Process
1. **Automated Review**: CI/CD pipeline validates all changes
2. **Peer Review**: At least one maintainer reviews the code
3. **Test Validation**: Review test results and artifacts
4. **Merge**: Approved PRs are merged to target branch

## Code Quality Standards

### WDL Best Practices
- Use WDL version 1.0
- Include comprehensive `meta` and `parameter_meta` blocks
- Specify resource requirements (`memory`, `cpu`, `disks`)
- Use appropriate Docker containers with version tags
- Include platform constraints for reproducibility

### Documentation Standards
- Clear task and workflow descriptions
- Complete parameter documentation
- Usage examples and expected outputs
- Error handling documentation

### Security Standards
- No hardcoded secrets or credentials
- Use secure base images
- Minimal container privileges
- Input validation and sanitization

## Security Guidelines

### Container Security
- Use specific version tags (not `latest`)
- Regularly update base images
- Scan for vulnerabilities weekly
- Follow principle of least privilege

### Data Security
- Maintain data sovereignty principles
- No data leakage outside containers
- Secure file handling practices
- Audit trail for all operations

### WDL Security
- Validate input parameters
- Use type-safe file inputs
- Avoid privileged container execution
- Sanitize user-provided values

## Issue Reporting

### Bug Reports
Use the bug report template and include:
- Environment details (local/cloud)
- Cromwell version and configuration
- Complete error logs and stack traces
- Steps to reproduce

### Security Issues
For security vulnerabilities:
1. **DO NOT** create public issues
2. Email maintainers directly
3. Allow 90 days for coordinated disclosure
4. Provide detailed vulnerability description

### Feature Requests
- Clear use case description
- Expected behavior and outputs
- Compatibility considerations
- Testing approach

## Development Tips

### Local Development Environment
```bash
# Required tools
- Docker (for container execution)
- Java 11+ (for Cromwell)
- WDL syntax checker
- Git (version control)

# Recommended tools
- VSCode with WDL extension
- Docker Desktop
- jq (for JSON processing)
```

### Testing Best Practices
- Test on different platforms (if possible)
- Use realistic test data sizes
- Validate all output files
- Test error conditions
- Document expected results

### Performance Considerations
- Monitor resource usage in tests
- Optimize Docker image layers
- Use appropriate computing resources
- Consider scalability implications

## Support

### Getting Help
- **Documentation**: Check `README.md` and `tests/README.md`
- **Issues**: Search existing issues before creating new ones
- **Discussions**: Use GitHub Discussions for questions
- **Direct Contact**: Maintainer email for urgent issues

### Maintainers
- **Primary**: GA4GH Hackathon 2025 Team
- **Email**: team@afrigenomics.org
- **Response Time**: 48-72 hours for issues
- **Availability**: UTC business hours primarily

---

## Ready to Contribute!

Thank you for contributing to the Federated Genotype Imputation Pipeline! Your contributions help advance African genomics research while maintaining the highest standards of quality and security.

**Remember**: Every contribution, no matter how small, makes a difference in advancing federated genomics research across African institutions! 