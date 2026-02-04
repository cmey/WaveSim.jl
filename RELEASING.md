# Release Process for WaveSim.jl

This document describes how to release a new version of WaveSim.jl to the Julia General registry.

## Overview

WaveSim.jl uses automated tooling for releases:
- **Registrator**: Handles registration of new versions to Julia's General registry
- **TagBot**: Automatically creates GitHub releases when a new version is registered

## Steps to Release a New Version

### 1. Update the Version Number

Edit `Project.toml` and update the version field following [semantic versioning](https://semver.org/):
- **Patch release** (0.1.x): Bug fixes and minor changes
- **Minor release** (0.x.0): New features, backwards compatible
- **Major release** (x.0.0): Breaking changes

Example:
```toml
version = "0.1.3"
```

Commit this change to the master branch.

### 2. Register the New Version

You can register a new version using one of these methods:

#### Option A: Comment on GitHub (Recommended)

1. Go to the commit or pull request that includes the version bump
2. Comment: `@JuliaRegistrator register`
3. Registrator will respond and create a pull request to the General registry

#### Option B: Registrator Web Interface

1. Visit https://juliahub.com/ui/Registrator
2. Enter the repository URL: `https://github.com/cmey/WaveSim.jl`
3. Follow the prompts

### 3. Wait for Registry PR to Merge

- Registrator creates a PR in the [General registry](https://github.com/JuliaRegistries/General)
- For patch/minor releases, this is usually auto-merged within minutes
- You can track the PR by following the link in Registrator's response

### 4. TagBot Creates the Release

Once the registry PR is merged:
1. JuliaTagBot will comment on an issue in this repository
2. This triggers the TagBot GitHub Actions workflow
3. TagBot creates:
   - A git tag (e.g., `v0.1.3`)
   - A GitHub release with auto-generated changelog

You don't need to do anything - this happens automatically!

## Verifying the Release

After TagBot completes:
- Check the [releases page](https://github.com/cmey/WaveSim.jl/releases) for the new release
- Verify the release notes look correct
- Test installation: `using Pkg; Pkg.add("WaveSim")`

## TagBot Configuration

The TagBot workflow is configured in `.github/workflows/TagBot.yml`. It:
- Runs when JuliaTagBot comments on an issue
- Can be manually triggered via workflow_dispatch
- Uses the latest version of TagBot (v1)

## Troubleshooting

### Release not created after registry merge

1. Check if TagBot workflow ran in the [Actions tab](https://github.com/cmey/WaveSim.jl/actions/workflows/TagBot.yml)
2. Manually trigger TagBot via the Actions tab → TagBot → Run workflow
3. Check TagBot logs for errors

### Registrator fails

- Ensure the version number is higher than all existing releases
- Verify `Project.toml` is valid TOML syntax
- Check that all dependencies have valid version constraints

## More Information

- [Registrator.jl](https://github.com/JuliaRegistries/Registrator.jl)
- [TagBot](https://github.com/JuliaRegistries/TagBot)
- [General Registry](https://github.com/JuliaRegistries/General)
