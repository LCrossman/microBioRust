# Changelog #
---

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project should adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## Unreleased ##
---

### Added ###
- docs
- windows install description for docs
- this changelog

### Changed ###
- We had two types of Records (gbk and embl) after parsing, they have been converted to one type of records after parsing - several structures and functions have been moved from gbk.rs and embl.rs into record.rs allowing a generic records type to be used for both
- Heatmap path fix
- Moved images folder to assets in docs windows install section

### Removed ###


