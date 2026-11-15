# Put the complete validation package on GitHub

Use GitHub Desktop or Git, rather than the browser uploader, for this complete repository. There is no need to extract, split or remove the nested scientific proof ZIP.

## Why the one large ZIP is included

`validation/ramp/r4_certificate/baseline/PAPER1_COMPLETE_DOMAIN_VALIDATION_R3.zip` is an input to the retained proof chain. It is **53,490,635 bytes** (about 51.0 MiB). GitHub's browser upload limit is 25 MiB per file; ordinary Git accepts files below its 100 MiB blocking threshold, with a warning possible above 50 MiB. Every file in this release is below 100 MiB. These files do not require Git LFS.

Keep this ZIP at its original internal path and with its original bytes. Deleting it, splitting it or uploading a Git LFS pointer in its place would make the checked release incomplete. See [GitHub's file-size documentation](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).

## Prepare a local repository

First extract into a short location and run the basic check in [README.md](README.md). The following steps put the validation files at the root of a dedicated Paper 1 validation repository.

1. Open GitHub Desktop. For an existing GitHub repository, choose **File → Clone repository**, select it or enter its URL, and choose a short local location, such as `C:\p1\repo`. If it is already cloned, select that checkout. To create a new repository instead, choose **File → New repository**, use a short name and local path, and do not generate another README, ignore file or license over this package's files.
2. In File Explorer, open the extracted `p1val` folder and copy **its contents** into the intended local checkout. Include `.github`, `.gitattributes` and `.gitignore`. Do not copy the outer `p1val` folder itself when installing at repository root. Keep the checkout's existing `.git` directory.
3. Open a terminal in that checkout and run `py -3.12 -B verify.py` (Windows), or `python3.12 -B verify.py` (Linux/macOS). Resolve any failure before uploading.

If the checkout already contains unrelated work, do not replace it indiscriminately. The complete `p1val` folder can instead be placed in a dedicated subfolder; run checks from that subfolder. Internal paths still work. The bundled GitHub Actions workflow is discovered automatically only when `.github/workflows/` is at the Git repository root. A subfolder installation needs a separate, deliberately adapted workflow if remote Actions checks are wanted.

GitHub's [cloning guide](https://docs.github.com/en/desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop) and [new-repository guide](https://docs.github.com/en/desktop/overview/creating-your-first-repository-using-github-desktop) explain the Desktop controls. Creating or committing locally does not publish the work.

## Check exactly what Git will upload

The ordinary check verifies the working files. The additional `--git-index` option checks that **every release file is present in Git's index with exactly the same bytes**. It needs the `git` command available in the terminal. Check with `git --version`; GitHub Desktop's installation does not by itself establish that this terminal can find `git`. If necessary, install the command-line Git described in [GitHub's setup instructions](https://docs.github.com/en/get-started/git-basics/set-up-git), then reopen the terminal.

Choose either route:

### GitHub Desktop route

1. In **Changes**, include all package files. Do not deselect the large nested ZIP, hidden dotfiles or retained history. Commit the files locally with a clear summary, such as `Add Paper 1 validation package, 15 November 2026`. Do not push yet. A local commit adds the selected files to the index.
2. From the validation root in the local checkout, run:

   ```powershell
   py -3.12 -B verify.py --git-index
   ```

   On Linux/macOS, use `python3.12` in place of `py -3.12`. A successful check includes `PASS. All ... release files have identical bytes in the Git index.` The number comes from the current release inventory.
3. If the check passes and these are the files you intend to share, use **Push origin** for an existing remote, or **Publish repository** for a new one. Choose the intended visibility when publishing. If the check fails, fix the cause and commit the corrected files before checking again.

Checkboxes in Desktop's Changes view are not a substitute for the index check. This procedure checks after the local commit so the selected content is actually in the index.

### Command-line Git route

In the validation root of the intended checkout, use:

```sh
git --version
git add .
```

Then run the platform's basic check **with `--git-index`**, as above. Inspect `git status` and the staged changes before committing and pushing. For an update, review any retained older files separately; the release manifest checks this release's required files, rather than deleting other repository content for you.

The supplied `.gitattributes` disables text conversion so that Git preserves the archived file bytes. Keep it in place. A missing-file index error usually means a file was omitted or ignored. A staged-byte error means Git's staged copy differs from the verified working copy; stage/commit the correct release bytes, then rerun the check.

## After pushing

Open the repository on GitHub and check that the expected root files and the nested proof ZIP are visible. For a stronger transfer check, clone the pushed repository to a new short location and run `verify.py --git-index` there. A successful local check before pushing does not itself verify the remote transfer.

The included Actions workflow is manual. When installed at repository root and enabled on GitHub, it can run the basic check or the optional compact suite. A completed local check does not mean an Actions run has occurred. This documentation update did not upload to your GitHub account or execute remote Actions.
