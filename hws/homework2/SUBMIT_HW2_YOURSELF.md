# Submit HW2 (GOOD_for_submission) yourself

Run these in **your terminal** (not Cursor). The good submission is already in the repo; this makes sure it’s tagged and pushed under your account.

## 1. Go to the repo

```bash
cd /Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232
```

## 2. If GOOD_for_submission is missing on disk (restore it from the good commit)

```bash
git checkout 67141e8 -- hws/homework2/GOOD_for_submission/
```

Then add and commit (you as author):

```bash
git add hws/homework2/GOOD_for_submission/
git commit -m "Submit HW2: GOOD_for_submission (solutions PDF, bvp, plots, run_q4.sh)"
```

## 3. Tag the submission and push

Tag the commit that contains GOOD_for_submission (use the current tip if you just committed in step 2, or the existing commit if you skipped step 2):

```bash
# If you did step 2, tag the commit you just made:
git tag -f -a hw2-submission -m "APMA 4302 HW2 submission (GOOD_for_submission)"

# If you did NOT do step 2, tag the existing good commit:
# git tag -f -a hw2-submission 67141e8 -m "APMA 4302 HW2 submission (GOOD_for_submission)"
```

Push your branch and the tag:

```bash
git push origin main
git push origin hw2-submission --force
```

Done. The submission is under your name on GitHub at tag `hw2-submission`.
