## Cristian
- A simple control rule based on MSY to recovery the stock (finished). Here I believe with Niels we'll have interesting overview points.
- Optimal allocation of catches by fleet (finished)(Operative Research)
- Impact of use different wt-at age (Peru and Chile) in the assessment and simulation analysis (finished)

Also, I'm preparing key information related to:  

- wt-at age by Chilean Fleet
- wt-at age for acoustic survey (F2)
- age comps of Chilean catches by fleet
- CPUE updated to 2013 (F2)
- length comps of acoutic survey (F1)

I expect to have all information at the end of September, while the documents have been send to text reviewing (translation).
In regards to the model, my main concern are related to the information and processes we are modeling: 

- How much informative are the length comps in Peru-Ecuador (what we are seeing there?). Availability vs abundance? Is there consistency between these comps and abundance indicators? 
- The growth between areas is too different, and we should do something about that (impact over wt-at age)
- Base model: natural mortality can not be the simple average, because the population size and it fishery off Chile is 10 times the Peruvian stock size.
- We should exclude non informative data, for example, Chilean acoustic survey in F1, Russian CPUE. Or at least to say something in regards to its quality for assessment purposes.
- Model 7 seems to be the best if more changes in selectivity are considered.
- We should establish the jack mackerel condition in explicit terms (overfished?), this in order to implement some recovery strategy. In this sense, what about the BRP that were shown last meeting?

## Niels
The work I’m doing involves:

- Management Strategy Evaluation to indicate recovery under different Harvest Control Rule designs
- Designing a template for advice (a 2-pager which mentions the core elements the commission needs to know to be able to adopt the SC advice)
- if necessary, update my R code to be able to read in all .dat files and output files from ADMB correctly

My wish list for the assessment would be:

- model 7, with extensions in selectivity patterns, especially in the later years where there is only 1 or 2 cohorts to be tracked
- model 6, with extensions in selectivity patterns, especially in the later years where there is only 1 / 2 cohorts to be tracked (so a blend of model 4 and 7)
- model 6/7 with further down weighting of those surveys that have not been standardized
- is there a way to update fisheries-wt-at-age for the 4 fleets?
- is there a way to update index-wt-at-age for the surveys?
- is there a way to make population wt-at-age variable over the years (e.g. informed by a combination of catches and surveys)

That would it be for me.
But let me ask you the same question: what would you like from me?
