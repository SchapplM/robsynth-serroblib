% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:46:08
% EndTime: 2019-03-09 16:46:19
% DurationCPUTime: 7.04s
% Computational Cost: add. (6306->493), mult. (14857->600), div. (0->0), fcn. (10162->6), ass. (0->201)
t470 = cos(qJ(2));
t568 = cos(qJ(3));
t517 = t568 * t470;
t501 = qJD(1) * t517;
t467 = sin(qJ(3));
t468 = sin(qJ(2));
t534 = qJD(1) * t468;
t516 = t467 * t534;
t417 = -t501 + t516;
t462 = qJD(2) + qJD(3);
t466 = sin(qJ(5));
t431 = t467 * t470 + t468 * t568;
t393 = t462 * t431;
t475 = qJD(1) * t393;
t469 = cos(qJ(5));
t530 = qJD(5) * t469;
t531 = qJD(5) * t466;
t346 = -t417 * t530 + t462 * t531 - t466 * t475;
t394 = -t469 * t417 + t462 * t466;
t535 = qJD(1) * t431;
t578 = qJD(5) + t535;
t575 = t394 * t578;
t581 = t346 - t575;
t582 = t581 * MDP(30);
t508 = t578 ^ 2;
t474 = t462 * t535;
t569 = -pkin(8) - pkin(7);
t442 = t569 * t470;
t437 = qJD(1) * t442;
t421 = t467 * t437;
t441 = t569 * t468;
t435 = qJD(1) * t441;
t391 = t435 * t568 + t421;
t513 = t568 * qJD(3);
t580 = -pkin(2) * t513 + t391;
t579 = t535 * t469 + t530;
t523 = qJD(1) * qJD(2);
t577 = -0.2e1 * t523;
t576 = MDP(5) * (t468 ^ 2 - t470 ^ 2);
t550 = t467 * t468;
t430 = -t517 + t550;
t455 = -pkin(2) * t470 - pkin(1);
t495 = -qJ(4) * t431 + t455;
t570 = pkin(3) + pkin(9);
t371 = t430 * t570 + t495;
t398 = -t441 * t568 - t467 * t442;
t379 = t431 * pkin(4) + t398;
t540 = t469 * t371 + t466 * t379;
t424 = t568 * t437;
t390 = t467 * t435 - t424;
t533 = qJD(3) * t467;
t521 = pkin(2) * t533;
t574 = t390 - t521;
t539 = -qJD(4) + t580;
t399 = -t467 * t441 + t568 * t442;
t494 = -pkin(5) * t469 - qJ(6) * t466 - pkin(4);
t573 = -pkin(5) * t530 - qJ(6) * t531 + t469 * qJD(6) + t494 * t535 - qJD(4);
t572 = t466 * pkin(5) - qJ(6) * t469;
t562 = qJD(2) * pkin(2);
t426 = t435 + t562;
t388 = -t568 * t426 - t421;
t526 = qJD(4) + t388;
t571 = -MDP(27) - MDP(29);
t522 = MDP(28) - MDP(31);
t567 = pkin(2) * t467;
t499 = t462 * t550;
t538 = t462 * t501;
t385 = qJD(1) * t499 - t538;
t566 = pkin(5) * t385;
t565 = pkin(5) * t417;
t564 = t417 * pkin(4);
t563 = t535 * pkin(4);
t561 = qJ(6) * t385;
t512 = t468 * t523;
t449 = pkin(2) * t512;
t335 = pkin(3) * t474 + t385 * qJ(4) - qJD(4) * t535 + t449;
t312 = pkin(9) * t475 + t335;
t518 = qJD(2) * t569;
t502 = qJD(1) * t518;
t427 = t468 * t502;
t428 = t470 * t502;
t345 = t426 * t533 + t467 * t427 - t568 * t428 - t437 * t513;
t332 = -pkin(4) * t385 + t345;
t527 = t563 + t526;
t351 = -t462 * t570 + t527;
t440 = t455 * qJD(1);
t479 = -qJ(4) * t535 + t440;
t354 = t417 * t570 + t479;
t504 = t466 * t312 - t469 * t332 + t351 * t531 + t354 * t530;
t298 = t504 + t566;
t297 = t298 * t469;
t503 = t426 * t513 + t568 * t427 + t467 * t428 + t437 * t533;
t342 = -t462 * qJD(4) - t503;
t327 = -pkin(4) * t475 - t342;
t384 = t469 * t475;
t396 = t417 * t466 + t462 * t469;
t532 = qJD(5) * t396;
t347 = -t384 + t532;
t303 = t347 * pkin(5) + t346 * qJ(6) - t396 * qJD(6) + t327;
t559 = t303 * t469;
t316 = t351 * t466 + t354 * t469;
t307 = qJ(6) * t578 + t316;
t558 = t307 * t535;
t557 = t346 * t469;
t436 = t468 * t518;
t438 = t470 * t518;
t355 = -t568 * t436 - t467 * t438 - t441 * t513 - t442 * t533;
t556 = t355 * t462;
t356 = -qJD(3) * t399 + t467 * t436 - t568 * t438;
t555 = t356 * t462;
t389 = t467 * t426 - t424;
t554 = t389 * t462;
t552 = t430 * t466;
t551 = t466 * t385;
t472 = qJD(2) ^ 2;
t549 = t468 * t472;
t368 = t389 - t564;
t548 = t469 * t368;
t383 = t469 * t385;
t547 = t470 * t472;
t473 = qJD(1) ^ 2;
t546 = t470 * t473;
t545 = t570 * t385;
t544 = -t388 + t573;
t543 = t573 + t580;
t386 = pkin(3) * t535 + qJ(4) * t417;
t378 = pkin(2) * t534 + t386;
t413 = t535 * pkin(9);
t357 = t378 + t413;
t372 = t390 - t564;
t542 = t469 * t357 + t466 * t372;
t364 = t386 + t413;
t541 = t469 * t364 + t466 * t368;
t537 = t563 - t539;
t529 = qJD(5) * t570;
t528 = qJD(6) * t578;
t315 = t351 * t469 - t354 * t466;
t524 = qJD(6) - t315;
t458 = t468 * t562;
t454 = -pkin(2) * t568 - pkin(3);
t450 = -pkin(9) + t454;
t520 = t450 * t383;
t519 = t469 * t545;
t515 = t450 * t530;
t514 = t469 * t529;
t511 = pkin(1) * t577;
t484 = -t469 * t312 - t466 * t332 - t351 * t530 + t354 * t531;
t296 = -t484 + t528 - t561;
t305 = -pkin(5) * t578 + t524;
t510 = -t305 * t535 - t296;
t509 = -t316 * t417 + t327 * t469;
t460 = t462 * qJ(4);
t361 = t368 + t460;
t507 = t578 * t361;
t506 = t578 * t396;
t439 = qJ(4) + t572;
t498 = t305 * t469 - t307 * t466;
t497 = t305 * t466 + t307 * t469;
t377 = -pkin(3) * t462 + t526;
t381 = -t460 - t389;
t496 = t377 * t417 - t381 * t535;
t330 = pkin(5) * t394 - qJ(6) * t396 + t361;
t493 = t303 * t466 - t305 * t417 + t330 * t579;
t492 = t315 * t417 + t327 * t466 + t361 * t579;
t491 = t393 * t466 + t430 * t530;
t490 = -t393 * t469 + t430 * t531;
t392 = -qJD(2) * t517 - t470 * t513 + t499;
t489 = qJ(4) * t392 - qJD(4) * t431 + t458;
t488 = t316 * t578 - t504;
t375 = pkin(3) * t417 + t479;
t487 = t375 * t535 + t345;
t486 = -t440 * t535 - t345;
t485 = t440 * t417 - t503;
t326 = t393 * t570 + t489;
t339 = -t392 * pkin(4) + t356;
t483 = t469 * t326 + t466 * t339 - t371 * t531 + t379 * t530;
t481 = t307 * t417 - t559 + (t466 * t535 + t531) * t330;
t480 = -t375 * t417 - t342;
t478 = t315 * t578 + t484;
t477 = qJD(5) * t497 + t296 * t466 - t297;
t476 = ((-t347 - t506) * t469 + (t346 + t575) * t466) * MDP(23) + (-t466 * t506 - t557) * MDP(22) + (t396 * t417 - t466 * t508 - t383) * MDP(24) + (-t394 * t417 - t469 * t508 + t551) * MDP(25) + (t538 + (t417 - t516) * t462) * MDP(13) + (-t417 ^ 2 + t535 ^ 2) * MDP(12) + (MDP(11) * t535 + MDP(26) * t578) * t417;
t451 = qJ(4) + t567;
t425 = t439 + t567;
t412 = t417 * qJ(6);
t387 = pkin(3) * t430 + t495;
t380 = -pkin(4) * t430 - t399;
t370 = t385 * t431;
t360 = pkin(5) * t396 + qJ(6) * t394;
t344 = t430 * t494 - t399;
t341 = pkin(3) * t393 + t489;
t338 = -pkin(4) * t393 - t355;
t334 = -pkin(5) * t431 + t371 * t466 - t379 * t469;
t333 = qJ(6) * t431 + t540;
t320 = t364 * t466 - t548 + t565;
t319 = -t412 + t541;
t318 = t357 * t466 - t372 * t469 + t565;
t317 = -t412 + t542;
t304 = (qJD(5) * t572 - qJD(6) * t466) * t430 + t494 * t393 - t355;
t300 = pkin(5) * t392 + qJD(5) * t540 + t326 * t466 - t339 * t469;
t299 = -qJ(6) * t392 + qJD(6) * t431 + t483;
t1 = [(-t504 * t431 - t315 * t392 + t338 * t394 + t380 * t347 + ((-qJD(5) * t379 - t326) * t578 + t371 * t385 + t361 * qJD(5) * t430) * t466 + ((-qJD(5) * t371 + t339) * t578 - t379 * t385 - t327 * t430 - t361 * t393) * t469) * MDP(27) + (-t392 * t578 - t370) * MDP(26) + (-t347 * t431 - t383 * t430 + t392 * t394 - t490 * t578) * MDP(25) + (-t346 * t431 - t392 * t396 - t430 * t551 + t491 * t578) * MDP(24) + (t296 * t431 + t299 * t578 - t303 * t552 - t304 * t396 - t307 * t392 - t330 * t491 - t333 * t385 + t344 * t346) * MDP(31) + (t316 * t392 + t327 * t552 + t338 * t396 - t380 * t346 + t361 * t491 + t385 * t540 + t431 * t484 - t483 * t578) * MDP(28) + (-t298 * t431 - t300 * t578 + t304 * t394 + t305 * t392 + t330 * t490 + t334 * t385 + t344 * t347 - t430 * t559) * MDP(29) + (-t299 * t394 + t300 * t396 - t333 * t347 - t334 * t346 + t497 * t393 + (qJD(5) * t498 + t296 * t469 + t298 * t466) * t430) * MDP(30) + (-MDP(13) * t392 - MDP(14) * t393) * t462 + MDP(6) * t547 + (t385 * t430 + t392 * t417 - t393 * t535 - t431 * t474) * MDP(12) + (-t335 * t430 - t341 * t417 - t375 * t393 - t387 * t474 + t555) * MDP(19) + (-t392 * t535 - t370) * MDP(11) + (t342 * t430 + t345 * t431 + t355 * t417 + t356 * t535 - t377 * t392 + t381 * t393 - t398 * t385 + t399 * t475) * MDP(18) + (-t335 * t431 - t341 * t535 + t375 * t392 + t385 * t387 - t556) * MDP(20) + (-t385 * t455 - t392 * t440 + 0.2e1 * t535 * t458 + t556) * MDP(17) + 0.2e1 * t470 * MDP(4) * t512 + t576 * t577 + (t296 * t333 + t298 * t334 + t299 * t307 + t300 * t305 + t303 * t344 + t304 * t330) * MDP(32) + (-pkin(7) * t547 + t468 * t511) * MDP(9) - MDP(7) * t549 + (pkin(7) * t549 + t470 * t511) * MDP(10) + (-t346 * t552 + t396 * t491) * MDP(22) + (t440 * t393 + t417 * t458 + t430 * t449 + t455 * t475 - t555) * MDP(16) + ((-t394 * t466 + t396 * t469) * t393 + (-t557 - t347 * t466 + (-t394 * t469 - t396 * t466) * qJD(5)) * t430) * MDP(23) + (t335 * t387 + t341 * t375 + t342 * t399 + t345 * t398 + t355 * t381 + t356 * t377) * MDP(21); t473 * t576 + (-t520 + t347 * t425 - t543 * t394 + (-t450 * t531 + t469 * t521 + t318) * t578 + t493) * MDP(29) + (-t520 + t451 * t347 + t537 * t394 + ((-t372 + t521) * t469 + (-qJD(5) * t450 + t357) * t466) * t578 + t492) * MDP(27) + t476 + (-t450 * t551 + t346 * t425 + t543 * t396 + (t466 * t521 - t317 + t515) * t578 + t481) * MDP(31) + (t303 * t425 - t305 * t318 - t307 * t317 - t330 * t543 + t450 * t477 - t498 * t521) * MDP(32) - t468 * MDP(4) * t546 + (-t454 * t385 + t417 * t539 - t451 * t475 - t535 * t574 + t496) * MDP(18) + (t390 * t462 + (-t417 * t534 - t462 * t533) * pkin(2) + t486) * MDP(16) + (t378 * t535 - t462 * t539 + t480) * MDP(20) + (t378 * t417 - t462 * t574 + t487) * MDP(19) + (t391 * t462 + (-t462 * t513 - t534 * t535) * pkin(2) + t485) * MDP(17) + (t317 * t394 - t318 * t396 + t297 + (-t396 * t521 - t558 + t346 * t450 + (-t394 * t450 - t307) * qJD(5)) * t469 + (-t394 * t521 - t347 * t450 + (t396 * t450 - t305) * qJD(5) + t510) * t466) * MDP(30) + (-t451 * t346 + (-t515 + t542) * t578 + t537 * t396 + (t450 * t385 - t521 * t578 - t507) * t466 + t509) * MDP(28) + (-t342 * t451 + t345 * t454 - t375 * t378 - t377 * t574 + t381 * t539) * MDP(21) + (MDP(9) * t468 * t473 + MDP(10) * t546) * pkin(1); (t303 * t439 - t305 * t320 - t307 * t319 - t330 * t544 - t477 * t570) * MDP(32) + (-pkin(3) * t345 - qJ(4) * t342 - t375 * t386 - t377 * t389 - t381 * t526) * MDP(21) + (t519 + t347 * t439 + (t466 * t529 + t320) * t578 - t544 * t394 + t493) * MDP(29) + (t466 * t545 + t346 * t439 + (-t319 - t514) * t578 + t544 * t396 + t481) * MDP(31) + t476 + (t319 * t394 - t320 * t396 + t297 + (-t558 - t346 * t570 + (t394 * t570 - t307) * qJD(5)) * t469 + (t347 * t570 + (-t396 * t570 - t305) * qJD(5) + t510) * t466) * MDP(30) + (t486 + t554) * MDP(16) + (t386 * t417 + t487 - t554) * MDP(19) + (t386 * t535 + t462 * t526 + t480) * MDP(20) + (-qJ(4) * t346 + (t514 + t541) * t578 + t527 * t396 + (-t507 - t545) * t466 + t509) * MDP(28) + (pkin(3) * t385 - qJ(4) * t475 - t389 * t535 - t417 * t526 + t496) * MDP(18) + (t519 + qJ(4) * t347 + (-t548 + (t364 + t529) * t466) * t578 + t527 * t394 + t492) * MDP(27) + (-t388 * t462 + t485) * MDP(17); (-qJD(3) * t516 - t467 * t512 + t538) * MDP(18) + t345 * MDP(21) - t383 * MDP(27) - t297 * MDP(32) + (-t417 * MDP(19) - MDP(20) * t535 + t375 * MDP(21)) * t535 + (t417 * MDP(18) - MDP(20) * t462 + t381 * MDP(21) - t330 * MDP(32) + t394 * t571 - t522 * t396) * t462 + (-t385 * MDP(29) + t582 + (MDP(32) * t307 - t522 * t578) * t578) * t469 + ((t396 * t535 - t347 + t532) * MDP(30) + (qJD(5) * t305 - t510) * MDP(32) + t522 * t385 + t571 * t508) * t466; -t581 * MDP(24) + (-t417 * t531 - t462 * t530 + t384) * MDP(25) - t385 * MDP(26) + t488 * MDP(27) + t478 * MDP(28) + (t488 - 0.2e1 * t566) * MDP(29) + (pkin(5) * t346 - qJ(6) * t347) * MDP(30) + (-t478 + 0.2e1 * t528 - 0.2e1 * t561) * MDP(31) + (-pkin(5) * t298 + qJ(6) * t296 - t305 * t316 + t307 * t524 - t330 * t360) * MDP(32) + (t578 * MDP(25) - t361 * MDP(27) - t330 * MDP(29) + (t307 - t316) * MDP(30) + t360 * MDP(31) + MDP(23) * t396) * t396 + (t396 * MDP(22) + t361 * MDP(28) - t360 * MDP(29) + (t305 - t524) * MDP(30) - t330 * MDP(31) - MDP(23) * t394) * t394; (t394 * t396 + t385) * MDP(29) - t582 + (-t396 ^ 2 - t508) * MDP(31) + (-t307 * t578 + t330 * t396 + t298) * MDP(32);];
tauc  = t1;
