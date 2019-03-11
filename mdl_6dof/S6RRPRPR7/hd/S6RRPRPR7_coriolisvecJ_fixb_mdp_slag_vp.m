% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:53
% EndTime: 2019-03-09 10:48:04
% DurationCPUTime: 6.30s
% Computational Cost: add. (4723->395), mult. (11191->528), div. (0->0), fcn. (7845->8), ass. (0->195)
t460 = sin(qJ(4));
t461 = sin(qJ(2));
t463 = cos(qJ(4));
t464 = cos(qJ(2));
t411 = t460 * t461 + t463 * t464;
t394 = t411 * qJD(1);
t535 = qJD(1) * t464;
t536 = qJD(1) * t461;
t396 = -t460 * t535 + t463 * t536;
t457 = sin(pkin(10));
t568 = cos(pkin(10));
t489 = -t457 * t394 + t396 * t568;
t459 = sin(qJ(6));
t530 = qJD(6) * t459;
t485 = t411 * qJD(4);
t526 = qJD(1) * qJD(2);
t518 = t464 * t526;
t519 = t461 * t526;
t542 = t460 * t519 + t463 * t518;
t361 = -qJD(1) * t485 + t542;
t532 = qJD(4) * t460;
t521 = t464 * t532;
t534 = qJD(2) * t461;
t484 = t463 * t534 + t521;
t531 = qJD(4) * t463;
t523 = t461 * t531;
t541 = qJD(1) * t523 + t460 * t518;
t475 = qJD(1) * t484 - t541;
t334 = t361 * t568 + t457 * t475;
t452 = qJD(2) - qJD(4);
t462 = cos(qJ(6));
t529 = qJD(6) * t462;
t548 = t462 * t334 - t452 * t529;
t314 = -t489 * t530 + t548;
t346 = -t452 * t459 + t462 * t489;
t562 = t334 * t459;
t315 = qJD(6) * t346 + t562;
t560 = t346 * t489;
t564 = t314 * t459;
t333 = t361 * t457 - t475 * t568;
t552 = t459 * t333;
t358 = t568 * t394 + t457 * t396;
t591 = qJD(6) + t358;
t596 = t591 * t462;
t582 = -t591 * t596 - t552;
t586 = t346 * t591;
t558 = t489 * t459;
t344 = t462 * t452 + t558;
t587 = t344 * t591;
t604 = -((t315 + t586) * t459 - (t314 - t587) * t462) * MDP(25) + (t346 * t596 + t564) * MDP(24) + (-t560 - t582) * MDP(26) + (-t396 * t452 + t475) * MDP(18) - t489 * MDP(28) * t591 - (t394 ^ 2 - t396 ^ 2) * MDP(16) + t396 * t394 * MDP(15);
t444 = pkin(7) * t536;
t598 = -pkin(8) * t536 + qJD(3) + t444;
t445 = pkin(7) * t535;
t417 = -pkin(8) * t535 + t445;
t465 = -pkin(2) - pkin(3);
t502 = -qJ(3) * t460 + t463 * t465;
t594 = qJD(4) * t502 - t460 * t417 + t598 * t463;
t420 = qJ(3) * t463 + t460 * t465;
t593 = qJD(4) * t420 + t463 * t417 + t598 * t460;
t590 = -pkin(5) * t489 - pkin(9) * t358;
t524 = t465 * qJD(2);
t379 = t524 + t598;
t454 = qJD(2) * qJ(3);
t397 = t417 + t454;
t514 = t463 * t379 - t397 * t460;
t566 = qJ(5) * t396;
t341 = t514 - t566;
t338 = -pkin(4) * t452 + t341;
t499 = t379 * t460 + t397 * t463;
t567 = qJ(5) * t394;
t342 = t499 - t567;
t339 = t568 * t342;
t317 = t457 * t338 + t339;
t313 = -pkin(9) * t452 + t317;
t398 = -qJD(1) * pkin(1) - pkin(2) * t535 - qJ(3) * t536;
t378 = pkin(3) * t535 - t398;
t354 = pkin(4) * t394 + qJD(5) + t378;
t320 = pkin(5) * t358 - pkin(9) * t489 + t354;
t304 = t313 * t462 + t320 * t459;
t588 = t304 * t489;
t561 = t344 * t489;
t500 = t313 * t459 - t320 * t462;
t585 = t489 * t500;
t584 = -t566 + t594;
t583 = t567 - t593;
t571 = pkin(7) - pkin(8);
t416 = t571 * t534;
t453 = qJD(2) * qJD(3);
t384 = -qJD(1) * t416 + t453;
t435 = pkin(7) * t518;
t405 = -pkin(8) * t518 + t435;
t473 = -t499 * qJD(4) - t460 * t384 + t463 * t405;
t581 = t378 * t396 - t473;
t579 = -0.2e1 * t526;
t455 = t461 ^ 2;
t578 = MDP(5) * (-t464 ^ 2 + t455);
t448 = t461 * qJD(3);
t540 = qJ(3) * t518 + qJD(1) * t448;
t468 = (-pkin(4) * t521 + (-pkin(4) * t463 + t465) * t534) * qJD(1) + pkin(4) * t541 + t540;
t487 = t379 * t531 + t463 * t384 - t397 * t532 + t460 * t405;
t311 = qJ(5) * t475 - t394 * qJD(5) + t487;
t469 = -qJ(5) * t361 - qJD(5) * t396 + t473;
t301 = t311 * t457 - t469 * t568;
t438 = pkin(4) * t457 + pkin(9);
t570 = pkin(4) * t396;
t575 = t591 * (qJD(6) * t438 + t570 - t590) + t301;
t414 = -pkin(4) + t502;
t368 = t457 * t414 + t568 * t420;
t366 = -pkin(9) + t368;
t442 = qJ(3) * t535;
t386 = t465 * t536 + t442;
t481 = t386 - t570;
t574 = t591 * (qJD(6) * t366 + t481 + t590) - t301;
t302 = t311 * t568 + t457 * t469;
t533 = qJD(2) * t464;
t369 = t460 * t533 - t484 + t523;
t425 = t571 * t464;
t418 = qJD(2) * t425;
t424 = t571 * t461;
t486 = -t463 * t416 + t460 * t418 + t424 * t531 - t425 * t532;
t323 = -qJ(5) * t369 - qJD(5) * t411 + t486;
t370 = qJD(2) * t411 - t485;
t412 = -t460 * t464 + t461 * t463;
t497 = -t424 * t460 - t425 * t463;
t476 = qJD(4) * t497 + t416 * t460 + t463 * t418;
t470 = -qJ(5) * t370 - qJD(5) * t412 + t476;
t308 = t323 * t568 + t457 * t470;
t554 = t457 * t342;
t316 = t338 * t568 - t554;
t312 = t452 * pkin(5) - t316;
t363 = t411 * t568 + t412 * t457;
t364 = -t457 * t411 + t412 * t568;
t422 = -t464 * pkin(2) - t461 * qJ(3) - pkin(1);
t407 = t464 * pkin(3) - t422;
t492 = pkin(4) * t411 + t407;
t327 = pkin(5) * t363 - pkin(9) * t364 + t492;
t336 = -t457 * t369 + t370 * t568;
t352 = -qJ(5) * t411 - t497;
t483 = -qJ(5) * t412 + t424 * t463 - t425 * t460;
t329 = t352 * t568 + t457 * t483;
t501 = t301 * t364 - t329 * t333;
t572 = t312 * t336 - (qJD(6) * t327 + t308) * t591 - (qJD(6) * t320 + t302) * t363 + t501;
t569 = qJD(2) * pkin(2);
t565 = t312 * t364;
t563 = t327 * t333;
t559 = t591 * t459;
t556 = t394 * t452;
t466 = qJD(2) ^ 2;
t551 = t461 * t466;
t330 = t462 * t333;
t550 = t464 * t466;
t467 = qJD(1) ^ 2;
t549 = t464 * t467;
t547 = t584 * t457 - t583 * t568;
t546 = t583 * t457 + t584 * t568;
t410 = t457 * t463 + t460 * t568;
t545 = t452 * t410;
t488 = -t457 * t460 + t463 * t568;
t544 = t452 * t488;
t539 = qJ(3) * t533 + t448;
t525 = t461 * t549;
t516 = pkin(1) * t579;
t515 = qJD(3) - t569;
t509 = qJD(1) * t422 + t398;
t504 = t452 ^ 2;
t503 = t461 * t524;
t496 = qJD(6) * t410 + t536;
t495 = -t358 * t559 - t530 * t591 + t330;
t377 = pkin(2) * t519 - t540;
t387 = pkin(2) * t534 - t539;
t493 = -pkin(7) * t466 - qJD(1) * t387 - t377;
t490 = t336 * t462 - t364 * t530;
t367 = t414 * t568 - t457 * t420;
t373 = t503 + t539;
t480 = pkin(4) * t369 + t373;
t479 = t378 * t394 - t487;
t319 = t341 * t568 - t554;
t478 = -t438 * t333 + (t312 + t319) * t591;
t477 = -t366 * t333 + (-t312 - t546) * t591;
t419 = -pkin(7) * t519 + t453;
t421 = t444 + t515;
t423 = t445 + t454;
t471 = t419 * t464 + (t421 * t464 + (-t423 + t445) * t461) * qJD(2);
t439 = -pkin(4) * t568 - pkin(5);
t413 = pkin(2) * t536 - t442;
t371 = qJD(1) * t503 + t540;
t365 = pkin(5) - t367;
t335 = t369 * t568 + t370 * t457;
t328 = t352 * t457 - t483 * t568;
t318 = t341 * t457 + t339;
t309 = pkin(5) * t335 - pkin(9) * t336 + t480;
t307 = t323 * t457 - t470 * t568;
t306 = t333 * pkin(5) - t334 * pkin(9) + t468;
t305 = t462 * t306;
t1 = [(-t361 * t411 - t396 * t369 - t370 * t394 + t412 * t475) * MDP(16) + (t378 * t369 + t371 * t411 + t373 * t394 - t407 * t475) * MDP(20) + (t461 * t493 - t509 * t533) * MDP(13) + (t464 * t493 + t509 * t534) * MDP(11) + (pkin(7) * t471 + t377 * t422 + t387 * t398) * MDP(14) + t471 * MDP(12) + (t314 * t364 * t462 + t346 * t490) * MDP(24) + ((-t344 * t462 - t346 * t459) * t336 + (-t564 - t315 * t462 + (t344 * t459 - t346 * t462) * qJD(6)) * t364) * MDP(25) + (-t304 * t335 + t307 * t346 + t328 * t314 + (-(-qJD(6) * t329 + t309) * t591 - t563 - (-qJD(6) * t313 + t306) * t363 - qJD(6) * t565) * t459 + t572 * t462) * MDP(30) + (-t500 * t335 + t305 * t363 + t307 * t344 + t328 * t315 + (t309 * t591 + t563 + (-t313 * t363 - t329 * t591 + t565) * qJD(6)) * t462 + t572 * t459) * MDP(29) + (t314 * t363 + t330 * t364 + t335 * t346 + t490 * t591) * MDP(26) - MDP(7) * t551 + (pkin(7) * t551 + t464 * t516) * MDP(10) + (-t364 * t552 - t315 * t363 - t335 * t344 + (-t336 * t459 - t364 * t529) * t591) * MDP(27) + (-t302 * t363 + t307 * t489 - t308 * t358 - t316 * t336 - t317 * t335 + t328 * t334 + t501) * MDP(22) + MDP(6) * t550 + (t301 * t328 + t302 * t329 - t316 * t307 + t317 * t308 + t354 * t480 + t468 * t492) * MDP(23) + (t407 * t361 + t378 * t370 + t371 * t412 + t373 * t396) * MDP(21) + t578 * t579 + (-pkin(7) * t550 + t461 * t516) * MDP(9) + 0.2e1 * t461 * MDP(4) * t518 + (t333 * t363 + t335 * t591) * MDP(28) + (t361 * t412 + t370 * t396) * MDP(15) + (-t370 * MDP(17) + t369 * MDP(18) - MDP(20) * t476 + MDP(21) * t486) * t452; -MDP(4) * t525 + t467 * t578 + 0.2e1 * t453 * MDP(13) + (qJ(3) * t419 + qJD(3) * t423 - t398 * t413) * MDP(14) + (-t542 + t556) * MDP(17) + (-t386 * t394 + t593 * t452 + t581) * MDP(20) + (-t386 * t396 + t594 * t452 - t479) * MDP(21) + (-t333 * t368 - t334 * t367 + (-t317 + t547) * t489 + (t316 - t546) * t358) * MDP(22) + (-t301 * t367 + t302 * t368 - t316 * t547 + t317 * t546 - t354 * t481) * MDP(23) + (t559 * t591 - t330 - t561) * MDP(27) + (t365 * t315 + t547 * t344 + t477 * t459 - t462 * t574 - t585) * MDP(29) + (t365 * t314 + t547 * t346 + t459 * t574 + t477 * t462 - t588) * MDP(30) + ((-t398 * t461 + t413 * t464) * MDP(11) + ((t423 - t454) * t461 + (-t421 + t515) * t464) * MDP(12) + (t398 * t464 + t413 * t461) * MDP(13) + (t423 * t461 + (-t421 - t569) * t464) * pkin(7) * MDP(14) + MDP(17) * t485) * qJD(1) + (MDP(9) * t461 * t467 + MDP(10) * t549) * pkin(1) - t604; -MDP(11) * t525 + (-t455 * t467 - t466) * MDP(13) + (-qJD(2) * t423 + t398 * t536 + t435) * MDP(14) + (-t394 * t536 - t460 * t504) * MDP(20) + (-t396 * t536 - t463 * t504) * MDP(21) + (-t333 * t410 - t334 * t488 + t358 * t544 - t489 * t545) * MDP(22) + (-t301 * t488 + t302 * t410 + t316 * t545 - t317 * t544 - t354 * t536) * MDP(23) + (-t410 * t552 - t488 * t315 - t545 * t344 + (t459 * t544 - t462 * t496) * t591) * MDP(29) + (-t410 * t330 - t488 * t314 - t545 * t346 + (t459 * t496 + t462 * t544) * t591) * MDP(30); (t361 - t556) * MDP(17) + (-t452 * t499 - t581) * MDP(20) + (-t452 * t514 + t479) * MDP(21) + ((-t333 * t457 - t334 * t568) * pkin(4) - (t316 - t319) * t358 + (t317 - t318) * t489) * MDP(22) + (t316 * t318 - t317 * t319 + (-t301 * t568 + t302 * t457 - t354 * t396) * pkin(4)) * MDP(23) + (t495 + t561) * MDP(27) + (t439 * t315 - t318 * t344 + t478 * t459 - t462 * t575 + t585) * MDP(29) + (t439 * t314 - t318 * t346 + t459 * t575 + t478 * t462 + t588) * MDP(30) + t604; (-t358 ^ 2 - t489 ^ 2) * MDP(22) + (t495 - t561) * MDP(29) + (-t560 + t582) * MDP(30) + (t316 * t489 + t317 * t358 + t468) * MDP(23); t346 * t344 * MDP(24) + (-t344 ^ 2 + t346 ^ 2) * MDP(25) + (t548 + t587) * MDP(26) + (-t562 + t586) * MDP(27) + t333 * MDP(28) + (-t302 * t459 + t304 * t591 - t312 * t346 + t305) * MDP(29) + (-t302 * t462 - t306 * t459 + t312 * t344 - t500 * t591) * MDP(30) + (-MDP(26) * t558 - MDP(27) * t346 - MDP(29) * t304 + MDP(30) * t500) * qJD(6);];
tauc  = t1;
