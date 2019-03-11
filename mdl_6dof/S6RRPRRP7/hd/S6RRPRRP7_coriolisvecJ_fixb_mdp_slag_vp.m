% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:16
% EndTime: 2019-03-09 12:20:27
% DurationCPUTime: 6.06s
% Computational Cost: add. (5658->474), mult. (12489->594), div. (0->0), fcn. (8041->6), ass. (0->196)
t438 = sin(qJ(4));
t441 = cos(qJ(4));
t442 = cos(qJ(2));
t515 = qJD(1) * t442;
t439 = sin(qJ(2));
t516 = qJD(1) * t439;
t388 = -t438 * t515 + t441 * t516;
t440 = cos(qJ(5));
t530 = t441 * t442;
t393 = t438 * t439 + t530;
t457 = t393 * qJD(4);
t354 = qJD(2) * t393 - t457;
t448 = t354 * qJD(1);
t502 = qJD(2) - qJD(4);
t482 = t440 * t502;
t437 = sin(qJ(5));
t508 = qJD(5) * t437;
t323 = qJD(5) * t482 + t388 * t508 - t440 * t448;
t359 = t440 * t388 - t437 * t502;
t509 = qJD(5) * t359;
t324 = t437 * t448 + t509;
t503 = qJD(1) * qJD(2);
t495 = t439 * t503;
t510 = qJD(4) * t441;
t511 = qJD(4) * t438;
t512 = qJD(2) * t442;
t584 = t438 * t512 + t439 * t510 - t442 * t511;
t346 = qJD(1) * t584 - t441 * t495;
t357 = t388 * t437 + t482;
t385 = t393 * qJD(1);
t531 = t440 * t346;
t533 = t437 * t346;
t545 = t323 * t437;
t575 = qJD(5) + t385;
t565 = t575 * t359;
t567 = t575 ^ 2;
t586 = t357 * t575;
t587 = -(t388 * t502 + t346) * MDP(18) - (t437 * (t324 + t565) + (t323 + t586) * t440) * MDP(23) + (t440 * t565 - t545) * MDP(22) + (-t359 * t388 + t440 * t567 + t533) * MDP(24) - (-t357 * t388 + t437 * t567 - t531) * MDP(25) - (t385 ^ 2 - t388 ^ 2) * MDP(16) + (MDP(15) * t385 - MDP(26) * t575) * t388;
t390 = -qJD(1) * pkin(1) - pkin(2) * t515 - qJ(3) * t516;
t370 = pkin(3) * t515 - t390;
t333 = pkin(4) * t385 - pkin(9) * t388 + t370;
t443 = -pkin(2) - pkin(3);
t499 = t443 * qJD(2);
t423 = pkin(7) * t516;
t399 = pkin(8) * t516 - t423;
t582 = qJD(3) - t399;
t373 = t499 + t582;
t424 = pkin(7) * t515;
t401 = -pkin(8) * t515 + t424;
t433 = qJD(2) * qJ(3);
t389 = t401 + t433;
t343 = t438 * t373 + t441 * t389;
t338 = -pkin(9) * t502 + t343;
t304 = t333 * t437 + t338 * t440;
t301 = qJ(6) * t575 + t304;
t585 = t575 * t301;
t476 = qJ(3) * t441 + t438 * t443;
t398 = -pkin(9) + t476;
t535 = t398 * t346;
t563 = -t438 * qJ(3) + t441 * t443;
t374 = qJD(3) * t441 + qJD(4) * t563;
t537 = t374 * t575;
t342 = t441 * t373 - t438 * t389;
t337 = pkin(4) * t502 - t342;
t309 = t357 * pkin(5) - t359 * qJ(6) + t337;
t566 = t575 * t309;
t577 = t566 + t535 + t537;
t552 = pkin(9) * t346;
t576 = t566 - t552;
t570 = -0.2e1 * t503;
t434 = t439 ^ 2;
t569 = MDP(5) * (-t442 ^ 2 + t434);
t568 = t337 * t575;
t351 = t399 * t441 + t401 * t438;
t564 = t351 - t374;
t406 = -t442 * pkin(2) - t439 * qJ(3) - pkin(1);
t392 = t442 * pkin(3) - t406;
t394 = -t438 * t442 + t439 * t441;
t341 = pkin(4) * t393 - pkin(9) * t394 + t392;
t555 = pkin(7) - pkin(8);
t408 = t555 * t439;
t409 = t555 * t442;
t362 = t408 * t438 + t409 * t441;
t522 = t437 * t341 + t440 * t362;
t521 = qJD(4) * t476 + t401 * t441 + t438 * t582;
t361 = -t408 * t441 + t438 * t409;
t558 = t502 ^ 2;
t477 = pkin(5) * t437 - qJ(6) * t440;
t557 = qJD(6) * t437 - t477 * t575;
t556 = t359 ^ 2;
t554 = pkin(5) * t346;
t553 = pkin(5) * t388;
t551 = qJD(2) * pkin(2);
t550 = qJ(6) * t346;
t549 = qJ(6) * t388;
t514 = qJD(2) * t439;
t400 = t555 * t514;
t432 = qJD(2) * qJD(3);
t378 = -qJD(1) * t400 + t432;
t494 = t442 * t503;
t417 = pkin(7) * t494;
t391 = -pkin(8) * t494 + t417;
t467 = -t373 * t511 - t438 * t378 - t389 * t510 + t441 * t391;
t298 = pkin(5) * t324 + qJ(6) * t323 - qJD(6) * t359 - t467;
t548 = t298 * t437;
t303 = t333 * t440 - t338 * t437;
t547 = t303 * t388;
t546 = t304 * t575;
t543 = t324 * t440;
t542 = t357 * t359;
t541 = t357 * t437;
t540 = t357 * t440;
t539 = t359 * t437;
t538 = t359 * t440;
t536 = t394 * t440;
t444 = qJD(2) ^ 2;
t532 = t439 * t444;
t528 = t442 * t444;
t445 = qJD(1) ^ 2;
t527 = t442 * t445;
t526 = t343 + t557;
t525 = -t557 - t521;
t348 = pkin(4) * t388 + pkin(9) * t385;
t524 = t440 * t342 + t437 * t348;
t421 = qJ(3) * t515;
t380 = t443 * t516 + t421;
t334 = -t348 + t380;
t523 = t437 * t334 + t440 * t351;
t427 = t439 * qJD(3);
t519 = qJ(3) * t494 + qJD(1) * t427;
t518 = qJ(3) * t512 + t427;
t513 = qJD(2) * t441;
t507 = qJD(5) * t440;
t504 = qJD(6) - t303;
t501 = pkin(9) * t508;
t500 = t439 * t527;
t498 = t398 * t508;
t493 = pkin(1) * t570;
t492 = qJD(3) - t551;
t485 = qJD(1) * t406 + t390;
t310 = t346 * pkin(4) + (pkin(9) * t457 + (-pkin(9) * t530 + (-pkin(9) * t438 + t443) * t439) * qJD(2)) * qJD(1) + t519;
t460 = -t373 * t510 - t441 * t378 + t389 * t511 - t438 * t391;
t481 = -t440 * t310 + t333 * t508 + t338 * t507 - t437 * t460;
t480 = t439 * t499;
t478 = pkin(5) * t440 + qJ(6) * t437;
t300 = -pkin(5) * t575 + t504;
t475 = -t298 * t440 + t300 * t388;
t474 = -t301 * t388 - t548;
t473 = t300 * t440 - t301 * t437;
t472 = t300 * t437 + t301 * t440;
t471 = t304 * t388 - t437 * t467;
t405 = -pkin(4) - t478;
t369 = pkin(2) * t495 - t519;
t382 = pkin(2) * t514 - t518;
t466 = -pkin(7) * t444 - qJD(1) * t382 - t369;
t464 = t354 * t437 + t394 * t507;
t463 = -t354 * t440 + t394 * t508;
t462 = -t552 + t568;
t461 = t309 * t359 + t481;
t365 = t480 + t518;
t459 = t437 * t310 + t333 * t507 - t338 * t508 - t440 * t460;
t353 = -t439 * t513 + t584;
t318 = pkin(4) * t353 - pkin(9) * t354 + t365;
t402 = qJD(2) * t409;
t329 = -qJD(4) * t361 - t400 * t441 + t402 * t438;
t458 = t437 * t318 + t440 * t329 + t341 * t507 - t362 * t508;
t456 = -t535 - t568;
t453 = -t370 * t388 + t467;
t452 = t370 * t385 + t460;
t294 = qJD(6) * t575 + t459 + t550;
t295 = t481 - t554;
t450 = qJD(5) * t473 + t294 * t440 + t295 * t437;
t330 = qJD(4) * t362 - t400 * t438 - t402 * t441;
t403 = -pkin(7) * t495 + t432;
t404 = t423 + t492;
t407 = t424 + t433;
t449 = t403 * t442 + (t404 * t442 + (-t407 + t424) * t439) * qJD(2);
t397 = pkin(4) - t563;
t396 = pkin(2) * t516 - t421;
t387 = t437 * t516 + t440 * t513;
t384 = t437 * t513 - t440 * t516;
t371 = -t405 - t563;
t356 = qJD(1) * t480 + t519;
t331 = pkin(5) * t359 + qJ(6) * t357;
t326 = t394 * t477 + t361;
t314 = -pkin(5) * t393 - t341 * t440 + t362 * t437;
t313 = qJ(6) * t393 + t522;
t312 = t342 * t437 - t348 * t440 - t553;
t311 = t524 + t549;
t306 = -t334 * t440 + t351 * t437 + t553;
t305 = t523 - t549;
t302 = -t323 + t586;
t299 = t477 * t354 + (qJD(5) * t478 - qJD(6) * t440) * t394 + t330;
t297 = -pkin(5) * t353 + qJD(5) * t522 - t318 * t440 + t329 * t437;
t296 = qJ(6) * t353 + qJD(6) * t393 + t458;
t1 = [(t370 * t354 + t356 * t394 + t365 * t388 + t392 * t448) * MDP(21) + ((-t539 - t540) * t354 + (t545 - t543 + (-t538 + t541) * qJD(5)) * t394) * MDP(23) + (-t323 * t536 - t359 * t463) * MDP(22) - MDP(7) * t532 + (pkin(7) * t532 + t442 * t493) * MDP(10) + (-pkin(7) * t528 + t439 * t493) * MDP(9) + (t439 * t466 - t485 * t512) * MDP(13) + (t442 * t466 + t485 * t514) * MDP(11) + t449 * MDP(12) + (pkin(7) * t449 + t369 * t406 + t382 * t390) * MDP(14) + (t388 * t354 + t394 * t448) * MDP(15) + (t392 * t346 + t370 * t353 + t356 * t393 + t365 * t385) * MDP(20) + (-t394 * t346 - t388 * t353 - t354 * t385 - t393 * t448) * MDP(16) + (-t354 * MDP(17) + t353 * MDP(18) + t330 * MDP(20) + t329 * MDP(21)) * t502 + 0.2e1 * t439 * MDP(4) * t494 + MDP(6) * t528 + (t294 * t313 + t295 * t314 + t296 * t301 + t297 * t300 + t298 * t326 + t299 * t309) * MDP(32) + (-t296 * t357 + t297 * t359 - t313 * t324 - t314 * t323 + t473 * t354 + (-qJD(5) * t472 - t294 * t437 + t295 * t440) * t394) * MDP(30) + (-t304 * t353 - t361 * t323 + t330 * t359 - t337 * t463 - t346 * t522 - t393 * t459 - t458 * t575 - t467 * t536) * MDP(28) + (-t481 * t393 + t303 * t353 + t330 * t357 + t361 * t324 + ((-qJD(5) * t362 + t318) * t575 + t341 * t346 + t337 * qJD(5) * t394) * t440 + ((-qJD(5) * t341 - t329) * t575 - t362 * t346 - t467 * t394 + t337 * t354) * t437) * MDP(27) + (-t295 * t393 - t297 * t575 + t299 * t357 - t300 * t353 + t309 * t464 - t314 * t346 + t324 * t326 + t394 * t548) * MDP(29) + (-t324 * t393 - t353 * t357 - t394 * t533 - t464 * t575) * MDP(25) + (t294 * t393 + t296 * t575 - t298 * t536 - t299 * t359 + t301 * t353 + t309 * t463 + t313 * t346 + t323 * t326) * MDP(31) + (-t323 * t393 + t353 * t359 + t394 * t531 - t463 * t575) * MDP(24) + (t346 * t393 + t353 * t575) * MDP(26) + t569 * t570; (t298 * t371 - t525 * t309 + (t374 * t440 - t305) * t301 + (t374 * t437 - t306) * t300 + t450 * t398) * MDP(32) + (t445 * t439 * MDP(9) + MDP(10) * t527) * pkin(1) + (qJ(3) * t403 + qJD(3) * t407 - t390 * t396) * MDP(14) + (t305 * t357 - t306 * t359 + (-t300 * t385 - t324 * t398 - t357 * t374 - t294 + (t359 * t398 - t300) * qJD(5)) * t440 + (t301 * t385 - t323 * t398 + t359 * t374 - t295 + (t357 * t398 + t301) * qJD(5)) * t437) * MDP(30) - MDP(4) * t500 + (-t380 * t388 - t502 * t564 - t452) * MDP(21) + (-t380 * t385 + t502 * t521 - t453) * MDP(20) + 0.2e1 * t432 * MDP(13) + (((t407 - t433) * t439 + (-t404 + t492) * t442) * MDP(12) + (t390 * t442 + t396 * t439) * MDP(13) + (-t390 * t439 + t396 * t442) * MDP(11) + (t407 * t439 + (-t404 - t551) * t442) * pkin(7) * MDP(14)) * qJD(1) + (-t397 * t323 + (t498 + t523) * t575 + t521 * t359 + (t456 - t537) * t440 - t471) * MDP(28) + (t323 * t371 + (-t305 - t498) * t575 + t525 * t359 + t577 * t440 - t474) * MDP(31) + (t324 * t371 + (-t398 * t507 + t306) * t575 - t525 * t357 - t577 * t437 - t475) * MDP(29) + (t547 + t397 * t324 + t521 * t357 + (-t467 + (-qJD(5) * t398 - t334) * t575) * t440 + (t564 * t575 + t456) * t437) * MDP(27) + t445 * t569 - t587; -MDP(11) * t500 + (-t434 * t445 - t444) * MDP(13) + (-qJD(2) * t407 + t390 * t516 + t417) * MDP(14) + (-t385 * t516 - t438 * t558) * MDP(20) + (-t388 * t516 - t441 * t558) * MDP(21) + (t357 * t387 - t359 * t384 + (t539 - t540) * t510 + (-t545 - t543 + (t538 + t541) * qJD(5)) * t438) * MDP(30) + (-t300 * t384 - t301 * t387 + (qJD(4) * t472 - t298) * t441 + (-t309 * t502 + t450) * t438) * MDP(32) + (MDP(27) + MDP(29)) * (-t324 * t441 + t357 * t511 + (-qJD(2) * t357 - t533) * t438 + (-t437 * t510 - t438 * t507 + t384) * t575) + (-MDP(28) + MDP(31)) * (t438 * (t359 * t502 - t508 * t575 + t531) + (t440 * t510 - t387) * t575 - t323 * t441); (-t385 * t502 + t448) * MDP(17) + (-t343 * t502 + t453) * MDP(20) + (-t342 * t502 + t452) * MDP(21) + (-pkin(4) * t324 - t547 - t343 * t357 + (t467 + (-pkin(9) * qJD(5) - t348) * t575) * t440 + (t342 * t575 + t462) * t437) * MDP(27) + (pkin(4) * t323 - t343 * t359 + (t501 + t524) * t575 + t462 * t440 + t471) * MDP(28) + (t324 * t405 + (-pkin(9) * t507 + t312) * t575 - t526 * t357 + t576 * t437 + t475) * MDP(29) + (t311 * t357 - t312 * t359 + (t294 + t575 * t300 + (-t324 + t509) * pkin(9)) * t440 + (t295 - t585 + (qJD(5) * t357 - t323) * pkin(9)) * t437) * MDP(30) + (t323 * t405 + (-t311 - t501) * t575 + t526 * t359 - t576 * t440 + t474) * MDP(31) + (pkin(9) * t450 + t298 * t405 - t300 * t312 - t301 * t311 - t309 * t526) * MDP(32) + t587; MDP(22) * t542 + (-t357 ^ 2 + t556) * MDP(23) + t302 * MDP(24) + (-t324 + t565) * MDP(25) + t346 * MDP(26) + (-t337 * t359 - t481 + t546) * MDP(27) + (t303 * t575 + t337 * t357 - t459) * MDP(28) + (-t331 * t357 - t461 + t546 + 0.2e1 * t554) * MDP(29) + (pkin(5) * t323 - qJ(6) * t324 + (t301 - t304) * t359 + (t300 - t504) * t357) * MDP(30) + (0.2e1 * t550 - t309 * t357 + t331 * t359 + (0.2e1 * qJD(6) - t303) * t575 + t459) * MDP(31) + (-pkin(5) * t295 + qJ(6) * t294 - t300 * t304 + t301 * t504 - t309 * t331) * MDP(32); (-t346 + t542) * MDP(29) + t302 * MDP(30) + (-t556 - t567) * MDP(31) + (t461 - t554 - t585) * MDP(32);];
tauc  = t1;
