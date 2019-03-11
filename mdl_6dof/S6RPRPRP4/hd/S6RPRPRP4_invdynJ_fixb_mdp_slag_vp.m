% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:12
% EndTime: 2019-03-09 03:13:19
% DurationCPUTime: 5.38s
% Computational Cost: add. (3228->492), mult. (6048->591), div. (0->0), fcn. (3487->10), ass. (0->208)
t429 = sin(pkin(9));
t404 = pkin(1) * t429 + pkin(7);
t382 = t404 * qJDD(1);
t569 = qJD(2) * qJD(3) + t382;
t432 = sin(qJ(3));
t508 = qJD(1) * qJD(3);
t492 = t432 * t508;
t435 = cos(qJ(3));
t505 = qJDD(1) * t435;
t567 = -t492 + t505;
t418 = t432 * qJ(4);
t568 = pkin(2) + t418;
t549 = pkin(4) + t404;
t384 = t404 * qJD(1);
t350 = -t435 * qJD(2) + t432 * t384;
t566 = qJD(4) + t350;
t524 = qJD(1) * t432;
t399 = qJD(5) + t524;
t565 = t399 ^ 2;
t491 = t435 * t508;
t506 = qJDD(1) * t432;
t460 = t491 + t506;
t368 = qJDD(5) + t460;
t501 = MDP(22) - MDP(25);
t564 = t501 * t368;
t424 = qJDD(3) * qJ(4);
t425 = qJD(3) * qJD(4);
t521 = qJD(3) * t432;
t485 = t432 * qJDD(2) - t384 * t521 + t569 * t435;
t312 = -t424 - t425 - t485;
t303 = t567 * pkin(4) - t312;
t431 = sin(qJ(5));
t434 = cos(qJ(5));
t522 = qJD(3) * t431;
t523 = qJD(1) * t435;
t373 = t434 * t523 + t522;
t515 = qJD(5) * t373;
t320 = -t434 * qJDD(3) + t567 * t431 + t515;
t493 = t431 * t523;
t321 = -qJD(5) * t493 + qJDD(3) * t431 + (qJD(3) * qJD(5) + t567) * t434;
t520 = qJD(3) * t434;
t375 = -t493 + t520;
t293 = pkin(5) * t321 + qJ(6) * t320 - qJD(6) * t375 + t303;
t550 = g(3) * t432;
t423 = qJ(1) + pkin(9);
t412 = sin(t423);
t413 = cos(t423);
t562 = -g(1) * t413 - g(2) * t412;
t449 = t562 * t435 - t550;
t437 = -pkin(3) - pkin(8);
t511 = qJD(5) * t437;
t442 = -t399 * t511 + t449;
t563 = t293 + t442;
t430 = cos(pkin(9));
t557 = pkin(1) * t430;
t405 = -pkin(2) - t557;
t420 = t435 * pkin(3);
t527 = t420 + t418;
t360 = t405 - t527;
t555 = pkin(8) * t435;
t347 = t360 - t555;
t361 = t549 * t432;
t531 = t434 * t347 + t431 * t361;
t351 = t432 * qJD(2) + t435 * t384;
t426 = qJD(3) * qJ(4);
t338 = -t426 - t351;
t335 = -qJD(3) * pkin(3) + t566;
t510 = pkin(4) * t524 + t566;
t502 = -MDP(21) - MDP(23);
t519 = qJD(3) * t435;
t560 = t384 * t519 + t569 * t432;
t536 = t431 * t432;
t344 = t412 * t434 + t413 * t536;
t346 = -t412 * t536 + t413 * t434;
t422 = g(3) * t435;
t504 = qJDD(2) * t435;
t457 = qJDD(4) - t504 + t560;
t302 = pkin(4) * t460 + qJDD(3) * t437 + t457;
t398 = pkin(3) * t492;
t548 = qJ(4) * t435;
t476 = pkin(8) * t432 - t548;
t518 = qJD(4) * t432;
t446 = qJD(3) * t476 - t518;
t448 = t435 * t437 - t557 - t568;
t309 = qJD(1) * t446 + qJDD(1) * t448 + t398;
t325 = qJD(3) * t437 + t510;
t513 = qJD(5) * t434;
t499 = t431 * t302 + t434 * t309 + t325 * t513;
t330 = t448 * qJD(1);
t516 = qJD(5) * t330;
t559 = -g(1) * t344 + g(2) * t346 + (-t516 + t422) * t431 + t499;
t558 = t375 ^ 2;
t556 = pkin(5) * t368;
t554 = g(1) * t412;
t551 = g(2) * t413;
t547 = qJ(6) * t368;
t546 = qJDD(3) * pkin(3);
t299 = t325 * t431 + t330 * t434;
t545 = t299 * t399;
t544 = t320 * t434;
t543 = t373 * t375;
t542 = t373 * t399;
t541 = t375 * t399;
t540 = t412 * t432;
t539 = t412 * t435;
t538 = t413 * t432;
t537 = t413 * t435;
t318 = t432 * t321;
t535 = t432 * t434;
t359 = t434 * t368;
t534 = t437 * t368;
t533 = t373 * t519 + t318;
t532 = -t432 * t320 + t375 * t519;
t334 = pkin(4) * t523 + t351;
t411 = pkin(3) * t524;
t352 = qJD(1) * t476 + t411;
t530 = t431 * t334 + t434 * t352;
t478 = pkin(5) * t434 + qJ(6) * t431;
t465 = -pkin(4) - t478;
t529 = qJD(5) * t478 - qJD(6) * t434 - t465 * t524 + t566;
t514 = qJD(5) * t431;
t495 = t399 * t514;
t496 = t432 * t520;
t528 = t399 * t496 + t435 * t495;
t362 = t549 * t435;
t427 = t432 ^ 2;
t428 = t435 ^ 2;
t526 = t427 - t428;
t377 = -qJ(4) * t523 + t411;
t525 = qJD(1) * t377;
t385 = qJD(1) * t405;
t298 = t325 * t434 - t330 * t431;
t509 = qJD(6) - t298;
t296 = -pkin(5) * t399 + t509;
t517 = qJD(5) * t296;
t512 = qJD(5) * t435;
t503 = qJDD(3) * t404;
t439 = qJD(1) ^ 2;
t500 = t432 * t435 * t439;
t498 = g(1) * t538 + g(2) * t540 - t422;
t497 = t431 * t521;
t433 = sin(qJ(1));
t490 = -pkin(1) * t433 + t413 * pkin(7);
t487 = -qJD(5) * t375 + t321;
t486 = t434 * t302 - t431 * t309 - t325 * t514 - t330 * t513;
t328 = t426 + t334;
t343 = t412 * t431 - t413 * t535;
t345 = t412 * t535 + t413 * t431;
t484 = -g(1) * t345 - g(2) * t343;
t483 = -g(1) * t346 - g(2) * t344;
t436 = cos(qJ(1));
t481 = g(1) * t433 - g(2) * t436;
t291 = qJD(6) * t399 - t330 * t514 + t499 + t547;
t480 = -t296 * t524 - t291;
t292 = qJDD(6) - t486 - t556;
t297 = qJ(6) * t399 + t299;
t479 = -t297 * t524 + t292;
t477 = pkin(5) * t431 - qJ(6) * t434;
t438 = qJD(3) ^ 2;
t475 = t404 * t438 + t551;
t473 = t436 * pkin(1) + pkin(3) * t537 + t412 * pkin(7) + t568 * t413;
t471 = -t296 * t434 + t297 * t431;
t470 = t296 * t431 + t297 * t434;
t308 = pkin(5) * t373 - qJ(6) * t375 + t328;
t468 = t338 * MDP(15) - t308 * MDP(26);
t466 = -t568 - t420;
t462 = -t368 * t431 - t399 * t513;
t461 = -qJ(4) * t519 - t518;
t410 = pkin(3) * t521;
t336 = t410 + t446;
t357 = t549 * t519;
t459 = t434 * t336 - t347 * t514 + t431 * t357 + t361 * t513;
t458 = -qJD(3) * t350 - t485;
t455 = t466 - t557;
t454 = t308 * t399 + t534;
t453 = -0.2e1 * qJDD(1) * t405 - t475;
t452 = -qJD(3) * t351 - t498 + t560;
t451 = 0.2e1 * qJD(3) * t385 - t503;
t339 = t455 * qJD(1);
t450 = t503 + (-qJD(1) * t360 - t339) * qJD(3);
t447 = g(1) * t343 - g(2) * t345 + t434 * t422 + t486;
t315 = t457 - t546;
t316 = qJD(1) * t461 + qJDD(1) * t455 + t398;
t358 = t410 + t461;
t444 = qJD(1) * t358 + qJDD(1) * t360 + t316 + t475;
t443 = -t312 * t435 + t315 * t432 + (t335 * t435 + t338 * t432) * qJD(3);
t441 = t308 * t375 + qJDD(6) - t447;
t392 = g(1) * t539;
t388 = qJ(4) * t537;
t386 = qJ(4) * t539;
t381 = qJDD(3) * t435 - t432 * t438;
t380 = qJDD(3) * t432 + t435 * t438;
t379 = qJ(4) + t477;
t365 = t399 * t497;
t356 = t549 * t521;
t342 = t373 * t497;
t331 = t339 * t524;
t327 = t435 * t478 + t362;
t326 = pkin(5) * t375 + qJ(6) * t373;
t314 = -pkin(5) * t432 + t347 * t431 - t361 * t434;
t313 = qJ(6) * t432 + t531;
t311 = -pkin(5) * t523 - t334 * t434 + t352 * t431;
t310 = qJ(6) * t523 + t530;
t307 = t542 - t320;
t306 = (-qJD(5) * t477 + qJD(6) * t431) * t435 + (-t404 + t465) * t521;
t295 = -pkin(5) * t519 + t531 * qJD(5) + t336 * t431 - t357 * t434;
t294 = qJ(6) * t519 + qJD(6) * t432 + t459;
t1 = [(t486 * t432 + t298 * t519 - t356 * t373 + t362 * t321 + ((-qJD(5) * t361 - t336) * t399 - t347 * t368 - t328 * t512) * t431 + ((-qJD(5) * t347 + t357) * t399 + t361 * t368 + t303 * t435 - t328 * t521) * t434 + t483) * MDP(21) + (-t459 * t399 - t531 * t368 - t356 * t375 - t362 * t320 + ((qJD(3) * t328 + t516) * t431 - t499) * t432 + (-qJD(3) * t299 - t303 * t431 - t328 * t513) * t435 - t484) * MDP(22) + (qJDD(1) * t427 + 0.2e1 * t432 * t491) * MDP(5) + (t432 * t451 + t435 * t453 + t392) * MDP(10) + (t432 * t450 + t435 * t444 - t392) * MDP(13) + qJDD(1) * MDP(1) + t380 * MDP(7) + (t481 + (t429 ^ 2 + t430 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t481 * MDP(2) + ((t427 + t428) * t382 + t443 + t562) * MDP(12) + t381 * MDP(8) + (t320 * t431 * t435 + (-t434 * t512 + t497) * t375) * MDP(16) + (t368 * t432 + t399 * t519) * MDP(20) + (-t295 * t399 + t306 * t373 - t314 * t368 + t321 * t327 + (-t308 * t520 - t292) * t432 + (-qJD(3) * t296 + t293 * t434 - t308 * t514) * t435 + t483) * MDP(23) + (t294 * t399 - t306 * t375 + t313 * t368 + t320 * t327 + (-t308 * t522 + t291) * t432 + (qJD(3) * t297 + t293 * t431 + t308 * t513) * t435 + t484) * MDP(25) + 0.2e1 * (t432 * t505 - t508 * t526) * MDP(6) + (t435 * t462 + t365 + t532) * MDP(18) + (-t318 + (-qJD(3) * t373 - t359) * t435 + t528) * MDP(19) + (g(1) * t436 + g(2) * t433) * MDP(3) + (t375 * t496 - t342 + (t544 + t321 * t431 + (t373 * t434 + t375 * t431) * qJD(5)) * t435) * MDP(17) + (-t294 * t373 + t295 * t375 - t313 * t321 - t314 * t320 + t392 + t470 * t521 + (qJD(5) * t471 - t291 * t434 - t292 * t431 - t551) * t435) * MDP(24) + (-g(1) * t490 - g(2) * t473 + t316 * t360 + t339 * t358 + t404 * t443 - t466 * t554) * MDP(15) + (t451 * t435 + (-t453 - t554) * t432) * MDP(11) + (t450 * t435 + (-t444 + t554) * t432) * MDP(14) + (t291 * t313 + t297 * t294 + t293 * t327 + t308 * t306 + t292 * t314 + t296 * t295 - g(1) * (pkin(4) * t413 + pkin(5) * t346 + qJ(6) * t345 + t490) - g(2) * (pkin(5) * t344 + pkin(8) * t537 + qJ(6) * t343 + t473) + (-g(1) * (t466 - t555) - g(2) * pkin(4)) * t412) * MDP(26); qJDD(2) * MDP(4) + (t528 + t533) * MDP(21) + t532 * MDP(22) + t533 * MDP(23) - t342 * MDP(24) + t365 * MDP(25) + (MDP(10) - MDP(13)) * t381 + (-MDP(11) + MDP(14)) * t380 + (-MDP(4) - MDP(15) - MDP(26)) * g(3) + (-t312 * MDP(15) + t320 * MDP(25) + t293 * MDP(26) + (t335 * MDP(15) - t375 * t434 * MDP(24) + t471 * MDP(26) + (-t431 * MDP(22) + t434 * MDP(23)) * t399) * qJD(3)) * t432 + (-t315 * MDP(15) + (-MDP(25) * t375 - t468) * qJD(3) + (qJD(5) * t399 * MDP(23) + t487 * MDP(24) + (-t291 - t517) * MDP(26) + t564) * t431 + (-t320 * MDP(24) + t292 * MDP(26) + t502 * t368 + (t373 * MDP(24) - t297 * MDP(26) + t399 * t501) * qJD(5)) * t434) * t435; -MDP(5) * t500 + t526 * MDP(6) * t439 + MDP(7) * t506 + MDP(8) * t505 + qJDD(3) * MDP(9) + (-t385 * t524 - t452 + t504) * MDP(10) + (t550 + (-qJD(1) * t385 - t562) * t435 + t458) * MDP(11) + (-pkin(3) * t432 + t548) * qJDD(1) * MDP(12) + (-0.2e1 * t546 + qJDD(4) + t331 + (-qJDD(2) - t525) * t435 + t452) * MDP(13) + (0.2e1 * t424 + 0.2e1 * t425 + (-g(3) + t525) * t432 + (qJD(1) * t339 + t562) * t435 - t458) * MDP(14) + (-t312 * qJ(4) - t315 * pkin(3) - t339 * t377 - t335 * t351 - g(1) * (-pkin(3) * t538 + t388) - g(2) * (-pkin(3) * t540 + t386) - g(3) * t527 - t566 * t338) * MDP(15) + (-t431 * t541 - t544) * MDP(16) + ((-t321 - t541) * t434 + (t320 + t542) * t431) * MDP(17) + (-t495 + t359 + (-t375 * t435 - t399 * t536) * qJD(1)) * MDP(18) + ((t373 * t435 - t399 * t535) * qJD(1) + t462) * MDP(19) - t399 * MDP(20) * t523 + (-t298 * t523 + qJ(4) * t321 + t510 * t373 + (t534 + (t328 - t334) * t399) * t434 + (t303 + (t352 - t511) * t399 + t449) * t431) * MDP(21) + (-qJ(4) * t320 + t530 * t399 + t299 * t523 + t510 * t375 + (-t328 * t399 - t534) * t431 + (t303 + t442) * t434) * MDP(22) + (t296 * t523 + t311 * t399 + t321 * t379 + t529 * t373 + t563 * t431 + t454 * t434) * MDP(23) + (t310 * t373 - t311 * t375 + (t320 * t437 + (-t373 * t437 - t297) * qJD(5) + t479) * t434 + (-t321 * t437 + (t375 * t437 - t296) * qJD(5) + t480) * t431 + t498) * MDP(24) + (-t297 * t523 - t310 * t399 + t320 * t379 - t529 * t375 + t454 * t431 - t563 * t434) * MDP(25) + (t293 * t379 - t297 * t310 - t296 * t311 - g(1) * t388 - g(2) * t386 - g(3) * (t432 * t477 + t527 + t555) + t529 * t308 + (qJD(5) * t470 + t291 * t431 - t292 * t434) * t437 + t562 * (t432 * t437 + t435 * t477)) * MDP(26); MDP(12) * t506 + (qJDD(3) + t500) * MDP(13) + (-t427 * t439 - t438) * MDP(14) + (t315 + t331 - t498) * MDP(15) + t359 * MDP(21) - t498 * MDP(26) + (t373 * t502 - t375 * t501 + t468) * qJD(3) + (t368 * MDP(23) + (-t373 * t524 + t320 - t515) * MDP(24) + (qJD(5) * t297 - t479) * MDP(26) - t501 * t565) * t434 + ((t375 * t524 - t487) * MDP(24) + (-t480 + t517) * MDP(26) - t564 + t502 * t565) * t431; MDP(16) * t543 + (-t373 ^ 2 + t558) * MDP(17) + t307 * MDP(18) + (-t321 + t541) * MDP(19) + t368 * MDP(20) + (-t328 * t375 + t447 + t545) * MDP(21) + (t298 * t399 + t328 * t373 - t559) * MDP(22) + (-t326 * t373 - t441 + t545 + 0.2e1 * t556) * MDP(23) + (pkin(5) * t320 - qJ(6) * t321 + (t297 - t299) * t375 + (t296 - t509) * t373) * MDP(24) + (0.2e1 * t547 - t308 * t373 + t326 * t375 + (0.2e1 * qJD(6) - t298) * t399 + t559) * MDP(25) + (t291 * qJ(6) - t292 * pkin(5) - t308 * t326 - t296 * t299 - g(1) * (-pkin(5) * t343 + qJ(6) * t344) - g(2) * (pkin(5) * t345 - qJ(6) * t346) + t478 * t422 + t509 * t297) * MDP(26); (-t368 + t543) * MDP(23) + t307 * MDP(24) + (-t558 - t565) * MDP(25) + (-t297 * t399 + t441 - t556) * MDP(26);];
tau  = t1;
