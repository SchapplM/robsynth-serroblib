% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:10
% EndTime: 2019-03-09 06:21:26
% DurationCPUTime: 8.93s
% Computational Cost: add. (8766->484), mult. (22162->603), div. (0->0), fcn. (16944->8), ass. (0->193)
t495 = cos(pkin(10));
t487 = -pkin(2) * t495 - pkin(1);
t472 = qJD(1) * t487 + qJD(2);
t630 = t472 * MDP(14);
t494 = sin(pkin(10));
t498 = sin(qJ(3));
t500 = cos(qJ(3));
t465 = t494 * t500 + t495 * t498;
t619 = -t494 * t498 + t500 * t495;
t629 = t619 * qJD(1);
t497 = sin(qJ(4));
t499 = cos(qJ(4));
t605 = t465 * qJD(1);
t432 = qJD(3) * t499 - t497 * t605;
t433 = qJD(3) * t497 + t499 * t605;
t496 = sin(qJ(5));
t595 = cos(qJ(5));
t375 = -t595 * t432 + t433 * t496;
t517 = t496 * t432 + t433 * t595;
t628 = t375 * t517;
t413 = pkin(3) * t605 - pkin(8) * t629;
t403 = t499 * t413;
t596 = -pkin(9) - pkin(8);
t541 = qJD(4) * t596;
t582 = t629 * t499;
t594 = pkin(7) + qJ(2);
t473 = t594 * t494;
t466 = qJD(1) * t473;
t474 = t594 * t495;
t467 = qJD(1) * t474;
t608 = -t466 * t500 - t498 * t467;
t627 = pkin(4) * t605 - pkin(9) * t582 - t497 * t608 - t499 * t541 + t403;
t559 = t497 * t413 + t499 * t608;
t583 = t629 * t497;
t626 = -pkin(9) * t583 - t497 * t541 + t559;
t504 = qJD(3) * (qJD(4) + t629);
t551 = qJD(4) * t497;
t503 = -t499 * t504 + t551 * t605;
t531 = t595 * qJD(5);
t459 = t619 * qJD(3);
t505 = qJD(1) * t459;
t550 = qJD(4) * t499;
t542 = qJD(3) * t551 + t497 * t505 + t550 * t605;
t549 = qJD(5) * t496;
t340 = -t432 * t531 + t433 * t549 + t496 * t542 + t595 * t503;
t452 = qJD(4) - t629;
t447 = qJD(5) + t452;
t327 = t375 * t447 - t340;
t341 = t432 * t549 + t433 * t531 - t496 * t503 + t595 * t542;
t460 = t465 * qJD(3);
t448 = qJD(1) * t460;
t597 = t517 ^ 2;
t625 = t448 * MDP(26) + (t447 * t517 - t341) * MDP(25) + MDP(22) * t628 + (-t375 ^ 2 + t597) * MDP(23) + t327 * MDP(24);
t411 = -qJD(3) * pkin(3) - t608;
t370 = -pkin(4) * t432 + t411;
t338 = pkin(5) * t375 - qJ(6) * t517 + t370;
t624 = t338 * t375;
t623 = t370 * t375;
t622 = (t494 ^ 2 + t495 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t572 = t496 * t497;
t515 = t595 * t499 - t572;
t398 = t515 * t629;
t604 = qJD(4) + qJD(5);
t606 = t595 * qJD(4) + t531;
t423 = -t499 * t606 + t572 * t604;
t556 = -t423 - t398;
t539 = t595 * t497;
t469 = t496 * t499 + t539;
t397 = t469 * t629;
t424 = t604 * t469;
t555 = t424 - t397;
t621 = t465 * qJD(2);
t620 = t459 * t497 + t465 * t550;
t618 = MDP(14) * t619;
t352 = pkin(5) * t517 + qJ(6) * t375;
t614 = t605 * MDP(8);
t570 = t497 * t448;
t415 = -pkin(3) * t619 - pkin(8) * t465 + t487;
t408 = t499 * t415;
t429 = -t473 * t498 + t474 * t500;
t579 = t465 * t499;
t360 = -pkin(4) * t619 - pkin(9) * t579 - t429 * t497 + t408;
t422 = t499 * t429;
t557 = t497 * t415 + t422;
t569 = t497 * t465;
t366 = -pkin(9) * t569 + t557;
t612 = t496 * t360 + t595 * t366;
t475 = t596 * t497;
t476 = t596 * t499;
t516 = t475 * t595 + t496 * t476;
t611 = -qJD(5) * t516 + t627 * t496 + t626 * t595;
t431 = t496 * t475 - t476 * t595;
t610 = -qJD(5) * t431 + t626 * t496 - t627 * t595;
t418 = -t498 * t466 + t500 * t467;
t528 = -t418 + (t551 - t583) * pkin(4);
t609 = t499 * t448 - t452 * t551;
t428 = t473 * t500 + t498 * t474;
t607 = -t433 * t551 - t503 * t499;
t443 = t448 * pkin(5);
t391 = -pkin(3) * t629 - pkin(8) * t605 + t472;
t399 = t448 * pkin(3) - pkin(8) * t505;
t396 = t499 * t399;
t412 = qJD(3) * pkin(8) + t418;
t525 = -t412 * t550 + t396;
t507 = t619 * qJD(2);
t380 = qJD(1) * t507 + qJD(3) * t608;
t571 = t497 * t380;
t321 = t448 * pkin(4) + pkin(9) * t503 - t391 * t551 + t525 - t571;
t512 = t499 * t380 + t391 * t550 + t497 * t399 - t412 * t551;
t326 = -pkin(9) * t542 + t512;
t364 = t499 * t391 - t412 * t497;
t354 = -pkin(9) * t433 + t364;
t346 = pkin(4) * t452 + t354;
t365 = t497 * t391 + t499 * t412;
t355 = pkin(9) * t432 + t365;
t529 = -t595 * t321 + t496 * t326 + t346 * t549 + t355 * t531;
t312 = -t443 + t529;
t506 = t338 * t517 + t312;
t603 = -t370 * t517 - t529;
t552 = qJD(3) * t500;
t553 = qJD(3) * t498;
t381 = qJD(1) * t621 - t466 * t553 + t467 * t552;
t602 = t381 * t497 + t411 * t550;
t601 = -t452 ^ 2 * t499 - t570;
t600 = -t340 * t515 - t517 * t555;
t392 = -qJD(3) * t428 + t507;
t414 = pkin(3) * t460 - pkin(8) * t459;
t404 = t499 * t414;
t580 = t459 * t499;
t333 = -pkin(9) * t580 + pkin(4) * t460 - t392 * t497 + t404 + (-t422 + (pkin(9) * t465 - t415) * t497) * qJD(4);
t511 = t499 * t392 + t497 * t414 + t415 * t550 - t429 * t551;
t337 = -pkin(9) * t620 + t511;
t599 = -qJD(5) * t612 + t333 * t595 - t496 * t337;
t540 = t595 * t355;
t323 = t496 * t346 + t540;
t592 = t323 * t447;
t590 = t375 * t605;
t589 = t517 * t605;
t588 = t516 * t448;
t587 = t431 * t448;
t586 = t432 * t605;
t585 = t433 * t452;
t584 = t433 * t605;
t420 = t515 * t448;
t421 = t469 * t448;
t573 = t496 * t355;
t568 = t499 * t432;
t329 = t354 * t595 - t573;
t566 = -pkin(4) * t531 - qJD(6) + t329;
t565 = qJ(6) * t605 + t611;
t564 = -pkin(5) * t605 + t610;
t563 = -pkin(5) * t555 + qJ(6) * t556 + qJD(6) * t469 - t528;
t560 = -t424 * t447 + t420;
t558 = -t423 * t447 + t421;
t322 = t346 * t595 - t573;
t547 = qJD(6) - t322;
t546 = qJD(1) * qJD(2);
t491 = -pkin(4) * t499 - pkin(3);
t537 = t465 * t551;
t530 = -t496 * t321 - t595 * t326 - t346 * t531 + t355 * t549;
t393 = -t473 * t553 + t474 * t552 + t621;
t328 = t496 * t354 + t540;
t527 = pkin(4) * t549 - t328;
t526 = -t469 * t341 - t375 * t556;
t394 = pkin(4) * t569 + t428;
t522 = t452 * t583 + t609;
t437 = t447 * qJD(6);
t440 = t448 * qJ(6);
t311 = t440 + t437 - t530;
t368 = pkin(4) * t620 + t393;
t520 = t360 * t595 - t496 * t366;
t513 = t322 * t447 + t530;
t510 = t496 * t333 + t595 * t337 + t360 * t531 - t366 * t549;
t359 = pkin(4) * t542 + t381;
t490 = -pkin(4) * t595 - pkin(5);
t486 = pkin(4) * t496 + qJ(6);
t419 = t448 * t619;
t416 = -pkin(5) * t515 - qJ(6) * t469 + t491;
t406 = t515 * t465;
t405 = t469 * t465;
t349 = pkin(5) * t405 - qJ(6) * t406 + t394;
t348 = t459 * t539 - t496 * t537 - t549 * t569 + (t459 * t496 + t465 * t606) * t499;
t347 = t424 * t465 - t459 * t515;
t343 = pkin(4) * t433 + t352;
t335 = pkin(5) * t619 - t520;
t334 = -qJ(6) * t619 + t612;
t318 = t447 * qJ(6) + t323;
t317 = -t447 * pkin(5) + t547;
t316 = pkin(5) * t348 + qJ(6) * t347 - qJD(6) * t406 + t368;
t315 = t341 * pkin(5) + t340 * qJ(6) - qJD(6) * t517 + t359;
t314 = -t460 * pkin(5) - t599;
t313 = qJ(6) * t460 - qJD(6) * t619 + t510;
t1 = [(t340 * t405 - t341 * t406 + t347 * t375 - t348 * t517) * MDP(23) + (-t311 * t405 + t312 * t406 - t313 * t375 + t314 * t517 - t317 * t347 - t318 * t348 - t334 * t341 - t335 * t340) * MDP(30) + (-t340 * t406 - t347 * t517) * MDP(22) + 0.2e1 * t546 * t622 + (t448 * t487 + t460 * t472) * MDP(13) + (t459 * t629 - t460 * t605 + t505 * t619) * MDP(9) + ((-t429 * t550 + t404) * t452 + t408 * t448 - t525 * t619 + t364 * t460 - t393 * t432 + t428 * t542 + ((-qJD(4) * t415 - t392) * t452 - t429 * t448 - (-qJD(4) * t391 - t380) * t619 + t411 * t459) * t497) * MDP(20) + (MDP(15) * t607 + MDP(17) * t609 + MDP(20) * t602 + t505 * MDP(8) - t448 * MDP(9)) * t465 + (t452 * t580 + t503 * t619) * MDP(17) + (-t365 * t460 + t381 * t579 - t428 * t503 - t448 * t557 - t452 * t511 + t619 * t512 + (-t537 + t580) * t411) * MDP(21) + (t580 * MDP(15) - MDP(16) * t620 + t460 * MDP(17) + t393 * MDP(21)) * t433 + (t341 * t619 - t348 * t447 - t375 * t460 - t405 * t448) * MDP(25) + (t312 * t619 - t314 * t447 + t315 * t405 + t316 * t375 - t317 * t460 - t335 * t448 + t338 * t348 + t341 * t349) * MDP(29) + (-t311 * t619 + t313 * t447 - t315 * t406 - t316 * t517 + t318 * t460 + t334 * t448 + t338 * t347 + t340 * t349) * MDP(31) + (t340 * t619 - t347 * t447 + t406 * t448 + t460 * t517) * MDP(24) + (-t323 * t460 - t394 * t340 - t370 * t347 + t359 * t406 + t368 * t517 - t447 * t510 - t448 * t612 - t530 * t619) * MDP(28) + (t322 * t460 + t394 * t341 + t370 * t348 + t359 * t405 + t368 * t375 + t447 * t599 + t520 * t448 + t529 * t619) * MDP(27) + (t432 * t460 - t448 * t569 - t452 * t620 + t542 * t619) * MDP(18) + (t452 * t460 - t419) * MDP(19) + (t447 * t460 - t419) * MDP(26) + (-t432 * t537 + t459 * t568 + t503 * t569 - t542 * t579) * MDP(16) + t459 * t614 + ((t487 * t629 - t392) * MDP(14) - t393 * MDP(13) + t459 * MDP(10) - t460 * MDP(11)) * qJD(3) + t459 * t630 + (t311 * t334 + t312 * t335 + t313 * t318 + t314 * t317 + t315 * t349 + t316 * t338) * MDP(32); (t522 + t586) * MDP(20) + (-t584 + t601) * MDP(21) + (t560 - t590) * MDP(27) + (-t421 - t589) * MDP(28) + (t420 - t590) * MDP(29) + (t526 - t600) * MDP(30) + (t558 + t589) * MDP(31) + (t311 * t469 - t312 * t515 + t317 * t555 + t318 * t556 - t338 * t605) * MDP(32) + (t397 * MDP(27) - MDP(28) * t556 - MDP(29) * t555 - t398 * MDP(31)) * t447 + (t605 * MDP(13) + t629 * MDP(14) + (MDP(13) * t465 + t618) * qJD(1)) * qJD(3) - qJD(1) ^ 2 * t622; (qJD(3) * t418 - t381) * MDP(13) - t546 * t618 + (t497 * t504 + t585) * t499 * MDP(15) + (t432 * t550 + t433 * t583 - t497 * t542 + t607) * MDP(16) + (-t584 - t601) * MDP(17) + (t522 - t586) * MDP(18) + (-pkin(8) * t570 - pkin(3) * t542 - t381 * t499 + t418 * t432 + (-pkin(8) * t550 - t403 + (t411 + t608) * t497) * t452) * MDP(20) + (pkin(3) * t503 - pkin(8) * t609 - t411 * t582 - t418 * t433 + t452 * t559 + t602) * MDP(21) + (-t340 * t469 + t517 * t556) * MDP(22) + (t526 + t600) * MDP(23) + (-t398 * t447 + t558 - t589) * MDP(24) + (t397 * t447 + t560 + t590) * MDP(25) + (t491 * t341 - t359 * t515 + t555 * t370 + t528 * t375 + t447 * t610 + t588) * MDP(27) + (-t491 * t340 + t359 * t469 + t556 * t370 + t447 * t611 + t528 * t517 - t587) * MDP(28) + (-t315 * t515 + t338 * t555 + t341 * t416 - t375 * t563 + t447 * t564 + t588) * MDP(29) + (t311 * t515 + t312 * t469 + t317 * t556 - t318 * t555 + t340 * t516 - t341 * t431 + t375 * t565 - t517 * t564) * MDP(30) + (-t315 * t469 - t338 * t556 + t340 * t416 - t447 * t565 + t517 * t563 + t587) * MDP(31) + (t311 * t431 - t312 * t516 + t315 * t416 - t317 * t564 - t318 * t565 - t338 * t563) * MDP(32) + (-MDP(16) * t568 - MDP(9) * t629 - t614 - t630) * t629 + (-qJD(4) * t497 ^ 2 * MDP(15) - t472 * MDP(13) - t452 * MDP(19) - t364 * MDP(20) + MDP(21) * t365 - t447 * MDP(26) - t322 * MDP(27) + t323 * MDP(28) + t317 * MDP(29) - t318 * MDP(31) + MDP(9) * t605) * t605; -t433 * t432 * MDP(15) + (-t432 ^ 2 + t433 ^ 2) * MDP(16) + (-t432 * t452 - t503) * MDP(17) + (-t542 + t585) * MDP(18) + t448 * MDP(19) + (-t411 * t433 + t396 - t571 + (-qJD(4) + t452) * t365) * MDP(20) + (t364 * t452 - t411 * t432 - t512) * MDP(21) + (t328 * t447 + (-t375 * t433 - t447 * t549 + t448 * t595) * pkin(4) + t603) * MDP(27) + (t329 * t447 + t623 + (-t433 * t517 - t447 * t531 - t496 * t448) * pkin(4) + t530) * MDP(28) + (-t343 * t375 - t447 * t527 - t448 * t490 - t506) * MDP(29) + (-t340 * t490 - t341 * t486 + (t318 + t527) * t517 + (t317 + t566) * t375) * MDP(30) + (t343 * t517 - t447 * t566 + t448 * t486 + t311 - t624) * MDP(31) + (t311 * t486 + t312 * t490 + t317 * t527 - t318 * t566 - t338 * t343) * MDP(32) + t625; (t592 + t603) * MDP(27) + (t513 + t623) * MDP(28) + (-t352 * t375 + t443 - t506 + t592) * MDP(29) + (pkin(5) * t340 - qJ(6) * t341 + (t318 - t323) * t517 + (t317 - t547) * t375) * MDP(30) + (t352 * t517 + 0.2e1 * t437 + 0.2e1 * t440 - t513 - t624) * MDP(31) + (-pkin(5) * t312 + qJ(6) * t311 - t317 * t323 + t318 * t547 - t338 * t352) * MDP(32) + t625; (-qJD(3) * t605 + t628) * MDP(29) + t327 * MDP(30) + (-t447 ^ 2 - t597) * MDP(31) + (-t318 * t447 + t506) * MDP(32);];
tauc  = t1;
