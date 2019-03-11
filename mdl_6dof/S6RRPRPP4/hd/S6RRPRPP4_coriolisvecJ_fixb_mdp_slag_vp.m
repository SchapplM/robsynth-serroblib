% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:48
% EndTime: 2019-03-09 10:01:57
% DurationCPUTime: 4.46s
% Computational Cost: add. (4895->445), mult. (10929->579), div. (0->0), fcn. (6574->6), ass. (0->205)
t552 = pkin(3) + pkin(7);
t453 = cos(qJ(4));
t451 = sin(qJ(4));
t519 = qJD(2) * t451;
t454 = cos(qJ(2));
t520 = qJD(1) * t454;
t397 = t453 * t520 + t519;
t493 = t451 * t520;
t517 = qJD(2) * t453;
t399 = -t493 + t517;
t449 = sin(pkin(9));
t450 = cos(pkin(9));
t348 = t450 * t397 + t399 * t449;
t452 = sin(qJ(2));
t521 = qJD(1) * t452;
t429 = qJD(4) + t521;
t568 = t348 * t429;
t513 = qJD(4) * t453;
t494 = t450 * t513;
t499 = t453 * t521;
t514 = qJD(4) * t451;
t539 = t449 * t451;
t526 = -t449 * t514 + t450 * t499 - t521 * t539 + t494;
t475 = t449 * t453 + t450 * t451;
t525 = t449 * t513 + t450 * t514 + t475 * t521;
t436 = pkin(7) * t521;
t567 = qJD(3) + t436;
t506 = qJD(1) * qJD(2);
t492 = t452 * t506;
t418 = t451 * t492;
t359 = -qJD(4) * t397 + t418;
t455 = -pkin(2) - pkin(8);
t490 = -qJ(3) * t452 - pkin(1);
t393 = t454 * t455 + t490;
t369 = t393 * qJD(1);
t509 = pkin(3) * t521 + t567;
t371 = qJD(2) * t455 + t509;
t335 = t369 * t453 + t371 * t451;
t428 = pkin(2) * t492;
t480 = pkin(8) * t452 - qJ(3) * t454;
t515 = qJD(3) * t452;
t461 = qJD(2) * t480 - t515;
t356 = qJD(1) * t461 + t428;
t491 = t454 * t506;
t427 = pkin(7) * t491;
t391 = pkin(3) * t491 + t427;
t484 = -t356 * t451 + t453 * t391;
t459 = -qJD(4) * t335 + t484;
t293 = pkin(4) * t491 - qJ(5) * t359 - qJD(5) * t399 + t459;
t505 = qJD(2) * qJD(4);
t360 = -qJD(4) * t493 + (-t492 + t505) * t453;
t502 = -t453 * t356 - t371 * t513 - t451 * t391;
t465 = -t369 * t514 - t502;
t296 = -qJ(5) * t360 - qJD(5) * t397 + t465;
t283 = t449 * t293 + t450 * t296;
t334 = -t369 * t451 + t453 * t371;
t324 = -qJ(5) * t399 + t334;
t320 = pkin(4) * t429 + t324;
t325 = -qJ(5) * t397 + t335;
t548 = t325 * t449;
t297 = t320 * t450 - t548;
t322 = t450 * t325;
t298 = t449 * t320 + t322;
t474 = -t450 * t453 + t539;
t531 = -t450 * t293 + t449 * t296;
t566 = t283 * t475 - t297 * t525 + t298 * t526 + t474 * t531;
t503 = qJ(6) * t491 + t283;
t510 = qJD(6) * t429;
t280 = t503 + t510;
t281 = -pkin(5) * t491 + t531;
t291 = -pkin(5) * t429 + qJD(6) - t297;
t292 = qJ(6) * t429 + t298;
t565 = t280 * t475 + t281 * t474 + t291 * t525 + t292 * t526;
t476 = -t397 * t449 + t450 * t399;
t564 = t476 ^ 2;
t563 = -0.2e1 * t506;
t447 = t452 ^ 2;
t448 = t454 ^ 2;
t562 = (t447 - t448) * MDP(5);
t437 = pkin(7) * t520;
t405 = pkin(3) * t520 + t437;
t446 = qJD(2) * qJ(3);
t385 = t446 + t405;
t353 = pkin(4) * t397 + qJD(5) + t385;
t309 = pkin(5) * t348 - qJ(6) * t476 + t353;
t557 = t309 * t476;
t441 = pkin(2) * t521;
t378 = qJD(1) * t480 + t441;
t483 = -t378 * t451 + t453 * t405;
t538 = t451 * t452;
t331 = (pkin(4) * t454 - qJ(5) * t538) * qJD(1) + t483;
t524 = t453 * t378 + t451 * t405;
t337 = qJ(5) * t499 + t524;
t532 = qJ(5) - t455;
t488 = t532 * t453;
t372 = -qJD(4) * t488 - qJD(5) * t451;
t462 = -qJD(5) * t453 + t514 * t532;
t527 = (-t331 + t462) * t450 + (t337 - t372) * t449;
t416 = t552 * t452;
t523 = t453 * t393 + t451 * t416;
t512 = qJD(4) * t454;
t496 = t451 * t512;
t500 = pkin(4) * t453 + pkin(3);
t518 = qJD(2) * t452;
t556 = (-pkin(7) - t500) * t518 - pkin(4) * t496;
t555 = t452 * t517 + t496;
t554 = MDP(22) + MDP(25);
t553 = pkin(4) * t513 + t500 * t521 + t567;
t551 = pkin(4) * t360;
t550 = qJD(2) * pkin(2);
t300 = t324 * t449 + t322;
t549 = t300 * t476;
t547 = t359 * t453;
t445 = qJD(2) * qJD(3);
t404 = t552 * t518;
t473 = qJD(1) * t404;
t375 = t445 - t473;
t546 = t375 * t451;
t545 = t375 * t453;
t544 = t397 * t429;
t543 = t399 * t429;
t542 = t399 * t454;
t541 = t429 * t453;
t540 = t429 * t455;
t456 = qJD(2) ^ 2;
t537 = t452 * t456;
t536 = t453 * t454;
t535 = t454 * t456;
t457 = qJD(1) ^ 2;
t534 = t454 * t457;
t533 = t451 * pkin(4) + qJ(3);
t307 = t449 * t331 + t450 * t337;
t302 = qJ(6) * t520 + t307;
t339 = t450 * t372 + t449 * t462;
t530 = -t302 + t339;
t529 = -pkin(5) * t520 + t527;
t440 = pkin(2) * t518;
t362 = t440 + t461;
t516 = qJD(2) * t454;
t406 = t552 * t516;
t390 = t453 * t406;
t485 = qJ(5) * t454 - t393;
t511 = qJD(5) * t454;
t305 = pkin(4) * t516 + t390 + t485 * t513 + (-qJ(5) * t518 - qJD(4) * t416 - t362 + t511) * t451;
t464 = t453 * t362 - t393 * t514 + t451 * t406 + t416 * t513;
t310 = qJ(5) * t555 - t453 * t511 + t464;
t287 = t449 * t305 + t450 * t310;
t528 = -t526 * pkin(5) - t525 * qJ(6) - qJD(6) * t474 - t553;
t401 = t453 * t416;
t341 = pkin(4) * t452 + t451 * t485 + t401;
t345 = -qJ(5) * t536 + t523;
t314 = t449 * t341 + t450 * t345;
t417 = t552 * t454;
t410 = -pkin(2) * t454 + t490;
t386 = qJD(1) * t410;
t301 = t324 * t450 - t548;
t507 = qJD(6) - t301;
t504 = t452 * t534;
t501 = pkin(4) * t536 + t417;
t498 = t451 * t518;
t495 = t453 * t512;
t489 = t445 + t551;
t487 = pkin(1) * t563;
t486 = qJD(3) - t550;
t327 = t359 * t449 + t450 * t360;
t328 = t359 * t450 - t360 * t449;
t408 = t532 * t451;
t354 = -t408 * t449 + t450 * t488;
t355 = -t450 * t408 - t449 * t488;
t481 = -t355 * t327 + t328 * t354 - t339 * t348;
t286 = t305 * t450 - t310 * t449;
t313 = t341 * t450 - t345 * t449;
t472 = -qJD(1) * t448 + t429 * t452;
t471 = -0.2e1 * qJD(2) * t386;
t470 = t429 * t451;
t466 = -qJ(3) * t516 - t515;
t367 = qJD(1) * t466 + t428;
t380 = t440 + t466;
t469 = pkin(7) * t456 + qJD(1) * t380 + t367;
t468 = t385 * t452 + t455 * t516;
t374 = t475 * t454;
t460 = pkin(5) * t327 - qJ(6) * t328 - qJD(6) * t476 + t489;
t407 = pkin(7) * t492 - t445;
t409 = t436 + t486;
t415 = -t437 - t446;
t458 = -t407 * t454 + (t409 * t454 + (t415 + t437) * t452) * qJD(2);
t435 = -pkin(4) * t450 - pkin(5);
t431 = pkin(4) * t449 + qJ(6);
t419 = t453 * t491;
t402 = -qJ(3) * t520 + t441;
t373 = t474 * t454;
t370 = t386 * t521;
t346 = pkin(5) * t475 + qJ(6) * t474 + t533;
t343 = -t449 * t555 - t450 * t498 + t454 * t494;
t342 = qJD(4) * t374 - t474 * t518;
t340 = t375 + t551;
t329 = -pkin(5) * t373 + qJ(6) * t374 + t501;
t316 = pkin(4) * t399 + pkin(5) * t476 + qJ(6) * t348;
t312 = -pkin(5) * t452 - t313;
t311 = qJ(6) * t452 + t314;
t299 = -pkin(5) * t342 + qJ(6) * t343 + qJD(6) * t374 + t556;
t288 = -t473 + t460;
t285 = -pkin(5) * t516 - t286;
t284 = qJ(6) * t516 + qJD(6) * t452 + t287;
t1 = [0.2e1 * t452 * MDP(4) * t491 + t562 * t563 + MDP(6) * t535 - MDP(7) * t537 + (-pkin(7) * t535 + t452 * t487) * MDP(9) + (pkin(7) * t537 + t454 * t487) * MDP(10) + t458 * MDP(11) + (t452 * t471 + t454 * t469) * MDP(12) + (-t452 * t469 + t454 * t471) * MDP(13) + (pkin(7) * t458 + t367 * t410 + t380 * t386) * MDP(14) + (-t359 * t451 * t454 + (-t495 + t498) * t399) * MDP(15) + ((-t397 * t451 + t399 * t453) * t518 + (-t547 + t360 * t451 + (t397 * t453 + t399 * t451) * qJD(4)) * t454) * MDP(16) + (-t429 * t495 + t359 * t452 + (t451 * t472 + t542) * qJD(2)) * MDP(17) + (t429 * t496 - t360 * t452 + (-t397 * t454 + t453 * t472) * qJD(2)) * MDP(18) + (t429 + t521) * MDP(19) * t516 + ((-t362 * t451 + t390) * t429 - t404 * t397 + t417 * t360 + (-t385 * t517 + t484) * t452 + (-t335 * t452 - t429 * t523) * qJD(4) + (-t385 * t514 + t545 + ((-t393 * t451 + t401) * qJD(1) + t334) * qJD(2)) * t454) * MDP(20) + (-t464 * t429 - t404 * t399 + t417 * t359 + ((qJD(2) * t385 + qJD(4) * t369) * t451 + t502) * t452 + (-t385 * t513 - t546 + (-qJD(1) * t523 - t335) * qJD(2)) * t454) * MDP(21) + (t283 * t373 - t286 * t476 - t287 * t348 + t297 * t343 + t298 * t342 - t313 * t328 - t314 * t327 - t374 * t531) * MDP(22) + (t283 * t314 + t297 * t286 + t298 * t287 - t313 * t531 + t340 * t501 + t353 * t556) * MDP(23) + (-t281 * t452 - t285 * t429 - t288 * t373 + t299 * t348 - t309 * t342 + t327 * t329 + (-qJD(1) * t312 - t291) * t516) * MDP(24) + (t280 * t373 - t281 * t374 - t284 * t348 + t285 * t476 - t291 * t343 + t292 * t342 - t311 * t327 + t312 * t328) * MDP(25) + (t280 * t452 + t284 * t429 + t288 * t374 - t299 * t476 + t309 * t343 - t328 * t329 + (qJD(1) * t311 + t292) * t516) * MDP(26) + (t280 * t311 + t281 * t312 + t284 * t292 + t285 * t291 + t288 * t329 + t299 * t309) * MDP(27); -MDP(4) * t504 + t457 * t562 + ((-t415 - t446) * t452 + (-t409 + t486) * t454) * qJD(1) * MDP(11) + (-t402 * t520 + t370) * MDP(12) + (0.2e1 * t445 + (t386 * t454 + t402 * t452) * qJD(1)) * MDP(13) + (-qJ(3) * t407 - qJD(3) * t415 - t386 * t402 + (-t415 * t452 + (-t409 - t550) * t454) * qJD(1) * pkin(7)) * MDP(14) + (-t399 * t470 + t547) * MDP(15) + ((-t360 - t543) * t453 + (-t359 + t544) * t451) * MDP(16) + (-t429 * t514 + t419 + (-t429 * t538 - t542) * qJD(1)) * MDP(17) + (-t429 * t513 + (-t452 * t541 + (t397 - t519) * t454) * qJD(1)) * MDP(18) - t429 * MDP(19) * t520 + (qJ(3) * t360 + t546 - t483 * t429 + t509 * t397 + (t385 * t453 - t451 * t540) * qJD(4) + (-t334 * t454 + t453 * t468) * qJD(1)) * MDP(20) + (qJ(3) * t359 + t545 + t524 * t429 + t509 * t399 + (-t385 * t451 - t453 * t540) * qJD(4) + (t335 * t454 - t451 * t468) * qJD(1)) * MDP(21) + (t307 * t348 - t476 * t527 + t481 - t566) * MDP(22) + (t283 * t355 + t531 * t354 + t340 * t533 + t553 * t353 + (t339 - t307) * t298 + t527 * t297) * MDP(23) + (t288 * t475 + t327 * t346 + t529 * t429 - t528 * t348 + t526 * t309 + (-qJD(2) * t354 + t291) * t520) * MDP(24) + (t302 * t348 - t476 * t529 + t481 - t565) * MDP(25) + (t288 * t474 - t328 * t346 + t530 * t429 + t528 * t476 + t525 * t309 + (qJD(2) * t355 - t292) * t520) * MDP(26) + (t280 * t355 + t281 * t354 + t288 * t346 - t291 * t529 + t292 * t530 - t309 * t528) * MDP(27) + (MDP(9) * t452 * t457 + MDP(10) * t534) * pkin(1); MDP(12) * t504 + (-t447 * t457 - t456) * MDP(13) + (t370 + t427) * MDP(14) + t419 * MDP(20) + t566 * MDP(23) + t565 * MDP(27) + (t415 * MDP(14) - t397 * MDP(20) + (-t399 - t493) * MDP(21) - t353 * MDP(23) + (-t474 * t520 - t348) * MDP(24) + (t475 * t520 + t476) * MDP(26) - t309 * MDP(27)) * qJD(2) + (-MDP(20) * t470 - MDP(21) * t541 + MDP(26) * t526) * t429 + (-t429 * MDP(24) + t476 * t554) * t525 + t554 * (-t475 * t327 + t328 * t474 - t526 * t348); t399 * t397 * MDP(15) + (-t397 ^ 2 + t399 ^ 2) * MDP(16) + (-t451 * t505 + t418 + t544) * MDP(17) + (-t360 + t543) * MDP(18) + (t335 * t429 - t385 * t399 + t459) * MDP(20) + (t334 * t429 + t385 * t397 - t465) * MDP(21) + (t298 * t476 - t549) * MDP(22) + (t297 * t300 - t298 * t301) * MDP(23) + (t300 * t429 - t531 - t557) * MDP(24) + (t292 * t476 - t327 * t431 + t328 * t435 - t549) * MDP(25) + (-t301 * t429 + t316 * t476 + t503 + 0.2e1 * t510) * MDP(26) + (t280 * t431 + t281 * t435 - t291 * t300 + t292 * t507 - t309 * t316) * MDP(27) + (-MDP(17) * t513 + (MDP(19) + (pkin(5) - t435) * MDP(24) + t431 * MDP(26)) * qJD(2)) * t520 + ((-t327 * t449 - t328 * t450) * MDP(22) + (t283 * t449 - t353 * t399 - t450 * t531) * MDP(23)) * pkin(4) + ((-t297 + t301) * MDP(22) - t316 * MDP(24) + (t291 - t507) * MDP(25) - t309 * MDP(26)) * t348; (t297 * t476 + t298 * t348 + t489) * MDP(23) + (t429 * t476 + t327) * MDP(24) + (-t328 + t568) * MDP(26) + (-t291 * t476 + t292 * t348 + t460) * MDP(27) + (-MDP(23) - MDP(27)) * t492 * t552 + t554 * (-t348 ^ 2 - t564); (t348 * t476 - t491) * MDP(24) + (t328 + t568) * MDP(25) + (-t429 ^ 2 - t564) * MDP(26) + (-t292 * t429 + t281 + t557) * MDP(27);];
tauc  = t1;
