% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:21:11
% EndTime: 2019-03-09 09:21:23
% DurationCPUTime: 7.12s
% Computational Cost: add. (3565->525), mult. (9093->703), div. (0->0), fcn. (6481->8), ass. (0->212)
t448 = sin(pkin(6));
t453 = sin(qJ(2));
t550 = qJD(1) * t453;
t521 = t448 * t550;
t456 = cos(qJ(2));
t449 = cos(pkin(6));
t537 = t449 * qJD(1);
t531 = pkin(1) * t537;
t554 = pkin(8) * t521 - t456 * t531;
t586 = qJD(3) + t554;
t589 = -qJ(4) * t521 + t586;
t409 = qJD(5) + t521;
t454 = cos(qJ(6));
t430 = qJD(2) + t537;
t452 = sin(qJ(5));
t455 = cos(qJ(5));
t551 = qJD(1) * t448;
t520 = t456 * t551;
t377 = t452 * t430 + t455 * t520;
t451 = sin(qJ(6));
t570 = t377 * t451;
t342 = -t454 * t409 - t570;
t376 = -t430 * t455 + t452 * t520;
t533 = qJD(6) - t376;
t588 = t342 * t533;
t482 = t377 * t454 - t409 * t451;
t587 = t482 * t533;
t530 = pkin(1) * qJD(2) * t449;
t426 = t456 * t530;
t439 = t449 * qJD(3);
t548 = qJD(2) * t453;
t519 = t448 * t548;
t583 = (pkin(8) * t548 + qJD(4) * t456) * t448;
t339 = -qJ(4) * t519 - t426 - t439 + t583;
t501 = qJD(1) * t530;
t407 = t456 * t501;
t410 = t430 * qJD(3);
t532 = qJD(1) * qJD(2);
t512 = t448 * t532;
t498 = t453 * t512;
t525 = qJ(4) * t498 + t407 + t410;
t326 = qJD(1) * t583 - t525;
t446 = t453 ^ 2;
t447 = t456 ^ 2;
t585 = MDP(5) * (t446 - t447);
t584 = t482 * t409;
t564 = t448 * t456;
t565 = t448 * t453;
t380 = -t448 * pkin(1) - pkin(2) * t564 - qJ(3) * t565;
t364 = pkin(3) * t564 - t380;
t477 = (pkin(4) * t453 + pkin(9) * t456) * t448;
t336 = t477 + t364;
t432 = pkin(8) * t565;
t494 = -qJ(4) * t565 + t432;
t523 = -pkin(1) * t456 - pkin(2);
t504 = -pkin(3) + t523;
t346 = (-pkin(9) + t504) * t449 + t494;
t582 = t452 * t336 + t455 * t346;
t578 = pkin(1) * t453;
t436 = t449 * t578;
t581 = pkin(8) * t564 + t436;
t412 = t456 * t512;
t457 = -pkin(2) - pkin(3);
t445 = -pkin(9) + t457;
t323 = t430 * t445 + t589;
t371 = -pkin(1) * t551 - pkin(2) * t520 - qJ(3) * t521;
t352 = pkin(3) * t520 + qJD(4) - t371;
t325 = qJD(1) * t477 + t352;
t306 = t323 * t455 + t325 * t452;
t464 = t448 * (pkin(4) * t456 + t445 * t453);
t460 = qJD(2) * t464;
t427 = qJD(3) * t565;
t556 = qJ(3) * t412 + qJD(1) * t427;
t320 = qJD(1) * t460 + t556;
t378 = pkin(8) * t412 + t453 * t501;
t546 = qJD(2) * t456;
t522 = qJ(4) * t546;
t545 = qJD(4) * t453;
t338 = (-t522 - t545) * t551 + t378;
t459 = -qJD(5) * t306 + t455 * t320 - t452 * t338;
t296 = -pkin(5) * t412 - t459;
t580 = t533 * (-pkin(5) * t377 + pkin(10) * t533) + t296;
t518 = t448 * t546;
t555 = qJ(3) * t518 + t427;
t324 = t460 + t555;
t351 = -t448 * t545 + (t436 + (pkin(8) - qJ(4)) * t564) * qJD(2);
t579 = -qJD(5) * t582 + t324 * t455 - t351 * t452;
t450 = qJ(3) + pkin(4);
t577 = pkin(8) * qJD(2);
t517 = qJD(5) * t564;
t497 = qJD(1) * t517;
t540 = qJD(5) * t455;
t347 = -t430 * t540 + t452 * t497 + t455 * t498;
t538 = qJD(6) * t454;
t526 = t454 * t347 + t409 * t538 + t451 * t412;
t539 = qJD(6) * t451;
t308 = t377 * t539 + t526;
t576 = t308 * t451;
t575 = t482 * t454;
t574 = t533 * t451;
t573 = t533 * t454;
t572 = t376 * t409;
t571 = t377 * t409;
t569 = t554 * t430;
t568 = t409 * t452;
t567 = t409 * t455;
t444 = t448 ^ 2;
t566 = t444 * qJD(1) ^ 2;
t542 = qJD(5) * t452;
t348 = t430 * t542 - t452 * t498 + t455 * t497;
t563 = t451 * t348;
t562 = t452 * t308;
t561 = t453 * t455;
t560 = t454 * t348;
t559 = -t589 + t409 * (pkin(5) * t452 - pkin(10) * t455);
t416 = qJ(3) * t520;
t334 = qJD(1) * t464 + t416;
t384 = pkin(8) * t520 + t453 * t531;
t368 = -qJ(4) * t520 + t384;
t557 = t452 * t334 + t455 * t368;
t552 = qJ(3) * qJD(2);
t549 = qJD(2) * t452;
t547 = qJD(2) * t455;
t544 = qJD(5) * t409;
t543 = qJD(5) * t451;
t541 = qJD(5) * t454;
t534 = -qJD(4) - t352;
t529 = t453 * t568;
t528 = t409 * t561;
t527 = t456 * t566;
t379 = t449 * qJ(3) + t581;
t516 = t409 * t542;
t515 = t409 * t540;
t513 = t444 * t532;
t468 = t452 * t320 - t323 * t542 + t325 * t540 + t455 * t338;
t295 = pkin(10) * t412 + t468;
t302 = -pkin(5) * t348 - pkin(10) * t347 - t326;
t510 = -t295 * t451 + t454 * t302;
t508 = t347 * t451 - t454 * t412;
t507 = t534 * t453;
t399 = pkin(5) * t455 + pkin(10) * t452 + t450;
t505 = -pkin(10) * t520 + qJD(6) * t399 - t557;
t502 = t457 * t565;
t500 = t453 * t527;
t496 = -0.2e1 * pkin(1) * t513;
t493 = t384 * t430 - t378;
t373 = (t451 * t456 + t454 * t561) * t551;
t492 = t454 * t540 + t373;
t491 = qJD(2) * t502;
t490 = t295 * t454 + t302 * t451;
t304 = pkin(10) * t409 + t306;
t413 = t430 * qJ(3);
t349 = -t368 - t413;
t333 = pkin(4) * t430 - t349;
t312 = -pkin(5) * t376 + pkin(10) * t377 + t333;
t298 = t304 * t454 + t312 * t451;
t489 = t304 * t451 - t312 * t454;
t311 = pkin(10) * t565 + t582;
t363 = qJ(4) * t564 - t379;
t353 = t449 * pkin(4) - t363;
t386 = -t449 * t455 + t452 * t564;
t387 = -t449 * t452 - t455 * t564;
t317 = -pkin(5) * t386 - pkin(10) * t387 + t353;
t488 = t311 * t454 + t317 * t451;
t487 = -t311 * t451 + t317 * t454;
t305 = -t323 * t452 + t325 * t455;
t485 = t334 * t455 - t368 * t452;
t484 = t336 * t455 - t346 * t452;
t385 = t581 * qJD(2);
t481 = -t378 * t449 - t385 * t430;
t479 = -pkin(8) * t519 + t426;
t478 = t430 * t520 - t412;
t475 = -t387 * t451 + t454 * t565;
t361 = t387 * t454 + t451 * t565;
t474 = -t333 * t453 - t445 * t546;
t473 = -t533 * t538 + t563;
t472 = t533 * t539 + t560;
t470 = t409 * t342;
t469 = -pkin(8) * t498 + t407;
t467 = t452 * t324 + t336 * t540 - t346 * t542 + t455 * t351;
t463 = qJD(5) * t342 + t473;
t462 = -qJD(5) * t482 + t472;
t303 = -pkin(5) * t409 - t305;
t461 = pkin(10) * t348 + (t303 + t305) * t533;
t388 = t430 * t521;
t382 = pkin(2) * t521 - t416;
t381 = t449 * t523 + t432;
t375 = t439 + t479;
t372 = t451 * t455 * t521 - t454 * t520;
t369 = pkin(2) * t519 - t555;
t367 = qJD(1) * t502 + t416;
t365 = t413 + t384;
t362 = -pkin(2) * t430 + t586;
t359 = qJD(5) * t386 + t455 * t519;
t358 = -t449 * t542 + t452 * t519 - t455 * t517;
t357 = t410 + t469;
t355 = pkin(2) * t498 - t556;
t354 = t449 * t504 + t494;
t350 = t491 + t555;
t337 = qJD(1) * t491 + t556;
t331 = t430 * t457 + t589;
t316 = qJD(6) * t475 + t359 * t454 + t451 * t518;
t315 = qJD(6) * t361 + t359 * t451 - t454 * t518;
t313 = -pkin(5) * t520 - t485;
t310 = -pkin(5) * t565 - t484;
t309 = -qJD(6) * t482 + t508;
t307 = pkin(5) * t358 - pkin(10) * t359 - t339;
t300 = -pkin(5) * t518 - t579;
t299 = pkin(10) * t518 + t467;
t294 = -qJD(6) * t298 + t510;
t293 = -qJD(6) * t489 + t490;
t1 = [0.2e1 * (t453 * t456 * MDP(4) - t585) * t513 + (t579 * t409 + t339 * t376 - t353 * t348 + t326 * t386 + t333 * t358 + (t459 * t453 + (qJD(1) * t484 + t305) * t546) * t448) * MDP(24) + (-t358 * t409 + (t348 * t453 + (qJD(1) * t386 + t376) * t546) * t448) * MDP(22) + (t357 * t449 + t375 * t430 + (-t371 * t546 - t355 * t453 + (-t369 * t453 - t380 * t546) * qJD(1)) * t448) * MDP(13) + (-t326 * t449 - t339 * t430 + (t352 * t546 + t337 * t453 + (t350 * t453 + t364 * t546) * qJD(1)) * t448) * MDP(15) + ((t371 * t548 - t355 * t456 + (-t369 * t456 + t380 * t548) * qJD(1)) * t448 + t481) * MDP(11) + (t338 * t449 + t351 * t430 + (t352 * t548 - t337 * t456 + (-t350 * t456 + t364 * t548) * qJD(1)) * t448) * MDP(16) + (t359 * t409 + (t347 * t453 + (qJD(1) * t387 - t377) * t546) * t448) * MDP(21) + ((-qJD(6) * t488 - t299 * t451 + t307 * t454) * t533 - t487 * t348 - t294 * t386 - t489 * t358 + t300 * t342 + t310 * t309 - t296 * t475 + t303 * t315) * MDP(31) + (t309 * t386 - t315 * t533 - t342 * t358 - t348 * t475) * MDP(29) + (t348 * t386 + t358 * t533) * MDP(30) + (-t467 * t409 + t339 * t377 + t353 * t347 - t326 * t387 + t333 * t359 + (-t468 * t453 + (-qJD(1) * t582 - t306) * t546) * t448) * MDP(25) + (t308 * t361 - t316 * t482) * MDP(26) + (t308 * t475 - t309 * t361 + t315 * t482 - t316 * t342) * MDP(27) + (-(qJD(6) * t487 + t299 * t454 + t307 * t451) * t533 + t488 * t348 + t293 * t386 - t298 * t358 - t300 * t482 + t310 * t308 + t296 * t361 + t303 * t316) * MDP(32) + (-t308 * t386 + t316 * t533 - t348 * t361 - t358 * t482) * MDP(28) + (MDP(6) * t518 - MDP(7) * t519) * (t430 + t537) + (t409 * t448 + t444 * t550) * MDP(23) * t546 + (t326 * t456 - t338 * t453 + (-t331 * t456 - t349 * t453) * qJD(2) + (t339 * t456 - t351 * t453 + (-t354 * t456 - t363 * t453) * qJD(2)) * qJD(1)) * t448 * MDP(17) + (-t430 * t479 - t449 * t469 + t456 * t496) * MDP(10) + (t453 * t496 + t481) * MDP(9) + (t326 * t363 + t331 * t351 + t337 * t364 + t338 * t354 + t339 * t349 + t350 * t352) * MDP(18) + (t355 * t380 + t357 * t379 + t362 * t385 + t365 * t375 + t369 * t371 + t378 * t381) * MDP(14) + (t347 * t387 - t359 * t377) * MDP(19) + (t347 * t386 + t348 * t387 + t358 * t377 + t359 * t376) * MDP(20) + (t357 * t456 + t378 * t453 + (t362 * t456 - t365 * t453) * qJD(2) + (t375 * t456 + t385 * t453 + (-t379 * t453 + t381 * t456) * qJD(2)) * qJD(1)) * t448 * MDP(12); (t388 - t498) * MDP(7) + t493 * MDP(11) + t516 * MDP(22) - t515 * MDP(21) + (-t368 * t430 + t378) * MDP(16) + (t566 * t578 + t493) * MDP(9) + ((-t347 - t572) * t455 + (-t348 - t571) * t452) * MDP(20) + (-t347 * t452 + t377 * t567) * MDP(19) - MDP(4) * t500 + (t342 * t373 - t482 * t372 + (t342 * t454 - t451 * t482) * t540 + (t576 + t309 * t454 + (-t342 * t451 - t575) * qJD(6)) * t452) * MDP(27) + (-t454 * t562 - (t452 * t539 - t492) * t482) * MDP(26) - t478 * MDP(6) + (t308 * t455 - t492 * t533 + (t472 + t584) * t452) * MDP(28) + (-t399 * t560 - t303 * t372 - t313 * t342 - (t451 * t505 + t454 * t559) * t533 + (-t303 * t543 + t445 * t463 + t294) * t455 + (t489 * t521 - t303 * t538 - t296 * t451 + t445 * t309 + (t445 * t574 + t489) * qJD(5)) * t452) * MDP(31) + (t399 * t563 - t303 * t373 + t313 * t482 - (-t451 * t559 + t454 * t505) * t533 + (-t303 * t541 + t445 * t462 - t293) * t455 + (t298 * t521 + t303 * t539 - t296 * t454 + t445 * t308 + (t445 * t573 + t298) * qJD(5)) * t452) * MDP(32) + (-t348 * t455 - t533 * t568) * MDP(30) + (-t309 * t455 - (-t451 * t540 - t372) * t533 + (t470 - t473) * t452) * MDP(29) + (t430 * t589 + t525) * MDP(15) + (t450 * t347 + t326 * t452 + t557 * t409 - t589 * t377 + (-t333 * t455 + t445 * t568) * qJD(5)) * MDP(25) + (-t450 * t348 - t326 * t455 - t485 * t409 - t589 * t376 + (-t333 * t452 - t445 * t567) * qJD(5)) * MDP(24) + (-qJ(3) * t326 - t331 * t368 + t338 * t457 - t349 * t589 - t352 * t367) * MDP(18) + ((t371 * t456 + (t382 - t577) * t453) * MDP(13) + (t534 * t456 + (-t367 - t577) * t453) * MDP(15) + (t306 * t456 + t455 * t474) * MDP(25) + (-t305 * t456 + t452 * t474) * MDP(24) + ((-qJ(4) * qJD(2) + t367) * t456 + t507) * MDP(16) + (-t371 * t453 + t382 * t456) * MDP(11) + (t529 + (-t376 - t547) * t456) * MDP(22) + (-t528 + (t377 - t549) * t456) * MDP(21) + ((t365 - t384 - t552) * t453 + (-pkin(2) * qJD(2) - t362 + t586) * t456) * MDP(12) + ((t349 + t368 + t552) * t453 + (-qJD(2) * t457 + t331 - t589) * t456) * MDP(17)) * t551 + (t569 + t407 + 0.2e1 * t410) * MDP(13) + (pkin(1) * t527 - t469 - t569) * MDP(10) - t409 * MDP(23) * t520 + t566 * t585 + (-pkin(2) * t378 + qJ(3) * t357 - t362 * t384 + t365 * t586 - t371 * t382) * MDP(14); (-t365 * t430 + t371 * t521 + t378) * MDP(14) + (t349 * t430 + (t507 - t522) * t551 + t378) * MDP(18) + (-t515 + t376 * t430 + (-t452 * t546 - t528) * t551) * MDP(24) + (t516 + t377 * t430 + (-t455 * t546 + t529) * t551) * MDP(25) + (t452 * t309 - (t454 * t430 - t451 * t568) * t533 + (t470 + t473) * t455) * MDP(31) + (t562 - (-t451 * t430 - t454 * t568) * t533 + (t472 - t584) * t455) * MDP(32) + (-MDP(11) + MDP(16)) * t500 + (MDP(13) + MDP(15)) * (-t430 ^ 2 - t446 * t566) + (-MDP(12) + MDP(17)) * t478; t412 * MDP(15) + t388 * MDP(16) + t556 * MDP(18) + (-t446 - t447) * MDP(17) * t566 - (MDP(31) * t372 + MDP(32) * t373) * t533 + (-MDP(25) * t544 + (-t533 * t543 - t309) * MDP(31) + (-t533 * t541 - t308) * MDP(32)) * t455 + (-MDP(24) * t544 + t463 * MDP(31) + t462 * MDP(32)) * t452 + ((t430 * MDP(15) - t349 * MDP(18) + (-t376 + t547) * MDP(24) + (-t377 - t549) * MDP(25)) * t456 + (-MDP(25) * t567 + t331 * MDP(18) + (MDP(18) * t457 + MDP(16)) * qJD(2) + (-MDP(24) * t409 + t342 * MDP(31) - MDP(32) * t482) * t452) * t453) * t551; -t376 ^ 2 * MDP(20) + (t347 - t572) * MDP(21) + (t348 - t571) * MDP(22) + MDP(23) * t412 + (t306 * t409 + t459) * MDP(24) + (t305 * t409 - t333 * t376 - t468) * MDP(25) + (-t533 * t575 + t576) * MDP(26) + ((t308 - t588) * t454 + (-t309 + t587) * t451) * MDP(27) + (t533 * t573 - t563) * MDP(28) + (-t533 * t574 - t560) * MDP(29) + (-pkin(5) * t309 - t306 * t342 + t461 * t451 - t454 * t580) * MDP(31) + (-pkin(5) * t308 + t306 * t482 + t451 * t580 + t461 * t454) * MDP(32) + (MDP(19) * t376 + t377 * MDP(20) + t333 * MDP(24) - MDP(28) * t482 - t342 * MDP(29) + MDP(30) * t533 - MDP(31) * t489 - MDP(32) * t298) * t377; -t482 * t342 * MDP(26) + (-t342 ^ 2 + t482 ^ 2) * MDP(27) + (t526 + t588) * MDP(28) + (-t508 - t587) * MDP(29) - t348 * MDP(30) + (t298 * t533 + t303 * t482 + t510) * MDP(31) + (t303 * t342 - t489 * t533 - t490) * MDP(32) + (MDP(28) * t570 + MDP(29) * t482 - MDP(31) * t298 + MDP(32) * t489) * qJD(6);];
tauc  = t1;
