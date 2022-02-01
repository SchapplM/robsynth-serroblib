% Calculate vector of inverse dynamics joint torques for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:32
% EndTime: 2022-01-20 11:49:38
% DurationCPUTime: 3.45s
% Computational Cost: add. (3354->357), mult. (5012->421), div. (0->0), fcn. (3237->12), ass. (0->192)
t568 = cos(qJ(4));
t507 = t568 * qJD(4);
t582 = t568 * qJD(3) + t507;
t455 = cos(qJ(3));
t446 = qJD(1) + qJD(2);
t453 = sin(qJ(2));
t567 = pkin(1) * t453;
t519 = qJD(1) * t567;
t392 = pkin(7) * t446 + t519;
t506 = pkin(8) * t446 + t392;
t351 = t506 * t455;
t451 = sin(qJ(4));
t343 = t451 * t351;
t452 = sin(qJ(3));
t350 = t506 * t452;
t346 = qJD(3) * pkin(3) - t350;
t502 = t568 * t346 - t343;
t382 = t451 * t455 + t568 * t452;
t361 = t382 * t446;
t553 = t361 * qJ(5);
t307 = -t553 + t502;
t527 = qJD(2) * t453;
t433 = pkin(1) * t527;
t456 = cos(qJ(2));
t566 = pkin(1) * t456;
t530 = -qJD(1) * t433 + qJDD(1) * t566;
t450 = qJ(1) + qJ(2);
t438 = cos(t450);
t562 = g(2) * t438;
t443 = qJDD(1) + qJDD(2);
t565 = pkin(2) * t443;
t581 = -t530 - t565 + t562;
t580 = t582 * t455;
t543 = t451 * t452;
t481 = t568 * t455 - t543;
t440 = t455 * pkin(3);
t560 = pkin(2) + t440;
t521 = qJDD(1) * t453;
t526 = qJD(2) * t456;
t365 = pkin(7) * t443 + (qJD(1) * t526 + t521) * pkin(1);
t524 = qJD(3) * t455;
t509 = t446 * t524;
t547 = t443 * t452;
t317 = -t392 * t524 + qJDD(3) * pkin(3) - t365 * t452 + (-t509 - t547) * pkin(8);
t525 = qJD(3) * t452;
t510 = t446 * t525;
t546 = t443 * t455;
t318 = -t392 * t525 + t365 * t455 + (-t510 + t546) * pkin(8);
t579 = t568 * t317 - t451 * t318;
t569 = -pkin(7) - pkin(8);
t514 = qJD(3) * t569;
t389 = t452 * t514;
t390 = t455 * t514;
t528 = qJD(1) * t456;
t518 = pkin(1) * t528;
t408 = t569 * t452;
t439 = t455 * pkin(8);
t409 = pkin(7) * t455 + t439;
t533 = t451 * t408 + t568 * t409;
t578 = -t533 * qJD(4) + t382 * t518 - t451 * t389 + t568 * t390;
t523 = qJD(4) * t451;
t577 = -t568 * t389 - t451 * t390 - t408 * t507 + t409 * t523 + t481 * t518;
t426 = pkin(7) + t567;
t559 = -pkin(8) - t426;
t378 = t559 * t452;
t379 = t426 * t455 + t439;
t535 = t451 * t378 + t568 * t379;
t442 = qJDD(3) + qJDD(4);
t445 = qJD(3) + qJD(4);
t563 = pkin(3) * t445;
t576 = -t451 * pkin(3) * t442 - t507 * t563;
t436 = sin(t450);
t575 = g(1) * t438 + g(2) * t436;
t432 = pkin(3) * t525;
t574 = t432 - t519;
t434 = t442 * pkin(4);
t486 = t445 * t543;
t494 = t382 * t443 + t580 * t446;
t313 = t446 * t486 - t494;
t557 = t313 * qJ(5);
t573 = t434 + t557;
t421 = g(1) * t436;
t572 = t562 - t421;
t449 = qJ(3) + qJ(4);
t435 = sin(t449);
t550 = t435 * t438;
t551 = t435 * t436;
t437 = cos(t449);
t561 = g(3) * t437;
t571 = g(1) * t550 + g(2) * t551 - t561;
t570 = t361 ^ 2;
t564 = pkin(2) * t446;
t558 = qJ(5) * t382;
t336 = t445 * t382;
t487 = t481 * t443;
t314 = t336 * t446 - t487;
t556 = t314 * qJ(5);
t359 = t481 * t446;
t555 = t359 * qJ(5);
t554 = t359 * t445;
t549 = t436 * t437;
t548 = t437 * t438;
t545 = t446 * t452;
t541 = t452 * t455;
t537 = -t336 * qJ(5) + qJD(5) * t481;
t540 = t537 - t577;
t335 = t486 - t580;
t485 = t335 * qJ(5) - t382 * qJD(5);
t539 = t485 + t578;
t306 = pkin(4) * t445 + t307;
t538 = t306 - t307;
t536 = -t568 * t350 - t343;
t393 = -t518 - t564;
t534 = t393 * t525 + t455 * t421;
t531 = pkin(4) * t437 + t440;
t447 = t452 ^ 2;
t529 = -t455 ^ 2 + t447;
t363 = -t446 * t560 - t518;
t503 = -pkin(4) * t359 + qJD(5);
t328 = t363 + t503;
t522 = qJD(5) + t328;
t520 = pkin(3) * t545;
t517 = pkin(1) * t526;
t516 = t393 * t524 + t581 * t452;
t428 = -pkin(2) - t566;
t345 = t568 * t351;
t511 = t446 * t527;
t331 = pkin(4) * t336 + t432;
t504 = qJD(3) * t559;
t501 = t350 * t451 - t345;
t500 = t568 * t378 - t379 * t451;
t388 = pkin(2) + t531;
t444 = qJ(5) - t569;
t499 = t438 * t388 + t436 * t444;
t498 = t568 * t408 - t409 * t451;
t497 = -t388 * t436 + t444 * t438;
t496 = t445 * t452;
t495 = t446 * t519;
t358 = -pkin(4) * t481 - t560;
t491 = -g(1) * t551 + g(2) * t550;
t490 = g(1) * t549 - g(2) * t548;
t489 = t331 - t519;
t484 = -t530 + t572;
t482 = -t451 * t346 - t345;
t467 = qJD(4) * t482 + t579;
t291 = -t361 * qJD(5) + t467 + t573;
t464 = t451 * t317 + t568 * t318 + t346 * t507 - t351 * t523;
t292 = t359 * qJD(5) + t464 - t556;
t308 = -t482 + t555;
t483 = -t291 * t382 + t292 * t481 + t306 * t335 - t308 * t336 - t575;
t332 = pkin(3) * t510 - t443 * t560 - t530;
t303 = pkin(4) * t314 + qJDD(5) + t332;
t480 = t303 * t382 - t328 * t335 + t491;
t479 = t332 * t382 - t363 * t335 + t491;
t478 = -t303 * t481 + t328 * t336 + t490;
t477 = -t332 * t481 + t363 * t336 + t490;
t347 = t452 * t504 + t455 * t517;
t348 = -t452 * t517 + t455 * t504;
t476 = t568 * t347 + t451 * t348 + t378 * t507 - t379 * t523;
t357 = t359 ^ 2;
t474 = -t361 * t359 * MDP(14) + (-t451 * t446 * t496 + t494 - t554) * MDP(16) + t487 * MDP(17) + (-t357 + t570) * MDP(15) + t442 * MDP(18);
t473 = -t393 * t446 - t365 + t575;
t458 = qJD(3) ^ 2;
t472 = (-t313 * t481 - t314 * t382 - t335 * t359 - t336 * t361) * MDP(15) + (-t313 * t382 - t335 * t361) * MDP(14) + (-t335 * t445 + t382 * t442) * MDP(16) + (-t336 * t445 + t442 * t481) * MDP(17) + 0.2e1 * (-t529 * t446 * qJD(3) + t443 * t541) * MDP(8) + (t443 * t447 + 0.2e1 * t452 * t509) * MDP(7) + (qJDD(3) * t455 - t452 * t458) * MDP(10) + (qJDD(3) * t452 + t455 * t458) * MDP(9) + t443 * MDP(4);
t471 = pkin(7) * t458 - t495 - t565;
t470 = pkin(1) * t511 + t426 * t458 + t428 * t443;
t469 = -pkin(7) * qJDD(3) + (t518 - t564) * qJD(3);
t468 = -qJDD(3) * t426 + (t428 * t446 - t517) * qJD(3);
t466 = -t535 * qJD(4) - t451 * t347 + t568 * t348;
t463 = g(1) * t548 + g(2) * t549 + g(3) * t435 - t464;
t462 = t467 + t571;
t461 = -t363 * t359 + t463;
t460 = -t363 * t361 + t462;
t459 = -t522 * t359 + t463 + t556;
t457 = cos(qJ(1));
t454 = sin(qJ(1));
t427 = t568 * pkin(3) + pkin(4);
t403 = t428 - t440;
t391 = t433 + t432;
t377 = t481 * qJ(5);
t349 = t358 - t566;
t337 = pkin(4) * t361 + t520;
t330 = t377 + t533;
t329 = t498 - t558;
t327 = t331 + t433;
t322 = t377 + t535;
t321 = t500 - t558;
t310 = -t553 + t536;
t309 = t501 - t555;
t295 = t466 + t485;
t294 = t476 + t537;
t1 = [t472 + ((t443 * t456 - t511) * pkin(1) - t484) * MDP(5) + (t468 * t452 + (-t470 - t581) * t455 + t534) * MDP(12) + (t468 * t455 + (t470 - t421) * t452 + t516) * MDP(13) + (t403 * t314 - t391 * t359 + t500 * t442 + t466 * t445 + t477) * MDP(19) + (((-qJDD(1) - t443) * t453 + (-qJD(1) - t446) * t526) * pkin(1) + t575) * MDP(6) + qJDD(1) * MDP(1) + (t295 * t445 + t314 * t349 + t321 * t442 - t327 * t359 + t478) * MDP(21) + (g(1) * t454 - g(2) * t457) * MDP(2) + (g(1) * t457 + g(2) * t454) * MDP(3) + (t292 * t322 + t308 * t294 + t291 * t321 + t306 * t295 + t303 * t349 + t328 * t327 - g(1) * (-pkin(1) * t454 + t497) - g(2) * (pkin(1) * t457 + t499)) * MDP(24) + (-t403 * t313 + t391 * t361 - t535 * t442 - t476 * t445 + t479) * MDP(20) + (-t294 * t445 - t313 * t349 - t322 * t442 + t327 * t361 + t480) * MDP(22) + (t294 * t359 - t295 * t361 + t313 * t321 - t314 * t322 + t483) * MDP(23); t472 + (-g(1) * t497 - g(2) * t499 + t291 * t329 + t292 * t330 + t303 * t358 + t539 * t306 + t540 * t308 + t489 * t328) * MDP(24) + ((-t521 + (-qJD(2) + t446) * t528) * pkin(1) + t575) * MDP(6) + (t469 * t452 + (-t471 - t581) * t455 + t534) * MDP(12) + (t469 * t455 + (t471 - t421) * t452 + t516) * MDP(13) + (-t314 * t560 - t574 * t359 + t498 * t442 + t578 * t445 + t477) * MDP(19) + (t314 * t358 + t329 * t442 - t489 * t359 + t539 * t445 + t478) * MDP(21) + (-t484 + t495) * MDP(5) + (t313 * t560 + t574 * t361 - t533 * t442 + t577 * t445 + t479) * MDP(20) + (-t313 * t358 - t330 * t442 + t489 * t361 - t540 * t445 + t480) * MDP(22) + (t313 * t329 - t314 * t330 + t540 * t359 - t539 * t361 + t483) * MDP(23); MDP(9) * t547 + MDP(10) * t546 + qJDD(3) * MDP(11) + (-g(3) * t455 + t452 * t473) * MDP(12) + (g(3) * t452 + t455 * t473) * MDP(13) + (-t501 * t445 + (t359 * t545 + t568 * t442 - t445 * t523) * pkin(3) + t460) * MDP(19) + (-t361 * t520 + t536 * t445 + t461 + t576) * MDP(20) + (-t309 * t445 + t337 * t359 + t427 * t442 - t522 * t361 + (-t345 + (-t346 - t563) * t451) * qJD(4) + t571 + t573 + t579) * MDP(21) + (t310 * t445 - t337 * t361 + t459 + t576) * MDP(22) + (t427 * t313 + (t308 + t309) * t361 - (-t306 + t310) * t359 + (-t314 * t451 + (t568 * t359 + t361 * t451) * qJD(4)) * pkin(3)) * MDP(23) + (t291 * t427 - t308 * t310 - t306 * t309 - t328 * t337 - g(3) * t531 - t575 * (-pkin(3) * t452 - pkin(4) * t435) + (t292 * t451 + (-t306 * t451 + t568 * t308) * qJD(4)) * pkin(3)) * MDP(24) + t474 + (-MDP(7) * t541 + t529 * MDP(8)) * t446 ^ 2; (-t445 * t482 + t460) * MDP(19) + (t502 * t445 + t461) * MDP(20) + (t557 + t308 * t445 + 0.2e1 * t434 + (-t328 - t503) * t361 + t462) * MDP(21) + (-t570 * pkin(4) + t307 * t445 + t459) * MDP(22) + (pkin(4) * t313 + t538 * t359) * MDP(23) + (t538 * t308 + (-t328 * t361 + t435 * t575 + t291 - t561) * pkin(4)) * MDP(24) + t474; (t361 * t445 - t487) * MDP(21) + (t494 + t554) * MDP(22) + (-t357 - t570) * MDP(23) + (t306 * t361 - t308 * t359 + t303 + t572) * MDP(24) + (t582 * t452 * MDP(21) + (t445 * MDP(21) * t455 - MDP(22) * t496) * t451) * t446;];
tau = t1;
