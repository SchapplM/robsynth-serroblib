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
% Datum: 2021-01-16 00:01
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 00:01:09
% EndTime: 2021-01-16 00:01:17
% DurationCPUTime: 3.46s
% Computational Cost: add. (3354->357), mult. (5012->421), div. (0->0), fcn. (3237->12), ass. (0->190)
t562 = cos(qJ(4));
t500 = t562 * qJD(4);
t577 = t562 * qJD(3) + t500;
t448 = sin(qJ(2));
t520 = qJD(2) * t448;
t428 = pkin(1) * t520;
t451 = cos(qJ(2));
t560 = pkin(1) * t451;
t523 = -qJD(1) * t428 + qJDD(1) * t560;
t438 = qJDD(1) + qJDD(2);
t559 = pkin(2) * t438;
t362 = -t523 - t559;
t445 = qJ(1) + qJ(2);
t431 = sin(t445);
t433 = cos(t445);
t481 = g(2) * t433 + g(3) * t431;
t576 = t362 + t481;
t450 = cos(qJ(3));
t441 = qJD(1) + qJD(2);
t561 = pkin(1) * t448;
t511 = qJD(1) * t561;
t390 = pkin(7) * t441 + t511;
t499 = pkin(8) * t441 + t390;
t349 = t499 * t450;
t446 = sin(qJ(4));
t341 = t446 * t349;
t447 = sin(qJ(3));
t348 = t499 * t447;
t344 = qJD(3) * pkin(3) - t348;
t496 = t562 * t344 - t341;
t380 = t446 * t450 + t562 * t447;
t359 = t380 * t441;
t546 = t359 * qJ(5);
t305 = -t546 + t496;
t575 = t577 * t450;
t536 = t446 * t447;
t474 = t562 * t450 - t536;
t435 = t450 * pkin(3);
t574 = -pkin(2) - t435;
t514 = qJDD(1) * t448;
t519 = qJD(2) * t451;
t363 = pkin(7) * t438 + (qJD(1) * t519 + t514) * pkin(1);
t517 = qJD(3) * t450;
t502 = t441 * t517;
t540 = t438 * t447;
t315 = -t390 * t517 + qJDD(3) * pkin(3) - t363 * t447 + (-t502 - t540) * pkin(8);
t518 = qJD(3) * t447;
t503 = t441 * t518;
t539 = t438 * t450;
t316 = -t390 * t518 + t363 * t450 + (-t503 + t539) * pkin(8);
t573 = t562 * t315 - t446 * t316;
t563 = -pkin(7) - pkin(8);
t507 = qJD(3) * t563;
t387 = t447 * t507;
t388 = t450 * t507;
t521 = qJD(1) * t451;
t510 = pkin(1) * t521;
t404 = t563 * t447;
t434 = t450 * pkin(8);
t405 = pkin(7) * t450 + t434;
t527 = t446 * t404 + t562 * t405;
t572 = -qJD(4) * t527 + t380 * t510 - t446 * t387 + t562 * t388;
t516 = qJD(4) * t446;
t571 = -t562 * t387 - t446 * t388 - t404 * t500 + t405 * t516 + t474 * t510;
t421 = pkin(7) + t561;
t552 = -pkin(8) - t421;
t376 = t552 * t447;
t377 = t421 * t450 + t434;
t528 = t446 * t376 + t562 * t377;
t437 = qJDD(3) + qJDD(4);
t440 = qJD(3) + qJD(4);
t557 = pkin(3) * t440;
t570 = -t446 * pkin(3) * t437 - t500 * t557;
t569 = g(2) * t431 - g(3) * t433;
t427 = pkin(3) * t518;
t568 = t427 - t511;
t429 = t437 * pkin(4);
t478 = t440 * t536;
t485 = t380 * t438 + t575 * t441;
t311 = t441 * t478 - t485;
t550 = t311 * qJ(5);
t567 = t429 + t550;
t453 = qJD(3) ^ 2;
t566 = pkin(7) * t453 - t559;
t444 = qJ(3) + qJ(4);
t430 = sin(t444);
t543 = t430 * t433;
t544 = t430 * t431;
t432 = cos(t444);
t555 = g(1) * t432;
t565 = g(2) * t544 - g(3) * t543 - t555;
t564 = t359 ^ 2;
t558 = pkin(2) * t441;
t551 = qJ(5) * t380;
t334 = t440 * t380;
t479 = t474 * t438;
t312 = t334 * t441 - t479;
t549 = t312 * qJ(5);
t357 = t474 * t441;
t548 = t357 * qJ(5);
t547 = t357 * t440;
t542 = t431 * t432;
t541 = t432 * t433;
t538 = t441 * t447;
t534 = t447 * t450;
t530 = -t334 * qJ(5) + qJD(5) * t474;
t533 = t530 - t571;
t333 = t478 - t575;
t477 = t333 * qJ(5) - t380 * qJD(5);
t532 = t477 + t572;
t304 = pkin(4) * t440 + t305;
t531 = t304 - t305;
t529 = -t562 * t348 - t341;
t526 = g(2) * t543 + g(3) * t544;
t524 = pkin(4) * t432 + t435;
t442 = t447 ^ 2;
t522 = -t450 ^ 2 + t442;
t361 = t441 * t574 - t510;
t497 = -pkin(4) * t357 + qJD(5);
t326 = t361 + t497;
t515 = qJD(5) + t326;
t513 = pkin(3) * t538;
t509 = pkin(1) * t519;
t423 = -pkin(2) - t560;
t343 = t562 * t349;
t504 = t441 * t520;
t329 = pkin(4) * t334 + t427;
t498 = qJD(3) * t552;
t495 = t348 * t446 - t343;
t494 = t562 * t376 - t377 * t446;
t386 = -pkin(2) - t524;
t439 = -qJ(5) + t563;
t493 = -t433 * t386 - t431 * t439;
t492 = t562 * t404 - t405 * t446;
t491 = -t386 * t431 + t433 * t439;
t490 = t440 * t447;
t489 = t441 * t511;
t330 = pkin(3) * t503 + t438 * t574 - t523;
t301 = pkin(4) * t312 + qJDD(5) + t330;
t488 = t301 * t380 - t326 * t333 + t526;
t487 = t330 * t380 - t361 * t333 + t526;
t391 = -t510 - t558;
t486 = t391 * t517 + t576 * t447;
t356 = -pkin(4) * t474 + t574;
t482 = t329 - t511;
t475 = -t446 * t344 - t343;
t462 = t475 * qJD(4) + t573;
t289 = -t359 * qJD(5) + t462 + t567;
t459 = t446 * t315 + t562 * t316 + t344 * t500 - t349 * t516;
t290 = t357 * qJD(5) + t459 - t549;
t306 = -t475 + t548;
t476 = -t289 * t380 + t290 * t474 + t304 * t333 - t306 * t334 - t569;
t345 = t447 * t498 + t450 * t509;
t346 = -t447 * t509 + t450 * t498;
t473 = t562 * t345 + t446 * t346 + t376 * t500 - t377 * t516;
t355 = t357 ^ 2;
t471 = -t359 * t357 * MDP(14) + (-t446 * t441 * t490 + t485 - t547) * MDP(16) + t479 * MDP(17) + (-t355 + t564) * MDP(15) + t437 * MDP(18);
t470 = -t391 * t441 - t363 + t569;
t469 = -g(2) * t541 - g(3) * t542 - t330 * t474 + t361 * t334;
t468 = -t301 * t474 + t326 * t334 - t481 * t432;
t467 = (-t311 * t474 - t312 * t380 - t333 * t357 - t334 * t359) * MDP(15) + (-t311 * t380 - t333 * t359) * MDP(14) + (-t333 * t440 + t380 * t437) * MDP(16) + (-t334 * t440 + t437 * t474) * MDP(17) + 0.2e1 * (-t522 * t441 * qJD(3) + t438 * t534) * MDP(8) + (t438 * t442 + 0.2e1 * t447 * t502) * MDP(7) + (qJDD(3) * t450 - t447 * t453) * MDP(10) + (qJDD(3) * t447 + t450 * t453) * MDP(9) + t438 * MDP(4);
t466 = -t481 + t489;
t465 = pkin(1) * t504 + t421 * t453 + t423 * t438;
t464 = -pkin(7) * qJDD(3) + (t510 - t558) * qJD(3);
t463 = -qJDD(3) * t421 + (t423 * t441 - t509) * qJD(3);
t461 = -qJD(4) * t528 - t446 * t345 + t562 * t346;
t458 = g(1) * t430 + g(2) * t542 - g(3) * t541 - t459;
t457 = t462 + t565;
t456 = -t361 * t357 + t458;
t455 = -t361 * t359 + t457;
t454 = -t515 * t357 + t458 + t549;
t452 = cos(qJ(1));
t449 = sin(qJ(1));
t422 = t562 * pkin(3) + pkin(4);
t399 = t423 - t435;
t389 = t428 + t427;
t375 = t474 * qJ(5);
t369 = t391 * t518;
t347 = t356 - t560;
t335 = pkin(4) * t359 + t513;
t328 = t375 + t527;
t327 = t492 - t551;
t325 = t329 + t428;
t320 = t375 + t528;
t319 = t494 - t551;
t308 = -t546 + t529;
t307 = t495 - t548;
t293 = t461 + t477;
t292 = t473 + t530;
t1 = [(t369 + t463 * t447 + (-t465 - t576) * t450) * MDP(12) + (t465 * t447 + t463 * t450 + t486) * MDP(13) + ((t438 * t451 - t504) * pkin(1) - t481 + t523) * MDP(5) + (t399 * t312 - t389 * t357 + t494 * t437 + t461 * t440 + t469) * MDP(19) + (t293 * t440 + t312 * t347 + t319 * t437 - t325 * t357 + t468) * MDP(21) + (-t399 * t311 + t389 * t359 - t528 * t437 - t473 * t440 + t487) * MDP(20) + (-t292 * t440 - t311 * t347 - t320 * t437 + t325 * t359 + t488) * MDP(22) + (t292 * t357 - t293 * t359 + t311 * t319 - t312 * t320 + t476) * MDP(23) + (t290 * t320 + t306 * t292 + t289 * t319 + t304 * t293 + t301 * t347 + t326 * t325 - g(2) * (pkin(1) * t452 + t493) - g(3) * (pkin(1) * t449 + t491)) * MDP(24) + (((-qJDD(1) - t438) * t448 + (-qJD(1) - t441) * t519) * pkin(1) + t569) * MDP(6) + (-g(2) * t452 - g(3) * t449) * MDP(2) + (g(2) * t449 - g(3) * t452) * MDP(3) + qJDD(1) * MDP(1) + t467; ((-t514 + (-qJD(2) + t441) * t521) * pkin(1) + t569) * MDP(6) + (t369 + t464 * t447 + (-t362 + t466 - t566) * t450) * MDP(12) + (t466 + t523) * MDP(5) + (t464 * t450 + (-t489 + t566) * t447 + t486) * MDP(13) + (t312 * t574 - t357 * t568 + t492 * t437 + t440 * t572 + t469) * MDP(19) + (t312 * t356 + t327 * t437 - t482 * t357 + t532 * t440 + t468) * MDP(21) + (-t311 * t574 + t359 * t568 - t527 * t437 + t440 * t571 + t487) * MDP(20) + (-t311 * t356 - t328 * t437 + t482 * t359 - t533 * t440 + t488) * MDP(22) + (t311 * t327 - t312 * t328 + t533 * t357 - t532 * t359 + t476) * MDP(23) + (-g(2) * t493 - g(3) * t491 + t289 * t327 + t290 * t328 + t301 * t356 + t532 * t304 + t533 * t306 + t482 * t326) * MDP(24) + t467; MDP(9) * t540 + MDP(10) * t539 + qJDD(3) * MDP(11) + (-g(1) * t450 + t470 * t447) * MDP(12) + (g(1) * t447 + t470 * t450) * MDP(13) + (-t495 * t440 + (t357 * t538 + t562 * t437 - t440 * t516) * pkin(3) + t455) * MDP(19) + (-t359 * t513 + t529 * t440 + t456 + t570) * MDP(20) + (-t307 * t440 + t335 * t357 + t422 * t437 - t515 * t359 + (-t343 + (-t344 - t557) * t446) * qJD(4) + t565 + t567 + t573) * MDP(21) + (t308 * t440 - t335 * t359 + t454 + t570) * MDP(22) + (t422 * t311 + (t306 + t307) * t359 - (-t304 + t308) * t357 + (-t312 * t446 + (t562 * t357 + t359 * t446) * qJD(4)) * pkin(3)) * MDP(23) + (t289 * t422 - t306 * t308 - t304 * t307 - t326 * t335 - g(1) * t524 + t569 * (pkin(3) * t447 + pkin(4) * t430) + (t290 * t446 + (-t304 * t446 + t562 * t306) * qJD(4)) * pkin(3)) * MDP(24) + t471 + (-MDP(7) * t534 + t522 * MDP(8)) * t441 ^ 2; (-t475 * t440 + t455) * MDP(19) + (t496 * t440 + t456) * MDP(20) + (t550 + t306 * t440 + 0.2e1 * t429 + (-t326 - t497) * t359 + t457) * MDP(21) + (-t564 * pkin(4) + t305 * t440 + t454) * MDP(22) + (pkin(4) * t311 + t531 * t357) * MDP(23) + (t531 * t306 + (-t326 * t359 + t430 * t569 + t289 - t555) * pkin(4)) * MDP(24) + t471; (t359 * t440 - t479) * MDP(21) + (t485 + t547) * MDP(22) + (-t355 - t564) * MDP(23) + (t304 * t359 - t306 * t357 + t301 + t481) * MDP(24) + (t577 * t447 * MDP(21) + (t440 * MDP(21) * t450 - MDP(22) * t490) * t446) * t441;];
tau = t1;
