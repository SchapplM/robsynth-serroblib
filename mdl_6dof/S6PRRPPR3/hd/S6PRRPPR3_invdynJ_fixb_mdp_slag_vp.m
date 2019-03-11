% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:11:10
% EndTime: 2019-03-08 21:11:16
% DurationCPUTime: 5.57s
% Computational Cost: add. (1935->479), mult. (4308->613), div. (0->0), fcn. (3127->10), ass. (0->217)
t433 = sin(qJ(2));
t426 = sin(pkin(6));
t542 = qJD(1) * t426;
t508 = t433 * t542;
t367 = qJD(2) * pkin(8) + t508;
t435 = cos(qJ(3));
t428 = cos(pkin(6));
t541 = qJD(1) * t428;
t390 = t435 * t541;
t432 = sin(qJ(3));
t548 = t432 * t367 - t390;
t579 = qJD(4) + t548;
t578 = t432 * qJ(4) + pkin(2);
t577 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t327 = -qJD(3) * pkin(3) + t579;
t522 = qJD(2) * qJD(3);
t502 = t435 * t522;
t518 = qJDD(2) * t432;
t463 = t502 + t518;
t503 = t432 * t522;
t517 = qJDD(2) * t435;
t576 = -t503 + t517;
t437 = -pkin(3) - pkin(4);
t509 = qJD(3) * t437;
t540 = qJD(2) * t432;
t323 = -qJ(5) * t540 + t548;
t524 = -qJD(4) - t323;
t309 = t509 - t524;
t538 = qJD(2) * t435;
t570 = pkin(8) - qJ(5);
t332 = t435 * t367 + t432 * t541;
t324 = -qJ(5) * t538 + t332;
t422 = qJD(3) * qJ(4);
t312 = -t324 - t422;
t575 = 0.2e1 * t577;
t328 = t422 + t332;
t574 = qJDD(3) * t437;
t573 = MDP(13) - MDP(18);
t421 = -pkin(9) + t437;
t436 = cos(qJ(2));
t537 = qJD(2) * t436;
t505 = qJD(1) * t537;
t519 = qJDD(1) * t433;
t566 = qJDD(2) * pkin(8);
t337 = t566 + (t505 + t519) * t426;
t532 = qJD(3) * t432;
t504 = qJD(1) * t532;
t520 = qJDD(1) * t428;
t531 = qJD(3) * t435;
t488 = t432 * t337 + t367 * t531 + t428 * t504 - t435 * t520;
t475 = -qJDD(4) - t488;
t521 = qJD(2) * qJD(5);
t441 = -qJ(5) * t463 - t432 * t521 - t475;
t293 = qJDD(3) * t421 + t441;
t487 = -qJD(3) * t390 - t435 * t337 + t367 * t532 - t432 * t520;
t302 = -t487 + t577;
t392 = qJ(5) * t503;
t298 = t435 * (qJ(5) * qJDD(2) + t521) - t302 - t392;
t296 = qJDD(3) * pkin(5) - t298;
t389 = t436 * t542;
t429 = qJD(2) * pkin(2);
t368 = -t389 - t429;
t334 = -pkin(3) * t538 - qJ(4) * t540 + t368;
t315 = pkin(4) * t538 + qJD(5) - t334;
t484 = pkin(5) * t432 + pkin(9) * t435;
t305 = qJD(2) * t484 + t315;
t310 = qJD(3) * pkin(5) - t312;
t369 = -t435 * pkin(3) - t578;
t359 = t435 * pkin(4) - t369;
t333 = t484 + t359;
t347 = -qJD(5) * t432 + t531 * t570;
t360 = qJDD(6) + t463;
t374 = t570 * t432;
t397 = qJD(6) + t540;
t569 = sin(pkin(10));
t498 = t569 * t436;
t427 = cos(pkin(10));
t553 = t427 * t433;
t350 = t428 * t553 + t498;
t499 = t569 * t433;
t552 = t427 * t436;
t352 = -t428 * t499 + t552;
t482 = g(1) * t352 + g(2) * t350;
t572 = -(qJD(6) * t333 + t347) * t397 + (qJD(3) * t310 - qJD(6) * t305 - t293) * t432 - t296 * t435 - t374 * t360 + t482;
t403 = qJ(4) * t538;
t317 = -t427 * t426 * t432 + t350 * t435;
t500 = t426 * t569;
t319 = t352 * t435 + t432 * t500;
t555 = t426 * t435;
t356 = t428 * t432 + t433 * t555;
t462 = g(1) * t319 + g(2) * t317 + g(3) * t356;
t464 = pkin(5) * t435 + t421 * t432;
t571 = t397 * (qJD(2) * t464 + qJD(6) * t421 + t403) - t296 + t462;
t568 = pkin(8) * qJDD(3);
t567 = qJ(4) * t435;
t425 = qJDD(2) * pkin(2);
t565 = qJDD(3) * pkin(3);
t434 = cos(qJ(6));
t431 = sin(qJ(6));
t528 = qJD(6) * t435;
t506 = t431 * t528;
t481 = qJD(2) * t506 - qJDD(3) * t431 + t434 * t503;
t307 = (-qJD(3) * qJD(6) - t517) * t434 + t481;
t564 = t307 * t431;
t349 = t428 * t552 - t499;
t563 = t349 * t435;
t351 = -t428 * t498 - t553;
t562 = t351 * t435;
t525 = t434 * qJD(3);
t361 = t431 * t538 - t525;
t561 = t361 * t397;
t560 = t361 * t435;
t533 = qJD(3) * t431;
t362 = t434 * t538 + t533;
t559 = t362 * t397;
t558 = t362 * t435;
t557 = t397 * t434;
t556 = t426 * t433;
t554 = t426 * t436;
t551 = t432 * t436;
t550 = t434 * t436;
t549 = t435 * t436;
t379 = qJD(2) * t508;
t547 = t435 * t379 + t504 * t554;
t411 = t432 * qJD(4);
t546 = qJ(4) * t531 + t411;
t423 = t432 ^ 2;
t424 = t435 ^ 2;
t545 = t423 - t424;
t544 = t423 + t424;
t543 = MDP(16) * t432;
t539 = qJD(2) * t433;
t536 = qJD(3) * t312;
t535 = qJD(3) * t361;
t534 = qJD(3) * t362;
t530 = qJD(6) * t431;
t529 = qJD(6) * t434;
t527 = qJD(6) * t436;
t526 = t334 * MDP(15);
t523 = -qJD(5) - t315;
t516 = MDP(12) - MDP(17);
t515 = t397 * t431 * t432;
t514 = t432 * t557;
t513 = t432 * t556;
t439 = qJD(2) ^ 2;
t512 = t432 * t435 * t439;
t511 = pkin(3) * t563 + t578 * t349;
t510 = pkin(3) * t562 + t578 * t351;
t386 = qJDD(1) * t554;
t336 = t379 - t386 - t425;
t507 = t426 * t537;
t501 = t436 * t522;
t497 = t368 - t429;
t316 = t350 * t432 + t427 * t555;
t495 = -t316 * pkin(3) + qJ(4) * t317;
t318 = t352 * t432 - t435 * t500;
t494 = -t318 * pkin(3) + qJ(4) * t319;
t355 = -t428 * t435 + t513;
t493 = -t355 * pkin(3) + qJ(4) * t356;
t492 = qJD(2) * t359 + t315;
t491 = qJD(2) * t369 + t334;
t486 = t432 * t509;
t483 = g(1) * t351 + g(2) * t349;
t306 = qJD(3) * t421 - t524;
t294 = t305 * t434 - t306 * t431;
t295 = t305 * t431 + t306 * t434;
t480 = t309 * t432 - t312 * t435;
t476 = g(3) * (pkin(2) * t554 + pkin(8) * t556 + (pkin(3) * t549 + qJ(4) * t551) * t426);
t474 = qJDD(2) * t436 - t433 * t439;
t321 = -qJD(3) * t513 + (qJD(3) * t428 + t507) * t435;
t473 = qJD(3) * t321 + qJDD(3) * t356;
t322 = qJD(3) * t356 + t432 * t507;
t472 = -qJD(3) * t322 - qJDD(3) * t355;
t396 = g(3) * t556;
t471 = t396 + t482;
t469 = -t434 * qJDD(3) + t576 * t431;
t468 = -t431 * t433 + t432 * t550;
t467 = t431 * t551 + t433 * t434;
t466 = -t360 * t431 - t397 * t529;
t465 = -t360 * t434 + t397 * t530;
t460 = pkin(3) * t517 + t463 * qJ(4) + qJD(2) * t411 - t336;
t459 = -qJD(6) * t306 - t483;
t458 = g(3) * t554 + t483;
t456 = t464 * qJD(3);
t455 = pkin(4) * t517 + qJDD(5) + t460;
t438 = qJD(3) ^ 2;
t454 = pkin(8) * t438 + t458;
t453 = g(1) * t318 + g(2) * t316 + g(3) * t355 - t488;
t452 = -t462 - t487;
t451 = -t421 * t360 + (-t310 + t324) * t397;
t450 = -(-qJD(6) * t374 + t456 + t546) * t397 - t333 * t360 + t310 * t528;
t449 = -qJDD(4) + t453;
t448 = qJD(3) * t332 + t453;
t447 = qJD(3) * t548 + t452;
t446 = t336 + t454 - t425;
t301 = qJD(2) * t486 + t455;
t335 = t486 + t546;
t445 = -qJD(2) * t335 - qJDD(2) * t359 - t301 + t458;
t444 = -t449 - t565;
t304 = pkin(3) * t503 - t460;
t348 = pkin(3) * t532 - t546;
t443 = -qJD(2) * t348 - qJDD(2) * t369 - t304 - t454;
t303 = -t475 - t565;
t442 = t302 * t435 + t303 * t432 + (t327 * t435 - t328 * t432) * qJD(3) - t482;
t430 = qJ(4) + pkin(5);
t375 = t570 * t435;
t365 = t432 * t379;
t363 = pkin(3) * t540 - t403;
t346 = qJD(5) * t435 + t532 * t570;
t342 = t437 * t540 + t403;
t326 = t355 * t434 + t431 * t554;
t325 = -t355 * t431 + t426 * t550;
t308 = qJD(6) * t362 + t469;
t297 = t441 + t574;
t292 = qJD(2) * t456 + qJDD(2) * t484 + t455;
t291 = t434 * t292;
t1 = [(qJDD(1) - g(3)) * MDP(1) + t472 * MDP(10) + (t302 * t356 + t303 * t355 + t321 * t328 + t322 * t327 - g(3)) * MDP(15) + t473 * MDP(16) + (t297 * t355 - t298 * t356 + t309 * t322 - t312 * t321 - g(3)) * MDP(19) + ((-t322 * t431 - t355 * t529) * t397 + t325 * t360 - t321 * t361 - t356 * t308) * MDP(25) + (-(t322 * t434 - t355 * t530) * t397 - t326 * t360 - t321 * t362 + t356 * t307) * MDP(26) + (-MDP(11) + MDP(14)) * ((t432 * t474 + t435 * t501) * t426 + t473) + t516 * ((-t432 * t501 + t435 * t474) * t426 + t472) + ((-qJDD(2) * MDP(4) + (-MDP(19) * t315 + t526) * qJD(2) + (-MDP(10) * t435 - MDP(3) - t543) * t439) * t433 + (t576 * MDP(10) - t304 * MDP(15) + t463 * MDP(16) + t301 * MDP(19) + qJDD(2) * MDP(3) - t439 * MDP(4)) * t436 + ((-t431 * t527 - t434 * t539) * MDP(25) - (-t431 * t539 + t434 * t527) * MDP(26)) * t397) * t426 + (qJDD(2) * (t355 * t432 + t356 * t435) + (t321 * t435 + t322 * t432 + (t355 * t435 - t356 * t432) * qJD(3)) * qJD(2)) * t573; (-t365 + (-t568 + (t497 + t389) * qJD(3)) * t435 + t446 * t432) * MDP(11) + ((qJD(3) * t491 - t568) * t432 + t443 * t435 + t547) * MDP(12) + (-t396 + t442 + (-t426 * t505 + t566) * t544) * MDP(13) + (t365 + (t568 + (-t491 - t389) * qJD(3)) * t435 + t443 * t432) * MDP(14) + (t304 * t369 + t334 * t348 - g(1) * t510 - g(2) * t511 - t476 + (-t334 * t433 + (-t327 * t432 - t328 * t435) * t436) * t542 + t442 * pkin(8)) * MDP(15) + (qJDD(3) * t375 + t365 + (-t346 + (t492 - t389) * t435) * qJD(3) - t445 * t432) * MDP(16) + (qJDD(3) * t374 + (t432 * t492 + t347) * qJD(3) + t445 * t435 - t547) * MDP(17) + ((-qJD(3) * t309 - qJDD(2) * t375 + t298) * t435 + (-qJDD(2) * t374 - t297 - t536) * t432 + (t346 * t435 - t347 * t432 + (-t374 * t435 + t375 * t432) * qJD(3) + t544 * t389) * qJD(2) + t471) * MDP(18) + (t297 * t374 + t309 * t347 - t298 * t375 + t312 * t346 + t301 * t359 + t315 * t335 - g(1) * (pkin(4) * t562 + t352 * t570 + t510) - g(2) * (pkin(4) * t563 + t350 * t570 + t511) - t476 + (-g(3) * (pkin(4) * t549 - qJ(5) * t433) + (t315 * t433 - t436 * t480) * qJD(1)) * t426) * MDP(19) + (-t307 * t434 * t435 + (-t432 * t525 - t506) * t362) * MDP(20) + ((t361 * t434 + t362 * t431) * t532 + (t564 - t308 * t434 + (t361 * t431 - t362 * t434) * qJD(6)) * t435) * MDP(21) + ((t397 * t525 + t307) * t432 + (t465 - t534) * t435) * MDP(22) + ((-t397 * t533 + t308) * t432 + (-t466 + t535) * t435) * MDP(23) + (t360 * t432 + t397 * t531) * MDP(24) + (t294 * t531 + t291 * t432 - t375 * t308 + t346 * t361 + (t432 * t459 - t450) * t434 + t572 * t431 + (-g(3) * t468 + (t361 * t549 + t397 * t467) * qJD(1)) * t426) * MDP(25) + (-t295 * t531 + t375 * t307 + t346 * t362 + ((-t292 - t459) * t432 + t450) * t431 + t572 * t434 + (g(3) * t467 + (t362 * t549 + t397 * t468) * qJD(1)) * t426) * MDP(26) + qJDD(2) * MDP(2) + (t386 - t458) * MDP(3) + (-t426 * t519 + t471) * MDP(4) + (qJDD(2) * t423 + 0.2e1 * t432 * t502) * MDP(5) + 0.2e1 * (t432 * t517 - t522 * t545) * MDP(6) + (qJDD(3) * t432 + t435 * t438) * MDP(7) + (qJDD(3) * t435 - t432 * t438) * MDP(8) + ((qJD(3) * t497 - t568) * t432 - t446 * t435 + t547) * MDP(10); -MDP(5) * t512 + t545 * t439 * MDP(6) + MDP(7) * t518 + MDP(8) * t517 + qJDD(3) * MDP(9) + (-t368 * t540 + t448) * MDP(10) + (-t368 * t538 - t447) * MDP(11) + (0.2e1 * t565 - qJDD(4) + (-t334 * t432 + t363 * t435) * qJD(2) + t448) * MDP(12) + ((t334 * t435 + t363 * t432) * qJD(2) + t447 + t575) * MDP(14) + (-t303 * pkin(3) - g(1) * t494 - g(2) * t495 - g(3) * t493 + t302 * qJ(4) - t327 * t332 + t579 * t328 - t334 * t363) * MDP(15) + (-qJ(5) * t517 + qJD(3) * t323 + t392 + (-t342 * t432 + t435 * t523) * qJD(2) + t452 + t575) * MDP(16) + (-qJ(5) * t518 - qJD(3) * t324 + 0.2e1 * t574 + ((-qJ(5) * qJD(3) + t342) * t435 + t523 * t432) * qJD(2) - t449) * MDP(17) + (t297 * t437 - t298 * qJ(4) - t309 * t324 - t315 * t342 - g(1) * (-pkin(4) * t318 + t494) - g(2) * (-pkin(4) * t316 + t495) - g(3) * (-pkin(4) * t355 + t493) + t524 * t312) * MDP(19) + (t362 * t557 - t564) * MDP(20) + ((-t307 - t561) * t434 + (-t308 - t559) * t431) * MDP(21) + ((-t514 + t558) * qJD(2) + t466) * MDP(22) + ((t515 - t560) * qJD(2) + t465) * MDP(23) - t397 * MDP(24) * t538 + (-t294 * t538 - t430 * t308 + t524 * t361 + t451 * t431 - t434 * t571) * MDP(25) + (t295 * t538 + t430 * t307 + t524 * t362 + t431 * t571 + t451 * t434) * MDP(26) + ((-pkin(3) * t432 + t567) * MDP(13) + (-t432 * t437 - t567) * MDP(18)) * qJDD(2); (-qJD(3) * t328 + t444) * MDP(15) + (-qJDD(3) * pkin(4) - qJ(5) * t502 + t444 + t536) * MDP(19) + (t466 + t535) * MDP(25) + (t465 + t534) * MDP(26) + (MDP(14) + MDP(16)) * (-t423 * t439 - t438) - t516 * (qJDD(3) + t512) + ((-MDP(19) * qJ(5) + t573) * qJDD(2) + (t526 + t523 * MDP(19) + (-MDP(25) * t434 + MDP(26) * t431) * t397) * qJD(2)) * t432; (t455 - t458) * MDP(19) - t465 * MDP(25) + t466 * MDP(26) - t544 * MDP(18) * t439 + (-MDP(17) * t435 + t543) * qJDD(2) + (t480 * MDP(19) + (-t515 - t560) * MDP(25) + (-t514 - t558) * MDP(26) + (0.2e1 * t435 * MDP(16) + (MDP(19) * t437 + 0.2e1 * MDP(17)) * t432) * qJD(3)) * qJD(2); t362 * t361 * MDP(20) + (-t361 ^ 2 + t362 ^ 2) * MDP(21) + (-t434 * t517 + t481 - t561) * MDP(22) + (t469 - t559) * MDP(23) + t360 * MDP(24) + (-t431 * t293 + t291 + t295 * t397 + t310 * t362 - g(1) * (-t318 * t431 + t351 * t434) - g(2) * (-t316 * t431 + t349 * t434) - g(3) * t325) * MDP(25) + (-t434 * t293 - t431 * t292 + t294 * t397 - t310 * t361 - g(1) * (-t318 * t434 - t351 * t431) - g(2) * (-t316 * t434 - t349 * t431) + g(3) * t326) * MDP(26) + (-MDP(22) * t525 + MDP(23) * t362 - MDP(25) * t295 - MDP(26) * t294) * qJD(6);];
tau  = t1;
