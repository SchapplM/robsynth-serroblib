% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:21
% EndTime: 2022-01-20 11:43:26
% DurationCPUTime: 3.46s
% Computational Cost: add. (3165->366), mult. (4849->460), div. (0->0), fcn. (3381->16), ass. (0->193)
t576 = qJ(4) + pkin(7);
t465 = qJD(1) + qJD(2);
t471 = cos(pkin(9));
t478 = cos(qJ(3));
t543 = t471 * t478;
t522 = t465 * t543;
t470 = sin(pkin(9));
t474 = sin(qJ(3));
t547 = t470 * t474;
t388 = t465 * t547 - t522;
t477 = cos(qJ(5));
t375 = t477 * t388;
t410 = t470 * t478 + t471 * t474;
t390 = t410 * t465;
t473 = sin(qJ(5));
t553 = t390 * t473;
t336 = -t375 - t553;
t464 = qJD(3) + qJD(5);
t554 = t336 * t464;
t502 = t388 * t473 - t477 * t390;
t555 = t502 * t464;
t475 = sin(qJ(2));
t532 = qJD(2) * t475;
t453 = pkin(1) * t532;
t479 = cos(qJ(2));
t567 = pkin(1) * t479;
t537 = qJD(1) * t453 - qJDD(1) * t567;
t469 = qJ(1) + qJ(2);
t459 = cos(t469);
t558 = g(2) * t459;
t463 = qJDD(1) + qJDD(2);
t566 = pkin(2) * t463;
t575 = t537 - t566 + t558;
t458 = sin(t469);
t444 = g(1) * t458;
t574 = t444 - t558;
t526 = qJDD(1) * t475;
t531 = qJD(2) * t479;
t491 = qJD(1) * t531 + t526;
t394 = t491 * pkin(1) + pkin(7) * t463;
t496 = qJ(4) * t463 + qJD(4) * t465 + t394;
t568 = pkin(1) * t475;
t525 = qJD(1) * t568;
t514 = t465 * t576 + t525;
t501 = qJD(3) * t514;
t319 = qJDD(3) * pkin(3) - t496 * t474 - t478 * t501;
t322 = -t474 * t501 + t496 * t478;
t296 = t471 * t319 - t322 * t470;
t529 = qJD(3) * t474;
t518 = t465 * t529;
t425 = t470 * t518;
t528 = qJD(3) * t478;
t517 = t465 * t528;
t342 = t410 * t463 + t471 * t517 - t425;
t290 = qJDD(3) * pkin(4) - pkin(8) * t342 + t296;
t297 = t470 * t319 + t471 * t322;
t402 = t410 * qJD(3);
t429 = t463 * t543;
t341 = t465 * t402 + t463 * t547 - t429;
t291 = -pkin(8) * t341 + t297;
t377 = t514 * t474;
t371 = qJD(3) * pkin(3) - t377;
t378 = t514 * t478;
t545 = t471 * t378;
t326 = t470 * t371 + t545;
t562 = pkin(8) * t388;
t309 = t326 - t562;
t448 = pkin(3) * t478 + pkin(2);
t533 = qJD(1) * t479;
t524 = pkin(1) * t533;
t387 = -t448 * t465 + qJD(4) - t524;
t343 = pkin(4) * t388 + t387;
t466 = qJ(3) + pkin(9);
t456 = qJ(5) + t466;
t440 = sin(t456);
t527 = qJD(5) * t473;
t441 = cos(t456);
t549 = t441 * t459;
t550 = t441 * t458;
t573 = g(1) * t549 + g(2) * t550 + g(3) * t440 - t473 * t290 - t477 * t291 + t309 * t527 - t343 * t336;
t551 = t440 * t459;
t552 = t440 * t458;
t572 = g(1) * t551 + g(2) * t552 - g(3) * t441 + t477 * t290 - t473 * t291 + t343 * t502;
t462 = qJDD(3) + qJDD(5);
t571 = t462 * MDP(22) + t336 * MDP(18) * t502 + (-t336 ^ 2 + t502 ^ 2) * MDP(19);
t457 = t478 * qJD(4);
t515 = qJD(3) * t576;
t399 = -t474 * t515 + t457;
t400 = -qJD(4) * t474 - t478 * t515;
t540 = -t399 * t470 + t471 * t400 + t410 * t524;
t409 = -t543 + t547;
t539 = t471 * t399 + t470 * t400 + t409 * t524;
t570 = g(1) * t459 + g(2) * t458;
t569 = pkin(3) * t518 + qJDD(4);
t513 = t477 * t341 + t342 * t473;
t295 = -t502 * qJD(5) + t513;
t565 = pkin(2) * t465;
t564 = pkin(3) * t470;
t563 = pkin(3) * t474;
t561 = pkin(8) * t390;
t403 = t409 * qJD(3);
t560 = pkin(8) * t403;
t559 = pkin(8) * t410;
t557 = g(3) * t478;
t365 = t470 * t378;
t325 = t471 * t371 - t365;
t307 = qJD(3) * pkin(4) + t325 - t561;
t556 = t307 * t477;
t548 = t465 * t474;
t542 = t474 * t478;
t447 = pkin(7) + t568;
t541 = -qJ(4) - t447;
t512 = qJD(3) * t541;
t523 = pkin(1) * t531;
t360 = t474 * t512 + t478 * t523 + t457;
t361 = (-qJD(4) - t523) * t474 + t478 * t512;
t321 = t471 * t360 + t470 * t361;
t328 = -t471 * t377 - t365;
t407 = t541 * t474;
t460 = t478 * qJ(4);
t408 = t447 * t478 + t460;
t352 = t470 * t407 + t471 * t408;
t424 = -t524 - t565;
t538 = t424 * t529 + t478 * t444;
t432 = t576 * t474;
t433 = pkin(7) * t478 + t460;
t370 = -t470 * t432 + t471 * t433;
t467 = t474 ^ 2;
t535 = -t478 ^ 2 + t467;
t530 = qJD(3) * t465;
t452 = pkin(3) * t529;
t521 = -qJD(5) * t375 - t473 * t341 + t477 * t342;
t520 = t424 * t528 + t474 * t575;
t519 = t465 * t532;
t374 = pkin(4) * t402 + t452;
t320 = -t360 * t470 + t471 * t361;
t327 = t377 * t470 - t545;
t351 = t471 * t407 - t408 * t470;
t369 = -t471 * t432 - t433 * t470;
t511 = t465 * t525;
t510 = t374 - t525;
t476 = sin(qJ(1));
t480 = cos(qJ(1));
t508 = g(1) * t476 - g(2) * t480;
t406 = t409 * pkin(8);
t347 = -t406 + t370;
t507 = qJD(5) * t347 - t540 - t560;
t346 = t369 - t559;
t398 = t402 * pkin(8);
t506 = -qJD(5) * t346 + t398 - t539;
t505 = -t307 * t473 - t309 * t477;
t329 = t351 - t559;
t330 = -t406 + t352;
t504 = t329 * t477 - t330 * t473;
t503 = t329 * t473 + t330 * t477;
t353 = t477 * t409 + t410 * t473;
t354 = -t409 * t473 + t410 * t477;
t382 = pkin(4) * t409 - t448;
t500 = -t537 + t574;
t442 = pkin(3) * t471 + pkin(4);
t499 = t442 * t473 + t477 * t564;
t498 = t442 * t477 - t473 * t564;
t497 = -t296 * t410 - t297 * t409 + t325 * t403 - t326 * t402 - t570;
t348 = -t448 * t463 + t537 + t569;
t308 = pkin(4) * t341 + t348;
t312 = -t353 * qJD(5) - t402 * t473 - t403 * t477;
t495 = -g(1) * t552 + g(2) * t551 + t308 * t354 + t343 * t312;
t454 = sin(t466);
t494 = t348 * t410 - t387 * t403 - t454 * t574;
t313 = t354 * qJD(5) + t477 * t402 - t403 * t473;
t493 = g(1) * t550 - g(2) * t549 + t308 * t353 + t343 * t313;
t455 = cos(t466);
t492 = t348 * t409 + t387 * t402 + t455 * t574;
t294 = -t390 * t527 + t521;
t489 = -t424 * t465 - t394 + t570;
t481 = qJD(3) ^ 2;
t488 = (-t294 * t353 - t295 * t354 + t312 * t336 + t313 * t502) * MDP(19) + (t294 * t354 - t312 * t502) * MDP(18) + (t312 * t464 + t354 * t462) * MDP(20) + (-t313 * t464 - t353 * t462) * MDP(21) + 0.2e1 * (t463 * t542 - t535 * t530) * MDP(8) + (t463 * t467 + 0.2e1 * t474 * t517) * MDP(7) + (qJDD(3) * t478 - t474 * t481) * MDP(10) + (qJDD(3) * t474 + t478 * t481) * MDP(9) + t463 * MDP(4);
t487 = pkin(7) * t481 - t511 - t566;
t449 = -pkin(2) - t567;
t485 = pkin(1) * t519 + t447 * t481 + t449 * t463;
t484 = -pkin(7) * qJDD(3) + (t524 - t565) * qJD(3);
t483 = -qJD(3) * t523 - qJDD(3) * t447 + t449 * t530;
t482 = -g(1) * (-t448 * t458 + t459 * t576) - g(2) * (t459 * t448 + t458 * t576);
t430 = -t448 - t567;
t422 = t453 + t452;
t373 = t382 - t567;
t364 = t374 + t453;
t359 = pkin(3) * t548 + pkin(4) * t390;
t311 = t328 - t561;
t310 = t327 + t562;
t305 = -t398 + t321;
t304 = t320 + t560;
t1 = [(t296 * t351 + t297 * t352 + t325 * t320 + t326 * t321 + t348 * t430 + t387 * t422 + t482) * MDP(17) + t488 + (-t320 * t390 - t321 * t388 - t341 * t352 - t342 * t351 + t497) * MDP(16) + (-qJD(3) * t321 - qJDD(3) * t352 + t342 * t430 + t390 * t422 + t494) * MDP(15) + t520 * MDP(13) + t508 * MDP(2) + (g(1) * t480 + g(2) * t476) * MDP(3) + t538 * MDP(12) + (-t364 * t502 + t373 * t294 - (t504 * qJD(5) + t304 * t473 + t305 * t477) * t464 - t503 * t462 + t495) * MDP(24) + (qJD(3) * t320 + qJDD(3) * t351 + t341 * t430 + t388 * t422 + t492) * MDP(14) + (-t364 * t336 + t373 * t295 + (-t503 * qJD(5) + t304 * t477 - t305 * t473) * t464 + t504 * t462 + t493) * MDP(23) + ((t463 * t479 - t519) * MDP(5) + (-t463 * t475 - t465 * t531 - t491) * MDP(6) + t508 * MDP(17)) * pkin(1) + qJDD(1) * MDP(1) + ((-t485 - t575) * MDP(12) + t483 * MDP(13)) * t478 + (t483 * MDP(12) + (t485 - t444) * MDP(13)) * t474 + t500 * MDP(5) + t570 * MDP(6); (t297 * t370 + t296 * t369 - t348 * t448 + (t452 - t525) * t387 + t539 * t326 + t540 * t325 + t482) * MDP(17) + t488 + (-t390 * t525 - qJDD(3) * t370 - t342 * t448 + (t390 * t563 - t539) * qJD(3) + t494) * MDP(15) + (-t341 * t370 - t342 * t369 - t539 * t388 - t540 * t390 + t497) * MDP(16) + (t382 * t294 - (t346 * t473 + t347 * t477) * t462 + (t507 * t473 + t506 * t477) * t464 - t510 * t502 + t495) * MDP(24) + (t382 * t295 + (t346 * t477 - t347 * t473) * t462 + (t506 * t473 - t507 * t477) * t464 - t510 * t336 + t493) * MDP(23) + (-t388 * t525 + qJDD(3) * t369 - t341 * t448 + (t388 * t563 + t540) * qJD(3) + t492) * MDP(14) + (t500 + t511) * MDP(5) + ((-t526 + (-qJD(2) + t465) * t533) * pkin(1) + t570) * MDP(6) + (t484 * t474 + (-t487 - t575) * t478 + t538) * MDP(12) + (t484 * t478 + (t487 - t444) * t474 + t520) * MDP(13); qJDD(3) * MDP(11) + (t489 * t474 - t557) * MDP(12) + (g(3) * t474 + t489 * t478) * MDP(13) + (-g(3) * t455 - qJD(3) * t327 - t387 * t390 + t570 * t454 + (qJDD(3) * t471 - t388 * t548) * pkin(3) + t296) * MDP(14) + (g(3) * t454 + qJD(3) * t328 + t387 * t388 + t570 * t455 + (-qJDD(3) * t470 - t390 * t548) * pkin(3) - t297) * MDP(15) + ((t326 + t327) * t390 + (-t325 + t328) * t388 + (-t341 * t470 - t342 * t471) * pkin(3)) * MDP(16) + (-t325 * t327 - t326 * t328 + (-t557 + t296 * t471 + t297 * t470 + (-t387 * t465 + t570) * t474) * pkin(3)) * MDP(17) + (t294 - t554) * MDP(20) + (-t295 - t555) * MDP(21) + (t498 * t462 + t359 * t336 - (t310 * t477 - t311 * t473) * t464 + (-t499 * t464 + t505) * qJD(5) + t572) * MDP(23) + (-t499 * t462 + t359 * t502 + (t310 * t473 + t311 * t477) * t464 + (-t498 * t464 - t556) * qJD(5) + t573) * MDP(24) + (t478 * MDP(10) + t474 * MDP(9)) * t463 + (-MDP(7) * t542 + t535 * MDP(8)) * t465 ^ 2 + t571; -t429 * MDP(14) - t425 * MDP(15) + (-t388 ^ 2 - t390 ^ 2) * MDP(16) + (t325 * t390 + t326 * t388 - t500 + t569) * MDP(17) + (t295 - t555) * MDP(23) + (t294 + t554) * MDP(24) + (MDP(14) * t547 + t410 * MDP(15) - t448 * MDP(17)) * t463 + (0.2e1 * t390 * MDP(14) + (-t388 + t522) * MDP(15)) * qJD(3); (t521 - t554) * MDP(20) + (-t513 - t555) * MDP(21) + (-t505 * t464 + t572) * MDP(23) + ((-t309 * t473 + t556) * t464 + t573) * MDP(24) + (-MDP(20) * t553 + t502 * MDP(21) + t505 * MDP(23) - MDP(24) * t556) * qJD(5) + t571;];
tau = t1;
