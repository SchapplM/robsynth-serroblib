% Calculate vector of inverse dynamics joint torques for
% S6PPRRRP2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:21
% EndTime: 2019-03-08 18:58:30
% DurationCPUTime: 6.28s
% Computational Cost: add. (4759->490), mult. (11935->660), div. (0->0), fcn. (11040->14), ass. (0->211)
t583 = cos(pkin(6));
t440 = qJD(1) * t583 + qJD(2);
t455 = sin(pkin(6));
t453 = sin(pkin(12));
t459 = sin(qJ(3));
t462 = cos(qJ(3));
t456 = cos(pkin(12));
t582 = cos(pkin(7));
t521 = t456 * t582;
t488 = t453 * t462 + t459 * t521;
t482 = t488 * t455;
t454 = sin(pkin(7));
t564 = t454 * t459;
t383 = qJD(1) * t482 + t440 * t564;
t458 = sin(qJ(4));
t461 = cos(qJ(4));
t512 = pkin(4) * t458 - pkin(10) * t461;
t430 = t512 * qJD(4);
t611 = t383 - t430;
t437 = t583 * qJDD(1) + qJDD(2);
t516 = t455 * t521;
t504 = qJD(1) * t516;
t551 = qJD(3) * t459;
t532 = t454 * t551;
t553 = qJD(1) * t455;
t533 = t453 * t553;
t565 = t453 * t459;
t535 = t455 * t565;
t549 = qJD(3) * t462;
t466 = -t462 * (qJDD(1) * t516 + t437 * t454) + qJDD(1) * t535 + t440 * t532 + t504 * t551 + t533 * t549;
t581 = cos(pkin(11));
t511 = t583 * t581;
t580 = sin(pkin(11));
t410 = t453 * t511 + t456 * t580;
t474 = t453 * t580 - t456 * t511;
t523 = t455 * t581;
t603 = t454 * t523 + t474 * t582;
t366 = t410 * t459 + t462 * t603;
t510 = t583 * t580;
t411 = -t453 * t510 + t456 * t581;
t475 = t453 * t581 + t456 * t510;
t522 = t455 * t580;
t602 = -t454 * t522 + t475 * t582;
t368 = t411 * t459 + t462 * t602;
t514 = t462 * t521;
t524 = t454 * t583;
t515 = t462 * t524;
t390 = -t455 * t514 - t515 + t535;
t491 = g(1) * t368 + g(2) * t366 + g(3) * t390;
t610 = t383 * qJD(3) - t466 + t491;
t537 = qJD(3) * qJD(4);
t526 = t461 * t537;
t536 = qJDD(3) * t458;
t609 = -t526 - t536;
t550 = qJD(3) * t461;
t608 = qJD(5) - t550;
t542 = qJD(5) * t458;
t607 = qJD(3) * t542 - qJDD(4);
t457 = sin(qJ(5));
t460 = cos(qJ(5));
t385 = t457 * (qJD(4) * (qJD(5) + t550) + t536) + t607 * t460;
t382 = -t459 * t533 + t462 * (t440 * t454 + t504);
t505 = pkin(4) * t461 + pkin(10) * t458 + pkin(3);
t543 = qJD(5) * t457;
t562 = t457 * t461;
t606 = -t382 * t562 + t460 * t611 - t505 * t543;
t541 = qJD(5) * t460;
t561 = t460 * t461;
t605 = -t382 * t561 - t457 * t611 - t505 * t541;
t381 = qJD(3) * pkin(9) + t383;
t534 = t455 * t456 * t454;
t406 = -qJD(1) * t534 + t440 * t582;
t604 = -t458 * t381 + t406 * t461;
t448 = t461 * qJDD(3);
t601 = -t458 * t537 + t448;
t600 = qJDD(1) * t488;
t599 = MDP(18) + MDP(20);
t598 = MDP(19) - MDP(22);
t367 = t410 * t462 - t459 * t603;
t369 = t411 * t462 - t459 * t602;
t391 = t459 * t524 + t482;
t483 = t582 * t583 - t534;
t370 = t391 * t458 - t461 * t483;
t392 = t454 * t474 - t523 * t582;
t393 = t454 * t475 + t522 * t582;
t597 = g(3) * t370 - g(2) * (-t367 * t458 + t392 * t461) - g(1) * (-t369 * t458 + t393 * t461);
t352 = -qJD(4) * pkin(4) - t604;
t539 = t460 * qJD(4);
t552 = qJD(3) * t458;
t421 = t457 * t552 - t539;
t546 = qJD(4) * t457;
t423 = t460 * t552 + t546;
t330 = pkin(5) * t421 - qJ(6) * t423 + t352;
t419 = qJDD(5) - t601;
t592 = pkin(10) * t419;
t596 = -t330 * t608 + t592;
t594 = t423 ^ 2;
t593 = pkin(5) * t419;
t585 = pkin(10) * qJD(5);
t584 = qJD(3) * pkin(3);
t579 = qJ(6) * t419;
t357 = t461 * t381 + t458 * t406;
t353 = qJD(4) * pkin(10) + t357;
t373 = -qJD(3) * t505 - t382;
t321 = t353 * t460 + t373 * t457;
t319 = qJ(6) * t608 + t321;
t577 = t319 * t608;
t576 = t321 * t608;
t575 = t382 * t421;
t574 = t382 * t423;
t384 = -qJD(5) * t539 + t457 * t607 + t460 * t609;
t573 = t384 * t457;
t405 = -qJDD(1) * t534 + t437 * t582;
t572 = t405 * t458;
t570 = t421 * t423;
t569 = t421 * t608;
t568 = t423 * t608;
t567 = t423 * t460;
t566 = t505 * t460;
t563 = t454 * t462;
t429 = t512 * qJD(3);
t560 = t457 * t429 + t460 * t604;
t540 = qJD(5) * t461;
t545 = qJD(4) * t458;
t559 = qJ(6) * t545 - qJD(6) * t461 + (-t457 * t540 - t458 * t539) * pkin(9) + t605;
t558 = -pkin(5) * t545 + (-t457 * t545 + t460 * t540) * pkin(9) + t606;
t508 = pkin(5) * t457 - qJ(6) * t460;
t557 = -qJD(6) * t457 + t508 * t608 - t357;
t555 = pkin(9) * t561 - t457 * t505;
t451 = t458 ^ 2;
t554 = -t461 ^ 2 + t451;
t548 = qJD(4) * t421;
t547 = qJD(4) * t423;
t544 = qJD(4) * t461;
t320 = -t353 * t457 + t373 * t460;
t538 = qJD(6) - t320;
t531 = t454 * t549;
t530 = t608 * t546;
t529 = t608 * t539;
t528 = t608 * t543;
t487 = t514 - t565;
t346 = qJDD(3) * pkin(9) + (t437 * t459 + t440 * t549) * t454 + (qJD(1) * qJD(3) * t487 + t600) * t455;
t312 = qJDD(4) * pkin(10) + qJD(4) * t604 + t346 * t461 + t572;
t329 = qJD(3) * t430 - qJDD(3) * t505 + t466;
t519 = t457 * t312 - t460 * t329 + t353 * t541 + t373 * t543;
t509 = pkin(5) * t460 + qJ(6) * t457;
t318 = -pkin(5) * t608 + t538;
t507 = t318 * t460 - t319 * t457;
t371 = t391 * t461 + t458 * t483;
t340 = t371 * t460 + t390 * t457;
t339 = t371 * t457 - t390 * t460;
t503 = pkin(4) + t509;
t502 = pkin(9) + t508;
t501 = t458 * t346 + t381 * t544 - t405 * t461 + t406 * t545;
t413 = t458 * t582 + t461 * t564;
t397 = t413 * t457 + t460 * t563;
t398 = t413 * t460 - t457 * t563;
t499 = -t419 * t457 - t541 * t608;
t498 = t419 * t460 - t528;
t496 = t460 * t312 + t457 * t329 - t353 * t543 + t373 * t541;
t322 = -t366 * t562 - t367 * t460;
t324 = -t368 * t562 - t369 * t460;
t354 = -t390 * t562 - t391 * t460;
t495 = g(1) * t324 + g(2) * t322 + g(3) * t354;
t323 = -t366 * t561 + t367 * t457;
t325 = -t368 * t561 + t369 * t457;
t355 = -t390 * t561 + t391 * t457;
t494 = -g(1) * t325 - g(2) * t323 - g(3) * t355;
t336 = t367 * t461 + t392 * t458;
t338 = t369 * t461 + t393 * t458;
t492 = g(1) * t338 + g(2) * t336 + g(3) * t371;
t412 = t458 * t564 - t461 * t582;
t313 = -qJDD(4) * pkin(4) + t501;
t486 = t352 * t608 - t592;
t380 = -t382 - t584;
t481 = -pkin(9) * qJDD(4) + (t380 + t382 - t584) * qJD(4);
t480 = -g(1) * t522 + g(2) * t523 - g(3) * t583;
t476 = -t585 * t608 + t597;
t314 = t336 * t457 - t366 * t460;
t316 = t338 * t457 - t368 * t460;
t473 = g(1) * t316 + g(2) * t314 + g(3) * t339 - t519;
t307 = pkin(5) * t385 + qJ(6) * t384 - qJD(6) * t423 + t313;
t472 = -t307 + t476;
t315 = t336 * t460 + t366 * t457;
t317 = t338 * t460 + t368 * t457;
t468 = -g(1) * t317 - g(2) * t315 - g(3) * t340 + t496;
t467 = t330 * t423 + qJDD(6) - t473;
t463 = qJD(4) ^ 2;
t465 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t463 + t610;
t464 = qJD(3) ^ 2;
t408 = t502 * t458;
t402 = t566 + (pkin(9) * t457 + pkin(5)) * t461;
t401 = -qJ(6) * t461 + t555;
t396 = qJD(4) * t413 + t458 * t531;
t395 = -qJD(4) * t412 + t461 * t531;
t394 = pkin(5) * t423 + qJ(6) * t421;
t387 = t391 * qJD(3);
t386 = (t455 * t487 + t515) * qJD(3);
t375 = (qJD(5) * t509 - qJD(6) * t460) * t458 + t502 * t544;
t365 = -t384 + t569;
t351 = qJD(5) * t398 + t395 * t457 - t460 * t532;
t350 = -qJD(5) * t397 + t395 * t460 + t457 * t532;
t334 = -pkin(5) * t552 - t429 * t460 + t457 * t604;
t333 = qJ(6) * t552 + t560;
t332 = -qJD(4) * t370 + t386 * t461;
t331 = qJD(4) * t371 + t386 * t458;
t309 = -qJD(5) * t339 + t332 * t460 + t387 * t457;
t308 = qJD(5) * t340 + t332 * t457 - t387 * t460;
t306 = qJDD(6) + t519 - t593;
t305 = qJD(6) * t608 + t496 + t579;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t437 * t583 - g(3) + (t453 ^ 2 + t456 ^ 2) * t455 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t387 - qJDD(3) * t390) * MDP(4) + (-qJD(3) * t386 - qJDD(3) * t391) * MDP(5) + (-t390 * t448 - qJD(4) * t331 - qJDD(4) * t370 + (-t387 * t461 + t390 * t545) * qJD(3)) * MDP(11) + (t390 * t536 - qJD(4) * t332 - qJDD(4) * t371 + (t387 * t458 + t390 * t544) * qJD(3)) * MDP(12) + (t308 * t423 - t309 * t421 - t339 * t384 - t340 * t385) * MDP(21) + (t305 * t340 + t306 * t339 + t307 * t370 + t308 * t318 + t309 * t319 + t330 * t331 - g(3)) * MDP(23) + t599 * (-t308 * t608 + t331 * t421 - t339 * t419 + t370 * t385) + t598 * (-t309 * t608 + t331 * t423 - t340 * t419 - t370 * t384); (t480 + t437) * MDP(2) + (-qJD(4) * t396 - qJDD(4) * t412) * MDP(11) + (-qJD(4) * t395 - qJDD(4) * t413) * MDP(12) + (-t350 * t421 + t351 * t423 - t384 * t397 - t385 * t398) * MDP(21) + (t305 * t398 + t306 * t397 + t307 * t412 + t318 * t351 + t319 * t350 + t330 * t396 + t480) * MDP(23) + ((-qJDD(3) * MDP(5) + (-MDP(11) * t461 + MDP(12) * t458 - MDP(4)) * t464) * t459 + (t601 * MDP(11) + MDP(12) * t609 + qJDD(3) * MDP(4) - t464 * MDP(5)) * t462) * t454 + t599 * (-t351 * t608 + t412 * t385 + t396 * t421 - t397 * t419) + t598 * (-t350 * t608 - t384 * t412 + t396 * t423 - t398 * t419); qJDD(3) * MDP(3) + t610 * MDP(4) + (-t437 * t564 + g(1) * t369 + g(2) * t367 + g(3) * t391 - t455 * t600 + (-t440 * t563 - t487 * t553 + t382) * qJD(3)) * MDP(5) + (qJDD(3) * t451 + 0.2e1 * t458 * t526) * MDP(6) + 0.2e1 * (t448 * t458 - t537 * t554) * MDP(7) + (qJDD(4) * t458 + t461 * t463) * MDP(8) + (qJDD(4) * t461 - t458 * t463) * MDP(9) + (t458 * t481 + t461 * t465) * MDP(11) + (-t458 * t465 + t461 * t481) * MDP(12) + (-t384 * t458 * t460 + (-t457 * t542 + t461 * t539) * t423) * MDP(13) + ((-t421 * t460 - t423 * t457) * t544 + (t573 - t385 * t460 + (t421 * t457 - t567) * qJD(5)) * t458) * MDP(14) + ((t384 + t529) * t461 + (t498 + t547) * t458) * MDP(15) + ((t385 - t530) * t461 + (t499 - t548) * t458) * MDP(16) + (-t419 * t461 + t545 * t608) * MDP(17) + (-t419 * t566 - t606 * t608 + (t352 * t546 + (t499 + t548) * pkin(9) + t519) * t461 + (t352 * t541 + t320 * qJD(4) + t313 * t457 - t575 + (t385 + t530) * pkin(9)) * t458 + t494) * MDP(18) + (-t555 * t419 - t605 * t608 + (t352 * t539 + (t528 + t547) * pkin(9) + t496) * t461 + (-t352 * t543 - t321 * qJD(4) + t313 * t460 - t574 + (-t384 + t529) * pkin(9)) * t458 + t495) * MDP(19) + (t375 * t421 + t385 * t408 - t402 * t419 + (t330 * t546 + t306) * t461 - t558 * t608 + (-qJD(4) * t318 + t307 * t457 + t330 * t541 - t575) * t458 + t494) * MDP(20) + (-t384 * t402 - t385 * t401 + t558 * t423 - t559 * t421 + t507 * t544 + (-t305 * t457 + t306 * t460 + (-t318 * t457 - t319 * t460) * qJD(5) + t491) * t458) * MDP(21) + (-t375 * t423 + t384 * t408 + t401 * t419 + (-t330 * t539 - t305) * t461 + t559 * t608 + (qJD(4) * t319 - t307 * t460 + t330 * t543 + t574) * t458 - t495) * MDP(22) + (t305 * t401 + t307 * t408 + t306 * t402 - g(1) * (pkin(5) * t325 + pkin(9) * t369 + qJ(6) * t324) - g(2) * (pkin(5) * t323 + pkin(9) * t367 + qJ(6) * t322) - g(3) * (pkin(5) * t355 + pkin(9) * t391 + qJ(6) * t354) + (-t382 * t458 + t375) * t330 + t559 * t319 + t558 * t318 + t491 * t505) * MDP(23); MDP(8) * t536 + MDP(9) * t448 + qJDD(4) * MDP(10) + (qJD(4) * t357 - t380 * t552 - t501 + t597) * MDP(11) + (-t572 + (-qJD(3) * t380 - t346) * t461 + t492) * MDP(12) + (t567 * t608 - t573) * MDP(13) + ((-t384 - t569) * t460 + (-t385 - t568) * t457) * MDP(14) + ((-t423 * t458 - t561 * t608) * qJD(3) - t499) * MDP(15) + ((t421 * t458 + t562 * t608) * qJD(3) + t498) * MDP(16) - t608 * MDP(17) * t552 + (-t320 * t552 - pkin(4) * t385 - t357 * t421 + (t604 * t608 + t486) * t457 + (-t313 - (t429 + t585) * t608 + t597) * t460) * MDP(18) + (pkin(4) * t384 + t560 * t608 + t321 * t552 - t357 * t423 + t486 * t460 + (t313 - t476) * t457) * MDP(19) + (t318 * t552 + t334 * t608 - t385 * t503 + t557 * t421 - t457 * t596 + t472 * t460) * MDP(20) + (t333 * t421 - t334 * t423 + (t305 + t608 * t318 + (qJD(5) * t423 - t385) * pkin(10)) * t460 + (t306 - t577 + (qJD(5) * t421 - t384) * pkin(10)) * t457 - t492) * MDP(21) + (-t319 * t552 - t333 * t608 - t384 * t503 - t557 * t423 + t472 * t457 + t460 * t596) * MDP(22) + (-t318 * t334 - t319 * t333 + t557 * t330 + (qJD(5) * t507 + t305 * t460 + t306 * t457 - t492) * pkin(10) + (-t307 + t597) * t503) * MDP(23) + (-MDP(6) * t458 * t461 + MDP(7) * t554) * t464; MDP(13) * t570 + (-t421 ^ 2 + t594) * MDP(14) + t365 * MDP(15) + (-t385 + t568) * MDP(16) + t419 * MDP(17) + (-t352 * t423 + t473 + t576) * MDP(18) + (t320 * t608 + t352 * t421 - t468) * MDP(19) + (-t394 * t421 - t467 + t576 + 0.2e1 * t593) * MDP(20) + (pkin(5) * t384 - qJ(6) * t385 + (t319 - t321) * t423 + (t318 - t538) * t421) * MDP(21) + (0.2e1 * t579 - t330 * t421 + t394 * t423 - (-0.2e1 * qJD(6) + t320) * t608 + t468) * MDP(22) + (t305 * qJ(6) - t306 * pkin(5) - t330 * t394 - t318 * t321 - g(1) * (-pkin(5) * t316 + qJ(6) * t317) - g(2) * (-pkin(5) * t314 + qJ(6) * t315) - g(3) * (-pkin(5) * t339 + qJ(6) * t340) + t538 * t319) * MDP(23); (-t419 + t570) * MDP(20) + t365 * MDP(21) + (-t608 ^ 2 - t594) * MDP(22) + (t467 - t577 - t593) * MDP(23);];
tau  = t1;
