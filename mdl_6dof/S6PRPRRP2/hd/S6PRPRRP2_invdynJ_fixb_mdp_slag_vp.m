% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:18
% EndTime: 2019-03-08 20:03:25
% DurationCPUTime: 5.91s
% Computational Cost: add. (3790->478), mult. (8667->629), div. (0->0), fcn. (7161->12), ass. (0->208)
t460 = cos(pkin(11));
t465 = sin(qJ(2));
t459 = sin(pkin(6));
t537 = qJD(1) * t459;
t514 = t465 * t537;
t429 = t460 * t514;
t457 = sin(pkin(11));
t468 = cos(qJ(2));
t513 = t468 * t537;
t391 = t457 * t513 + t429;
t464 = sin(qJ(4));
t467 = cos(qJ(4));
t502 = pkin(4) * t464 - pkin(9) * t467;
t423 = t502 * qJD(4);
t599 = t391 - t423;
t536 = qJD(2) * t459;
t510 = qJD(1) * t536;
t521 = qJDD(1) * t459;
t598 = t465 * t521 + t468 * t510;
t534 = qJD(2) * t467;
t597 = qJD(5) - t534;
t527 = qJD(5) * t464;
t596 = qJD(2) * t527 - qJDD(4);
t463 = sin(qJ(5));
t466 = cos(qJ(5));
t520 = qJDD(2) * t464;
t374 = t463 * (qJD(4) * (qJD(5) + t534) + t520) + t596 * t466;
t428 = t457 * t514;
t394 = t460 * t513 - t428;
t553 = t463 * t467;
t595 = -t394 * t553 + t466 * t599;
t494 = pkin(4) * t467 + pkin(9) * t464 + pkin(3);
t586 = pkin(2) * t460;
t410 = -t494 - t586;
t526 = qJD(5) * t466;
t550 = t466 * t467;
t594 = -t394 * t550 + t410 * t526 - t463 * t599;
t425 = qJD(2) * pkin(2) + t513;
t386 = t457 * t425 + t429;
t384 = qJD(2) * pkin(8) + t386;
t462 = cos(pkin(6));
t442 = qJD(1) * t462 + qJD(3);
t593 = -t464 * t384 + t442 * t467;
t447 = pkin(2) * t457 + pkin(8);
t540 = t463 * t410 + t447 * t550;
t458 = sin(pkin(10));
t461 = cos(pkin(10));
t554 = t462 * t468;
t592 = -t458 * t554 - t461 * t465;
t519 = MDP(18) + MDP(20);
t518 = MDP(19) - MDP(22);
t549 = t468 * t460;
t413 = t457 * t465 - t549;
t399 = t413 * t462;
t495 = t457 * t468 + t460 * t465;
t368 = -t399 * t461 - t458 * t495;
t371 = t399 * t458 - t461 * t495;
t559 = t459 * t465;
t397 = t457 * t559 - t459 * t549;
t482 = -g(1) * t371 - g(2) * t368 + g(3) * t397;
t400 = t495 * t462;
t369 = t400 * t461 - t413 * t458;
t370 = t400 * t458 + t413 * t461;
t398 = t495 * t459;
t381 = t398 * t464 - t462 * t467;
t558 = t459 * t467;
t591 = g(3) * t381 - g(2) * (-t369 * t464 - t461 * t558) - g(1) * (t370 * t464 + t458 * t558);
t355 = -qJD(4) * pkin(4) - t593;
t524 = t466 * qJD(4);
t535 = qJD(2) * t464;
t415 = t463 * t535 - t524;
t532 = qJD(4) * t463;
t417 = t466 * t535 + t532;
t326 = pkin(5) * t415 - qJ(6) * t417 + t355;
t453 = t467 * qJDD(2);
t522 = qJD(2) * qJD(4);
t412 = t464 * t522 + qJDD(5) - t453;
t584 = pkin(9) * t412;
t590 = -t326 * t597 + t584;
t589 = t466 * MDP(19) + t463 * t519;
t587 = t417 ^ 2;
t585 = pkin(5) * t412;
t577 = pkin(9) * qJD(5);
t576 = qJ(6) * t412;
t360 = t467 * t384 + t464 * t442;
t356 = qJD(4) * pkin(9) + t360;
t385 = t425 * t460 - t428;
t366 = -qJD(2) * t494 - t385;
t318 = t356 * t466 + t366 * t463;
t312 = qJ(6) * t597 + t318;
t575 = t312 * t597;
t574 = t318 * t597;
t573 = t394 * t415;
t572 = t394 * t417;
t571 = t415 * t417;
t570 = t415 * t597;
t569 = t417 * t597;
t568 = t417 * t464;
t567 = t417 * t466;
t439 = t462 * qJDD(1) + qJDD(3);
t565 = t439 * t464;
t563 = t447 * t463;
t562 = t447 * t466;
t561 = t458 * t465;
t560 = t459 * t464;
t557 = t459 * t468;
t555 = t462 * t465;
t552 = t464 * t466;
t551 = t466 * t412;
t548 = qJDD(1) - g(3);
t528 = qJD(5) * t463;
t531 = qJD(4) * t464;
t547 = (-t447 * t528 - qJD(6)) * t467 + (qJ(6) - t562) * t531 + t594;
t506 = pkin(5) + t563;
t546 = t540 * qJD(5) - t506 * t531 + t595;
t422 = t502 * qJD(2);
t545 = t463 * t422 + t466 * t593;
t512 = t467 * t524;
t544 = -t374 * t552 - t415 * t512;
t543 = t464 * t551 + t512 * t597;
t500 = pkin(5) * t463 - qJ(6) * t466;
t541 = -qJD(6) * t463 + t500 * t597 - t360;
t455 = t464 ^ 2;
t539 = -t467 ^ 2 + t455;
t538 = MDP(21) * t463;
t533 = qJD(4) * t415;
t530 = qJD(4) * t467;
t529 = qJD(5) * t415;
t317 = -t356 * t463 + t366 * t466;
t523 = qJD(6) - t317;
t516 = t461 * t554;
t438 = t468 * t521;
t395 = qJDD(2) * pkin(2) - t465 * t510 + t438;
t353 = t457 * t395 + t460 * t598;
t349 = qJDD(2) * pkin(8) + t353;
t309 = qJDD(4) * pkin(9) + qJD(4) * t593 + t349 * t467 + t565;
t352 = t395 * t460 - t457 * t598;
t325 = qJD(2) * t423 - qJDD(2) * t494 - t352;
t515 = t466 * t309 + t463 * t325 + t366 * t526;
t511 = t463 * t527;
t509 = t467 * t522;
t504 = t463 * t309 - t466 * t325 + t356 * t526 + t366 * t528;
t501 = pkin(5) * t466 + qJ(6) * t463;
t487 = -t356 * t528 + t515;
t302 = qJD(6) * t597 + t487 + t576;
t303 = qJDD(6) + t504 - t585;
t499 = t302 * t466 + t303 * t463;
t311 = -pkin(5) * t597 + t523;
t498 = t311 * t466 - t312 * t463;
t497 = t311 * t463 + t312 * t466;
t382 = t398 * t467 + t462 * t464;
t338 = t382 * t466 + t397 * t463;
t337 = t382 * t463 - t397 * t466;
t493 = pkin(4) + t501;
t492 = t464 * t349 + t384 * t530 - t439 * t467 + t442 * t531;
t490 = t447 + t500;
t489 = -t412 * t463 - t526 * t597;
t488 = -g(3) * t462 + (-g(1) * t458 + g(2) * t461) * t459;
t319 = t368 * t553 - t369 * t466;
t321 = t370 * t466 + t371 * t553;
t350 = -t397 * t553 - t398 * t466;
t486 = g(1) * t321 + g(2) * t319 + g(3) * t350;
t320 = t368 * t550 + t369 * t463;
t322 = -t370 * t463 + t371 * t550;
t351 = -t397 * t550 + t398 * t463;
t485 = -g(1) * t322 - g(2) * t320 - g(3) * t351;
t343 = t369 * t467 - t461 * t560;
t345 = -t370 * t467 + t458 * t560;
t483 = g(1) * t345 + g(2) * t343 + g(3) * t382;
t310 = -qJDD(4) * pkin(4) + t492;
t480 = t355 * t597 - t584;
t373 = -qJD(5) * t524 + (-t509 - t520) * t466 + t596 * t463;
t383 = -qJD(2) * pkin(3) - t385;
t448 = -pkin(3) - t586;
t479 = -qJDD(4) * t447 + (qJD(2) * t448 + t383 + t394) * qJD(4);
t477 = -t577 * t597 + t591;
t313 = t343 * t463 + t368 * t466;
t315 = t345 * t463 + t371 * t466;
t476 = g(1) * t315 + g(2) * t313 + g(3) * t337 - t504;
t475 = -g(1) * t592 - g(3) * t557;
t304 = pkin(5) * t374 + qJ(6) * t373 - qJD(6) * t417 + t310;
t474 = -t304 + t477;
t314 = t343 * t466 - t368 * t463;
t316 = t345 * t466 - t371 * t463;
t473 = -g(1) * t316 - g(2) * t314 - g(3) * t338 + t487;
t472 = t326 * t417 + qJDD(6) - t476;
t469 = qJD(4) ^ 2;
t471 = -qJD(2) * t391 + t447 * t469 - t352 - t482 + (-pkin(3) + t448) * qJDD(2);
t470 = qJD(2) ^ 2;
t434 = pkin(2) * t516;
t433 = qJDD(4) * t467 - t464 * t469;
t432 = qJDD(4) * t464 + t467 * t469;
t404 = t417 * t531;
t393 = t413 * t536;
t392 = qJD(2) * t398;
t388 = t490 * t464;
t380 = pkin(5) * t417 + qJ(6) * t415;
t376 = -t410 * t466 + t467 * t506;
t375 = -qJ(6) * t467 + t540;
t346 = (qJD(5) * t501 - qJD(6) * t466) * t464 + t490 * t530;
t336 = -t373 + t570;
t335 = -qJD(4) * t381 - t393 * t467;
t334 = qJD(4) * t382 - t393 * t464;
t330 = -pkin(5) * t535 - t422 * t466 + t463 * t593;
t329 = qJ(6) * t535 + t545;
t306 = -qJD(5) * t337 + t335 * t466 + t392 * t463;
t305 = qJD(5) * t338 + t335 * t463 - t392 * t466;
t1 = [t548 * MDP(1) + (-t352 * t397 + t353 * t398 - t385 * t392 - t386 * t393 + t439 * t462 - g(3)) * MDP(5) + (-qJD(4) * t334 - qJDD(4) * t381 - t397 * t453) * MDP(11) + (-qJD(4) * t335 - qJDD(4) * t382 + t397 * t520) * MDP(12) + (t305 * t417 - t306 * t415 - t337 * t373 - t338 * t374) * MDP(21) + (t302 * t338 + t303 * t337 + t304 * t381 + t305 * t311 + t306 * t312 + t326 * t334 - g(3)) * MDP(23) + ((-t392 * t467 + t397 * t531) * MDP(11) + (t392 * t464 + t397 * t530) * MDP(12)) * qJD(2) + ((qJDD(2) * t468 - t465 * t470) * MDP(3) + (-qJDD(2) * t465 - t468 * t470) * MDP(4)) * t459 + t519 * (-t305 * t597 + t334 * t415 - t337 * t412 + t381 * t374) + t518 * (-t306 * t597 + t334 * t417 - t338 * t412 - t373 * t381); qJDD(2) * MDP(2) + (t438 - g(2) * (t516 - t561) + t475) * MDP(3) + (-g(1) * (t458 * t555 - t461 * t468) - g(2) * (-t458 * t468 - t461 * t555) - t548 * t559) * MDP(4) + (-g(2) * t434 + t385 * t391 - t386 * t394 + (g(2) * t561 + t352 * t460 + t353 * t457 + t475) * pkin(2)) * MDP(5) + (qJDD(2) * t455 + 0.2e1 * t464 * t509) * MDP(6) + 0.2e1 * (t453 * t464 - t522 * t539) * MDP(7) + t432 * MDP(8) + t433 * MDP(9) + (t464 * t479 - t467 * t471) * MDP(11) + (t464 * t471 + t467 * t479) * MDP(12) + (-t373 * t552 + (-t511 + t512) * t417) * MDP(13) + (-t526 * t568 + (-t417 * t530 + (t373 + t529) * t464) * t463 + t544) * MDP(14) + (t373 * t467 - t511 * t597 + t404 + t543) * MDP(15) + ((-t532 * t597 + t374) * t467 + (t489 - t533) * t464) * MDP(16) + (-t412 * t467 + t531 * t597) * MDP(17) + (t410 * t551 - (t410 * t528 + t595) * t597 + (t355 * t532 + (t489 + t533) * t447 + t504) * t467 + (t355 * t526 + t310 * t463 + t447 * t374 - t573 + (t563 * t597 + t317) * qJD(4)) * t464 + t485) * MDP(18) + (-t540 * t412 - t594 * t597 + ((t447 * t597 - t356) * t528 + (t355 * t466 + t417 * t447) * qJD(4) + t515) * t467 + (-t355 * t528 + t310 * t466 - t447 * t373 - t572 + (t562 * t597 - t318) * qJD(4)) * t464 + t486) * MDP(19) + (t346 * t415 + t374 * t388 - t376 * t412 + (t326 * t532 + t303) * t467 - t546 * t597 + (-qJD(4) * t311 + t304 * t463 + t326 * t526 - t573) * t464 + t485) * MDP(20) + (-t373 * t376 - t374 * t375 + t546 * t417 - t547 * t415 + t498 * t530 + (-qJD(5) * t497 - t302 * t463 + t303 * t466 + t482) * t464) * MDP(21) + (-t346 * t417 + t373 * t388 + t375 * t412 + (-t326 * t524 - t302) * t467 + t547 * t597 + (qJD(4) * t312 - t304 * t466 + t326 * t528 + t572) * t464 - t486) * MDP(22) + (t302 * t375 + t304 * t388 + t303 * t376 - g(1) * (pkin(2) * t592 + pkin(5) * t322 - pkin(8) * t370 + qJ(6) * t321) - g(2) * (-pkin(2) * t561 + pkin(5) * t320 + pkin(8) * t369 + qJ(6) * t319 + t434) - g(3) * (pkin(2) * t557 + pkin(5) * t351 + pkin(8) * t398 + qJ(6) * t350) + (-t394 * t464 + t346) * t326 + t547 * t312 + t546 * t311 + t482 * t494) * MDP(23); (t488 + t439) * MDP(5) + t433 * MDP(11) - t432 * MDP(12) + t404 * MDP(19) + t544 * MDP(21) + t543 * MDP(22) + t488 * MDP(23) + t519 * t415 * t531 + (-t304 * MDP(23) - t519 * t374 + t518 * t373 + (t497 * MDP(23) + t417 * t538 - t589 * t597) * qJD(4)) * t467 + (-t373 * t538 - qJD(4) * t417 * MDP(22) + (qJD(4) * t326 + t499) * MDP(23) - t589 * t412 + ((t415 * t463 + t567) * MDP(21) + t498 * MDP(23) - (-t463 * t518 + t466 * t519) * t597) * qJD(5)) * t464; MDP(8) * t520 + MDP(9) * t453 + qJDD(4) * MDP(10) + (qJD(4) * t360 - t383 * t535 - t492 + t591) * MDP(11) + (-t565 + (-qJD(2) * t383 - t349) * t467 + t483) * MDP(12) + (-t373 * t463 + t567 * t597) * MDP(13) + ((-t373 - t570) * t466 + (-t374 - t569) * t463) * MDP(14) + ((-t550 * t597 - t568) * qJD(2) - t489) * MDP(15) + (-t597 * t528 + t551 + (t415 * t464 + t553 * t597) * qJD(2)) * MDP(16) - t597 * MDP(17) * t535 + (-t317 * t535 - pkin(4) * t374 - t360 * t415 + (t593 * t597 + t480) * t463 + (-t310 - (t422 + t577) * t597 + t591) * t466) * MDP(18) + (pkin(4) * t373 + t545 * t597 + t318 * t535 - t360 * t417 + t480 * t466 + (t310 - t477) * t463) * MDP(19) + (t311 * t535 + t330 * t597 - t374 * t493 + t541 * t415 - t463 * t590 + t474 * t466) * MDP(20) + (t329 * t415 - t330 * t417 + (t302 + t597 * t311 + (qJD(5) * t417 - t374) * pkin(9)) * t466 + (t303 - t575 + (-t373 + t529) * pkin(9)) * t463 - t483) * MDP(21) + (-t312 * t535 - t329 * t597 - t373 * t493 - t541 * t417 + t474 * t463 + t466 * t590) * MDP(22) + (-t311 * t330 - t312 * t329 + t541 * t326 + (qJD(5) * t498 - t483 + t499) * pkin(9) + (-t304 + t591) * t493) * MDP(23) + (-MDP(6) * t464 * t467 + MDP(7) * t539) * t470; MDP(13) * t571 + (-t415 ^ 2 + t587) * MDP(14) + t336 * MDP(15) + (-t374 + t569) * MDP(16) + t412 * MDP(17) + (-t355 * t417 + t476 + t574) * MDP(18) + (t317 * t597 + t355 * t415 - t473) * MDP(19) + (-t380 * t415 - t472 + t574 + 0.2e1 * t585) * MDP(20) + (pkin(5) * t373 - qJ(6) * t374 + (t312 - t318) * t417 + (t311 - t523) * t415) * MDP(21) + (0.2e1 * t576 - t326 * t415 + t380 * t417 - (-0.2e1 * qJD(6) + t317) * t597 + t473) * MDP(22) + (t302 * qJ(6) - t303 * pkin(5) - t326 * t380 - t311 * t318 - g(1) * (-pkin(5) * t315 + qJ(6) * t316) - g(2) * (-pkin(5) * t313 + qJ(6) * t314) - g(3) * (-pkin(5) * t337 + qJ(6) * t338) + t523 * t312) * MDP(23); (-t412 + t571) * MDP(20) + t336 * MDP(21) + (-t597 ^ 2 - t587) * MDP(22) + (t472 - t575 - t585) * MDP(23);];
tau  = t1;
