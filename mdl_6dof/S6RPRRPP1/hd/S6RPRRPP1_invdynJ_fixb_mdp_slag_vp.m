% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:30:05
% EndTime: 2019-03-09 04:30:13
% DurationCPUTime: 6.03s
% Computational Cost: add. (5890->535), mult. (12228->675), div. (0->0), fcn. (7987->14), ass. (0->222)
t536 = sin(qJ(4));
t539 = cos(qJ(4));
t598 = t539 * qJD(3);
t537 = sin(qJ(3));
t609 = qJD(1) * t537;
t467 = t536 * t609 - t598;
t605 = qJD(3) * t536;
t469 = t539 * t609 + t605;
t531 = sin(pkin(10));
t533 = cos(pkin(10));
t415 = t533 * t467 + t469 * t531;
t540 = cos(qJ(3));
t608 = qJD(1) * t540;
t502 = -qJD(4) + t608;
t663 = t415 * t502;
t462 = t531 * t539 + t533 * t536;
t449 = t462 * qJD(4);
t616 = t462 * t608 - t449;
t602 = qJD(4) * t536;
t581 = t531 * t602;
t586 = t536 * t608;
t600 = qJD(4) * t539;
t629 = t533 * t539;
t615 = -t531 * t586 - t533 * t600 + t608 * t629 + t581;
t532 = sin(pkin(9));
t510 = pkin(1) * t532 + pkin(7);
t484 = t510 * qJDD(1);
t662 = -qJD(2) * qJD(3) - t484;
t515 = pkin(4) * t539 + pkin(3);
t643 = qJ(5) + pkin(8);
t574 = t540 * t515 + t643 * t537;
t603 = qJD(3) * t540;
t584 = t536 * t603;
t661 = t537 * t600 + t584;
t486 = t510 * qJD(1);
t441 = qJD(2) * t540 - t537 * t486;
t660 = qJD(3) * t441;
t601 = qJD(4) * t537;
t659 = -qJD(1) * t601 + qJDD(3);
t658 = -pkin(2) - t574;
t594 = qJDD(1) * t537;
t411 = t536 * (qJD(3) * (qJD(4) + t608) + t594) - t659 * t539;
t565 = -t467 * t531 + t533 * t469;
t657 = t565 ^ 2;
t528 = qJ(1) + pkin(9);
t517 = sin(t528);
t519 = cos(t528);
t569 = g(1) * t519 + g(2) * t517;
t656 = t537 * t569;
t570 = pkin(3) * t537 - pkin(8) * t540;
t471 = t570 * qJD(1);
t453 = t539 * t471;
t623 = t539 * t540;
t563 = pkin(4) * t537 - qJ(5) * t623;
t390 = qJD(1) * t563 - t441 * t536 + t453;
t617 = t539 * t441 + t536 * t471;
t394 = -qJ(5) * t586 + t617;
t576 = qJD(4) * t643;
t599 = qJD(5) * t539;
t446 = -t536 * t576 + t599;
t554 = -qJD(5) * t536 - t539 * t576;
t620 = (-t390 + t554) * t533 + (t394 - t446) * t531;
t442 = t537 * qJD(2) + t540 * t486;
t655 = -t442 + (-t586 + t602) * pkin(4);
t654 = MDP(19) + MDP(22);
t561 = qJDD(2) * t540 - t486 * t603 + t537 * t662;
t398 = -qJDD(3) * pkin(3) - t561;
t526 = g(3) * t540;
t549 = -t526 + t656;
t653 = qJD(4) * pkin(8) * t502 - t398 + t549;
t428 = -qJD(3) * pkin(3) - t441;
t409 = pkin(4) * t467 + qJD(5) + t428;
t362 = pkin(5) * t415 - qJ(6) * t565 + t409;
t527 = qJ(4) + pkin(10);
t516 = sin(t527);
t518 = cos(t527);
t633 = t517 * t540;
t422 = t516 * t633 + t518 * t519;
t630 = t519 * t540;
t424 = t516 * t630 - t517 * t518;
t429 = qJD(3) * pkin(8) + t442;
t534 = cos(pkin(9));
t512 = -pkin(1) * t534 - pkin(2);
t456 = -pkin(3) * t540 - pkin(8) * t537 + t512;
t432 = t456 * qJD(1);
t393 = t429 * t539 + t432 * t536;
t397 = qJDD(3) * pkin(8) + qJDD(2) * t537 + t484 * t540 + t660;
t472 = t570 * qJD(3);
t413 = qJD(1) * t472 + qJDD(1) * t456;
t407 = t539 * t413;
t596 = qJD(1) * qJD(3);
t579 = t540 * t596;
t410 = qJD(4) * t598 + (t579 + t594) * t539 + t659 * t536;
t522 = t540 * qJDD(1);
t460 = t537 * t596 + qJDD(4) - t522;
t344 = pkin(4) * t460 - qJ(5) * t410 - qJD(4) * t393 - qJD(5) * t469 - t397 * t536 + t407;
t589 = t539 * t397 + t536 * t413 + t432 * t600;
t558 = -t429 * t602 + t589;
t349 = -qJ(5) * t411 - qJD(5) * t467 + t558;
t339 = t533 * t344 - t531 * t349;
t577 = -qJDD(6) + t339;
t644 = g(3) * t537;
t652 = g(1) * t424 + g(2) * t422 - t362 * t565 + t516 * t644 + t577;
t651 = t540 * t569 + t644;
t649 = pkin(5) * t460;
t648 = g(1) * t517;
t645 = g(2) * t519;
t386 = -qJ(5) * t467 + t393;
t383 = t533 * t386;
t392 = -t429 * t536 + t539 * t432;
t385 = -qJ(5) * t469 + t392;
t356 = t385 * t531 + t383;
t642 = t356 * t565;
t641 = t386 * t531;
t640 = t410 * t536;
t639 = t467 * t502;
t638 = t469 * t502;
t637 = t502 * t539;
t636 = t510 * t536;
t635 = t517 * t536;
t634 = t517 * t539;
t632 = t519 * t536;
t631 = t519 * t539;
t627 = t643 * t540;
t626 = t536 * t537;
t625 = t536 * t540;
t624 = t537 * t539;
t622 = qJDD(2) - g(3);
t340 = t531 * t344 + t533 * t349;
t470 = t510 * t623;
t604 = qJD(3) * t537;
t613 = t539 * t472 + t604 * t636;
t368 = -t537 * t599 + t563 * qJD(3) + (-t470 + (qJ(5) * t537 - t456) * t536) * qJD(4) + t613;
t606 = qJD(3) * t510;
t614 = t456 * t600 + t536 * t472;
t373 = (-qJ(5) * qJD(4) - t606) * t624 + (-qJD(5) * t537 + (-qJ(5) * qJD(3) - qJD(4) * t510) * t540) * t536 + t614;
t351 = t531 * t368 + t533 * t373;
t382 = -pkin(4) * t502 + t385;
t355 = t531 * t382 + t383;
t621 = -pkin(5) * t616 + qJ(6) * t615 - qJD(6) * t462 + t655;
t364 = t531 * t390 + t533 * t394;
t619 = pkin(5) * t609 - t620;
t359 = qJ(6) * t609 + t364;
t403 = t533 * t446 + t531 * t554;
t618 = t403 - t359;
t444 = t539 * t456;
t401 = -qJ(5) * t624 + t444 + (-pkin(4) - t636) * t540;
t612 = t536 * t456 + t470;
t408 = -qJ(5) * t626 + t612;
t372 = t531 * t401 + t533 * t408;
t504 = pkin(4) * t626;
t611 = t537 * t510 + t504;
t529 = t537 ^ 2;
t610 = -t540 ^ 2 + t529;
t487 = qJD(1) * t512;
t357 = t385 * t533 - t641;
t597 = qJD(6) - t357;
t591 = t519 * t625;
t590 = t460 * qJ(6) + t340;
t587 = pkin(4) * t661 + t510 * t603;
t585 = t502 * t605;
t583 = t540 * t598;
t580 = t643 * t536;
t376 = t410 * t531 + t533 * t411;
t575 = t502 * t510 + t429;
t573 = -qJD(4) * t432 - t397;
t493 = t537 * t648;
t571 = -t537 * t645 + t493;
t538 = sin(qJ(1));
t541 = cos(qJ(1));
t568 = g(1) * t538 - g(2) * t541;
t567 = pkin(5) * t518 + qJ(6) * t516;
t350 = t368 * t533 - t373 * t531;
t354 = t382 * t533 - t641;
t371 = t401 * t533 - t408 * t531;
t377 = t410 * t533 - t411 * t531;
t434 = t517 * t625 + t631;
t560 = -t460 * t536 + t502 * t600;
t559 = -t460 * t539 - t502 * t602;
t556 = -qJD(1) * t487 + t569;
t555 = -t536 * t601 + t583;
t553 = -pkin(8) * t460 - t428 * t502;
t552 = t541 * pkin(1) + pkin(4) * t635 + t517 * pkin(7) - t519 * t658;
t542 = qJD(3) ^ 2;
t551 = 0.2e1 * qJDD(1) * t512 + t510 * t542 + t645;
t550 = 0.2e1 * qJD(3) * t487 - qJDD(3) * t510;
t547 = -pkin(1) * t538 + pkin(4) * t632 + t519 * pkin(7) + t517 * t658;
t375 = pkin(4) * t411 + qJDD(5) + t398;
t488 = t643 * t539;
t419 = t488 * t531 + t533 * t580;
t420 = t533 * t488 - t531 * t580;
t545 = -t420 * t376 + t377 * t419 - t403 * t415 - t651;
t341 = pkin(5) * t376 - qJ(6) * t377 - qJD(6) * t565 + t375;
t511 = -pkin(4) * t533 - pkin(5);
t505 = pkin(4) * t531 + qJ(6);
t491 = pkin(4) * t634;
t478 = qJDD(3) * t540 - t537 * t542;
t477 = qJDD(3) * t537 + t540 * t542;
t461 = t531 * t536 - t629;
t445 = t469 * t604;
t440 = -t531 * t626 + t533 * t624;
t439 = t462 * t537;
t437 = t519 * t623 + t635;
t436 = -t591 + t634;
t435 = -t517 * t623 + t632;
t425 = t516 * t517 + t518 * t630;
t423 = -t519 * t516 + t518 * t633;
t412 = pkin(5) * t461 - qJ(6) * t462 - t515;
t400 = t449 * t537 + t531 * t584 - t533 * t583;
t399 = -t531 * t583 - t533 * t661 + t537 * t581;
t387 = pkin(5) * t439 - qJ(6) * t440 + t611;
t374 = pkin(4) * t469 + pkin(5) * t565 + qJ(6) * t415;
t369 = pkin(5) * t540 - t371;
t367 = -qJ(6) * t540 + t372;
t358 = -pkin(5) * t399 + qJ(6) * t400 - qJD(6) * t440 + t587;
t353 = -qJ(6) * t502 + t355;
t352 = pkin(5) * t502 + qJD(6) - t354;
t348 = -pkin(5) * t604 - t350;
t345 = qJ(6) * t604 - qJD(6) * t540 + t351;
t338 = -t577 - t649;
t337 = -qJD(6) * t502 + t590;
t1 = [qJDD(1) * MDP(1) + t568 * MDP(2) + (g(1) * t541 + g(2) * t538) * MDP(3) + (t568 + (t532 ^ 2 + t534 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t529 + 0.2e1 * t537 * t579) * MDP(5) + 0.2e1 * (t522 * t537 - t596 * t610) * MDP(6) + t477 * MDP(7) + t478 * MDP(8) + (t550 * t537 + (-t551 + t648) * t540) * MDP(10) + (t537 * t551 + t540 * t550 - t493) * MDP(11) + (t410 * t624 + t469 * t555) * MDP(12) + ((-t467 * t539 - t469 * t536) * t603 + (-t640 - t411 * t539 + (t467 * t536 - t469 * t539) * qJD(4)) * t537) * MDP(13) + (-t410 * t540 + t460 * t624 - t502 * t555 + t445) * MDP(14) + ((t411 + t585) * t540 + (-qJD(3) * t467 + t560) * t537) * MDP(15) + (-t460 * t540 - t502 * t604) * MDP(16) + (-(-t456 * t602 + t613) * t502 + t444 * t460 - g(1) * t435 - g(2) * t437 + (t467 * t606 - t407 + t575 * t600 + (qJD(3) * t428 - t460 * t510 - t573) * t536) * t540 + (qJD(3) * t392 + t398 * t536 + t411 * t510 + t428 * t600) * t537) * MDP(17) + (t614 * t502 - t612 * t460 - g(1) * t434 - g(2) * t436 + (-t575 * t602 + (t428 * t539 + t469 * t510) * qJD(3) + t589) * t540 + (-t428 * t602 + t398 * t539 + t510 * t410 + (-t510 * t637 - t393) * qJD(3)) * t537) * MDP(18) + (-t339 * t440 - t340 * t439 - t350 * t565 - t351 * t415 + t354 * t400 + t355 * t399 - t371 * t377 - t372 * t376 + t571) * MDP(19) + (-g(1) * t547 - g(2) * t552 + t339 * t371 + t340 * t372 + t354 * t350 + t355 * t351 + t375 * t611 + t409 * t587) * MDP(20) + (g(1) * t423 - g(2) * t425 + t338 * t540 + t341 * t439 + t348 * t502 - t352 * t604 + t358 * t415 - t362 * t399 - t369 * t460 + t376 * t387) * MDP(21) + (-t337 * t439 + t338 * t440 - t345 * t415 + t348 * t565 - t352 * t400 + t353 * t399 - t367 * t376 + t369 * t377 + t571) * MDP(22) + (g(1) * t422 - g(2) * t424 - t337 * t540 - t341 * t440 - t345 * t502 + t353 * t604 - t358 * t565 + t362 * t400 + t367 * t460 - t377 * t387) * MDP(23) + (t337 * t367 + t353 * t345 + t341 * t387 + t362 * t358 + t338 * t369 + t352 * t348 - g(1) * (-pkin(5) * t423 - qJ(6) * t422 + t547) - g(2) * (pkin(5) * t425 + qJ(6) * t424 + t552)) * MDP(24); t622 * MDP(4) + t478 * MDP(10) - t477 * MDP(11) + t445 * MDP(18) + (-t339 * t439 + t340 * t440 + t354 * t399 - t355 * t400 - g(3)) * MDP(20) + (-t399 * t502 - t439 * t460) * MDP(21) + (t400 * t502 + t440 * t460) * MDP(23) + (t337 * t440 + t338 * t439 - t352 * t399 - t353 * t400 - g(3)) * MDP(24) + ((-t411 + t585) * MDP(17) + (t502 * t598 - t410) * MDP(18) - t375 * MDP(20) - t376 * MDP(21) + t377 * MDP(23) - t341 * MDP(24)) * t540 + (t560 * MDP(17) + t559 * MDP(18) + (MDP(17) * t467 + MDP(20) * t409 + MDP(21) * t415 - MDP(23) * t565 + MDP(24) * t362) * qJD(3)) * t537 + t654 * (-t440 * t376 + t377 * t439 - t399 * t565 + t400 * t415); MDP(7) * t594 + MDP(8) * t522 + qJDD(3) * MDP(9) + (qJD(3) * t442 + t537 * t556 - t526 + t561) * MDP(10) + (t660 + (qJD(3) * t486 - t622) * t537 + (t556 + t662) * t540) * MDP(11) + (-t469 * t637 + t640) * MDP(12) + ((t410 + t639) * t539 + (-t411 + t638) * t536) * MDP(13) + ((-t469 * t537 + t502 * t623) * qJD(1) - t560) * MDP(14) + ((t467 * t537 - t502 * t625) * qJD(1) - t559) * MDP(15) + t502 * MDP(16) * t609 + (-t392 * t609 - pkin(3) * t411 - t442 * t467 + t453 * t502 + (-t441 * t502 + t553) * t536 + t653 * t539) * MDP(17) + (-pkin(3) * t410 + t393 * t609 - t442 * t469 - t617 * t502 - t536 * t653 + t553 * t539) * MDP(18) + (-t339 * t462 - t340 * t461 + t354 * t615 + t355 * t616 + t364 * t415 - t565 * t620 + t545) * MDP(19) + (t340 * t420 - t339 * t419 - t375 * t515 - g(3) * t574 + t655 * t409 + (t403 - t364) * t355 + t620 * t354 + t569 * (t515 * t537 - t627)) * MDP(20) + (t341 * t461 + t352 * t609 - t362 * t616 + t376 * t412 + t415 * t621 - t419 * t460 + t502 * t619 + t518 * t549) * MDP(21) + (-t337 * t461 + t338 * t462 - t352 * t615 + t353 * t616 + t359 * t415 + t565 * t619 + t545) * MDP(22) + (-t341 * t462 - t353 * t609 + t362 * t615 - t377 * t412 + t420 * t460 - t502 * t618 + t516 * t549 - t565 * t621) * MDP(23) + (t337 * t420 + t341 * t412 + t338 * t419 - g(3) * (t540 * t567 + t574) + t621 * t362 + t618 * t353 + t619 * t352 + t569 * (-t627 - (-t515 - t567) * t537)) * MDP(24) + (-MDP(5) * t537 * t540 + MDP(6) * t610) * qJD(1) ^ 2; t469 * t467 * MDP(12) + (-t467 ^ 2 + t469 ^ 2) * MDP(13) + (t410 - t639) * MDP(14) + (-t411 - t638) * MDP(15) + t460 * MDP(16) + (-t429 * t600 - g(1) * t436 + g(2) * t434 - t393 * t502 - t428 * t469 + t407 + (t573 + t644) * t536) * MDP(17) + (g(1) * t437 - g(2) * t435 + g(3) * t624 - t392 * t502 + t428 * t467 - t558) * MDP(18) + (t355 * t565 - t642 + (-t376 * t531 - t377 * t533) * pkin(4) + (-t354 + t357) * t415) * MDP(19) + (-g(1) * t491 + t354 * t356 - t355 * t357 + (g(2) * t631 + t339 * t533 + t340 * t531 - t409 * t469 + t536 * t651) * pkin(4)) * MDP(20) + (-t356 * t502 - t374 * t415 + (pkin(5) - t511) * t460 + t652) * MDP(21) + (t353 * t565 - t376 * t505 + t377 * t511 - t642 + (t352 - t597) * t415) * MDP(22) + (-t518 * t644 - g(1) * t425 - g(2) * t423 - t362 * t415 + t374 * t565 + t460 * t505 + (-0.2e1 * qJD(6) + t357) * t502 + t590) * MDP(23) + (t337 * t505 + t338 * t511 - t362 * t374 - t352 * t356 - g(1) * (-pkin(4) * t591 - pkin(5) * t424 + qJ(6) * t425 + t491) - g(2) * (-pkin(4) * t434 - pkin(5) * t422 + qJ(6) * t423) - g(3) * (-t504 + (-pkin(5) * t516 + qJ(6) * t518) * t537) + t597 * t353) * MDP(24); (t354 * t565 + t355 * t415 + t375 + t526) * MDP(20) + (-t502 * t565 + t376) * MDP(21) + (-t377 - t663) * MDP(23) + (-t352 * t565 + t353 * t415 + t341 + t526) * MDP(24) + (-MDP(20) - MDP(24)) * t656 + t654 * (-t415 ^ 2 - t657); (t415 * t565 - t460) * MDP(21) + (t377 - t663) * MDP(22) + (-t502 ^ 2 - t657) * MDP(23) + (t353 * t502 - t649 - t652) * MDP(24);];
tau  = t1;
