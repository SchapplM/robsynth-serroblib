% Calculate vector of inverse dynamics joint torques for
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:38
% EndTime: 2019-03-08 19:05:49
% DurationCPUTime: 7.25s
% Computational Cost: add. (4166->529), mult. (10328->765), div. (0->0), fcn. (9480->18), ass. (0->235)
t649 = cos(pkin(6));
t487 = qJD(1) * t649 + qJD(2);
t506 = sin(pkin(6));
t504 = sin(pkin(13));
t512 = sin(qJ(3));
t516 = cos(qJ(3));
t507 = cos(pkin(13));
t508 = cos(pkin(7));
t631 = t507 * t508;
t545 = t504 * t516 + t512 * t631;
t533 = t545 * t506;
t505 = sin(pkin(7));
t633 = t505 * t512;
t408 = qJD(1) * t533 + t487 * t633;
t511 = sin(qJ(4));
t515 = cos(qJ(4));
t563 = pkin(4) * t511 - pkin(10) * t515;
t472 = t563 * qJD(4);
t669 = t408 - t472;
t510 = sin(qJ(5));
t514 = cos(qJ(5));
t600 = t514 * qJD(4);
t614 = qJD(3) * t511;
t461 = t510 * t614 - t600;
t513 = cos(qJ(6));
t609 = qJD(4) * t510;
t463 = t514 * t614 + t609;
t509 = sin(qJ(6));
t638 = t463 * t509;
t414 = t461 * t513 + t638;
t612 = qJD(3) * t515;
t488 = -qJD(5) + t612;
t483 = -qJD(6) + t488;
t668 = t414 * t483;
t552 = t461 * t509 - t463 * t513;
t667 = t483 * t552;
t485 = qJDD(1) * t649 + qJDD(2);
t615 = qJD(1) * t506;
t590 = t507 * t615;
t568 = t508 * t590;
t598 = qJDD(1) * t506;
t579 = t507 * t598;
t613 = qJD(3) * t512;
t589 = t505 * t613;
t591 = t504 * t615;
t634 = t504 * t512;
t596 = t506 * t634;
t611 = qJD(3) * t516;
t520 = -t516 * (t485 * t505 + t508 * t579) + qJDD(1) * t596 + t487 * t589 + t568 * t613 + t591 * t611;
t648 = cos(pkin(12));
t560 = t649 * t648;
t647 = sin(pkin(12));
t443 = t504 * t560 + t507 * t647;
t526 = t504 * t647 - t507 * t560;
t576 = t506 * t648;
t658 = t505 * t576 + t508 * t526;
t391 = t443 * t512 + t516 * t658;
t559 = t649 * t647;
t444 = -t504 * t559 + t507 * t648;
t527 = t504 * t648 + t507 * t559;
t575 = t506 * t647;
t657 = -t505 * t575 + t508 * t527;
t393 = t444 * t512 + t516 * t657;
t577 = t505 * t649;
t567 = t516 * t577;
t595 = t516 * t631;
t419 = -t506 * t595 - t567 + t596;
t539 = g(1) * t393 + g(2) * t391 + g(3) * t419;
t666 = qJD(3) * t408 - t520 + t539;
t604 = qJD(5) * t511;
t665 = -qJD(3) * t604 + qJDD(4);
t403 = qJD(3) * pkin(9) + t408;
t436 = t487 * t508 - t505 * t590;
t385 = t403 * t515 + t436 * t511;
t383 = qJD(4) * pkin(10) + t385;
t407 = -t512 * t591 + t516 * (t487 * t505 + t568);
t476 = -pkin(4) * t515 - pkin(10) * t511 - pkin(3);
t398 = qJD(3) * t476 - t407;
t354 = t383 * t514 + t398 * t510;
t352 = -pkin(11) * t461 + t354;
t602 = qJD(6) * t509;
t350 = t352 * t602;
t392 = t443 * t516 - t512 * t658;
t421 = t505 * t526 - t508 * t576;
t367 = t392 * t515 + t421 * t511;
t394 = t444 * t516 - t512 * t657;
t422 = t505 * t527 + t508 * t575;
t369 = t394 * t515 + t422 * t511;
t659 = -t403 * t511 + t436 * t515;
t382 = -qJD(4) * pkin(4) - t659;
t372 = pkin(5) * t461 + t382;
t420 = t512 * t577 + t533;
t534 = -t505 * t506 * t507 + t508 * t649;
t397 = t420 * t515 + t511 * t534;
t503 = qJ(5) + qJ(6);
t498 = sin(t503);
t499 = cos(t503);
t664 = t372 * t414 - g(1) * (-t369 * t499 - t393 * t498) - g(2) * (-t367 * t499 - t391 * t498) - g(3) * (-t397 * t499 - t419 * t498) + t350;
t599 = qJD(3) * qJD(4);
t582 = t515 * t599;
t597 = qJDD(3) * t511;
t409 = qJD(5) * t600 + (t582 + t597) * t514 + t665 * t510;
t496 = t515 * qJDD(3);
t458 = t511 * t599 + qJDD(5) - t496;
t546 = t595 - t634;
t376 = qJDD(3) * pkin(9) + (t485 * t512 + t487 * t611) * t505 + (qJD(1) * qJD(3) * t546 + qJDD(1) * t545) * t506;
t434 = t485 * t508 - t505 * t579;
t642 = t434 * t511;
t343 = qJDD(4) * pkin(10) + qJD(4) * t659 + t376 * t515 + t642;
t359 = qJD(3) * t472 + qJDD(3) * t476 + t520;
t358 = t514 * t359;
t524 = -qJD(5) * t354 - t510 * t343 + t358;
t335 = pkin(5) * t458 - pkin(11) * t409 + t524;
t410 = t510 * (qJD(4) * (qJD(5) + t612) + t597) - t665 * t514;
t603 = qJD(5) * t514;
t594 = t343 * t514 + t359 * t510 + t398 * t603;
t605 = qJD(5) * t510;
t541 = t383 * t605 - t594;
t336 = -pkin(11) * t410 - t541;
t574 = t335 * t513 - t509 * t336;
t663 = t372 * t552 - g(1) * (-t369 * t498 + t393 * t499) - g(2) * (-t367 * t498 + t391 * t499) - g(3) * (-t397 * t498 + t419 * t499) + t574;
t455 = qJDD(6) + t458;
t662 = t455 * MDP(24) + (-t414 ^ 2 + t552 ^ 2) * MDP(21) - t414 * MDP(20) * t552;
t465 = t509 * t514 + t510 * t513;
t440 = t465 * t511;
t608 = qJD(4) * t511;
t628 = t510 * t515;
t651 = pkin(9) * t510;
t661 = -t407 * t628 + t514 * t669 - t608 * t651;
t626 = t514 * t515;
t660 = -t407 * t626 + t476 * t603 - t510 * t669;
t656 = -t510 * t604 + t515 * t600;
t655 = qJD(5) + qJD(6);
t653 = qJD(5) * (pkin(9) * t488 + t383) + t539;
t573 = t409 * t509 + t410 * t513;
t356 = -qJD(6) * t552 + t573;
t652 = pkin(10) + pkin(11);
t650 = qJD(3) * pkin(3);
t353 = -t383 * t510 + t398 * t514;
t351 = -pkin(11) * t463 + t353;
t349 = -pkin(5) * t488 + t351;
t645 = t349 * t513;
t644 = t352 * t513;
t643 = t409 * t510;
t640 = t461 * t488;
t639 = t463 * t488;
t637 = t463 * t514;
t636 = t498 * t515;
t635 = t499 * t515;
t632 = t505 * t516;
t630 = t509 * t335;
t629 = t510 * t511;
t627 = t511 * t514;
t489 = pkin(9) * t626;
t550 = pkin(5) * t511 - pkin(11) * t626;
t625 = -t550 * qJD(4) - (-t489 + (pkin(11) * t511 - t476) * t510) * qJD(5) + t661;
t607 = qJD(4) * t515;
t585 = t510 * t607;
t536 = t511 * t603 + t585;
t624 = -t536 * pkin(11) + (-t511 * t600 - t515 * t605) * pkin(9) + t660;
t469 = t563 * qJD(3);
t623 = t469 * t510 + t514 * t659;
t464 = t509 * t510 - t513 * t514;
t542 = t515 * t464;
t622 = qJD(3) * t542 - t464 * t655;
t621 = (-t612 + t655) * t465;
t618 = t476 * t510 + t489;
t501 = t511 ^ 2;
t617 = -t515 ^ 2 + t501;
t610 = qJD(4) * t461;
t606 = qJD(5) * t488;
t601 = qJD(6) * t513;
t593 = t409 * t513 - t410 * t509 - t461 * t601;
t592 = qJD(5) * t652;
t588 = t505 * t611;
t587 = t510 * t612;
t586 = t488 * t600;
t581 = t516 * t599;
t572 = pkin(5) * t536 + pkin(9) * t607 - t407 * t511;
t571 = qJD(6) * t349 + t336;
t564 = -t385 + (-t587 + t605) * pkin(5);
t452 = t514 * t469;
t481 = t652 * t514;
t562 = qJD(3) * t550 + qJD(6) * t481 - t510 * t659 + t514 * t592 + t452;
t480 = t652 * t510;
t561 = pkin(11) * t587 - qJD(6) * t480 - t510 * t592 - t623;
t338 = t349 * t509 + t644;
t370 = -t397 * t510 + t419 * t514;
t371 = t397 * t514 + t419 * t510;
t558 = t370 * t513 - t371 * t509;
t557 = t370 * t509 + t371 * t513;
t460 = t514 * t476;
t418 = -pkin(11) * t627 + t460 + (-pkin(5) - t651) * t515;
t429 = -pkin(11) * t629 + t618;
t555 = t418 * t509 + t429 * t513;
t449 = t508 * t511 + t515 * t633;
t427 = -t449 * t510 - t514 * t632;
t547 = -t449 * t514 + t510 * t632;
t554 = t427 * t513 + t509 * t547;
t553 = t427 * t509 - t513 * t547;
t518 = qJD(3) ^ 2;
t551 = qJDD(3) * t516 - t512 * t518;
t549 = t376 * t511 + t403 * t607 - t434 * t515 + t436 * t608;
t448 = -t508 * t515 + t511 * t633;
t544 = t458 * t510 - t488 * t603;
t543 = t458 * t514 + t488 * t605;
t355 = -t463 * t602 + t593;
t396 = t420 * t511 - t515 * t534;
t540 = g(1) * (-t394 * t511 + t422 * t515) + g(2) * (-t392 * t511 + t421 * t515) - g(3) * t396;
t538 = g(1) * t394 + g(2) * t392 + g(3) * t420;
t344 = -qJDD(4) * pkin(4) + t549;
t531 = -pkin(10) * t458 - t382 * t488;
t402 = -t407 - t650;
t530 = -pkin(9) * qJDD(4) + (t402 + t407 - t650) * qJD(4);
t523 = pkin(10) * t606 - t344 - t540;
t517 = qJD(4) ^ 2;
t519 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t517 + t666;
t492 = -pkin(5) * t514 - pkin(4);
t473 = (pkin(5) * t510 + pkin(9)) * t511;
t441 = t464 * t511;
t426 = qJD(4) * t449 + t511 * t588;
t425 = -qJD(4) * t448 + t515 * t588;
t412 = t420 * qJD(3);
t411 = (t506 * t546 + t567) * qJD(3);
t389 = -t602 * t629 + (t627 * t655 + t585) * t513 + t656 * t509;
t388 = -qJD(4) * t542 - t440 * t655;
t381 = qJD(5) * t547 - t425 * t510 + t514 * t589;
t380 = qJD(5) * t427 + t425 * t514 + t510 * t589;
t365 = -qJD(4) * t396 + t411 * t515;
t364 = qJD(4) * t397 + t411 * t511;
t341 = pkin(5) * t410 + t344;
t340 = qJD(5) * t370 + t365 * t514 + t412 * t510;
t339 = -qJD(5) * t371 - t365 * t510 + t412 * t514;
t337 = -t352 * t509 + t645;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t485 * t649 - g(3) + (t504 ^ 2 + t507 ^ 2) * t506 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t412 - qJDD(3) * t419) * MDP(4) + (-qJD(3) * t411 - qJDD(3) * t420) * MDP(5) + (-t419 * t496 - qJD(4) * t364 - qJDD(4) * t396 + (-t412 * t515 + t419 * t608) * qJD(3)) * MDP(11) + (t419 * t597 - qJD(4) * t365 - qJDD(4) * t397 + (t412 * t511 + t419 * t607) * qJD(3)) * MDP(12) + (-t339 * t488 + t364 * t461 + t370 * t458 + t396 * t410) * MDP(18) + (t340 * t488 + t364 * t463 - t371 * t458 + t396 * t409) * MDP(19) + (-(-qJD(6) * t557 + t339 * t513 - t340 * t509) * t483 + t558 * t455 + t364 * t414 + t396 * t356) * MDP(25) + ((qJD(6) * t558 + t339 * t509 + t340 * t513) * t483 - t557 * t455 - t364 * t552 + t396 * t355) * MDP(26); (-g(3) * t649 + (-g(1) * t647 + g(2) * t648) * t506 + t485) * MDP(2) + (-qJD(4) * t426 - qJDD(4) * t448) * MDP(11) + (-qJD(4) * t425 - qJDD(4) * t449) * MDP(12) + (-t381 * t488 + t410 * t448 + t426 * t461 + t427 * t458) * MDP(18) + (t380 * t488 + t409 * t448 + t426 * t463 + t458 * t547) * MDP(19) + (-(-qJD(6) * t553 - t380 * t509 + t381 * t513) * t483 + t554 * t455 + t426 * t414 + t448 * t356) * MDP(25) + ((qJD(6) * t554 + t380 * t513 + t381 * t509) * t483 - t553 * t455 - t426 * t552 + t448 * t355) * MDP(26) + (t551 * MDP(4) + (-qJDD(3) * t512 - t516 * t518) * MDP(5) + (-t511 * t581 + t515 * t551) * MDP(11) + (-t511 * t551 - t515 * t581) * MDP(12)) * t505; qJDD(3) * MDP(3) + t666 * MDP(4) + (-t485 * t633 - t545 * t598 + (-t487 * t632 - t546 * t615 + t407) * qJD(3) + t538) * MDP(5) + (qJDD(3) * t501 + 0.2e1 * t511 * t582) * MDP(6) + 0.2e1 * (t496 * t511 - t599 * t617) * MDP(7) + (qJDD(4) * t511 + t515 * t517) * MDP(8) + (qJDD(4) * t515 - t511 * t517) * MDP(9) + (t511 * t530 + t515 * t519) * MDP(11) + (-t511 * t519 + t515 * t530) * MDP(12) + (t409 * t627 + t463 * t656) * MDP(13) + ((-t461 * t514 - t463 * t510) * t607 + (-t643 - t410 * t514 + (t461 * t510 - t637) * qJD(5)) * t511) * MDP(14) + ((-t409 - t586) * t515 + (qJD(4) * t463 + t543) * t511) * MDP(15) + ((t488 * t609 + t410) * t515 + (-t544 - t610) * t511) * MDP(16) + (-t458 * t515 - t488 * t608) * MDP(17) + (t460 * t458 + t661 * t488 + (t476 * t606 - t538) * t510 + (pkin(9) * t610 - t358 + (-pkin(9) * t458 + qJD(4) * t382 + qJD(5) * t398 + t343) * t510 + t653 * t514) * t515 + (pkin(9) * t410 + qJD(4) * t353 + t344 * t510 + t382 * t603 - t407 * t461) * t511) * MDP(18) + (-t618 * t458 + t660 * t488 - t538 * t514 + ((pkin(9) * t463 + t382 * t514) * qJD(4) - t653 * t510 + t594) * t515 + (-t382 * t605 - t354 * qJD(4) + t344 * t514 - t407 * t463 + (t409 - t586) * pkin(9)) * t511) * MDP(19) + (-t355 * t441 - t388 * t552) * MDP(20) + (-t355 * t440 + t356 * t441 - t388 * t414 + t389 * t552) * MDP(21) + (-t355 * t515 - t388 * t483 - t441 * t455 - t552 * t608) * MDP(22) + (t356 * t515 + t389 * t483 - t414 * t608 - t440 * t455) * MDP(23) + (-t455 * t515 - t483 * t608) * MDP(24) + ((t418 * t513 - t429 * t509) * t455 - t574 * t515 + t337 * t608 + t473 * t356 + t341 * t440 + t372 * t389 - g(1) * (-t393 * t635 + t394 * t498) - g(2) * (-t391 * t635 + t392 * t498) - g(3) * (-t419 * t635 + t420 * t498) + (t509 * t624 + t513 * t625) * t483 + t572 * t414 + (t338 * t515 + t483 * t555) * qJD(6)) * MDP(25) + (-t555 * t455 + (t571 * t513 - t350 + t630) * t515 - t338 * t608 + t473 * t355 - t341 * t441 + t372 * t388 - g(1) * (t393 * t636 + t394 * t499) - g(2) * (t391 * t636 + t392 * t499) - g(3) * (t419 * t636 + t420 * t499) + ((qJD(6) * t418 + t624) * t513 + (-qJD(6) * t429 - t625) * t509) * t483 - t572 * t552) * MDP(26); MDP(8) * t597 + MDP(9) * t496 + qJDD(4) * MDP(10) + (qJD(4) * t385 - t540 - t549) * MDP(11) + (g(1) * t369 + g(2) * t367 + g(3) * t397 - t642 + (-qJD(3) * t402 - t376) * t515) * MDP(12) + (-t488 * t637 + t643) * MDP(13) + ((t409 + t640) * t514 + (-t410 + t639) * t510) * MDP(14) + ((-t463 * t511 + t488 * t626) * qJD(3) + t544) * MDP(15) + ((t461 * t511 - t488 * t628) * qJD(3) + t543) * MDP(16) + (-pkin(4) * t410 - t385 * t461 + t452 * t488 + (-t488 * t659 + t531) * t510 + t523 * t514) * MDP(18) + (-pkin(4) * t409 - t385 * t463 - t488 * t623 - t510 * t523 + t514 * t531) * MDP(19) + (t355 * t465 - t552 * t622) * MDP(20) + (-t355 * t464 - t356 * t465 - t414 * t622 + t552 * t621) * MDP(21) + (t455 * t465 - t483 * t622) * MDP(22) + (-t455 * t464 + t483 * t621) * MDP(23) + ((-t480 * t513 - t481 * t509) * t455 + t492 * t356 + t341 * t464 + (t509 * t561 + t513 * t562) * t483 + t564 * t414 + t621 * t372 - t540 * t499) * MDP(25) + (-(-t480 * t509 + t481 * t513) * t455 + t492 * t355 + t341 * t465 + (-t509 * t562 + t513 * t561) * t483 - t564 * t552 + t622 * t372 + t540 * t498) * MDP(26) + (-MDP(11) * t402 + t488 * MDP(17) - t353 * MDP(18) + MDP(19) * t354 + MDP(22) * t552 + t414 * MDP(23) + t483 * MDP(24) - t337 * MDP(25) + t338 * MDP(26)) * t614 + (-MDP(6) * t511 * t515 + MDP(7) * t617) * t518; t463 * t461 * MDP(13) + (-t461 ^ 2 + t463 ^ 2) * MDP(14) + (t409 - t640) * MDP(15) + (-t410 - t639) * MDP(16) + t458 * MDP(17) + (-t354 * t488 - t382 * t463 - g(1) * (-t369 * t510 + t393 * t514) - g(2) * (-t367 * t510 + t391 * t514) - g(3) * t370 + t524) * MDP(18) + (-t353 * t488 + t382 * t461 - g(1) * (-t369 * t514 - t393 * t510) - g(2) * (-t367 * t514 - t391 * t510) + g(3) * t371 + t541) * MDP(19) + (t355 - t668) * MDP(22) + (-t356 + t667) * MDP(23) + ((-t351 * t509 - t644) * t483 - t338 * qJD(6) + (-t414 * t463 + t455 * t513 + t483 * t602) * pkin(5) + t663) * MDP(25) + ((t352 * t483 - t335) * t509 + (-t351 * t483 - t571) * t513 + (-t455 * t509 + t463 * t552 + t483 * t601) * pkin(5) + t664) * MDP(26) + t662; (t593 - t668) * MDP(22) + (-t573 + t667) * MDP(23) + (-t338 * t483 + t663) * MDP(25) + (-t513 * t336 - t337 * t483 - t630 + t664) * MDP(26) + (-MDP(22) * t638 + MDP(23) * t552 - MDP(25) * t338 - MDP(26) * t645) * qJD(6) + t662;];
tau  = t1;
