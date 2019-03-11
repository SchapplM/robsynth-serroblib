% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:21
% EndTime: 2019-03-09 15:30:30
% DurationCPUTime: 7.14s
% Computational Cost: add. (4476->550), mult. (9642->628), div. (0->0), fcn. (6650->10), ass. (0->230)
t521 = sin(qJ(3));
t525 = cos(qJ(2));
t648 = cos(qJ(3));
t593 = t648 * t525;
t570 = qJD(1) * t593;
t522 = sin(qJ(2));
t608 = qJD(1) * t522;
t434 = t521 * t608 - t570;
t511 = qJD(2) + qJD(3);
t520 = sin(qJ(6));
t524 = cos(qJ(6));
t409 = t434 * t520 + t524 * t511;
t449 = t521 * t525 + t522 * t648;
t436 = t449 * qJD(1);
t662 = qJD(6) + t436;
t666 = t409 * t662;
t411 = t434 * t524 - t511 * t520;
t665 = t411 * t662;
t506 = t525 * pkin(2);
t497 = t506 + pkin(1);
t600 = qJDD(1) * t525;
t602 = qJD(1) * qJD(2);
t664 = t522 * t602 - t600;
t527 = -pkin(8) - pkin(7);
t458 = t527 * t525;
t452 = qJD(1) * t458;
t438 = t521 * t452;
t457 = t527 * t522;
t450 = qJD(1) * t457;
t643 = qJD(2) * pkin(2);
t442 = t450 + t643;
t614 = -t648 * t442 - t438;
t663 = qJD(4) + t614;
t516 = qJ(2) + qJ(3);
t504 = sin(t516);
t492 = g(3) * t504;
t505 = cos(t516);
t526 = cos(qJ(1));
t623 = t505 * t526;
t523 = sin(qJ(1));
t624 = t505 * t523;
t661 = g(1) * t623 + g(2) * t624 + t492;
t392 = t436 * pkin(3) + t434 * qJ(4);
t649 = -pkin(4) - pkin(9);
t358 = -pkin(5) * t434 + t436 * t649 - t392;
t513 = -pkin(3) + t649;
t660 = t662 * (qJD(6) * t513 + t358);
t496 = -pkin(2) * t648 - pkin(3);
t487 = -pkin(4) + t496;
t475 = -pkin(9) + t487;
t599 = pkin(2) * t608;
t659 = t662 * (qJD(6) * t475 + t358 - t599);
t441 = t648 * t452;
t401 = t521 * t450 - t441;
t607 = qJD(3) * t521;
t566 = pkin(2) * t607 - t401;
t402 = t450 * t648 + t438;
t591 = t648 * qJD(3);
t612 = pkin(2) * t591 + qJD(4) - t402;
t611 = t505 * pkin(3) + t504 * qJ(4);
t510 = qJDD(2) + qJDD(3);
t499 = t510 * qJ(4);
t501 = t511 * qJD(4);
t657 = -t499 - t501;
t588 = g(1) * t523 - g(2) * t526;
t509 = g(1) * t526;
t610 = g(2) * t523 + t509;
t656 = t648 * qJD(2) + t591;
t502 = t510 * pkin(3);
t654 = qJDD(4) - t502;
t404 = t511 * t449;
t587 = qJDD(1) * t648;
t601 = qJDD(1) * t522;
t562 = t521 * t601 - t525 * t587;
t381 = qJD(1) * t404 + t562;
t589 = t525 * t602;
t407 = qJDD(2) * pkin(2) - t527 * (-t589 - t601);
t408 = t527 * t664;
t573 = -t521 * t407 - t648 * t408 - t442 * t591 - t452 * t607;
t556 = t573 + t657;
t341 = -t381 * qJ(5) - t434 * qJD(5) + t556;
t339 = pkin(5) * t510 - t341;
t413 = t521 * t457 - t458 * t648;
t594 = qJD(2) * t527;
t451 = t522 * t594;
t453 = t525 * t594;
t369 = qJD(3) * t413 + t521 * t451 - t453 * t648;
t619 = t521 * t522;
t563 = t511 * t619;
t403 = -t656 * t525 + t563;
t353 = t403 * qJ(5) - t449 * qJD(5) + t369;
t448 = -t593 + t619;
t582 = -t448 * pkin(3) + t449 * qJ(4) + t497;
t360 = pkin(5) * t449 + t448 * t649 + t582;
t395 = t521 * t442 - t441;
t632 = t434 * qJ(5);
t379 = t395 + t632;
t503 = t511 * qJ(4);
t371 = -t379 - t503;
t367 = pkin(5) * t511 - t371;
t572 = t511 * t570 + t521 * t600 + t522 * t587;
t380 = qJD(1) * t563 - t572;
t375 = -qJDD(6) + t380;
t412 = -t457 * t648 - t521 * t458;
t389 = -t449 * qJ(5) + t412;
t574 = -t648 * t407 + t521 * t408 + t442 * t607 - t452 * t591;
t351 = t574 + t654;
t641 = qJ(5) * t380;
t536 = -qJD(5) * t436 + t351 + t641;
t338 = t510 * t649 + t536;
t456 = t497 * qJD(1);
t386 = t434 * pkin(3) - t436 * qJ(4) - t456;
t561 = qJD(5) - t386;
t354 = pkin(5) * t436 + t434 * t649 + t561;
t580 = qJD(6) * t354 + t338;
t652 = t339 * t448 + t367 * t404 + t389 * t375 - (qJD(6) * t360 + t353) * t662 - t449 * t580;
t433 = t436 ^ 2;
t650 = pkin(3) + pkin(4);
t647 = pkin(2) * t522;
t646 = pkin(4) * t436;
t645 = pkin(4) * t510;
t493 = g(3) * t505;
t490 = t505 * pkin(4);
t642 = qJ(4) * t381;
t517 = qJDD(1) * pkin(1);
t605 = qJD(6) * t524;
t596 = t524 * t381 - t520 * t510 - t511 * t605;
t606 = qJD(6) * t520;
t356 = -t434 * t606 + t596;
t640 = t356 * t520;
t639 = t360 * t375;
t638 = t367 * t448;
t489 = pkin(2) * t521 + qJ(4);
t637 = t381 * t489;
t636 = t614 * t511;
t635 = t395 * t511;
t634 = t409 * t434;
t633 = t411 * t434;
t631 = t434 * t436;
t630 = t434 * t511;
t629 = t436 * qJ(5);
t627 = t489 * t510;
t626 = t504 * t523;
t625 = t504 * t526;
t622 = t520 * t375;
t621 = t520 * t523;
t620 = t520 * t526;
t618 = t523 * t524;
t617 = t524 * t375;
t616 = t524 * t526;
t615 = qJ(5) + t527;
t613 = -t629 + t612;
t514 = t522 ^ 2;
t609 = -t525 ^ 2 + t514;
t558 = -t629 + t663;
t598 = t522 * t643;
t429 = pkin(2) * t664 - t517;
t595 = t506 + t611;
t364 = t511 * t513 + t558;
t344 = t354 * t524 - t364 * t520;
t586 = t339 * t524 + t344 * t434;
t585 = t520 * t662;
t584 = t524 * t662;
t583 = t662 * t367;
t348 = t381 * pkin(3) + t380 * qJ(4) - t436 * qJD(4) + t429;
t549 = qJDD(5) - t348;
t335 = -pkin(5) * t380 + t381 * t649 + t549;
t581 = qJD(6) * t364 - t335;
t576 = MDP(31) * t662;
t575 = t511 * t522;
t571 = g(2) * (pkin(3) * t623 + qJ(4) * t625 + t526 * t497);
t569 = g(1) * t626 - g(2) * t625;
t568 = -g(1) * t624 + g(2) * t623;
t567 = -t632 + t566;
t565 = -pkin(3) * t504 - t647;
t388 = t599 + t392;
t345 = t354 * t520 + t364 * t524;
t559 = -t345 * t434 + t520 * t661;
t557 = -t497 - t611;
t555 = -0.2e1 * pkin(1) * t602 - pkin(7) * qJDD(2);
t554 = t404 * t524 - t448 * t606;
t359 = t404 * pkin(3) + t403 * qJ(4) - t449 * qJD(4) + t598;
t368 = t648 * t451 + t521 * t453 + t457 * t591 + t458 * t607;
t551 = -t573 - t661;
t550 = -g(1) * t625 - g(2) * t626 + t493 + t574;
t547 = -t505 * t610 - t492;
t546 = -t550 - t654;
t545 = t513 * t375 + t379 * t662 - t583;
t529 = qJD(2) ^ 2;
t544 = -pkin(7) * t529 + 0.2e1 * t517 + t588;
t530 = qJD(1) ^ 2;
t543 = pkin(1) * t530 - pkin(7) * qJDD(1) + t610;
t542 = t368 * t511 + t413 * t510 + t569;
t541 = -t369 * t511 - t412 * t510 - t568;
t336 = -pkin(4) * t381 + t549;
t540 = -t386 * t434 + t551;
t539 = -t456 * t434 - t551;
t538 = t456 * t436 - t550;
t537 = t610 * t650 * t504;
t535 = -t386 * t436 + t546;
t534 = t475 * t375 - t567 * t662 - t583;
t366 = -pkin(4) * t434 + t561;
t533 = t366 * t434 - t341 - t661;
t483 = t524 * t510;
t357 = qJD(6) * t411 + t381 * t520 + t483;
t363 = -qJD(1) * t521 * t575 + t572 + t630;
t432 = t434 ^ 2;
t532 = t662 * t434 * MDP(30) + MDP(11) * t631 + ((-t356 + t666) * t524 + (t357 + t665) * t520) * MDP(27) + (-t411 * t584 - t640) * MDP(26) + (t585 * t662 + t617 - t634) * MDP(29) + (-t584 * t662 + t622 + t633) * MDP(28) + t363 * MDP(13) - t562 * MDP(14) + (-t432 + t433) * MDP(12) + t510 * MDP(15);
t531 = t641 + (-qJD(5) - t366) * t436 - t546;
t519 = qJ(4) + pkin(5);
t484 = pkin(5) + t489;
t461 = qJ(4) * t623;
t459 = qJ(4) * t624;
t425 = t504 * t616 - t621;
t424 = -t504 * t620 - t618;
t423 = -t504 * t618 - t620;
t422 = t504 * t621 - t616;
t391 = t503 + t395;
t390 = t448 * qJ(5) + t413;
t387 = -pkin(3) * t511 + t663;
t383 = -pkin(4) * t448 + t582;
t377 = -t392 - t646;
t370 = -t388 - t646;
t365 = -t511 * t650 + t558;
t352 = -t404 * qJ(5) - t448 * qJD(5) - t368;
t349 = -pkin(4) * t404 - t359;
t342 = -pkin(5) * t403 + t404 * t649 - t359;
t340 = t536 - t645;
t334 = t524 * t335;
t1 = [(qJDD(1) * t514 + 0.2e1 * t522 * t589) * MDP(4) + t588 * MDP(2) + (t351 * t449 - t368 * t434 + t369 * t436 - t380 * t412 - t381 * t413 - t387 * t403 - t391 * t404 + t448 * t556 - t610) * MDP(19) + (-t556 * t413 + t391 * t368 - t348 * t582 + t386 * t359 + t351 * t412 + t387 * t369 + t527 * t509 - t571 + (-g(1) * t557 + g(2) * t527) * t523) * MDP(21) + (t348 * t448 + t359 * t434 - t381 * t582 + t386 * t404 + t541) * MDP(18) + (-t348 * t449 - t359 * t436 - t380 * t582 + t386 * t403 + t542) * MDP(20) + (-t403 * t511 + t449 * t510) * MDP(13) + (-t404 * t511 - t448 * t510) * MDP(14) + (t380 * t448 - t381 * t449 + t403 * t434 - t404 * t436) * MDP(12) + (-t380 * t449 - t403 * t436) * MDP(11) + (t336 * t448 + t349 * t434 + t353 * t511 + t366 * t404 + t381 * t383 + t389 * t510 + t568) * MDP(23) + (t336 * t449 + t349 * t436 - t352 * t511 - t366 * t403 - t380 * t383 + t390 * t510 + t569) * MDP(22) + (-t340 * t449 - t341 * t448 - t352 * t434 - t353 * t436 + t365 * t403 - t371 * t404 + t380 * t389 + t381 * t390 + t610) * MDP(24) + (t522 * t555 + t525 * t544) * MDP(9) + (-t522 * t544 + t525 * t555) * MDP(10) + (t356 * t448 * t524 + t411 * t554) * MDP(26) + (t380 * t497 + t403 * t456 + t429 * t449 + t436 * t598 - t542) * MDP(17) + (-t381 * t497 - t404 * t456 + t429 * t448 + t434 * t598 + t541) * MDP(16) + (-t375 * t449 - t403 * t662) * MDP(30) + (-g(1) * t423 - g(2) * t425 + t334 * t449 - t344 * t403 - t352 * t409 + t390 * t357 + (t342 * t662 - t639 + (-t364 * t449 - t389 * t662 + t638) * qJD(6)) * t524 + t652 * t520) * MDP(31) + (-g(1) * t422 - g(2) * t424 + t345 * t403 - t352 * t411 + t390 * t356 + (-(-qJD(6) * t389 + t342) * t662 + t639 + t581 * t449 - qJD(6) * t638) * t520 + t652 * t524) * MDP(32) + (t356 * t449 - t403 * t411 - t448 * t617 + t554 * t662) * MDP(28) + (t448 * t622 - t357 * t449 + t403 * t409 + (-t404 * t520 - t448 * t605) * t662) * MDP(29) + (qJDD(2) * t522 + t525 * t529) * MDP(6) + (qJDD(2) * t525 - t522 * t529) * MDP(7) + qJDD(1) * MDP(1) + (t340 * t389 + t365 * t353 - t341 * t390 + t371 * t352 + t336 * t383 + t366 * t349 - t571 + (g(1) * t615 - g(2) * t490) * t526 + (-g(1) * (t557 - t490) + g(2) * t615) * t523) * MDP(25) + 0.2e1 * (t522 * t600 - t602 * t609) * MDP(5) + t610 * MDP(3) + ((-t409 * t524 - t411 * t520) * t404 + (-t640 - t357 * t524 + (t409 * t520 - t411 * t524) * qJD(6)) * t448) * MDP(27); (t484 * t356 + t613 * t411 + (-t339 + t659) * t520 + t534 * t524 + t559) * MDP(32) + (t484 * t357 + t613 * t409 + (t547 - t659) * t524 + t534 * t520 + t586) * MDP(31) + (t402 * t511 + (-t436 * t608 - t510 * t521 - t511 * t591) * pkin(2) + t539) * MDP(17) + (t401 * t511 + (-t434 * t608 + t510 * t648 - t511 * t607) * pkin(2) + t538) * MDP(16) + (-t370 * t436 + t511 * t613 + t533 + t627) * MDP(22) + (-t388 * t434 - t496 * t510 - t511 * t566 + t535) * MDP(18) + (t380 * t487 + t637 + (t371 - t567) * t436 + (-t365 + t613) * t434) * MDP(24) + (-t380 * t496 - t637 + (t391 + t566) * t436 + (t387 - t612) * t434) * MDP(19) + (t388 * t436 + t511 * t612 + t540 + t627 - t657) * MDP(20) + t532 + (-g(3) * t525 + t522 * t543) * MDP(9) + (g(3) * t522 + t525 * t543) * MDP(10) + ((-pkin(4) + t487) * t510 + t567 * t511 + t531 - t370 * t434) * MDP(23) + (-t556 * t489 + t351 * t496 - t386 * t388 - g(1) * (t526 * t565 + t461) - g(2) * (t523 * t565 + t459) - g(3) * t595 + t612 * t391 + t566 * t387) * MDP(21) + MDP(7) * t600 + MDP(6) * t601 + qJDD(2) * MDP(8) + (t340 * t487 - t341 * t489 - t366 * t370 - g(1) * (-t526 * t647 + t461) - g(2) * (-t523 * t647 + t459) - g(3) * (t490 + t595) + t537 - t613 * t371 + t567 * t365) * MDP(25) + (-MDP(4) * t522 * t525 + MDP(5) * t609) * t530; (t538 + t635) * MDP(16) + (-t392 * t434 + t502 + t535 + t635) * MDP(18) + (pkin(3) * t380 - t642 + (t391 - t395) * t436 + (t387 - t663) * t434) * MDP(19) + (t392 * t436 + 0.2e1 * t499 + 0.2e1 * t501 + t540 + t636) * MDP(20) + (-t377 * t436 + t511 * t558 + t499 + t533) * MDP(22) + (t642 - t380 * t650 + (t371 + t379) * t436 + (-t365 + t558) * t434) * MDP(24) + (t519 * t357 + t558 * t409 + t545 * t520 + (t547 - t660) * t524 + t586) * MDP(31) + (t519 * t356 + t558 * t411 + (-t339 + t660) * t520 + t545 * t524 + t559) * MDP(32) + t532 + (t539 - t636) * MDP(17) + (t531 + (-pkin(4) - t650) * t510 - t377 * t434 - t379 * t511) * MDP(23) + (-t556 * qJ(4) - t351 * pkin(3) - t386 * t392 - t387 * t395 - g(1) * (-pkin(3) * t625 + t461) - g(2) * (-pkin(3) * t626 + t459) - g(3) * t611 + t663 * t391) * MDP(21) + (-t340 * t650 - t341 * qJ(4) - t365 * t379 - t366 * t377 - g(1) * t461 - g(2) * t459 - g(3) * (t490 + t611) + t537 - t558 * t371) * MDP(25); t363 * MDP(19) + (-t391 * t511 - t535) * MDP(21) + (t380 - t630) * MDP(24) + (t371 * t511 + t531 - t645) * MDP(25) + (-t409 * t511 + t622) * MDP(31) + (-t411 * t511 + t617) * MDP(32) + (MDP(32) * t585 - t524 * t576) * t662 + (-MDP(18) + MDP(23)) * (t510 - t631) + (MDP(22) + MDP(20)) * (-t511 ^ 2 - t433); (t572 - t630) * MDP(22) + (t436 * t511 + t562) * MDP(23) + (-t433 - t432) * MDP(24) + (t365 * t436 + t371 * t434 + t336 + t588) * MDP(25) + (-t617 - t634) * MDP(31) + (t622 - t633) * MDP(32) + (-MDP(32) * t584 - t520 * t576) * t662 + (t656 * t522 * MDP(23) + (MDP(23) * t511 * t525 - MDP(22) * t575) * t521) * qJD(1); t411 * t409 * MDP(26) + (-t409 ^ 2 + t411 ^ 2) * MDP(27) + (t596 + t666) * MDP(28) + (-t483 + t665) * MDP(29) - t375 * MDP(30) + (-g(1) * t424 + g(2) * t422 + t345 * t662 - t367 * t411 + t334) * MDP(31) + (g(1) * t425 - g(2) * t423 + t344 * t662 + t367 * t409) * MDP(32) + ((-t338 - t493) * MDP(32) + (-MDP(29) * t434 - MDP(31) * t364 - MDP(32) * t354) * qJD(6)) * t524 + (-qJD(6) * t434 * MDP(28) + (qJD(6) * t511 - t381) * MDP(29) + (-t580 - t493) * MDP(31) + t581 * MDP(32)) * t520;];
tau  = t1;
