% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:33
% EndTime: 2019-03-08 22:00:45
% DurationCPUTime: 9.66s
% Computational Cost: add. (5175->557), mult. (12193->767), div. (0->0), fcn. (9948->18), ass. (0->251)
t538 = sin(pkin(12));
t547 = sin(qJ(3));
t619 = qJD(2) * t547;
t541 = cos(pkin(12));
t551 = cos(qJ(3));
t641 = t541 * t551;
t493 = qJD(2) * t641 - t538 * t619;
t676 = qJD(5) + qJD(6);
t687 = t493 - t676;
t685 = qJ(4) + pkin(8);
t595 = qJD(3) * t685;
t487 = qJD(4) * t551 - t547 * t595;
t488 = -qJD(4) * t547 - t551 * t595;
t423 = t487 * t541 + t488 * t538;
t502 = t538 * t547 - t641;
t540 = sin(pkin(6));
t552 = cos(qJ(2));
t643 = t540 * t552;
t603 = qJD(1) * t643;
t468 = t502 * t603;
t628 = t423 + t468;
t503 = t538 * t551 + t541 * t547;
t494 = t503 * qJD(3);
t497 = t502 * qJD(3);
t548 = sin(qJ(2));
t621 = qJD(1) * t548;
t604 = t540 * t621;
t669 = qJD(3) * pkin(3);
t608 = t547 * t669;
t686 = pkin(4) * t494 + pkin(9) * t497 - t604 + t608;
t495 = t503 * qJD(2);
t546 = sin(qJ(5));
t550 = cos(qJ(5));
t614 = t550 * qJD(3);
t471 = t495 * t546 - t614;
t549 = cos(qJ(6));
t473 = qJD(3) * t546 + t495 * t550;
t545 = sin(qJ(6));
t657 = t473 * t545;
t404 = t549 * t471 + t657;
t484 = qJD(5) - t493;
t482 = qJD(6) + t484;
t684 = t404 * t482;
t575 = t471 * t545 - t549 * t473;
t683 = t482 * t575;
t637 = t545 * t550;
t506 = t546 * t549 + t637;
t625 = t687 * t506;
t618 = qJD(5) * t546;
t656 = t493 * t546;
t682 = t618 - t656;
t588 = qJD(2) * t685 + t604;
t543 = cos(pkin(6));
t622 = qJD(1) * t543;
t463 = -t547 * t588 + t551 * t622;
t456 = t463 + t669;
t464 = t547 * t622 + t551 * t588;
t642 = t541 * t464;
t393 = t538 * t456 + t642;
t388 = qJD(3) * pkin(9) + t393;
t526 = pkin(3) * t551 + pkin(2);
t483 = -qJD(2) * t526 + qJD(4) - t603;
t409 = -pkin(4) * t493 - pkin(9) * t495 + t483;
t368 = t388 * t550 + t409 * t546;
t362 = -pkin(10) * t471 + t368;
t616 = qJD(6) * t545;
t360 = t362 * t616;
t452 = t538 * t464;
t392 = t456 * t541 - t452;
t387 = -qJD(3) * pkin(4) - t392;
t375 = pkin(5) * t471 + t387;
t539 = sin(pkin(11));
t542 = cos(pkin(11));
t640 = t543 * t548;
t490 = t539 * t552 + t542 * t640;
t534 = qJ(3) + pkin(12);
t527 = sin(t534);
t528 = cos(t534);
t646 = t540 * t542;
t449 = t490 * t528 - t527 * t646;
t492 = -t539 * t640 + t542 * t552;
t647 = t539 * t540;
t451 = t492 * t528 + t527 * t647;
t645 = t540 * t548;
t480 = t527 * t543 + t528 * t645;
t639 = t543 * t552;
t489 = t539 * t548 - t542 * t639;
t491 = t539 * t639 + t542 * t548;
t537 = qJ(5) + qJ(6);
t532 = sin(t537);
t533 = cos(t537);
t681 = t375 * t404 - g(1) * (-t451 * t533 - t491 * t532) - g(2) * (-t449 * t533 - t489 * t532) - g(3) * (-t480 * t533 + t532 * t643) + t360;
t612 = qJD(2) * qJD(3);
t597 = t551 * t612;
t598 = t547 * t612;
t447 = qJDD(2) * t503 - t538 * t598 + t541 * t597;
t389 = qJD(5) * t614 + t546 * qJDD(3) + t550 * t447 - t495 * t618;
t610 = qJDD(2) * t551;
t611 = qJDD(2) * t547;
t579 = -t538 * t611 + t541 * t610;
t442 = qJD(2) * t494 + qJDD(5) - t579;
t609 = t543 * qJDD(1);
t517 = t551 * t609;
t613 = qJD(1) * qJD(2);
t475 = qJDD(2) * pkin(8) + (qJDD(1) * t548 + t552 * t613) * t540;
t558 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t622 + t475;
t574 = t588 * qJD(3);
t382 = qJDD(3) * pkin(3) - t547 * t558 - t551 * t574 + t517;
t383 = (-t574 + t609) * t547 + t558 * t551;
t359 = t538 * t382 + t541 * t383;
t357 = qJDD(3) * pkin(9) + t359;
t599 = t548 * t613;
t515 = t540 * t599;
t561 = pkin(3) * t598 - qJDD(2) * t526 + qJDD(4) + t515;
t596 = qJDD(1) * t643;
t445 = t561 - t596;
t446 = -qJD(3) * t495 + t579;
t373 = -pkin(4) * t446 - pkin(9) * t447 + t445;
t372 = t550 * t373;
t557 = -qJD(5) * t368 - t546 * t357 + t372;
t345 = pkin(5) * t442 - pkin(10) * t389 + t557;
t390 = qJD(5) * t473 - t550 * qJDD(3) + t447 * t546;
t617 = qJD(5) * t550;
t568 = -t550 * t357 - t546 * t373 + t388 * t618 - t409 * t617;
t346 = -pkin(10) * t390 - t568;
t592 = t549 * t345 - t545 * t346;
t680 = t375 * t575 - g(1) * (-t451 * t532 + t491 * t533) - g(2) * (-t449 * t532 + t489 * t533) - g(3) * (-t480 * t532 - t533 * t643) + t592;
t438 = qJDD(6) + t442;
t679 = t438 * MDP(25) + (-t404 ^ 2 + t575 ^ 2) * MDP(22) - t404 * MDP(21) * t575;
t433 = t506 * t503;
t678 = -t468 * t546 + t550 * t686;
t629 = t487 * t538 - t541 * t488 - t503 * t603;
t444 = pkin(4) * t502 - pkin(9) * t503 - t526;
t510 = t685 * t547;
t511 = t685 * t551;
t470 = -t510 * t538 + t511 * t541;
t677 = -t444 * t617 + t470 * t618 - t686 * t546 - t550 * t628;
t505 = t545 * t546 - t549 * t550;
t626 = t687 * t505;
t675 = -t438 * t506 - t482 * t626;
t553 = qJD(3) ^ 2;
t585 = g(1) * t491 + g(2) * t489;
t674 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t553 + t540 * (-g(3) * t552 + t599) - t515 + t585 + t596;
t591 = t389 * t545 + t549 * t390;
t352 = -qJD(6) * t575 + t591;
t498 = t543 * t551 - t547 * t645;
t673 = g(3) * t498;
t672 = g(3) * t540;
t522 = pkin(3) * t538 + pkin(9);
t671 = pkin(10) + t522;
t670 = qJD(2) * pkin(2);
t367 = -t388 * t546 + t550 * t409;
t361 = -pkin(10) * t473 + t367;
t354 = pkin(5) * t484 + t361;
t667 = t354 * t549;
t666 = t362 * t549;
t665 = t389 * t546;
t664 = t404 * t495;
t663 = t575 * t495;
t661 = t471 * t484;
t660 = t471 * t495;
t659 = t473 * t484;
t658 = t473 * t495;
t655 = t497 * t546;
t654 = t497 * t550;
t653 = t503 * t546;
t652 = t503 * t550;
t651 = t528 * t532;
t650 = t528 * t533;
t649 = t528 * t546;
t648 = t528 * t552;
t644 = t540 * t551;
t638 = t545 * t345;
t636 = t546 * t442;
t635 = t546 * t552;
t432 = t550 * t442;
t457 = t550 * t470;
t634 = qJDD(1) - g(3);
t633 = pkin(10) * t654 + pkin(5) * t494 - t423 * t546 + (-t457 + (pkin(10) * t503 - t444) * t546) * qJD(5) + t678;
t570 = t503 * t617 - t655;
t632 = pkin(10) * t570 + t677;
t631 = pkin(5) * t570 + t629;
t397 = t463 * t541 - t452;
t424 = pkin(3) * t619 + pkin(4) * t495 - pkin(9) * t493;
t630 = t550 * t397 + t546 * t424;
t627 = t546 * t444 + t457;
t535 = t547 ^ 2;
t624 = -t551 ^ 2 + t535;
t620 = qJD(2) * t540;
t615 = qJD(6) * t549;
t607 = t540 * t635;
t606 = t550 * t643;
t605 = t549 * t389 - t545 * t390 - t471 * t615;
t523 = -pkin(3) * t541 - pkin(4);
t602 = t548 * t620;
t601 = t552 * t620;
t600 = t503 * t618;
t594 = qJD(5) * t671;
t593 = t540 * t634;
t358 = t382 * t541 - t538 * t383;
t395 = t463 * t538 + t642;
t469 = t541 * t510 + t511 * t538;
t589 = t484 * t550;
t587 = qJD(6) * t354 + t346;
t586 = pkin(5) * t682 - t395;
t584 = g(1) * t492 + g(2) * t490;
t583 = g(1) * t539 - g(2) * t542;
t582 = -t505 * t438 + t482 * t625;
t417 = t550 * t424;
t501 = t671 * t550;
t581 = pkin(5) * t495 + qJD(6) * t501 - t397 * t546 + t417 + (-pkin(10) * t493 + t594) * t550;
t500 = t671 * t546;
t580 = -pkin(10) * t656 + qJD(6) * t500 + t546 * t594 + t630;
t348 = t354 * t545 + t666;
t436 = t550 * t444;
t374 = pkin(5) * t502 - pkin(10) * t652 - t470 * t546 + t436;
t376 = -pkin(10) * t653 + t627;
t578 = t374 * t545 + t376 * t549;
t499 = t543 * t547 + t548 * t644;
t430 = t498 * t538 + t499 * t541;
t410 = -t430 * t546 - t606;
t571 = -t430 * t550 + t607;
t577 = t410 * t549 + t545 * t571;
t576 = t410 * t545 - t549 * t571;
t573 = -t484 * t682 + t432;
t356 = -qJDD(3) * pkin(4) - t358;
t569 = -t600 - t654;
t351 = -t473 * t616 + t605;
t566 = t387 * t484 - t522 * t442;
t565 = g(1) * (-t492 * t527 + t528 * t647) + g(2) * (-t490 * t527 - t528 * t646) + g(3) * (-t527 * t645 + t528 * t543);
t564 = -g(3) * t645 - t584;
t508 = -t603 - t670;
t562 = -qJD(2) * t508 - t475 + t584;
t556 = qJD(5) * t484 * t522 + t356 + t565;
t555 = -pkin(8) * qJDD(3) + (t508 + t603 - t670) * qJD(3);
t554 = qJD(2) ^ 2;
t509 = -pkin(5) * t550 + t523;
t462 = -qJD(3) * t499 - t547 * t601;
t461 = qJD(3) * t498 + t551 * t601;
t434 = t505 * t503;
t429 = -t541 * t498 + t499 * t538;
t421 = pkin(5) * t653 + t469;
t396 = t461 * t541 + t462 * t538;
t394 = t461 * t538 - t541 * t462;
t370 = -t497 * t637 - t545 * t600 - t616 * t653 + (t652 * t676 - t655) * t549;
t369 = -t433 * t676 + t505 * t497;
t365 = qJD(5) * t571 - t396 * t546 + t550 * t602;
t364 = qJD(5) * t410 + t396 * t550 + t546 * t602;
t349 = pkin(5) * t390 + t356;
t347 = -t362 * t545 + t667;
t1 = [t634 * MDP(1) + (qJD(3) * t462 + qJDD(3) * t498) * MDP(10) + (-qJD(3) * t461 - qJDD(3) * t499) * MDP(11) + (t394 * t495 + t396 * t493 + t429 * t447 + t430 * t446) * MDP(12) + (-t358 * t429 + t359 * t430 - t392 * t394 + t393 * t396 - g(3)) * MDP(13) + (t365 * t484 + t390 * t429 + t394 * t471 + t410 * t442) * MDP(19) + (-t364 * t484 + t389 * t429 + t394 * t473 + t442 * t571) * MDP(20) + ((-qJD(6) * t576 - t364 * t545 + t365 * t549) * t482 + t577 * t438 + t394 * t404 + t429 * t352) * MDP(26) + (-(qJD(6) * t577 + t364 * t549 + t365 * t545) * t482 - t576 * t438 - t394 * t575 + t429 * t351) * MDP(27) + ((t483 * qJD(2) * MDP(13) - qJDD(2) * MDP(4) + (-t551 * MDP(10) + t547 * MDP(11) - MDP(3)) * t554) * t548 + (qJDD(2) * MDP(3) - t554 * MDP(4) + (-t598 + t610) * MDP(10) + (-t597 - t611) * MDP(11) - t445 * MDP(13)) * t552) * t540; (t359 * t470 - t358 * t469 - t445 * t526 + t483 * t608 - g(1) * (-t491 * t526 + t492 * t685) - g(2) * (-t489 * t526 + t490 * t685) + t628 * t393 - t629 * t392 + (-t483 * t621 - g(3) * (t526 * t552 + t548 * t685)) * t540) * MDP(13) + (-t358 * t503 - t359 * t502 + t392 * t497 - t393 * t494 + t446 * t470 + t447 * t469 + t493 * t628 + t495 * t629 + t564) * MDP(12) + (-(-t471 * t550 - t473 * t546) * t497 + (-t665 - t390 * t550 + (t471 * t546 - t473 * t550) * qJD(5)) * t503) * MDP(15) + (-t627 * t442 + t568 * t502 - t368 * t494 + t469 * t389 + t356 * t652 - g(1) * (t491 * t649 + t492 * t550) - g(2) * (t489 * t649 + t490 * t550) - (-t528 * t635 + t548 * t550) * t672 + t677 * t484 + t629 * t473 + t569 * t387) * MDP(20) + (t367 * t494 + t372 * t502 + t469 * t390 + t436 * t442 + t678 * t484 + t629 * t471 + ((-g(3) * t643 + t585) * t528 + (t387 * t503 - t388 * t502 - t470 * t484) * qJD(5)) * t550 + ((-qJD(5) * t444 - t423) * t484 - t470 * t442 + (-qJD(5) * t409 - t357) * t502 + t356 * t503 - t387 * t497 + t564) * t546) * MDP(19) + (t555 * t547 + t551 * t674) * MDP(10) + (-t547 * t674 + t555 * t551) * MDP(11) + (t438 * t502 + t482 * t494) * MDP(25) + (t442 * t502 + t484 * t494) * MDP(18) + (-t352 * t502 - t370 * t482 - t404 * t494 - t433 * t438) * MDP(24) + (qJDD(3) * t547 + t551 * t553) * MDP(7) + (qJDD(3) * t551 - t547 * t553) * MDP(8) + (t351 * t502 + t369 * t482 - t434 * t438 - t494 * t575) * MDP(23) + (-t351 * t433 + t352 * t434 - t369 * t404 + t370 * t575) * MDP(22) + (-t351 * t434 - t369 * t575) * MDP(21) + (-t578 * t438 - (t587 * t549 - t360 + t638) * t502 - t348 * t494 + t421 * t351 - t349 * t434 + t375 * t369 - g(1) * (t491 * t651 + t492 * t533) - g(2) * (t489 * t651 + t490 * t533) - (-t532 * t648 + t533 * t548) * t672 + ((-qJD(6) * t374 + t632) * t549 + (qJD(6) * t376 - t633) * t545) * t482 - t631 * t575) * MDP(27) + (-t548 * t593 + t584) * MDP(4) + (qJDD(2) * t535 + 0.2e1 * t547 * t597) * MDP(5) + 0.2e1 * (t547 * t610 - t612 * t624) * MDP(6) + (t389 * t502 + t432 * t503 + t473 * t494 + t484 * t569) * MDP(16) + (-t390 * t502 - t471 * t494 - t484 * t570 - t503 * t636) * MDP(17) + (t634 * t643 + t585) * MDP(3) + (t389 * t652 + t473 * t569) * MDP(14) + ((t374 * t549 - t376 * t545) * t438 + t592 * t502 + t347 * t494 + t421 * t352 + t349 * t433 + t375 * t370 - g(1) * (-t491 * t650 + t492 * t532) - g(2) * (-t489 * t650 + t490 * t532) - (t532 * t548 + t533 * t648) * t672 + (t545 * t632 + t549 * t633) * t482 + t631 * t404 + (-t348 * t502 - t482 * t578) * qJD(6)) * MDP(26) + qJDD(2) * MDP(2); MDP(7) * t611 + MDP(8) * t610 + qJDD(3) * MDP(9) + (t547 * t562 - t583 * t644 + t517 - t673) * MDP(10) + (g(3) * t499 + (t540 * t583 - t609) * t547 + t562 * t551) * MDP(11) + ((t392 - t397) * t493 + (t446 * t538 - t447 * t541) * pkin(3)) * MDP(12) + (t392 * t395 - t393 * t397 + (t359 * t538 + t358 * t541 - t483 * t619 - g(1) * (-t492 * t547 + t539 * t644) - g(2) * (-t490 * t547 - t542 * t644) - t673) * pkin(3)) * MDP(13) + (t473 * t589 + t665) * MDP(14) + ((t389 - t661) * t550 + (-t390 - t659) * t546) * MDP(15) + (t484 * t589 + t636 - t658) * MDP(16) + (t573 + t660) * MDP(17) + (t523 * t390 - t395 * t471 - t417 * t484 + (t397 * t484 + t566) * t546 - t556 * t550) * MDP(19) + (t523 * t389 - t395 * t473 + t484 * t630 + t546 * t556 + t550 * t566) * MDP(20) + (t351 * t506 - t575 * t626) * MDP(21) + (-t351 * t505 - t352 * t506 - t404 * t626 - t575 * t625) * MDP(22) + (t663 - t675) * MDP(23) + (t582 + t664) * MDP(24) + ((-t500 * t549 - t501 * t545) * t438 + t509 * t352 + t349 * t505 + (t545 * t580 - t549 * t581) * t482 + t586 * t404 - t625 * t375 - t565 * t533) * MDP(26) + (-(-t500 * t545 + t501 * t549) * t438 + t509 * t351 + t349 * t506 + (t545 * t581 + t549 * t580) * t482 - t586 * t575 + t626 * t375 + t565 * t532) * MDP(27) + (-MDP(5) * t547 * t551 + MDP(6) * t624) * t554 - ((-t393 + t395) * MDP(12) + t484 * MDP(18) + t367 * MDP(19) - t368 * MDP(20) + t482 * MDP(25) + t347 * MDP(26) - t348 * MDP(27)) * t495; (-t493 ^ 2 - t495 ^ 2) * MDP(12) + (t392 * t495 - t393 * t493 - t552 * t593 + t561 - t585) * MDP(13) + (t573 - t660) * MDP(19) + (-t484 ^ 2 * t550 - t636 - t658) * MDP(20) + (t582 - t664) * MDP(26) + (t663 + t675) * MDP(27); t473 * t471 * MDP(14) + (-t471 ^ 2 + t473 ^ 2) * MDP(15) + (t389 + t661) * MDP(16) + (-t390 + t659) * MDP(17) + t442 * MDP(18) + (t368 * t484 - t387 * t473 - g(1) * (-t451 * t546 + t491 * t550) - g(2) * (-t449 * t546 + t489 * t550) - g(3) * (-t480 * t546 - t606) + t557) * MDP(19) + (t367 * t484 + t387 * t471 - g(1) * (-t451 * t550 - t491 * t546) - g(2) * (-t449 * t550 - t489 * t546) - g(3) * (-t480 * t550 + t607) + t568) * MDP(20) + (t351 + t684) * MDP(23) + (-t352 - t683) * MDP(24) + (-(-t361 * t545 - t666) * t482 - t348 * qJD(6) + (-t404 * t473 + t549 * t438 - t482 * t616) * pkin(5) + t680) * MDP(26) + ((-t362 * t482 - t345) * t545 + (t361 * t482 - t587) * t549 + (-t545 * t438 + t473 * t575 - t482 * t615) * pkin(5) + t681) * MDP(27) + t679; (t605 + t684) * MDP(23) + (-t591 - t683) * MDP(24) + (t348 * t482 + t680) * MDP(26) + (-t549 * t346 + t347 * t482 - t638 + t681) * MDP(27) + (-MDP(23) * t657 + MDP(24) * t575 - MDP(26) * t348 - MDP(27) * t667) * qJD(6) + t679;];
tau  = t1;
