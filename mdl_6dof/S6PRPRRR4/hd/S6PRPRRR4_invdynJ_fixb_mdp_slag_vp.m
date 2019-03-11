% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:39:15
% EndTime: 2019-03-08 20:39:25
% DurationCPUTime: 9.13s
% Computational Cost: add. (4822->540), mult. (11339->728), div. (0->0), fcn. (9580->18), ass. (0->249)
t535 = cos(pkin(12));
t672 = cos(qJ(4));
t602 = t672 * t535;
t515 = qJD(2) * t602;
t539 = sin(qJ(4));
t532 = sin(pkin(12));
t623 = qJD(2) * t532;
t599 = t539 * t623;
t486 = -t515 + t599;
t676 = qJD(5) + qJD(6);
t692 = t486 + t676;
t563 = -t539 * t532 + t602;
t670 = pkin(8) + qJ(3);
t500 = t670 * t532;
t501 = t670 * t535;
t564 = -t500 * t672 - t539 * t501;
t408 = qJD(3) * t563 + qJD(4) * t564;
t534 = sin(pkin(6));
t543 = cos(qJ(2));
t642 = t534 * t543;
t550 = t563 * t642;
t465 = qJD(1) * t550;
t691 = -t408 + t465;
t489 = t563 * qJD(4);
t493 = t532 * t672 + t539 * t535;
t490 = t493 * qJD(4);
t540 = sin(qJ(2));
t625 = qJD(1) * t534;
t601 = t540 * t625;
t690 = pkin(4) * t490 - pkin(9) * t489 - t601;
t499 = qJD(2) * qJ(3) + t601;
t536 = cos(pkin(6));
t624 = qJD(1) * t536;
t511 = t535 * t624;
t455 = t511 + (-pkin(8) * qJD(2) - t499) * t532;
t471 = t535 * t499 + t532 * t624;
t622 = qJD(2) * t535;
t456 = pkin(8) * t622 + t471;
t389 = t539 * t455 + t672 * t456;
t689 = t389 * qJD(4);
t488 = t493 * qJD(2);
t538 = sin(qJ(5));
t542 = cos(qJ(5));
t616 = t542 * qJD(4);
t466 = t488 * t538 - t616;
t541 = cos(qJ(6));
t468 = qJD(4) * t538 + t488 * t542;
t537 = sin(qJ(6));
t657 = t468 * t537;
t397 = t541 * t466 + t657;
t479 = qJD(5) + t486;
t477 = qJD(6) + t479;
t688 = t397 * t477;
t574 = t466 * t537 - t541 * t468;
t687 = t477 * t574;
t639 = t537 * t542;
t496 = t538 * t541 + t639;
t628 = t692 * t496;
t595 = qJDD(2) * t672;
t611 = qJDD(2) * t539;
t578 = t532 * t611 - t535 * t595;
t444 = qJD(2) * t490 + t578;
t620 = qJD(5) * t538;
t656 = t486 * t538;
t686 = t620 + t656;
t647 = t532 * MDP(6);
t685 = t535 * MDP(5) - t647;
t385 = qJD(4) * pkin(9) + t389;
t517 = -pkin(3) * t535 - pkin(2);
t600 = t543 * t625;
t579 = qJD(3) - t600;
t478 = qJD(2) * t517 + t579;
t401 = pkin(4) * t486 - pkin(9) * t488 + t478;
t366 = t385 * t542 + t401 * t538;
t361 = -pkin(10) * t466 + t366;
t618 = qJD(6) * t537;
t359 = t361 * t618;
t678 = t672 * t455 - t539 * t456;
t384 = -qJD(4) * pkin(4) - t678;
t374 = t466 * pkin(5) + t384;
t669 = cos(pkin(11));
t593 = t669 * t540;
t533 = sin(pkin(11));
t644 = t533 * t543;
t483 = t536 * t593 + t644;
t530 = pkin(12) + qJ(4);
t521 = sin(t530);
t522 = cos(t530);
t594 = t534 * t669;
t446 = t483 * t522 - t521 * t594;
t592 = t669 * t543;
t645 = t533 * t540;
t485 = -t536 * t645 + t592;
t646 = t533 * t534;
t448 = t485 * t522 + t521 * t646;
t643 = t534 * t540;
t476 = t521 * t536 + t522 * t643;
t482 = -t536 * t592 + t645;
t484 = t536 * t644 + t593;
t531 = qJ(5) + qJ(6);
t526 = sin(t531);
t527 = cos(t531);
t684 = t374 * t397 - g(1) * (-t448 * t527 - t484 * t526) - g(2) * (-t446 * t527 - t482 * t526) - g(3) * (-t476 * t527 + t526 * t642) + t359;
t605 = qJD(4) * t515 + t532 * t595 + t535 * t611;
t443 = -qJD(4) * t599 + t605;
t382 = qJD(5) * t616 + t538 * qJDD(4) + t542 * t443 - t488 * t620;
t437 = qJDD(5) + t444;
t613 = qJDD(2) * qJ(3);
t615 = qJDD(1) * t534;
t469 = t540 * t615 + t613 + (qJD(3) + t600) * qJD(2);
t614 = qJDD(1) * t536;
t509 = t535 * t614;
t413 = t509 + (-pkin(8) * qJDD(2) - t469) * t532;
t433 = t535 * t469 + t532 * t614;
t612 = qJDD(2) * t535;
t414 = pkin(8) * t612 + t433;
t567 = t539 * t413 + t414 * t672;
t357 = qJDD(4) * pkin(9) + qJD(4) * t678 + t567;
t621 = qJD(2) * t540;
t596 = qJD(1) * t621;
t610 = t534 * t596 + qJDD(3);
t568 = -t543 * t615 + t610;
t458 = qJDD(2) * t517 + t568;
t373 = pkin(4) * t444 - pkin(9) * t443 + t458;
t371 = t542 * t373;
t546 = -qJD(5) * t366 - t538 * t357 + t371;
t346 = pkin(5) * t437 - pkin(10) * t382 + t546;
t383 = qJD(5) * t468 - t542 * qJDD(4) + t443 * t538;
t619 = qJD(5) * t542;
t558 = -t542 * t357 - t538 * t373 + t385 * t620 - t401 * t619;
t347 = -pkin(10) * t383 - t558;
t590 = t541 * t346 - t537 * t347;
t683 = t374 * t574 - g(1) * (-t448 * t526 + t484 * t527) - g(2) * (-t446 * t526 + t482 * t527) - g(3) * (-t476 * t526 - t527 * t642) + t590;
t431 = qJDD(6) + t437;
t682 = t431 * MDP(27) + (-t397 ^ 2 + t574 ^ 2) * MDP(24) - t397 * MDP(23) * t574;
t627 = t532 ^ 2 + t535 ^ 2;
t681 = MDP(7) * t627;
t426 = t496 * t493;
t463 = -t539 * t500 + t501 * t672;
t551 = t493 * t642;
t631 = -qJD(1) * t551 + qJD(3) * t493 + qJD(4) * t463;
t680 = t465 * t538 + t542 * t690;
t440 = -pkin(4) * t563 - pkin(9) * t493 + t517;
t679 = -t440 * t619 + t463 * t620 - t690 * t538 + t542 * t691;
t495 = t537 * t538 - t541 * t542;
t629 = t692 * t495;
t675 = -t431 * t496 + t477 * t629;
t589 = t382 * t537 + t541 * t383;
t353 = -qJD(6) * t574 + t589;
t673 = pkin(9) + pkin(10);
t671 = g(3) * t534;
t668 = qJDD(2) * pkin(2);
t365 = -t385 * t538 + t542 * t401;
t360 = -pkin(10) * t468 + t365;
t355 = pkin(5) * t479 + t360;
t667 = t355 * t541;
t666 = t361 * t541;
t665 = t382 * t538;
t664 = t397 * t488;
t663 = t574 * t488;
t661 = t466 * t479;
t660 = t466 * t488;
t659 = t468 * t479;
t658 = t468 * t488;
t655 = t489 * t538;
t654 = t489 * t542;
t653 = t493 * t538;
t652 = t493 * t542;
t651 = t522 * t526;
t650 = t522 * t527;
t649 = t522 * t538;
t648 = t522 * t543;
t640 = t537 * t346;
t638 = t538 * t437;
t637 = t538 * t543;
t421 = t542 * t437;
t452 = t542 * t463;
t636 = qJDD(1) - g(3);
t635 = -pkin(10) * t654 + pkin(5) * t490 - t408 * t538 + (-t452 + (pkin(10) * t493 - t440) * t538) * qJD(5) + t680;
t562 = t493 * t619 + t655;
t634 = pkin(10) * t562 + t679;
t633 = pkin(5) * t562 + t631;
t438 = pkin(4) * t488 + pkin(9) * t486;
t632 = t538 * t438 + t542 * t678;
t630 = t538 * t440 + t452;
t617 = qJD(6) * t541;
t609 = g(3) * t643;
t608 = t534 * t637;
t607 = t542 * t642;
t606 = t541 * t382 - t537 * t383 - t466 * t617;
t604 = qJD(5) * t673;
t598 = t534 * t621;
t597 = t493 * t620;
t591 = t534 * t636;
t587 = t479 * t542;
t586 = qJD(6) * t355 + t347;
t585 = pkin(5) * t686 - t389;
t584 = g(1) * t484 + g(2) * t482;
t583 = g(1) * t485 + g(2) * t483;
t582 = -t495 * t431 - t477 * t628;
t424 = t542 * t438;
t505 = t673 * t542;
t581 = pkin(5) * t488 + qJD(6) * t505 - t678 * t538 + t424 + (pkin(10) * t486 + t604) * t542;
t504 = t673 * t538;
t580 = pkin(10) * t656 + qJD(6) * t504 + t538 * t604 + t632;
t349 = t355 * t537 + t666;
t429 = t542 * t440;
t372 = -pkin(5) * t563 - pkin(10) * t652 - t463 * t538 + t429;
t375 = -pkin(10) * t653 + t630;
t577 = t372 * t537 + t375 * t541;
t480 = -t532 * t643 + t535 * t536;
t481 = t532 * t536 + t535 * t643;
t418 = t539 * t480 + t481 * t672;
t402 = -t418 * t538 - t607;
t570 = -t418 * t542 + t608;
t576 = t402 * t541 + t537 * t570;
t575 = t402 * t537 - t541 * t570;
t470 = -t499 * t532 + t511;
t573 = t470 * t532 - t471 * t535;
t572 = -t479 * t686 + t421;
t569 = MDP(3) + t685;
t565 = t480 * t672 - t539 * t481;
t561 = -t597 + t654;
t560 = -pkin(9) * t437 + t384 * t479;
t559 = -t413 * t672 + t539 * t414 + t689;
t352 = -t468 * t618 + t606;
t556 = g(1) * (-t485 * t521 + t522 * t646) + g(2) * (-t483 * t521 - t522 * t594) + g(3) * (-t521 * t643 + t522 * t536);
t555 = t584 + t668;
t554 = g(3) * t642 - t584;
t358 = -qJDD(4) * pkin(4) + t559;
t549 = t554 * t522;
t432 = -t469 * t532 + t509;
t548 = -t432 * t532 + t433 * t535 - t583;
t545 = pkin(9) * qJD(5) * t479 + t358 + t556;
t544 = qJD(2) ^ 2;
t520 = -pkin(5) * t542 - pkin(4);
t494 = -qJD(2) * pkin(2) + t579;
t472 = t568 - t668;
t427 = t495 * t493;
t410 = pkin(5) * t653 - t564;
t387 = qJD(2) * t551 + qJD(4) * t418;
t386 = qJD(2) * t550 + qJD(4) * t565;
t369 = t489 * t639 - t537 * t597 - t618 * t653 + (t652 * t676 + t655) * t541;
t368 = -t426 * t676 - t495 * t489;
t364 = qJD(5) * t570 - t386 * t538 + t542 * t598;
t363 = qJD(5) * t402 + t386 * t542 + t538 * t598;
t350 = t383 * pkin(5) + t358;
t348 = -t361 * t537 + t667;
t1 = [t636 * MDP(1) + (t432 * t480 + t433 * t481 - g(3)) * MDP(8) + (-qJD(4) * t387 + qJDD(4) * t565) * MDP(14) + (-qJD(4) * t386 - qJDD(4) * t418) * MDP(15) + (t364 * t479 - t383 * t565 + t387 * t466 + t402 * t437) * MDP(21) + (-t363 * t479 - t382 * t565 + t387 * t468 + t437 * t570) * MDP(22) + ((-qJD(6) * t575 - t363 * t537 + t364 * t541) * t477 + t576 * t431 + t387 * t397 - t565 * t353) * MDP(28) + (-(qJD(6) * t576 + t363 * t541 + t364 * t537) * t477 - t575 * t431 - t387 * t574 - t565 * t352) * MDP(29) + (-t480 * t532 + t481 * t535) * MDP(7) * qJDD(2) + ((-qJDD(2) * MDP(4) - t569 * t544 + (MDP(14) * t486 + MDP(15) * t488 + MDP(8) * t494) * qJD(2)) * t540 + ((-t470 * t623 + t471 * t622 - t472) * MDP(8) - t444 * MDP(14) - t443 * MDP(15) + (-MDP(4) + t681) * t544 + t569 * qJDD(2)) * t543) * t534; (-t540 * t591 + t583) * MDP(4) + (-t352 * t427 - t368 * t574) * MDP(23) + (-t352 * t426 + t353 * t427 - t368 * t397 + t369 * t574) * MDP(24) + (-qJD(4) * t490 + qJDD(4) * t563) * MDP(12) + (t443 * t563 - t444 * t493 - t486 * t489 - t488 * t490) * MDP(10) + (-t382 * t563 + t421 * t493 + t468 * t490 + t479 * t561) * MDP(18) + ((t372 * t541 - t375 * t537) * t431 - t590 * t563 + t348 * t490 + t410 * t353 + t350 * t426 + t374 * t369 - g(1) * (-t484 * t650 + t485 * t526) - g(2) * (-t482 * t650 + t483 * t526) - (t526 * t540 + t527 * t648) * t671 + (t537 * t634 + t541 * t635) * t477 + t633 * t397 + (t349 * t563 - t477 * t577) * qJD(6)) * MDP(28) + (t383 * t563 - t466 * t490 - t479 * t562 - t493 * t638) * MDP(19) + (-t431 * t563 + t477 * t490) * MDP(27) + (t353 * t563 - t369 * t477 - t397 * t490 - t426 * t431) * MDP(26) + (-t437 * t563 + t479 * t490) * MDP(20) + (-t577 * t431 + (t541 * t586 - t359 + t640) * t563 - t349 * t490 + t410 * t352 - t350 * t427 + t374 * t368 - g(1) * (t484 * t651 + t485 * t527) - g(2) * (t482 * t651 + t483 * t527) - (-t526 * t648 + t527 * t540) * t671 + ((-qJD(6) * t372 + t634) * t541 + (qJD(6) * t375 - t635) * t537) * t477 - t633 * t574) * MDP(29) + (-t352 * t563 + t368 * t477 - t427 * t431 - t490 * t574) * MDP(25) + (-qJD(4) * t631 + qJDD(4) * t564 + t444 * t517 - t458 * t563 + t478 * t490 - t486 * t601 - t549) * MDP(14) + (t365 * t490 - t371 * t563 - t564 * t383 + t429 * t437 + t680 * t479 + t631 * t466 + (-t549 + (t384 * t493 + t385 * t563 - t463 * t479) * qJD(5)) * t542 + ((-qJD(5) * t440 - t408) * t479 - t463 * t437 - (-qJD(5) * t401 - t357) * t563 + t358 * t493 + t384 * t489 - t609 - t583) * t538) * MDP(21) + (-t630 * t437 - t558 * t563 - t366 * t490 - t564 * t382 + t358 * t652 - g(1) * (t484 * t649 + t485 * t542) - g(2) * (t482 * t649 + t483 * t542) - (-t522 * t637 + t540 * t542) * t671 + t679 * t479 + t631 * t468 + t561 * t384) * MDP(22) + (-t609 + t548 + (qJD(2) * t579 + t613) * t627) * MDP(7) + (qJD(4) * t691 - qJDD(4) * t463 + t443 * t517 + t458 * t493 + t478 * t489 - t488 * t601 + t554 * t521) * MDP(15) + t685 * (t534 * (-g(3) * t543 + t596) - t472 + t555) + qJDD(2) * MDP(2) + (t382 * t652 + t468 * t561) * MDP(16) + (-t573 * qJD(3) + (-t472 + t584) * pkin(2) + t548 * qJ(3) + (-g(3) * (pkin(2) * t543 + qJ(3) * t540) + (-t494 * t540 + t543 * t573) * qJD(1)) * t534) * MDP(8) + ((-t466 * t542 - t468 * t538) * t489 + (-t665 - t383 * t542 + (t466 * t538 - t468 * t542) * qJD(5)) * t493) * MDP(17) + (t636 * t642 + t584) * MDP(3) + (qJD(4) * t489 + qJDD(4) * t493) * MDP(11) + (t443 * t493 + t488 * t489) * MDP(9); -MDP(5) * t612 + qJDD(2) * t647 - t544 * t681 + (qJD(2) * t573 - t543 * t591 - t555 + t610) * MDP(8) + (0.2e1 * qJD(4) * t488 + t578) * MDP(14) + ((-t486 - t599) * qJD(4) + t605) * MDP(15) + (t572 - t660) * MDP(21) + (-t479 ^ 2 * t542 - t638 - t658) * MDP(22) + (t582 - t664) * MDP(28) + (t663 + t675) * MDP(29); -t486 ^ 2 * MDP(10) + ((t486 - t599) * qJD(4) + t605) * MDP(11) - t444 * MDP(12) + qJDD(4) * MDP(13) + (-t556 - t559 + t689) * MDP(14) + (g(1) * t448 + g(2) * t446 + g(3) * t476 + t478 * t486 - t567) * MDP(15) + (t468 * t587 + t665) * MDP(16) + ((t382 - t661) * t542 + (-t383 - t659) * t538) * MDP(17) + (t479 * t587 + t638 - t658) * MDP(18) + (t572 + t660) * MDP(19) + (-pkin(4) * t383 - t389 * t466 - t424 * t479 + (t479 * t678 + t560) * t538 - t545 * t542) * MDP(21) + (-pkin(4) * t382 - t389 * t468 + t479 * t632 + t538 * t545 + t542 * t560) * MDP(22) + (t352 * t496 + t574 * t629) * MDP(23) + (-t352 * t495 - t353 * t496 + t397 * t629 + t574 * t628) * MDP(24) + (t663 - t675) * MDP(25) + (t582 + t664) * MDP(26) + ((-t504 * t541 - t505 * t537) * t431 + t520 * t353 + t350 * t495 + (t537 * t580 - t541 * t581) * t477 + t585 * t397 + t628 * t374 - t556 * t527) * MDP(28) + (-(-t504 * t537 + t505 * t541) * t431 + t520 * t352 + t350 * t496 + (t537 * t581 + t541 * t580) * t477 - t585 * t574 - t629 * t374 + t556 * t526) * MDP(29) + (MDP(10) * t488 + MDP(12) * qJD(4) - MDP(14) * t478 - MDP(20) * t479 - MDP(21) * t365 + MDP(22) * t366 - MDP(27) * t477 - MDP(28) * t348 + MDP(29) * t349 + MDP(9) * t486) * t488; t468 * t466 * MDP(16) + (-t466 ^ 2 + t468 ^ 2) * MDP(17) + (t382 + t661) * MDP(18) + (-t383 + t659) * MDP(19) + t437 * MDP(20) + (t366 * t479 - t384 * t468 - g(1) * (-t448 * t538 + t484 * t542) - g(2) * (-t446 * t538 + t482 * t542) - g(3) * (-t476 * t538 - t607) + t546) * MDP(21) + (t365 * t479 + t384 * t466 - g(1) * (-t448 * t542 - t484 * t538) - g(2) * (-t446 * t542 - t482 * t538) - g(3) * (-t476 * t542 + t608) + t558) * MDP(22) + (t352 + t688) * MDP(25) + (-t353 - t687) * MDP(26) + (-(-t360 * t537 - t666) * t477 - t349 * qJD(6) + (-t397 * t468 + t431 * t541 - t477 * t618) * pkin(5) + t683) * MDP(28) + ((-t361 * t477 - t346) * t537 + (t360 * t477 - t586) * t541 + (-t431 * t537 + t468 * t574 - t477 * t617) * pkin(5) + t684) * MDP(29) + t682; (t606 + t688) * MDP(25) + (-t589 - t687) * MDP(26) + (t349 * t477 + t683) * MDP(28) + (-t541 * t347 + t348 * t477 - t640 + t684) * MDP(29) + (-MDP(25) * t657 + MDP(26) * t574 - MDP(28) * t349 - MDP(29) * t667) * qJD(6) + t682;];
tau  = t1;
