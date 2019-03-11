% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:25
% EndTime: 2019-03-09 04:52:36
% DurationCPUTime: 8.23s
% Computational Cost: add. (3963->572), mult. (7336->664), div. (0->0), fcn. (4332->6), ass. (0->233)
t666 = MDP(22) - MDP(27);
t501 = sin(qJ(4));
t504 = cos(qJ(4));
t502 = sin(qJ(3));
t595 = t504 * qJD(3);
t505 = cos(qJ(3));
t600 = qJD(4) * t505;
t528 = t501 * t600 + t502 * t595;
t588 = qJDD(1) * t505;
t376 = qJD(1) * t528 - qJD(4) * t595 - t501 * qJDD(3) - t504 * t588;
t607 = qJD(1) * t505;
t438 = t501 * t607 - t595;
t608 = qJD(1) * t502;
t469 = qJD(4) + t608;
t641 = t438 * t469;
t679 = -t376 - t641;
t678 = t376 - t641;
t577 = t501 * t608;
t606 = qJD(3) * t501;
t440 = t504 * t607 + t606;
t603 = qJD(4) * t440;
t377 = -qJD(3) * t577 - t504 * qJDD(3) + t501 * t588 + t603;
t638 = t440 * t469;
t677 = -t377 - t638;
t591 = qJD(1) * qJD(3);
t570 = t505 * t591;
t589 = qJDD(1) * t502;
t436 = qJDD(4) + t570 + t589;
t676 = t436 * qJ(5) + t469 * qJD(5);
t661 = pkin(4) + pkin(5);
t584 = t661 * t501;
t650 = qJ(5) * t504;
t535 = t584 - t650;
t675 = MDP(23) + MDP(26);
t634 = t501 * qJ(5);
t534 = -t504 * t661 - t634;
t674 = -pkin(3) + t534;
t673 = t377 * qJ(6) + t438 * qJD(6);
t672 = 0.2e1 * t676;
t508 = -pkin(1) - pkin(7);
t465 = qJD(1) * t508 + qJD(2);
t446 = t502 * t465;
t671 = qJD(5) * t501 + t446;
t492 = t505 * pkin(8);
t669 = pkin(3) * t502 + qJ(2) - t492;
t506 = cos(qJ(1));
t623 = t505 * t506;
t477 = g(2) * t623;
t668 = g(3) * t502 + t477;
t414 = t669 * qJD(1);
t421 = qJD(3) * pkin(8) + t446;
t373 = t504 * t414 - t501 * t421;
t593 = qJD(5) - t373;
t667 = MDP(21) + MDP(25);
t425 = t436 * pkin(4);
t665 = t425 - qJDD(5);
t422 = -qJD(3) * pkin(3) - t465 * t505;
t527 = qJ(5) * t440 - t422;
t367 = pkin(4) * t438 - t527;
t658 = pkin(8) * t436;
t664 = t367 * t469 - t658;
t357 = -t438 * t661 + qJD(6) + t527;
t622 = t506 * t504;
t503 = sin(qJ(1));
t627 = t503 * t501;
t415 = t502 * t627 - t622;
t626 = t503 * t504;
t629 = t502 * t506;
t417 = t501 * t629 + t626;
t557 = pkin(3) * t505 + pkin(8) * t502;
t437 = qJD(3) * t557 + qJD(2);
t383 = qJD(1) * t437 + qJDD(1) * t669;
t462 = qJDD(1) * t508 + qJDD(2);
t604 = qJD(3) * t505;
t389 = qJDD(3) * pkin(8) + t462 * t502 + t465 * t604;
t601 = qJD(4) * t504;
t602 = qJD(4) * t501;
t561 = t504 * t383 - t501 * t389 - t414 * t602 - t421 * t601;
t632 = t501 * t505;
t522 = g(1) * t415 - g(2) * t417 + g(3) * t632 + t561;
t519 = t522 + t665;
t649 = qJ(6) * t376;
t663 = (qJD(6) + t357) * t440 + t519 - t649;
t662 = t438 ^ 2;
t435 = t440 ^ 2;
t660 = pkin(4) * t501;
t659 = pkin(5) * t436;
t657 = g(1) * t503;
t497 = g(2) * t506;
t656 = pkin(8) - qJ(6);
t655 = pkin(8) * qJD(4);
t654 = pkin(1) * qJDD(1);
t510 = qJD(1) ^ 2;
t653 = qJ(2) * t510;
t652 = qJ(5) * t377;
t651 = qJ(5) * t438;
t648 = qJ(6) * t505;
t374 = t501 * t414 + t504 * t421;
t361 = qJ(6) * t438 + t374;
t461 = t469 * qJ(5);
t356 = t361 + t461;
t647 = t356 * t469;
t365 = t461 + t374;
t646 = t365 * t469;
t645 = t374 * t469;
t644 = t376 * t501;
t643 = t377 * t504;
t642 = t438 * t440;
t640 = t438 * t501;
t639 = t438 * t504;
t637 = t440 * t501;
t636 = t440 * t504;
t635 = t469 * t504;
t633 = t501 * t502;
t631 = t502 * t503;
t630 = t502 * t504;
t628 = t502 * t508;
t625 = t503 * t505;
t624 = t504 * t505;
t621 = -t469 * t535 + t671;
t443 = t557 * qJD(1);
t616 = t501 * t443 + t465 * t624;
t380 = qJ(5) * t607 + t616;
t596 = qJD(6) * t504;
t619 = qJ(6) * t577 - t602 * t656 - t380 - t596;
t552 = -t650 + t660;
t618 = t469 * t552 - t671;
t455 = t656 * t504;
t432 = t465 * t632;
t564 = -t443 * t504 + t432;
t583 = t661 * t505;
t617 = qJD(4) * t455 - qJD(6) * t501 - (qJ(6) * t630 - t583) * qJD(1) - t564;
t615 = g(3) * t633 + t501 * t477;
t614 = g(3) * t630 + t504 * t477;
t613 = t501 * t669 + t504 * t628;
t612 = g(1) * t631 + g(3) * t505;
t611 = t506 * pkin(1) + t503 * qJ(2);
t499 = t505 ^ 2;
t610 = t502 ^ 2 - t499;
t509 = qJD(3) ^ 2;
t609 = -t509 - t510;
t605 = qJD(3) * t502;
t599 = qJD(4) * t508;
t597 = qJD(5) * t504;
t360 = qJ(6) * t440 + t373;
t594 = qJD(5) - t360;
t590 = qJDD(1) * qJ(2);
t587 = g(1) * t625;
t586 = g(2) * t629;
t585 = 0.2e1 * qJD(1) * qJD(2);
t582 = t503 * t624;
t581 = -t501 * t383 - t504 * t389 - t414 * t601;
t575 = t508 * t604;
t580 = t501 * t437 + t504 * t575 + t601 * t669;
t571 = t502 * t599;
t579 = t501 * t575 + t504 * t571 + t602 * t669;
t444 = t505 * t462;
t388 = -qJDD(3) * pkin(3) + t465 * t605 - t444;
t391 = t502 * qJ(5) + t613;
t576 = t469 * t604;
t574 = t357 * t602;
t573 = t357 * t601;
t572 = t469 * t601;
t568 = -t497 + t657;
t567 = MDP(20) - t675;
t416 = t501 * t506 + t502 * t626;
t566 = -t415 * pkin(4) + qJ(5) * t416;
t418 = t502 * t622 - t627;
t565 = t417 * pkin(4) - qJ(5) * t418;
t463 = t501 * t628;
t563 = t504 * t669 - t463;
t562 = qJDD(2) - t654;
t560 = -qJD(1) * t635 - t436 * t633 - t501 * t576 - t502 * t572;
t559 = pkin(4) * t582 + pkin(8) * t631 + (pkin(3) + t634) * t625;
t350 = t377 * pkin(4) + t376 * qJ(5) - t440 * qJD(5) + t388;
t347 = -pkin(5) * t377 + qJDD(6) - t350;
t558 = t347 - t587;
t556 = g(1) * t417 + g(2) * t415;
t555 = -g(1) * t418 - g(2) * t416;
t554 = g(1) * t506 + g(2) * t503;
t553 = pkin(4) * t504 + t634;
t551 = -t653 - t657;
t550 = qJ(5) * t604 + t502 * qJD(5) + t580;
t529 = -t421 * t602 - t581;
t348 = t529 + t676;
t349 = -t561 - t665;
t549 = t348 * t504 + t349 * t501;
t353 = -t469 * t661 + t594;
t548 = t353 * t504 - t356 * t501;
t547 = t353 * t501 + t356 * t504;
t364 = -pkin(4) * t469 + t593;
t546 = t364 * t504 - t365 * t501;
t545 = t364 * t501 + t365 * t504;
t541 = t586 - t612;
t540 = -MDP(19) * t501 - MDP(20) * t504;
t539 = pkin(3) + t553;
t537 = t437 * t504 - t579;
t533 = t436 * t501 + t572;
t532 = t436 * t504 - t469 * t602;
t531 = 0.2e1 * qJ(2) * t591 + qJDD(3) * t508;
t530 = -t469 * t655 - t587;
t526 = pkin(3) * t631 + t416 * pkin(4) + t506 * pkin(7) + qJ(5) * t415 + t611;
t525 = -t350 + t530;
t524 = t422 * t469 - t658;
t491 = t506 * qJ(2);
t523 = pkin(3) * t629 + t418 * pkin(4) - pkin(8) * t623 + qJ(5) * t417 + t491;
t521 = -t554 + t585 + 0.2e1 * t590;
t517 = -t508 * t509 + t521;
t515 = t367 * t440 - t519;
t513 = g(1) * t416 - g(2) * t418 + g(3) * t624 - t529;
t512 = t373 * t469 + t513;
t511 = t546 * MDP(24) + t548 * MDP(28) + (-t504 * MDP(19) + t501 * t567) * t469 + t666 * (t636 + t640);
t486 = qJDD(3) * t505;
t468 = qJ(5) * t624;
t454 = t656 * t501;
t404 = t438 * t605;
t394 = -t468 + (-t508 + t660) * t505;
t392 = -pkin(4) * t502 - t563;
t390 = t468 + (t508 - t584) * t505;
t386 = pkin(4) * t440 + t651;
t382 = qJ(6) * t632 + t391;
t381 = -pkin(4) * t607 + t564;
t371 = t463 + (-t669 - t648) * t504 - t661 * t502;
t368 = -t440 * t661 - t651;
t363 = (qJD(4) * t553 - t597) * t505 + (t508 - t552) * t605;
t358 = -pkin(4) * t604 - t537;
t355 = -t501 * t571 + t550;
t354 = (qJD(4) * t534 + t597) * t505 + (-t508 + t535) * t605;
t352 = qJ(6) * t504 * t600 + (qJD(6) * t505 + (-qJ(6) * qJD(3) - t599) * t502) * t501 + t550;
t351 = (qJ(6) * t605 - t437) * t504 + (qJ(6) * t602 - qJD(3) * t661 - t596) * t505 + t579;
t346 = t348 + t673;
t345 = -qJD(6) * t440 + t349 + t649 - t659;
t1 = [t521 * MDP(5) + (-qJDD(3) * t502 - t505 * t509) * MDP(10) + (-t502 * t509 + t486) * MDP(9) + (t502 * t517 + t505 * t531) * MDP(12) + (-t502 * t531 + t505 * t517) * MDP(13) + qJDD(1) * MDP(1) + (t346 * t382 + t356 * t352 + t345 * t371 + t353 * t351 + t347 * t390 + t357 * t354 - g(1) * (pkin(5) * t418 + qJ(6) * t623 + t523) - g(2) * (pkin(5) * t416 + t526) + (g(2) * t505 * t656 - g(1) * t508) * t503) * MDP(28) + (t348 * t391 + t365 * t355 + t350 * t394 + t367 * t363 + t349 * t392 + t364 * t358 - g(1) * (t503 * t508 + t523) - g(2) * (-pkin(8) * t625 + t526)) * MDP(24) + (qJDD(2) - t568 - 0.2e1 * t654) * MDP(4) + ((t637 + t639) * t605 + (t644 - t643 + (-t636 + t640) * qJD(4)) * t505) * MDP(15) + (-t376 * t624 - t440 * t528) * MDP(14) + 0.2e1 * (-t502 * t588 + t591 * t610) * MDP(8) + (-t562 * pkin(1) - g(1) * (-pkin(1) * t503 + t491) - g(2) * t611 + (t585 + t590) * qJ(2)) * MDP(6) + (-t580 * t469 - t613 * t436 + ((t469 * t508 + t421) * t602 + (-t422 * t504 + t440 * t508) * qJD(3) + t581) * t502 + (-qJD(3) * t374 + t376 * t508 + t388 * t504 - t422 * t602) * t505 + t556) * MDP(20) + (-t351 * t440 + t352 * t438 + t371 * t376 + t377 * t382 + t548 * t605 + (qJD(4) * t547 - t345 * t504 + t346 * t501 - t554) * t505) * MDP(27) + (-t358 * t469 + t363 * t438 + t377 * t394 - t392 * t436 + (-t367 * t606 - t349) * t502 + (-qJD(3) * t364 + t350 * t501 + t367 * t601) * t505 + t555) * MDP(21) + ((t469 * t606 - t377) * t502 + (-qJD(3) * t438 - t533) * t505) * MDP(17) + (-t351 * t469 - t354 * t438 - t371 * t436 - t377 * t390 + (t357 * t606 - t345) * t502 + (-qJD(3) * t353 - t347 * t501 - t573) * t505 + t555) * MDP(25) + (t355 * t469 - t363 * t440 + t376 * t394 + t391 * t436 + (t367 * t595 + t348) * t502 + (qJD(3) * t365 - t350 * t504 + t367 * t602) * t505 - t556) * MDP(23) + (-t355 * t438 + t358 * t440 - t376 * t392 - t377 * t391 - t546 * t605 + (-qJD(4) * t545 - t348 * t501 + t349 * t504 + t554) * t505) * MDP(22) + (t537 * t469 + t563 * t436 + ((-t422 * t501 + t438 * t508) * qJD(3) + t561) * t502 + (qJD(3) * t373 - t377 * t508 + t388 * t501 + t422 * t601) * t505 + t555) * MDP(19) + ((-t469 * t595 - t376) * t502 + (qJD(3) * t440 + t532) * t505) * MDP(16) + (t352 * t469 + t354 * t440 - t376 * t390 + t382 * t436 + (-t357 * t595 + t346) * t502 + (qJD(3) * t356 + t347 * t504 - t574) * t505 - t556) * MDP(26) + (t436 * t502 + t576) * MDP(18) + t554 * MDP(3) + t568 * MDP(2) + (qJDD(1) * t499 - 0.2e1 * t502 * t570) * MDP(7); qJDD(1) * MDP(4) - t510 * MDP(5) + (t497 + t551 + t562) * MDP(6) + t486 * MDP(12) + t404 * MDP(19) + (t404 + t560) * MDP(21) + t560 * MDP(25) + t511 * qJD(1) + (t609 * MDP(13) - t350 * MDP(24) + t347 * MDP(28) + (-MDP(19) - t667) * t377 + t567 * t376 + (MDP(24) * t545 + MDP(28) * t547 + t469 * t540 - t666 * (-t637 + t639)) * qJD(3)) * t505 + (t609 * MDP(12) - qJDD(3) * MDP(13) + t549 * MDP(24) + (t345 * t501 + t346 * t504) * MDP(28) + t540 * t436 + (t367 * MDP(24) + t438 * MDP(25) - t357 * MDP(28) + t440 * t567) * qJD(3) + t511 * qJD(4) - t666 * (t643 + t644)) * t502 + t675 * (t436 * t630 + t504 * t576) + (-MDP(24) - MDP(28)) * t568; MDP(9) * t588 - MDP(10) * t589 + qJDD(3) * MDP(11) + (t505 * t551 + t444 + t668) * MDP(12) + ((-t462 + t653 - t497) * t502 + t612) * MDP(13) + (t440 * t635 - t644) * MDP(14) + (t677 * t501 + t504 * t679) * MDP(15) + ((-t440 * t505 + t469 * t630) * qJD(1) + t533) * MDP(16) + ((t438 * t505 - t469 * t633) * qJD(1) + t532) * MDP(17) - t469 * MDP(18) * t607 + (-t373 * t607 - t438 * t446 - pkin(3) * t377 + t432 * t469 + (-t587 - t388 + (-t443 - t655) * t469) * t504 + t524 * t501 + t614) * MDP(19) + (pkin(3) * t376 + t616 * t469 + t374 * t607 - t440 * t446 + t524 * t504 + (t388 - t530) * t501 - t615) * MDP(20) + (t364 * t607 - t377 * t539 + t381 * t469 + t618 * t438 + t501 * t664 + t525 * t504 + t614) * MDP(21) + (t380 * t438 - t381 * t440 + (t348 + t469 * t364 + (-t377 + t603) * pkin(8)) * t504 + (t349 - t646 + (qJD(4) * t438 - t376) * pkin(8)) * t501 + t541) * MDP(22) + (-t365 * t607 - t376 * t539 - t380 * t469 - t618 * t440 + t525 * t501 - t504 * t664 + t615) * MDP(23) + (-t365 * t380 - t364 * t381 - g(1) * t559 - g(3) * t492 + t618 * t367 + (qJD(4) * t546 + t549 + t586) * pkin(8) + (-t350 + t668) * t539) * MDP(24) + (-t574 + t377 * t674 - t436 * t454 + t558 * t504 - t617 * t469 - t621 * t438 + (t353 * t505 - t357 * t633) * qJD(1) + t614) * MDP(25) + (t573 + t376 * t674 + t436 * t455 + t558 * t501 + t619 * t469 + t621 * t440 + (-t356 * t505 + t357 * t630) * qJD(1) + t615) * MDP(26) + (t376 * t454 + t377 * t455 - t617 * t440 + t619 * t438 + (-t353 * t469 - t346) * t504 + (-t345 + t647) * t501 - t541) * MDP(27) + (t346 * t455 + t345 * t454 - t347 * t674 - g(1) * (pkin(5) * t582 + t559) - g(3) * (t492 - t648) + t621 * t357 + t619 * t356 + t617 * t353 + (qJ(6) * t657 - g(3) * (-pkin(5) * t504 - t539)) * t502 - (-t656 * t502 + t505 * t674) * t497) * MDP(28) + (t505 * t502 * MDP(7) - t610 * MDP(8)) * t510; MDP(14) * t642 + (t435 - t662) * MDP(15) - t678 * MDP(16) + (-t377 + t638) * MDP(17) + t436 * MDP(18) + (-t422 * t440 + t522 + t645) * MDP(19) + (t422 * t438 + t512) * MDP(20) + (-t386 * t438 + t425 - t515 + t645) * MDP(21) + (pkin(4) * t376 - t652 + (t365 - t374) * t440 + (t364 - t593) * t438) * MDP(22) + (-t367 * t438 + t386 * t440 - t512 + t672) * MDP(23) + (t348 * qJ(5) - t349 * pkin(4) - t367 * t386 - t364 * t374 - g(1) * t566 - g(2) * t565 - g(3) * (-pkin(4) * t632 + t468) + t593 * t365) * MDP(24) + (t361 * t469 + t368 * t438 + (pkin(5) + t661) * t436 + t663) * MDP(25) + (t357 * t438 - t360 * t469 - t368 * t440 - t513 + t672 + t673) * MDP(26) + (t652 - t376 * t661 + (-t356 + t361) * t440 + (-t353 + t594) * t438) * MDP(27) + (t346 * qJ(5) - t345 * t661 - t353 * t361 - t357 * t368 - g(1) * (-pkin(5) * t415 + t566) - g(2) * (pkin(5) * t417 + t565) - g(3) * (-t501 * t583 + t468) + t594 * t356) * MDP(28); (t515 - t646) * MDP(24) + (-t647 - t659 - t663) * MDP(28) + t667 * (-t436 + t642) + t675 * (-t469 ^ 2 - t435) - t666 * t678; t677 * MDP(25) + t679 * MDP(26) + (-t435 - t662) * MDP(27) + (t353 * t440 - t356 * t438 + t558 + t668) * MDP(28);];
tau  = t1;
