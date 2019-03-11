% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:33
% EndTime: 2019-03-09 06:55:41
% DurationCPUTime: 6.29s
% Computational Cost: add. (7185->444), mult. (15647->581), div. (0->0), fcn. (11763->18), ass. (0->220)
t545 = sin(qJ(4));
t546 = sin(qJ(3));
t636 = qJD(1) * t546;
t615 = t545 * t636;
t550 = cos(qJ(4));
t551 = cos(qJ(3));
t635 = qJD(1) * t551;
t616 = t550 * t635;
t472 = -t615 + t616;
t473 = -t545 * t635 - t550 * t636;
t544 = sin(qJ(5));
t549 = cos(qJ(5));
t430 = t549 * t472 + t473 * t544;
t548 = cos(qJ(6));
t628 = qJD(6) * t548;
t691 = -t430 * t548 + t628;
t578 = t472 * t544 - t549 * t473;
t536 = qJD(3) + qJD(4);
t625 = qJD(1) * qJD(3);
t614 = t551 * t625;
t623 = qJDD(1) * t551;
t624 = qJDD(1) * t546;
t417 = qJD(4) * t616 - t536 * t615 + t545 * t623 + (t614 + t624) * t550;
t480 = t545 * t551 + t546 * t550;
t441 = t536 * t480;
t587 = t545 * t624 - t550 * t623;
t418 = qJD(1) * t441 + t587;
t630 = qJD(5) * t549;
t631 = qJD(5) * t544;
t385 = t549 * t417 - t544 * t418 + t472 * t630 + t473 * t631;
t535 = qJDD(3) + qJDD(4);
t526 = qJDD(5) + t535;
t529 = qJD(5) + t536;
t543 = sin(qJ(6));
t617 = t548 * t385 + t543 * t526 + t529 * t628;
t629 = qJD(6) * t543;
t365 = -t578 * t629 + t617;
t363 = t365 * t543;
t364 = t365 * t548;
t421 = t529 * t543 + t548 * t578;
t502 = t548 * t526;
t366 = qJD(6) * t421 + t385 * t543 - t502;
t386 = qJD(5) * t578 + t417 * t544 + t549 * t418;
t384 = qJDD(6) + t386;
t378 = t543 * t384;
t379 = t548 * t384;
t419 = -t548 * t529 + t543 * t578;
t681 = qJD(6) - t430;
t689 = t681 * t543;
t690 = t526 * MDP(23) - t386 * MDP(22) - t430 ^ 2 * MDP(20) + (-t430 * t529 + t385) * MDP(21) + (-MDP(19) * t430 + MDP(20) * t578 + MDP(22) * t529 - MDP(30) * t681) * t578 + (t691 * t421 + t363) * MDP(26) + (-t421 * t578 + t691 * t681 + t378) * MDP(28) + (-t543 * t366 - t691 * t419 - t421 * t689 + t364) * MDP(27) + (t419 * t578 - t681 * t689 + t379) * MDP(29);
t540 = qJ(3) + qJ(4);
t534 = qJ(5) + t540;
t519 = sin(t534);
t537 = qJ(1) + pkin(11);
t527 = sin(t537);
t528 = cos(t537);
t682 = g(1) * t528 + g(2) * t527;
t688 = t682 * t519;
t541 = sin(pkin(11));
t514 = pkin(1) * t541 + pkin(7);
t662 = pkin(8) + t514;
t470 = t473 * pkin(9);
t611 = t662 * qJD(1);
t449 = qJD(2) * t546 + t551 * t611;
t443 = t545 * t449;
t448 = t551 * qJD(2) - t611 * t546;
t661 = qJD(3) * pkin(3);
t446 = t448 + t661;
t604 = t550 * t446 - t443;
t401 = t470 + t604;
t398 = pkin(4) * t536 + t401;
t445 = t550 * t449;
t579 = -t446 * t545 - t445;
t666 = pkin(9) * t472;
t402 = -t579 + t666;
t657 = t402 * t544;
t373 = t398 * t549 - t657;
t371 = -pkin(5) * t529 - t373;
t660 = t371 * t430;
t530 = t551 * qJDD(2);
t492 = t514 * qJDD(1);
t609 = pkin(8) * qJDD(1) + t492;
t409 = qJDD(3) * pkin(3) - qJD(3) * t449 - t609 * t546 + t530;
t415 = qJD(3) * t448 + t546 * qJDD(2) + t609 * t551;
t568 = qJD(4) * t579 + t550 * t409 - t545 * t415;
t360 = pkin(4) * t535 - pkin(9) * t417 + t568;
t633 = qJD(4) * t545;
t671 = (qJD(4) * t446 + t415) * t550 + t545 * t409 - t449 * t633;
t361 = -pkin(9) * t418 + t671;
t656 = t402 * t549;
t374 = t398 * t544 + t656;
t670 = qJD(5) * t374 - t549 * t360 + t544 * t361;
t347 = -pkin(5) * t526 + t670;
t520 = cos(t534);
t663 = g(3) * t520;
t684 = t347 + t663;
t542 = cos(pkin(11));
t515 = -pkin(1) * t542 - pkin(2);
t489 = -pkin(3) * t551 + t515;
t474 = t489 * qJD(1);
t438 = -pkin(4) * t472 + t474;
t511 = g(3) * t519;
t672 = (qJD(5) * t398 + t361) * t549 + t544 * t360 - t402 * t631;
t561 = -t438 * t430 + t520 * t682 + t511 - t672;
t678 = pkin(5) * t578;
t677 = t681 * (pkin(10) * t681 + t678);
t477 = t662 * t546;
t478 = t662 * t551;
t640 = -t545 * t477 + t550 * t478;
t372 = pkin(10) * t529 + t374;
t388 = -pkin(5) * t430 - pkin(10) * t578 + t438;
t355 = -t372 * t543 + t388 * t548;
t676 = -t355 * t578 + t371 * t629 + t548 * t688;
t356 = t372 * t548 + t388 * t543;
t675 = t356 * t578 + t371 * t628 + t543 * t684;
t558 = -t438 * t578 - t663 - t670 + t688;
t612 = qJD(3) * t662;
t468 = t546 * t612;
t469 = t551 * t612;
t632 = qJD(4) * t550;
t573 = -t550 * t468 - t545 * t469 - t477 * t632 - t478 * t633;
t389 = -pkin(9) * t441 + t573;
t479 = t545 * t546 - t550 * t551;
t440 = t536 * t479;
t567 = -qJD(4) * t640 + t468 * t545 - t550 * t469;
t390 = pkin(9) * t440 + t567;
t602 = -t550 * t477 - t478 * t545;
t423 = -pkin(9) * t480 + t602;
t424 = -pkin(9) * t479 + t640;
t581 = t423 * t549 - t424 * t544;
t353 = qJD(5) * t581 + t389 * t549 + t390 * t544;
t392 = t423 * t544 + t424 * t549;
t436 = t549 * t479 + t480 * t544;
t394 = -qJD(5) * t436 - t440 * t549 - t441 * t544;
t437 = -t479 * t544 + t480 * t549;
t447 = pkin(4) * t479 + t489;
t397 = pkin(5) * t436 - pkin(10) * t437 + t447;
t346 = pkin(10) * t526 + t672;
t597 = qJD(6) * t388 + t346;
t669 = t347 * t437 + t371 * t394 - t392 * t384 - (qJD(6) * t397 + t353) * t681 - t436 * t597;
t667 = pkin(4) * t473;
t659 = t371 * t437;
t658 = t397 * t384;
t653 = t527 * t543;
t652 = t527 * t548;
t651 = t528 * t543;
t650 = t528 * t548;
t649 = t544 * t545;
t648 = t545 * t549;
t647 = qJDD(2) - g(3);
t395 = qJD(5) * t437 - t440 * t544 + t549 * t441;
t646 = t365 * t436 + t421 * t395;
t603 = -t448 * t545 - t445;
t406 = t603 - t666;
t641 = t550 * t448 - t443;
t407 = t470 + t641;
t523 = pkin(3) * t550 + pkin(4);
t643 = t406 * t544 + t407 * t549 - t523 * t630 - (-t545 * t631 + (t549 * t550 - t649) * qJD(4)) * pkin(3);
t642 = t406 * t549 - t407 * t544 + t523 * t631 + (t545 * t630 + (t544 * t550 + t648) * qJD(4)) * pkin(3);
t639 = pkin(3) * t648 + t544 * t523;
t538 = t546 ^ 2;
t638 = -t551 ^ 2 + t538;
t496 = qJD(1) * t515;
t525 = t546 * t661;
t524 = pkin(3) * t636;
t622 = t437 * t378;
t621 = t437 * t379;
t432 = pkin(4) * t441 + t525;
t455 = qJD(3) * t524 + qJDD(1) * t489;
t403 = pkin(4) * t418 + t455;
t349 = pkin(5) * t386 - pkin(10) * t385 + t403;
t596 = qJD(6) * t372 - t349;
t396 = -pkin(10) * t430 - t667 + t678;
t467 = pkin(10) + t639;
t594 = qJD(6) * t467 + t396 + t524;
t521 = pkin(4) * t544 + pkin(10);
t593 = qJD(6) * t521 + t396;
t375 = t401 * t544 + t656;
t592 = pkin(4) * t631 - t375;
t376 = t401 * t549 - t657;
t591 = -pkin(4) * t630 + t376;
t589 = g(1) * t527 - g(2) * t528;
t547 = sin(qJ(1));
t552 = cos(qJ(1));
t588 = g(1) * t547 - g(2) * t552;
t585 = -t366 * t436 - t395 * t419;
t584 = -t384 * t467 - t660;
t583 = -t384 * t521 - t660;
t582 = t394 * t529 + t437 * t526;
t580 = -t440 * t536 + t480 * t535;
t577 = -pkin(3) * t649 + t523 * t549;
t575 = -t394 * t543 - t437 * t628;
t574 = -t394 * t548 + t437 * t629;
t572 = -pkin(10) * t384 + t373 * t681 - t660;
t570 = -qJD(1) * t496 - t492 + t682;
t569 = 0.2e1 * qJD(3) * t496 - qJDD(3) * t514;
t565 = -t548 * t684 + t676;
t553 = qJD(3) ^ 2;
t564 = -0.2e1 * qJDD(1) * t515 - t514 * t553 + t589;
t562 = -t543 * t688 + t675;
t532 = sin(t540);
t533 = cos(t540);
t560 = g(3) * t532 - t474 * t472 + t533 * t682 - t671;
t557 = -g(3) * t533 + t474 * t473 + t532 * t682 + t568;
t556 = t473 * t472 * MDP(12) + (-t472 * t536 + t417) * MDP(14) + (-t587 + (-qJD(1) * t480 - t473) * t536) * MDP(15) + (-t472 ^ 2 + t473 ^ 2) * MDP(13) + t535 * MDP(16) + t690;
t522 = -pkin(4) * t549 - pkin(5);
t491 = qJDD(3) * t551 - t546 * t553;
t490 = qJDD(3) * t546 + t551 * t553;
t466 = -pkin(5) - t577;
t454 = t520 * t650 + t653;
t453 = -t520 * t651 + t652;
t452 = -t520 * t652 + t651;
t451 = t520 * t653 + t650;
t450 = t524 - t667;
t416 = -t441 * t536 - t479 * t535;
t387 = -t395 * t529 - t436 * t526;
t358 = pkin(5) * t395 - pkin(10) * t394 + t432;
t354 = qJD(5) * t392 + t389 * t544 - t390 * t549;
t348 = t548 * t349;
t1 = [(-t574 * t681 + t621 + t646) * MDP(28) + (t384 * t436 + t395 * t681) * MDP(30) + (t575 * t681 + t585 - t622) * MDP(29) + (-g(1) * t451 - g(2) * t453 + t354 * t421 - t356 * t395 - t581 * t365 + (-(-qJD(6) * t392 + t358) * t681 - t658 + t596 * t436 - qJD(6) * t659) * t543 + t669 * t548) * MDP(32) + (-g(1) * t452 - g(2) * t454 + t348 * t436 + t354 * t419 + t355 * t395 - t581 * t366 + (t358 * t681 + t658 + (-t372 * t436 - t392 * t681 + t659) * qJD(6)) * t548 + t669 * t543) * MDP(31) + (-t353 * t529 + t385 * t447 - t392 * t526 + t394 * t438 + t403 * t437 + t432 * t578 - t519 * t589) * MDP(25) + (t385 * t437 + t394 * t578) * MDP(19) + qJDD(1) * MDP(1) + (t546 * t569 + t551 * t564) * MDP(10) + (-t546 * t564 + t551 * t569) * MDP(11) + (-t385 * t436 - t386 * t437 + t394 * t430 - t395 * t578) * MDP(20) + (-t354 * t529 + t386 * t447 + t395 * t438 + t403 * t436 - t430 * t432 + t520 * t589 + t526 * t581) * MDP(24) + (t417 * t480 + t440 * t473) * MDP(12) + (-t417 * t479 - t418 * t480 - t440 * t472 + t441 * t473) * MDP(13) + (t489 * t417 - t474 * t440 + t455 * t480 - t473 * t525 - t532 * t589 - t535 * t640 - t536 * t573) * MDP(18) + t387 * MDP(22) + (t588 + (t541 ^ 2 + t542 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t580 * MDP(14) + t582 * MDP(21) + t416 * MDP(15) + t588 * MDP(2) + (qJDD(1) * t538 + 0.2e1 * t546 * t614) * MDP(5) + t490 * MDP(7) + t491 * MDP(8) + (t489 * t418 + t474 * t441 + t455 * t479 - t472 * t525 + t533 * t589 + t535 * t602 + t536 * t567) * MDP(17) + 0.2e1 * (t546 * t623 - t625 * t638) * MDP(6) + (g(1) * t552 + g(2) * t547) * MDP(3) + (t364 * t437 - t421 * t574) * MDP(26) + ((-t419 * t548 - t421 * t543) * t394 + (-t363 - t366 * t548 + (t419 * t543 - t421 * t548) * qJD(6)) * t437) * MDP(27); t647 * MDP(4) + t491 * MDP(10) - t490 * MDP(11) + t416 * MDP(17) - t580 * MDP(18) + t387 * MDP(24) - t582 * MDP(25) + (-t585 - t622) * MDP(31) + (-t621 + t646) * MDP(32) + (MDP(31) * t575 + MDP(32) * t574) * t681; (t466 * t365 + t584 * t548 + t642 * t421 + (t543 * t594 + t548 * t643) * t681 + t562) * MDP(32) + (t641 * t536 + (t473 * t636 - t535 * t545 - t536 * t632) * pkin(3) + t560) * MDP(18) + (-t603 * t536 + (t472 * t636 + t535 * t550 - t536 * t633) * pkin(3) + t557) * MDP(17) + (t466 * t366 + t584 * t543 + t642 * t419 + (t543 * t643 - t548 * t594) * t681 + t565) * MDP(31) + MDP(8) * t623 + MDP(7) * t624 + qJDD(3) * MDP(9) + t556 + (-g(3) * t551 + t546 * t570 + t530) * MDP(10) + (t430 * t450 + t526 * t577 - t529 * t642 + t558) * MDP(24) + (-t546 * t647 + t551 * t570) * MDP(11) + (-t450 * t578 - t526 * t639 + t529 * t643 + t561) * MDP(25) + (-MDP(5) * t546 * t551 + MDP(6) * t638) * qJD(1) ^ 2; (t375 * t529 + (-t430 * t473 + t526 * t549 - t529 * t631) * pkin(4) + t558) * MDP(24) + (t522 * t365 + t583 * t548 + t592 * t421 + (t543 * t593 + t548 * t591) * t681 + t562) * MDP(32) + (t376 * t529 + (t473 * t578 - t526 * t544 - t529 * t630) * pkin(4) + t561) * MDP(25) + (t522 * t366 + t583 * t543 + t592 * t419 + (t543 * t591 - t548 * t593) * t681 + t565) * MDP(31) + t556 + (t536 * t604 + t560) * MDP(18) + (-t536 * t579 + t557) * MDP(17); (t374 * t529 + t558) * MDP(24) + (t373 * t529 + t561) * MDP(25) + (-pkin(5) * t366 - t374 * t419 + t572 * t543 + (-t684 - t677) * t548 + t676) * MDP(31) + (-pkin(5) * t365 - t374 * t421 + t572 * t548 + (-t688 + t677) * t543 + t675) * MDP(32) + t690; t421 * t419 * MDP(26) + (-t419 ^ 2 + t421 ^ 2) * MDP(27) + (t419 * t681 + t617) * MDP(28) + (t421 * t681 + t502) * MDP(29) + t384 * MDP(30) + (-g(1) * t453 + g(2) * t451 + t356 * t681 - t371 * t421 + t348) * MDP(31) + (g(1) * t454 - g(2) * t452 + t355 * t681 + t371 * t419) * MDP(32) + ((-t346 + t511) * MDP(32) + (-MDP(29) * t578 - MDP(31) * t372 - MDP(32) * t388) * qJD(6)) * t548 + (-qJD(6) * t578 * MDP(28) + (-qJD(6) * t529 - t385) * MDP(29) + (-t597 + t511) * MDP(31) + t596 * MDP(32)) * t543;];
tau  = t1;
