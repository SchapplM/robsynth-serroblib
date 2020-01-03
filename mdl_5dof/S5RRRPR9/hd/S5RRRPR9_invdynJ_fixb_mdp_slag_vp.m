% Calculate vector of inverse dynamics joint torques for
% S5RRRPR9
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:52
% EndTime: 2019-12-31 21:25:05
% DurationCPUTime: 8.03s
% Computational Cost: add. (4581->483), mult. (10491->652), div. (0->0), fcn. (7519->12), ass. (0->217)
t547 = cos(qJ(5));
t544 = sin(qJ(3));
t545 = sin(qJ(2));
t616 = qJD(1) * t545;
t595 = t544 * t616;
t548 = cos(qJ(3));
t603 = t548 * qJD(2);
t498 = t595 - t603;
t612 = qJD(2) * t544;
t500 = t548 * t616 + t612;
t540 = sin(pkin(9));
t541 = cos(pkin(9));
t572 = -t498 * t541 - t500 * t540;
t631 = t547 * t572;
t435 = t498 * t540 - t500 * t541;
t543 = sin(qJ(5));
t641 = t435 * t543;
t390 = t631 + t641;
t549 = cos(qJ(2));
t615 = qJD(1) * t549;
t520 = -qJD(3) + t615;
t512 = -qJD(5) + t520;
t643 = t390 * t512;
t602 = qJD(1) * qJD(2);
t588 = t549 * t602;
t601 = qJDD(1) * t545;
t608 = qJD(3) * t545;
t664 = -qJD(1) * t608 + qJDD(2);
t428 = qJD(3) * t603 + (t588 + t601) * t548 + t664 * t544;
t581 = pkin(2) * t545 - pkin(7) * t549;
t502 = t581 * qJD(2);
t506 = -pkin(2) * t549 - pkin(7) * t545 - pkin(1);
t447 = qJD(1) * t502 + t506 * qJDD(1);
t441 = t548 * t447;
t488 = t506 * qJD(1);
t530 = pkin(6) * t615;
t510 = qJD(2) * pkin(7) + t530;
t446 = t488 * t544 + t510 * t548;
t534 = t549 * qJDD(1);
t659 = -t545 * t602 + t534;
t474 = pkin(6) * t659 + qJDD(2) * pkin(7);
t490 = qJDD(3) - t659;
t355 = pkin(3) * t490 - qJ(4) * t428 - t446 * qJD(3) - qJD(4) * t500 - t474 * t544 + t441;
t429 = t544 * ((qJD(3) + t615) * qJD(2) + t601) - t664 * t548;
t607 = qJD(3) * t548;
t609 = qJD(3) * t544;
t559 = t544 * t447 + t548 * t474 + t488 * t607 - t510 * t609;
t360 = -qJ(4) * t429 - qJD(4) * t498 + t559;
t343 = t541 * t355 - t360 * t540;
t378 = t428 * t541 - t429 * t540;
t341 = pkin(4) * t490 - pkin(8) * t378 + t343;
t344 = t540 * t355 + t541 * t360;
t377 = -t428 * t540 - t429 * t541;
t342 = pkin(8) * t377 + t344;
t445 = t548 * t488 - t510 * t544;
t410 = -qJ(4) * t500 + t445;
t402 = -pkin(3) * t520 + t410;
t411 = -qJ(4) * t498 + t446;
t637 = t541 * t411;
t367 = t540 * t402 + t637;
t662 = pkin(8) * t572;
t359 = t367 + t662;
t605 = qJD(5) * t543;
t357 = t359 * t605;
t509 = -qJD(2) * pkin(2) + pkin(6) * t616;
t451 = pkin(3) * t498 + qJD(4) + t509;
t396 = -pkin(4) * t572 + t451;
t535 = qJ(3) + pkin(9) + qJ(5);
t523 = sin(t535);
t524 = cos(t535);
t550 = cos(qJ(1));
t546 = sin(qJ(1));
t633 = t546 * t549;
t455 = t523 * t550 - t524 * t633;
t629 = t549 * t550;
t457 = t523 * t546 + t524 * t629;
t647 = g(3) * t545;
t674 = g(1) * t457 - g(2) * t455 - t543 * t341 - t547 * t342 - t396 * t390 + t524 * t647 + t357;
t487 = qJDD(5) + t490;
t663 = -t547 * t435 + t543 * t572;
t673 = t487 * MDP(24) + (-t390 ^ 2 + t663 ^ 2) * MDP(21) - t390 * t663 * MDP(20);
t644 = t663 * t512;
t630 = t548 * t549;
t653 = pkin(3) * t545;
t568 = -qJ(4) * t630 + t653;
t645 = qJ(4) + pkin(7);
t586 = qJD(3) * t645;
t501 = t581 * qJD(1);
t621 = pkin(6) * t595 + t548 * t501;
t671 = -t568 * qJD(1) - qJD(4) * t544 - t548 * t586 - t621;
t483 = t544 * t501;
t606 = qJD(4) * t548;
t634 = t545 * t548;
t635 = t544 * t549;
t670 = t483 + (-pkin(6) * t634 - qJ(4) * t635) * qJD(1) + t544 * t586 - t606;
t528 = pkin(6) * t601;
t475 = -qJDD(2) * pkin(2) + pkin(6) * t588 + t528;
t580 = g(1) * t550 + g(2) * t546;
t646 = g(3) * t549;
t557 = t580 * t545 - t646;
t669 = qJD(3) * pkin(7) * t520 - t475 + t557;
t454 = t523 * t633 + t524 * t550;
t456 = -t523 * t629 + t524 * t546;
t585 = t547 * t341 - t543 * t342;
t668 = -g(1) * t456 + g(2) * t454 - t396 * t663 + t523 * t647 + t585;
t667 = pkin(8) * t435;
t492 = t540 * t548 + t541 * t544;
t560 = t492 * t549;
t666 = qJD(1) * t560 - t492 * qJD(3);
t569 = t540 * t544 - t541 * t548;
t665 = t520 * t569;
t626 = t540 * t670 + t541 * t671;
t625 = t540 * t671 - t541 * t670;
t654 = pkin(3) * t544;
t658 = pkin(3) * t609 - t615 * t654 - t530;
t479 = t544 * t633 + t548 * t550;
t481 = -t544 * t629 + t546 * t548;
t657 = -g(1) * t481 + g(2) * t479;
t584 = -t547 * t377 + t378 * t543;
t348 = qJD(5) * t663 + t584;
t655 = pkin(3) * t540;
t652 = pkin(6) * t544;
t642 = t428 * t544;
t640 = t498 * t520;
t639 = t500 * t520;
t638 = t500 * t548;
t406 = t540 * t411;
t636 = t544 * t545;
t366 = t541 * t402 - t406;
t356 = -pkin(4) * t520 + t366 + t667;
t632 = t547 * t356;
t573 = -t492 * t543 - t547 * t569;
t628 = t573 * qJD(5) + t543 * t666 + t547 * t665;
t434 = t492 * t547 - t543 * t569;
t627 = t434 * qJD(5) + t543 * t665 - t547 * t666;
t522 = pkin(6) * t630;
t611 = qJD(2) * t545;
t622 = t548 * t502 + t611 * t652;
t383 = -t545 * t606 + t568 * qJD(2) + (-t522 + (qJ(4) * t545 - t506) * t544) * qJD(3) + t622;
t623 = t544 * t502 + t506 * t607;
t392 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t634 + (-qJD(4) * t545 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t549) * t544 + t623;
t354 = t540 * t383 + t541 * t392;
t371 = t541 * t410 - t406;
t494 = t548 * t506;
t442 = -qJ(4) * t634 + t494 + (-pkin(3) - t652) * t549;
t620 = t544 * t506 + t522;
t448 = -qJ(4) * t636 + t620;
t394 = t540 * t442 + t541 * t448;
t624 = -pkin(4) * t666 + t658;
t507 = t645 * t544;
t508 = t645 * t548;
t450 = -t540 * t507 + t541 * t508;
t618 = pkin(3) * t636 + t545 * pkin(6);
t538 = t545 ^ 2;
t617 = -t549 ^ 2 + t538;
t614 = qJD(2) * t498;
t613 = qJD(2) * t500;
t610 = qJD(2) * t549;
t598 = qJD(5) * t631 + t543 * t377 + t547 * t378;
t593 = t544 * t610;
t597 = pkin(3) * t593 + pkin(6) * t610 + t607 * t653;
t527 = pkin(3) * t548 + pkin(2);
t596 = pkin(6) + t654;
t594 = t520 * t603;
t592 = t549 * t603;
t591 = t520 * t609;
t590 = t520 * t607;
t353 = t541 * t383 - t392 * t540;
t370 = -t410 * t540 - t637;
t393 = t541 * t442 - t448 * t540;
t449 = -t541 * t507 - t508 * t540;
t583 = -qJD(3) * t488 - t474;
t579 = g(1) * t546 - g(2) * t550;
t578 = t510 * t607 - t441;
t422 = -pkin(8) * t569 + t450;
t577 = pkin(4) * t616 + t665 * pkin(8) + qJD(5) * t422 - t626;
t421 = -pkin(8) * t492 + t449;
t576 = t666 * pkin(8) + qJD(5) * t421 + t625;
t575 = -pkin(7) * t490 + qJD(3) * t509;
t346 = t543 * t356 + t547 * t359;
t468 = t492 * t545;
t469 = t569 * t545;
t574 = -t547 * t468 + t469 * t543;
t414 = -t468 * t543 - t469 * t547;
t571 = t527 * t549 + t545 * t645;
t525 = pkin(3) * t541 + pkin(4);
t567 = t525 * t543 + t547 * t655;
t566 = t525 * t547 - t543 * t655;
t564 = pkin(1) + t571;
t563 = -0.2e1 * pkin(1) * t602 - pkin(6) * qJDD(2);
t562 = t490 * t544 - t590;
t561 = t490 * t548 + t591;
t347 = t435 * t605 + t598;
t552 = qJD(1) ^ 2;
t558 = pkin(1) * t552 + t580;
t401 = pkin(3) * t429 + qJDD(4) + t475;
t551 = qJD(2) ^ 2;
t554 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t551 + t579;
t482 = t544 * t546 + t548 * t629;
t480 = t544 * t550 - t546 * t630;
t461 = pkin(4) * t569 - t527;
t443 = pkin(4) * t468 + t618;
t416 = t492 * t608 + t540 * t593 - t541 * t592;
t415 = -qJD(2) * t560 + t569 * t608;
t403 = pkin(3) * t500 - pkin(4) * t435;
t395 = -pkin(4) * t415 + t597;
t373 = -pkin(8) * t468 + t394;
t372 = -pkin(4) * t549 + pkin(8) * t469 + t393;
t365 = t371 + t667;
t364 = t370 - t662;
t363 = t414 * qJD(5) - t547 * t415 - t416 * t543;
t362 = t574 * qJD(5) + t415 * t543 - t416 * t547;
t361 = -pkin(4) * t377 + t401;
t350 = pkin(8) * t415 + t354;
t349 = pkin(4) * t611 + pkin(8) * t416 + t353;
t345 = -t359 * t543 + t632;
t1 = [((-t498 * t548 - t500 * t544) * t610 + (-t642 - t429 * t548 + (t498 * t544 - t638) * qJD(3)) * t545) * MDP(12) + (t563 * t545 + t554 * t549) * MDP(9) + (-t554 * t545 + t563 * t549) * MDP(10) + (t347 * t574 - t348 * t414 + t362 * t390 - t363 * t663) * MDP(21) + (t348 * t549 + t363 * t512 + t390 * t611 + t487 * t574) * MDP(23) + (-(t349 * t547 - t350 * t543) * t512 + (t372 * t547 - t373 * t543) * t487 - t585 * t549 + t345 * t611 - t395 * t390 + t443 * t348 - t361 * t574 + t396 * t363 - g(1) * t455 - g(2) * t457 + (-(-t372 * t543 - t373 * t547) * t512 + t346 * t549) * qJD(5)) * MDP(25) + (-t347 * t549 - t362 * t512 + t414 * t487 + t611 * t663) * MDP(22) + (-t346 * t611 - g(1) * t454 - g(2) * t456 + t443 * t347 - t357 * t549 + t361 * t414 + t396 * t362 + t395 * t663 + ((-qJD(5) * t373 + t349) * t512 - t372 * t487 + t341 * t549) * t543 + ((qJD(5) * t372 + t350) * t512 - t373 * t487 + (qJD(5) * t356 + t342) * t549) * t547) * MDP(26) + (t347 * t414 + t362 * t663) * MDP(20) + 0.2e1 * (t545 * t534 - t617 * t602) * MDP(5) + (t344 * t394 + t367 * t354 + t343 * t393 + t366 * t353 + t401 * t618 + t451 * t597 + (-g(1) * t596 - g(2) * t564) * t550 + (g(1) * t564 - g(2) * t596) * t546) * MDP(19) + (-(-t506 * t609 + t622) * t520 + t494 * t490 - g(1) * t480 - g(2) * t482 + ((t590 + t614) * pkin(6) + (-pkin(6) * t490 + qJD(2) * t509 - t583) * t544 + t578) * t549 + (pkin(6) * t429 + qJD(2) * t445 + t475 * t544 + t509 * t607) * t545) * MDP(16) + (t623 * t520 - t620 * t490 - g(1) * t479 - g(2) * t481 + (t509 * t603 + (-t591 + t613) * pkin(6) + t559) * t549 + (-t509 * t609 - t446 * qJD(2) + t475 * t548 + (t428 - t594) * pkin(6)) * t545) * MDP(17) + (t343 * t469 - t344 * t468 + t353 * t435 + t354 * t572 + t366 * t416 + t367 * t415 + t377 * t394 - t378 * t393 + t579 * t545) * MDP(18) + (-t490 * t549 - t520 * t611) * MDP(15) + (-t487 * t549 - t512 * t611) * MDP(24) + ((-t428 - t594) * t549 + (t561 + t613) * t545) * MDP(13) + ((t520 * t612 + t429) * t549 + (-t562 - t614) * t545) * MDP(14) + (t428 * t634 + (-t544 * t608 + t592) * t500) * MDP(11) + (qJDD(1) * t538 + 0.2e1 * t545 * t588) * MDP(4) + (qJDD(2) * t545 + t549 * t551) * MDP(6) + (qJDD(2) * t549 - t545 * t551) * MDP(7) + t579 * MDP(2) + t580 * MDP(3) + qJDD(1) * MDP(1); MDP(6) * t601 + MDP(7) * t534 + qJDD(2) * MDP(8) + (t558 * t545 - t528 - t646) * MDP(9) + (t647 + (-pkin(6) * qJDD(1) + t558) * t549) * MDP(10) + (-t520 * t638 + t642) * MDP(11) + ((t428 + t640) * t548 + (-t429 + t639) * t544) * MDP(12) + ((-t500 * t545 + t520 * t630) * qJD(1) + t562) * MDP(13) + ((t498 * t545 - t520 * t635) * qJD(1) + t561) * MDP(14) + (-pkin(2) * t429 + t621 * t520 + t575 * t544 + (-t445 * t545 + (-pkin(6) * t498 - t509 * t544) * t549) * qJD(1) + t669 * t548) * MDP(16) + (-pkin(2) * t428 - t483 * t520 + t575 * t548 + (-t509 * t630 + t446 * t545 + (-t500 * t549 + t520 * t634) * pkin(6)) * qJD(1) - t669 * t544) * MDP(17) + (-t343 * t492 - t344 * t569 - t366 * t665 + t367 * t666 + t377 * t450 - t378 * t449 + t626 * t435 - t580 * t549 + t625 * t572 - t647) * MDP(18) + (t344 * t450 + t343 * t449 - t401 * t527 - g(3) * t571 + t658 * t451 + t625 * t367 + t626 * t366 + t580 * (t527 * t545 - t549 * t645)) * MDP(19) + (t347 * t434 + t628 * t663) * MDP(20) + (t347 * t573 - t348 * t434 + t390 * t628 - t627 * t663) * MDP(21) + (t434 * t487 - t628 * t512) * MDP(22) + (t487 * t573 + t627 * t512) * MDP(23) + ((t421 * t547 - t422 * t543) * t487 + t461 * t348 - t361 * t573 + (t576 * t543 + t577 * t547) * t512 + t627 * t396 - t624 * t390 + t557 * t524) * MDP(25) + (-(t421 * t543 + t422 * t547) * t487 + t461 * t347 + t361 * t434 + (-t577 * t543 + t576 * t547) * t512 + t628 * t396 + t624 * t663 - t557 * t523) * MDP(26) + (t520 * MDP(15) - MDP(22) * t663 - MDP(23) * t390 + t512 * MDP(24) - t345 * MDP(25) + t346 * MDP(26)) * t616 + (-MDP(4) * t545 * t549 + MDP(5) * t617) * t552; t500 * t498 * MDP(11) + (-t498 ^ 2 + t500 ^ 2) * MDP(12) + (t428 - t640) * MDP(13) + (-t429 - t639) * MDP(14) + t490 * MDP(15) + (-t446 * t520 - t500 * t509 + (t583 + t647) * t544 - t578 + t657) * MDP(16) + (g(1) * t482 - g(2) * t480 + g(3) * t634 - t445 * t520 + t498 * t509 - t559) * MDP(17) + ((t377 * t540 - t378 * t541) * pkin(3) + (t366 - t371) * t572 + (-t367 - t370) * t435) * MDP(18) + (-t366 * t370 - t367 * t371 + (g(3) * t636 + t343 * t541 + t344 * t540 - t451 * t500 + t657) * pkin(3)) * MDP(19) + (t347 + t643) * MDP(22) + (-t348 - t644) * MDP(23) + (t566 * t487 + (t364 * t547 - t365 * t543) * t512 + t403 * t390 + (t567 * t512 - t346) * qJD(5) + t668) * MDP(25) + (-t567 * t487 - (t364 * t543 + t365 * t547) * t512 - t403 * t663 + (t566 * t512 - t632) * qJD(5) + t674) * MDP(26) + t673; (-t435 ^ 2 - t572 ^ 2) * MDP(18) + (-t366 * t435 - t367 * t572 + t401 - t557) * MDP(19) + (t348 - t644) * MDP(25) + (t347 - t643) * MDP(26); (t598 + t643) * MDP(22) + (-t584 - t644) * MDP(23) + (-t346 * t512 + t668) * MDP(25) + (-t345 * t512 + t674) * MDP(26) + (MDP(22) * t641 - MDP(23) * t663 - MDP(25) * t346 - MDP(26) * t632) * qJD(5) + t673;];
tau = t1;
