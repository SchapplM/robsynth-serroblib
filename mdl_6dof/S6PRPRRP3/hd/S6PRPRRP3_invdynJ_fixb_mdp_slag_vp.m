% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP3
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:48
% EndTime: 2021-01-16 01:39:06
% DurationCPUTime: 9.27s
% Computational Cost: add. (5067->517), mult. (11805->668), div. (0->0), fcn. (9645->14), ass. (0->232)
t673 = 2 * qJD(4);
t505 = sin(pkin(11));
t508 = cos(pkin(11));
t513 = sin(qJ(4));
t651 = cos(qJ(4));
t538 = -t513 * t505 + t651 * t508;
t511 = qJ(3) + pkin(8);
t478 = t511 * t505;
t479 = t511 * t508;
t539 = -t478 * t651 - t513 * t479;
t388 = qJD(3) * t538 + qJD(4) * t539;
t507 = sin(pkin(6));
t516 = cos(qJ(2));
t609 = t507 * t516;
t524 = t538 * t609;
t441 = qJD(1) * t524;
t672 = t388 - t441;
t469 = t538 * qJD(4);
t473 = t505 * t651 + t513 * t508;
t470 = t473 * qJD(4);
t514 = sin(qJ(2));
t590 = qJD(1) * t507;
t572 = t514 * t590;
t671 = pkin(4) * t470 - pkin(9) * t469 - t572;
t475 = qJD(2) * qJ(3) + t572;
t509 = cos(pkin(6));
t589 = qJD(1) * t509;
t488 = t508 * t589;
t640 = pkin(8) * qJD(2);
t433 = t488 + (-t475 - t640) * t505;
t450 = t508 * t475 + t505 * t589;
t434 = t508 * t640 + t450;
t376 = t513 * t433 + t651 * t434;
t670 = t376 * qJD(4);
t669 = t538 * qJD(2);
t591 = t505 ^ 2 + t508 ^ 2;
t668 = MDP(6) * t591;
t468 = t473 * qJD(2);
t515 = cos(qJ(5));
t512 = sin(qJ(5));
t587 = qJD(4) * t512;
t445 = t468 * t515 + t587;
t626 = t445 * t512;
t439 = -t513 * t478 + t479 * t651;
t525 = t473 * t609;
t595 = -qJD(1) * t525 + qJD(3) * t473 + qJD(4) * t439;
t667 = t441 * t512 + t671 * t515;
t526 = t473 * qJDD(2);
t519 = qJD(2) * t469 + t526;
t666 = -(qJD(4) * qJD(5)) - t519;
t494 = t508 * pkin(3) + pkin(2);
t418 = -pkin(4) * t538 - pkin(9) * t473 - t494;
t585 = qJD(5) * t515;
t665 = t418 * t585 + t671 * t512 + t672 * t515;
t639 = cos(pkin(10));
t561 = t639 * t514;
t506 = sin(pkin(10));
t611 = t506 * t516;
t464 = t509 * t561 + t611;
t504 = pkin(11) + qJ(4);
t497 = sin(t504);
t498 = cos(t504);
t565 = t507 * t639;
t423 = t464 * t497 + t498 * t565;
t560 = t639 * t516;
t612 = t506 * t514;
t462 = t509 * t612 - t560;
t615 = t506 * t507;
t425 = t462 * t497 + t498 * t615;
t610 = t507 * t514;
t454 = t497 * t610 - t509 * t498;
t530 = -g(1) * t425 + g(2) * t423 + g(3) * t454;
t664 = g(3) * t507;
t457 = qJD(5) - t669;
t663 = t457 * t626;
t463 = -t509 * t560 + t612;
t465 = t509 * t611 + t561;
t555 = g(1) * t465 + g(2) * t463;
t528 = -g(3) * t609 + t555;
t662 = t528 * t497;
t571 = t516 * t590;
t581 = qJDD(2) * qJ(3);
t583 = qJDD(1) * t507;
t447 = t514 * t583 + t581 + (qJD(3) + t571) * qJD(2);
t582 = qJDD(1) * t509;
t486 = t508 * t582;
t395 = t486 + (-pkin(8) * qJDD(2) - t447) * t505;
t413 = t508 * t447 + t505 * t582;
t580 = qJDD(2) * t508;
t396 = pkin(8) * t580 + t413;
t542 = t513 * t395 + t396 * t651;
t658 = t651 * t433 - t513 * t434;
t346 = qJDD(4) * pkin(9) + qJD(4) * t658 + t542;
t372 = qJD(4) * pkin(9) + t376;
t552 = qJD(3) - t571;
t456 = -qJD(2) * t494 + t552;
t383 = -pkin(4) * t669 - pkin(9) * t468 + t456;
t356 = t372 * t515 + t383 * t512;
t551 = t538 * qJDD(2);
t421 = qJD(2) * t470 - t551;
t567 = t516 * t583;
t588 = qJD(2) * t514;
t568 = qJD(1) * t588;
t655 = t507 * t568 + qJDD(3);
t543 = -t567 + t655;
t636 = qJDD(2) * pkin(2);
t453 = t543 - t636;
t361 = -pkin(3) * t580 + t421 * pkin(4) - pkin(9) * t519 + t453;
t359 = t515 * t361;
t521 = -qJD(5) * t356 - t346 * t512 + t359;
t586 = qJD(5) * t512;
t369 = -t512 * qJDD(4) + t468 * t586 + t666 * t515;
t638 = qJ(6) * t369;
t415 = qJDD(5) + t421;
t650 = pkin(5) * t415;
t337 = -qJD(6) * t445 + t521 + t638 + t650;
t443 = -t515 * qJD(4) + t468 * t512;
t349 = -qJ(6) * t443 + t356;
t635 = t349 * t457;
t661 = t337 + t635;
t534 = -t395 * t651 + t513 * t396 + t670;
t347 = -qJDD(4) * pkin(4) + t534;
t553 = -t515 * qJDD(4) + t468 * t585;
t370 = t512 * t526 + (qJD(5) + t669) * t587 + t553;
t341 = t370 * pkin(5) + qJDD(6) + t347;
t660 = -t341 + t530;
t432 = t515 * t439;
t547 = -qJ(6) * t469 - qJD(6) * t473;
t601 = pkin(5) * t470 - t388 * t512 + t547 * t515 + (-t432 + (qJ(6) * t473 - t418) * t512) * qJD(5) + t667;
t569 = t473 * t585;
t659 = -qJ(6) * t569 + (-qJD(5) * t439 + t547) * t512 + t665;
t537 = t469 * t512 + t569;
t598 = pkin(5) * t537 + t595;
t657 = MDP(20) + MDP(22);
t656 = MDP(21) + MDP(23);
t422 = t462 * t498 - t497 * t615;
t424 = t464 * t498 - t497 * t565;
t455 = t497 * t509 + t498 * t610;
t604 = t515 * t516;
t578 = t507 * t604;
t654 = -g(3) * (-t455 * t512 - t578) - g(2) * (-t424 * t512 + t463 * t515) - g(1) * (t422 * t512 + t465 * t515);
t550 = (-t475 * t505 + t488) * t505 - t450 * t508;
t653 = t516 * t550 - (-qJD(2) * pkin(2) + t552) * t514;
t652 = t445 ^ 2;
t641 = t443 * pkin(5);
t510 = -qJ(6) - pkin(9);
t637 = qJ(6) * t370;
t634 = t369 * t512;
t412 = -t447 * t505 + t486;
t633 = t412 * t505;
t632 = t413 * t508;
t631 = t443 * t457;
t630 = t443 * t669;
t629 = t443 * t468;
t628 = t445 * t457;
t627 = t445 * t468;
t624 = t462 * t512;
t623 = t464 * t512;
t622 = t669 * t512;
t621 = t669 * t515;
t620 = t473 * t512;
t619 = t473 * t515;
t617 = t498 * t512;
t616 = t498 * t515;
t614 = t506 * t509;
t613 = t506 * t511;
t608 = t508 * MDP(5);
t607 = t512 * t415;
t606 = t512 * t516;
t402 = t515 * t415;
t517 = qJD(2) ^ 2;
t603 = t516 * t517;
t602 = qJDD(1) - g(3);
t355 = -t372 * t512 + t515 * t383;
t348 = -qJ(6) * t445 + t355;
t344 = pkin(5) * t457 + t348;
t599 = -t348 + t344;
t597 = -t512 * t370 - t443 * t585;
t416 = pkin(4) * t468 - pkin(9) * t669;
t596 = t512 * t416 + t515 * t658;
t594 = t512 * t418 + t432;
t566 = qJD(5) * t510;
t593 = qJ(6) * t622 + qJD(6) * t515 + t512 * t566 - t596;
t405 = t515 * t416;
t592 = -pkin(5) * t468 + qJ(6) * t621 + t515 * t566 - t405 + (-qJD(6) + t658) * t512;
t579 = t507 * t606;
t576 = t530 * t512;
t573 = t639 * pkin(2);
t570 = t507 * t588;
t564 = t639 * qJ(3);
t563 = t639 * t494;
t562 = t639 * t511;
t559 = qJD(6) + t641;
t558 = t457 * t515;
t557 = t515 * t346 + t512 * t361 - t372 * t586 + t383 * t585;
t556 = -g(1) * t462 + g(2) * t464;
t496 = pkin(5) * t515 + pkin(4);
t549 = -t496 * t498 + t497 * t510;
t338 = -qJD(6) * t443 + t557 - t637;
t548 = -t344 * t457 + t338;
t545 = t402 + (-t586 + t622) * t457;
t460 = -t505 * t610 + t508 * t509;
t461 = t505 * t509 + t508 * t610;
t400 = t513 * t460 + t461 * t651;
t384 = -t400 * t512 - t578;
t544 = -t400 * t515 + t579;
t540 = t460 * t651 - t513 * t461;
t536 = t469 * t515 - t473 * t586;
t371 = -qJD(4) * pkin(4) - t658;
t535 = -pkin(9) * t415 + t371 * t457;
t533 = -g(1) * (-t462 * t515 + t465 * t617) - g(2) * (t463 * t617 + t464 * t515) - (-t498 * t606 + t514 * t515) * t664;
t532 = -g(1) * (-t465 * t616 - t624) - g(2) * (-t463 * t616 + t623) - (t498 * t604 + t512 * t514) * t664;
t531 = -g(1) * t422 + g(2) * t424 + g(3) * t455;
t523 = t528 + t567;
t522 = -g(1) * (t422 * t515 - t465 * t512) - g(2) * (-t424 * t515 - t463 * t512) - g(3) * (-t455 * t515 + t579) - t557;
t520 = t521 + t654;
t518 = t666 * t512 - t553;
t481 = t510 * t515;
t480 = t510 * t512;
t442 = t443 ^ 2;
t436 = -qJDD(2) * t494 + t543;
t408 = t515 * t418;
t390 = pkin(5) * t620 - t539;
t374 = qJD(2) * t525 + qJD(4) * t400;
t373 = qJD(2) * t524 + qJD(4) * t540;
t364 = pkin(5) * t622 + t376;
t363 = -qJ(6) * t620 + t594;
t362 = t371 + t559;
t360 = -pkin(5) * t538 - qJ(6) * t619 - t439 * t512 + t408;
t354 = qJD(5) * t544 - t373 * t512 + t515 * t570;
t353 = qJD(5) * t384 + t373 * t515 + t512 * t570;
t1 = [t602 * MDP(1) + (-t460 * t505 + t461 * t508) * qJDD(2) * MDP(6) + (t412 * t460 + t413 * t461 - g(3)) * MDP(7) + (-qJD(4) * t374 + qJDD(4) * t540) * MDP(13) + (-t373 * qJD(4) - t400 * qJDD(4)) * MDP(14) + (-t353 * t443 - t354 * t445 + t369 * t384 + t370 * t544) * MDP(24) + (t337 * t384 - t338 * t544 - t341 * t540 + t344 * t354 + t349 * t353 + t362 * t374 - g(3)) * MDP(25) + t657 * (t354 * t457 - t370 * t540 + t374 * t443 + t384 * t415) + t656 * (-t353 * t457 + t369 * t540 + t374 * t445 + t415 * t544) + ((-qJDD(2) * t514 - t603) * MDP(4) + t603 * t668 + (-qJD(2) * t653 - t453 * t516) * MDP(7) + (-t421 * t516 - t588 * t669) * MDP(13) + (-t516 * t526 + (t514 * t468 - t469 * t516) * qJD(2)) * MDP(14) + (MDP(3) + t608) * (qJDD(2) * t516 - t514 * t517)) * t507; qJDD(2) * MDP(2) + t523 * MDP(3) + (-t602 * t610 + t556) * MDP(4) + (t636 - t453 + (-g(3) * t516 + t568) * t507 + t555) * t608 + (-g(3) * t610 - t556 + t632 - t633 + (qJD(2) * t552 + t581) * t591) * MDP(6) + (qJ(3) * t632 - qJ(3) * t633 - t453 * pkin(2) - g(1) * (-(qJ(3) * t614 + t573) * t514 + (-pkin(2) * t614 + t564) * t516) - g(2) * (-(t506 * pkin(2) - t509 * t564) * t514 + (t506 * qJ(3) + t509 * t573) * t516) - t550 * qJD(3) + (-g(3) * (pkin(2) * t516 + qJ(3) * t514) + t653 * qJD(1)) * t507) * MDP(7) + (t468 * t469 + t473 * t519) * MDP(8) + (-t473 * t421 - t468 * t470 + t469 * t669 + t519 * t538) * MDP(9) + (qJD(4) * t469 + qJDD(4) * t473) * MDP(10) + (-qJD(4) * t470 + qJDD(4) * t538) * MDP(11) + (-qJD(4) * t595 + qJDD(4) * t539 - t421 * t494 - t436 * t538 + t456 * t470 + t498 * t528 + t572 * t669) * MDP(13) + (-t672 * qJD(4) - t439 * qJDD(4) + t436 * t473 + t456 * t469 - t468 * t572 - t494 * t519 - t662) * MDP(14) + (-t369 * t619 + t445 * t536) * MDP(15) + ((-t443 * t515 - t626) * t469 + (t634 - t370 * t515 + (t443 * t512 - t445 * t515) * qJD(5)) * t473) * MDP(16) + (t369 * t538 + t402 * t473 + t445 * t470 + t457 * t536) * MDP(17) + (t370 * t538 - t443 * t470 - t457 * t537 - t473 * t607) * MDP(18) + (-t415 * t538 + t457 * t470) * MDP(19) + (t408 * t415 - (-t372 * t585 + t359) * t538 + t355 * t470 - t539 * t370 + t371 * t569 + (-t439 * t585 + t667) * t457 + t595 * t443 + ((-qJD(5) * t418 - t388) * t457 - t439 * t415 - (-qJD(5) * t383 - t346) * t538 + t347 * t473 + t371 * t469) * t512 + t532) * MDP(20) + (-t594 * t415 + t557 * t538 - t356 * t470 + t539 * t369 + t347 * t619 + (t439 * t586 - t665) * t457 + t595 * t445 + t536 * t371 + t533) * MDP(21) + (-t337 * t538 + t341 * t620 + t344 * t470 + t360 * t415 + t362 * t537 + t370 * t390 + t443 * t598 + t457 * t601 + t532) * MDP(22) + (t338 * t538 + t341 * t619 - t349 * t470 + t362 * t536 - t363 * t415 - t369 * t390 + t445 * t598 - t457 * t659 + t533) * MDP(23) + (t360 * t369 - t363 * t370 + (-t344 * t515 - t349 * t512) * t469 - t601 * t445 - t659 * t443 + t662 + (-t337 * t515 - t338 * t512 + (t344 * t512 - t349 * t515) * qJD(5)) * t473) * MDP(24) + (t338 * t363 + t337 * t360 + t341 * t390 - g(1) * (-pkin(5) * t624 - (t509 * t613 + t563) * t514 + (-t494 * t614 + t562) * t516 + t549 * t465) - g(2) * (pkin(5) * t623 - (t506 * t494 - t509 * t562) * t514 + (t509 * t563 + t613) * t516 + t549 * t463) - ((pkin(5) * t512 + t511) * t514 + (t494 - t549) * t516) * t664 + t598 * t362 + t659 * t349 + t601 * t344) * MDP(25); -MDP(5) * t580 - t517 * t668 + (qJD(2) * t550 - t523 - t636 + t655) * MDP(7) + (t468 * t673 - t551) * MDP(13) + (t669 * t673 + t526) * MDP(14) + ((t369 + t630) * t515 + t663 + t597) * MDP(24) + (-t362 * t468 + t548 * t512 + t515 * t661 - t528) * MDP(25) + t656 * (-t457 ^ 2 * t515 - t607 - t627) + t657 * (t545 - t629); -t669 ^ 2 * MDP(9) + t526 * MDP(10) + t551 * MDP(11) + qJDD(4) * MDP(12) + (t530 - t534 + t670) * MDP(13) + (-t456 * t669 + t531 - t542) * MDP(14) + (t445 * t558 - t634) * MDP(15) + ((-t369 + t630) * t515 - t663 + t597) * MDP(16) + (t457 * t558 + t607 - t627) * MDP(17) + (t545 + t629) * MDP(18) + (-pkin(4) * t370 - t376 * t443 - t405 * t457 + (t457 * t658 + t535) * t512 + (-pkin(9) * qJD(5) * t457 - t347 + t530) * t515) * MDP(20) + (pkin(4) * t369 + t347 * t512 - t376 * t445 + (pkin(9) * t586 + t596) * t457 + t535 * t515 - t576) * MDP(21) + (-t364 * t443 - t370 * t496 + t415 * t480 + t592 * t457 + (-t362 * t669 + (t362 + t641) * qJD(5)) * t512 + t660 * t515) * MDP(22) + (-t362 * t621 + t341 * t512 - t364 * t445 + t369 * t496 + t415 * t481 - t593 * t457 + (pkin(5) * t626 + t362 * t515) * qJD(5) - t576) * MDP(23) + (t369 * t480 + t370 * t481 - t593 * t443 - t592 * t445 - t512 * t661 + t548 * t515 - t531) * MDP(24) + (-t338 * t481 + t337 * t480 - t341 * t496 - g(1) * (t422 * t510 + t425 * t496) - g(2) * (-t423 * t496 - t424 * t510) - g(3) * (-t454 * t496 - t455 * t510) + (pkin(5) * t586 - t364) * t362 + t593 * t349 + t592 * t344) * MDP(25) + (-MDP(13) * t456 - MDP(19) * t457 - MDP(20) * t355 + MDP(21) * t356 - MDP(22) * t344 + MDP(23) * t349 - MDP(8) * t669 + MDP(9) * t468) * t468; t445 * t443 * MDP(15) + (-t442 + t652) * MDP(16) + (-t369 + t631) * MDP(17) + (t518 + t628) * MDP(18) + t415 * MDP(19) + (t356 * t457 - t371 * t445 + t520) * MDP(20) + (t355 * t457 + t371 * t443 + t522) * MDP(21) + (0.2e1 * t650 + t638 + t635 + (-t362 - t559) * t445 + t520) * MDP(22) + (-pkin(5) * t652 + t637 + t348 * t457 + (qJD(6) + t362) * t443 + t522) * MDP(23) + (pkin(5) * t369 - t443 * t599) * MDP(24) + (t599 * t349 + (-t362 * t445 + t337 + t654) * pkin(5)) * MDP(25); (-t518 + t628) * MDP(22) + (-t369 - t631) * MDP(23) + (-t442 - t652) * MDP(24) + (t344 * t445 + t349 * t443 - t660) * MDP(25);];
tau = t1;
