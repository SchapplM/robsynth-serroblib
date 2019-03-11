% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:12
% EndTime: 2019-03-09 05:09:31
% DurationCPUTime: 11.61s
% Computational Cost: add. (41735->587), mult. (81880->791), div. (0->0), fcn. (100209->10), ass. (0->332)
t384 = cos(pkin(11));
t564 = t384 * pkin(5);
t373 = -pkin(4) - t564;
t383 = sin(pkin(10));
t385 = cos(pkin(10));
t568 = sin(qJ(3));
t570 = cos(qJ(3));
t360 = -t383 * t568 + t385 * t570;
t362 = t383 * t570 + t385 * t568;
t387 = sin(qJ(4));
t569 = cos(qJ(4));
t325 = -t569 * t360 + t362 * t387;
t382 = sin(pkin(11));
t386 = sin(qJ(6));
t388 = cos(qJ(6));
t361 = t382 * t388 + t386 * t384;
t638 = t361 * t325;
t665 = t638 * mrSges(7,1);
t426 = t386 * t382 - t384 * t388;
t635 = t426 * t325;
t666 = t635 * mrSges(7,2);
t674 = t666 - t665;
t689 = t373 * t674;
t600 = -m(7) / 0.2e1;
t611 = t665 / 0.2e1 - t666 / 0.2e1;
t602 = -m(6) / 0.2e1;
t563 = pkin(7) + qJ(2);
t442 = t563 * t385;
t443 = t563 * t383;
t338 = -t442 * t568 - t443 * t570;
t291 = -t362 * pkin(8) + t338;
t340 = t442 * t570 - t443 * t568;
t292 = t360 * pkin(8) + t340;
t609 = t291 * t387 + t569 * t292;
t613 = t609 * t602;
t637 = t382 * t325;
t655 = -pkin(5) * t637 + t609;
t408 = t655 * t600 + t611 + t613;
t528 = t382 * mrSges(6,1);
t455 = t528 / 0.2e1;
t526 = t384 * mrSges(6,2);
t688 = t408 + t325 * (t526 / 0.2e1 + t455);
t636 = t384 * t325;
t448 = -t636 / 0.2e1;
t687 = mrSges(6,2) * t448 - t408;
t418 = t387 * t360 + t362 * t569;
t229 = t426 * t418;
t415 = t361 * t418;
t115 = mrSges(7,1) * t415 - mrSges(7,2) * t229;
t198 = t569 * t291 - t387 * t292;
t498 = t418 * t382;
t143 = pkin(5) * t498 - t198;
t434 = t526 + t528;
t247 = t434 * t418;
t497 = t418 * t384;
t454 = -pkin(2) * t385 - pkin(1);
t343 = -pkin(3) * t360 + t454;
t195 = pkin(4) * t325 - qJ(5) * t418 + t343;
t93 = t384 * t195 - t382 * t609;
t67 = pkin(5) * t325 - pkin(9) * t497 + t93;
t94 = t382 * t195 + t384 * t609;
t71 = -pkin(9) * t498 + t94;
t40 = -t386 * t71 + t388 * t67;
t41 = t386 * t67 + t388 * t71;
t590 = -t229 / 0.2e1;
t592 = -t415 / 0.2e1;
t634 = t434 * t325;
t645 = Ifges(7,4) * t635 + Ifges(7,2) * t638 + Ifges(7,6) * t418;
t646 = Ifges(7,1) * t635 + Ifges(7,4) * t638 + Ifges(7,5) * t418;
t650 = mrSges(7,1) * t418 - mrSges(7,3) * t635;
t651 = -mrSges(7,2) * t418 + mrSges(7,3) * t638;
t652 = mrSges(6,1) * t418 + mrSges(6,3) * t636;
t653 = -mrSges(6,2) * t418 + mrSges(6,3) * t637;
t686 = t655 * t115 + t143 * t674 + t198 * t634 + t609 * t247 + t40 * t650 + t41 * t651 + t646 * t590 + t645 * t592 + t93 * t652 + t94 * t653;
t365 = -mrSges(6,1) * t384 + mrSges(6,2) * t382;
t631 = t609 * t365;
t640 = t609 * mrSges(5,1);
t556 = Ifges(7,4) * t361;
t334 = -Ifges(7,2) * t426 + t556;
t659 = t638 * t334;
t356 = Ifges(7,4) * t426;
t336 = Ifges(7,1) * t361 - t356;
t660 = t635 * t336;
t667 = t198 * mrSges(5,2);
t331 = mrSges(7,1) * t426 + mrSges(7,2) * t361;
t676 = t655 * t331;
t677 = t426 * t645;
t680 = t361 * t646;
t685 = t631 - t640 + t659 / 0.2e1 + t660 / 0.2e1 - t667 + t676 - t677 / 0.2e1 + t680 / 0.2e1;
t558 = Ifges(6,4) * t384;
t369 = Ifges(6,1) * t382 + t558;
t484 = t384 * t369;
t432 = -Ifges(6,2) * t382 + t558;
t630 = Ifges(6,6) * t418 - t325 * t432;
t489 = t384 * t630;
t559 = Ifges(6,4) * t382;
t490 = t382 * (Ifges(6,2) * t384 + t559);
t433 = Ifges(6,1) * t384 - t559;
t629 = Ifges(6,5) * t418 - t325 * t433;
t495 = t382 * t629;
t684 = -t631 / 0.2e1 - (-t484 / 0.4e1 + t490 / 0.4e1) * t325 + t640 / 0.2e1 - t489 / 0.4e1 - t495 / 0.4e1 - t659 / 0.4e1 - t660 / 0.4e1 + t667 / 0.2e1 - t676 / 0.2e1 + t677 / 0.4e1 - t680 / 0.4e1;
t585 = -t325 / 0.2e1;
t683 = t653 / 0.2e1;
t473 = t382 ^ 2 + t384 ^ 2;
t642 = mrSges(6,3) * t473;
t681 = t143 * t655;
t467 = t569 * pkin(3);
t374 = -t467 - pkin(4);
t363 = t374 - t564;
t679 = t363 * t655;
t678 = t373 * t655;
t551 = Ifges(6,6) * t382;
t554 = Ifges(6,5) * t384;
t431 = -t551 + t554;
t571 = t384 / 0.2e1;
t572 = -t382 / 0.2e1;
t624 = t418 / 0.2e1;
t625 = -t418 / 0.2e1;
t643 = t325 / 0.2e1;
t675 = 0.2e1 * Ifges(5,4) * t625 + t629 * t571 + t630 * t572 + Ifges(7,5) * t590 + Ifges(7,6) * t592 + t431 * t624 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3)) * t643;
t670 = t362 ^ 2;
t669 = -t634 / 0.2e1;
t567 = pkin(3) * t362;
t668 = pkin(9) * t637;
t664 = t198 * t382;
t663 = t198 * t387;
t511 = t198 * t609;
t575 = -t361 / 0.2e1;
t662 = t229 * t575;
t661 = t384 * t198;
t658 = t331 + t365;
t657 = -Ifges(5,1) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,3) / 0.2e1;
t654 = pkin(5) * t418 + pkin(9) * t636;
t648 = t585 * t642 + t658 * t624;
t641 = Ifges(7,3) * t625;
t565 = pkin(4) * t418;
t257 = qJ(5) * t325 + t565;
t557 = Ifges(7,4) * t229;
t85 = -Ifges(7,2) * t415 + t325 * Ifges(7,6) - t557;
t627 = -t85 / 0.2e1;
t612 = t115 + t247;
t566 = pkin(3) * t387;
t372 = qJ(5) + t566;
t349 = (-pkin(9) - t372) * t382;
t377 = t384 * pkin(9);
t496 = t372 * t384;
t350 = t377 + t496;
t309 = t349 * t388 - t386 * t350;
t310 = t386 * t349 + t350 * t388;
t476 = -t309 * t361 - t310 * t426;
t364 = (-pkin(9) - qJ(5)) * t382;
t366 = qJ(5) * t384 + t377;
t337 = t364 * t388 - t386 * t366;
t339 = t386 * t364 + t366 * t388;
t475 = -t337 * t361 - t339 * t426;
t253 = -mrSges(6,2) * t325 - mrSges(6,3) * t498;
t487 = t384 * t253;
t256 = t325 * mrSges(6,1) - mrSges(6,3) * t497;
t493 = t382 * t256;
t607 = -t493 / 0.2e1 + t487 / 0.2e1;
t430 = -t382 * t93 + t384 * t94;
t606 = -Ifges(5,4) * t325 + t343 * mrSges(5,2) + Ifges(5,1) * t624 + (Ifges(6,6) * t325 + t418 * t432) * t572 + (Ifges(6,5) * t325 + t418 * t433) * t571;
t149 = -mrSges(7,2) * t325 - mrSges(7,3) * t415;
t152 = mrSges(7,1) * t325 + mrSges(7,3) * t229;
t576 = -t426 / 0.2e1;
t419 = t149 * t576 + t152 * t575;
t577 = t426 / 0.2e1;
t401 = (-t415 * t577 - t662) * mrSges(7,3) + t419;
t18 = t401 - t611;
t604 = t18 * qJD(1);
t330 = t361 * mrSges(7,1) - mrSges(7,2) * t426;
t60 = 0.2e1 * mrSges(7,1) * t590 + 0.2e1 * t592 * mrSges(7,2);
t603 = qJD(1) * t60 + (qJD(3) + qJD(4)) * t330;
t601 = m(6) / 0.2e1;
t599 = m(7) / 0.2e1;
t598 = m(5) * pkin(3);
t597 = m(6) * pkin(3);
t596 = -mrSges(7,1) / 0.2e1;
t595 = mrSges(7,2) / 0.2e1;
t225 = Ifges(7,4) * t415;
t88 = -Ifges(7,1) * t229 + t325 * Ifges(7,5) - t225;
t594 = t88 / 0.2e1;
t593 = t638 / 0.2e1;
t591 = t635 / 0.2e1;
t589 = t309 / 0.2e1;
t588 = t310 / 0.2e1;
t580 = -t337 / 0.2e1;
t579 = -t339 / 0.2e1;
t452 = t388 * t569;
t453 = t386 * t569;
t345 = (-t382 * t453 + t384 * t452) * pkin(3);
t578 = t345 / 0.2e1;
t574 = t361 / 0.2e1;
t573 = t363 / 0.2e1;
t562 = -t361 * t40 - t41 * t426;
t561 = m(7) * qJD(2);
t555 = Ifges(5,5) * t325;
t553 = Ifges(7,5) * t635;
t552 = Ifges(5,6) * t418;
t550 = Ifges(7,6) * t638;
t321 = t418 * mrSges(5,1);
t353 = t362 * mrSges(4,1);
t202 = t257 + t567;
t96 = t384 * t202 - t664;
t69 = t654 + t96;
t97 = t382 * t202 + t661;
t72 = t97 + t668;
t44 = -t386 * t72 + t388 * t69;
t45 = t386 * t69 + t388 * t72;
t3 = (m(5) * t567 + t321) * t343 + m(6) * (t93 * t96 + t94 * t97 - t511) + t454 * t353 + (t454 * mrSges(4,2) + Ifges(4,4) * t360 + (Ifges(4,1) - Ifges(4,2)) * t362) * t360 + t96 * t256 + t97 * t253 + t44 * t152 + t45 * t149 - t638 * t627 + t635 * t594 + (t550 + t553) * t643 + m(7) * (t40 * t44 + t41 * t45 + t681) + (mrSges(5,2) * t567 + t675) * t418 + (mrSges(5,1) * t567 + t657 * t418 - t431 * t643 - t606) * t325 - t670 * Ifges(4,4) + t686;
t539 = t3 * qJD(1);
t538 = t310 * mrSges(7,3);
t537 = t325 * mrSges(5,2);
t535 = t339 * mrSges(7,3);
t534 = t426 * mrSges(7,3);
t531 = t361 * mrSges(7,3);
t108 = t384 * t257 - t664;
t109 = t382 * t257 + t661;
t423 = Ifges(7,5) * t591 + Ifges(7,6) * t593;
t70 = t108 + t654;
t80 = t109 + t668;
t52 = -t386 * t80 + t388 * t70;
t53 = t386 * t70 + t388 * t80;
t4 = ((-t554 / 0.2e1 + t551 / 0.2e1) * t325 + t423 - t606) * t325 + t88 * t591 + t85 * t593 + m(6) * (t108 * t93 + t109 * t94 - t511) + t108 * t256 + t109 * t253 + t52 * t152 + t53 * t149 + m(7) * (t40 * t52 + t41 * t53 + t681) + (t343 * mrSges(5,1) + t325 * t657 + t675) * t418 + t686;
t524 = t4 * qJD(1);
t112 = -t229 * mrSges(7,1) - t415 * mrSges(7,2);
t116 = Ifges(7,2) * t229 - t225;
t117 = -Ifges(7,1) * t415 + t557;
t478 = -Ifges(7,5) * t415 + Ifges(7,6) * t229;
t9 = t143 * t112 + t478 * t643 + t40 * t149 - t41 * t152 - (t627 + t117 / 0.2e1 - t41 * mrSges(7,3)) * t229 - (t116 / 0.2e1 + t594 - t40 * mrSges(7,3)) * t415;
t523 = t9 * qJD(1);
t522 = t96 * t382;
t521 = t97 * t384;
t470 = t598 / 0.2e1;
t396 = (t382 * t97 + t384 * t96) * t601 + (t361 * t45 - t426 * t44) * t599 + t650 * t576 + t651 * t574 + t382 * t683 + t652 * t571 + t362 * t470;
t440 = t473 * t372;
t398 = (-t325 * t440 + t374 * t418) * t601 + (t309 * t638 + t310 * t635 + t363 * t418) * t599 + (-t325 * t387 - t418 * t569) * t470;
t456 = t531 / 0.2e1;
t459 = t534 / 0.2e1;
t11 = t360 * mrSges(4,2) + t638 * t456 + t635 * t459 + t321 + t353 + t396 - t398 - t537 - t648;
t520 = qJD(1) * t11;
t22 = m(7) * (t229 * t40 - t41 * t415) - t415 * t149 + t229 * t152 + (-t382 * t253 - t384 * t256 + m(6) * (-t382 * t94 - t384 * t93)) * t418;
t519 = qJD(1) * t22;
t509 = t198 * t418;
t10 = t635 * t149 + t638 * t152 + (t360 ^ 2 + t670) * mrSges(4,3) - (-mrSges(5,3) * t325 + t487 - t493) * t325 + (mrSges(5,3) * t418 + t612) * t418 + m(7) * (t143 * t418 + t40 * t638 + t41 * t635) + m(6) * (-t325 * t430 - t509) + m(5) * (-t325 * t609 - t509) + m(4) * (-t338 * t362 + t340 * t360) + (m(3) * qJ(2) + mrSges(3,3)) * (t383 ^ 2 + t385 ^ 2);
t517 = t10 * qJD(1);
t516 = t108 * t382;
t515 = t109 * t384;
t332 = Ifges(7,5) * t361 - Ifges(7,6) * t426;
t501 = t418 * t332;
t367 = Ifges(6,5) * t382 + Ifges(6,6) * t384;
t499 = t418 * t367;
t494 = t382 * t652;
t488 = t384 * t653;
t407 = (-t229 * t426 - t361 * t415) * t599 + m(6) * t473 * t625;
t59 = 0.2e1 * (-m(6) / 0.4e1 - m(7) / 0.4e1) * t418 + t407;
t480 = t59 * qJD(1);
t474 = -Ifges(7,5) * t426 - Ifges(7,6) * t361;
t328 = t330 * qJD(6);
t469 = mrSges(6,3) * t522;
t468 = mrSges(6,3) * t521;
t344 = (-t382 * t452 - t384 * t453) * pkin(3);
t260 = -t344 * t426 + t345 * t361;
t466 = t260 * t599;
t464 = qJD(4) * t466;
t458 = -t534 / 0.2e1;
t457 = -t531 / 0.2e1;
t333 = -Ifges(7,2) * t361 - t356;
t446 = t333 / 0.4e1 + t336 / 0.4e1;
t335 = -Ifges(7,1) * t426 - t556;
t445 = -t334 / 0.4e1 + t335 / 0.4e1;
t444 = t334 / 0.2e1 - t335 / 0.2e1;
t441 = t473 * qJ(5);
t436 = (-t333 / 0.2e1 - t336 / 0.2e1) * t426;
t435 = (t599 + t601) * t566;
t429 = t521 - t522;
t399 = (-t415 * t576 + t662) * mrSges(7,3) + t430 * t601 + t419 + t607;
t394 = (t229 * t309 - t310 * t415 + t562) * t599 + t399;
t13 = t394 + t688;
t425 = (t361 ^ 2 + t426 ^ 2) * mrSges(7,3) + t642;
t98 = m(6) * t440 + m(7) * t476 + t425;
t428 = qJD(1) * t13 + qJD(3) * t98;
t427 = t515 - t516;
t422 = mrSges(7,2) * t578 + t344 * t596;
t420 = t638 * t457 + t635 * t458 + t648;
t397 = (t108 * t384 + t109 * t382) * t602 + (t361 * t53 - t426 * t52) * t600 + mrSges(5,1) * t625 + t650 * t577 + t651 * t575 + t653 * t572 - t384 * t652 / 0.2e1;
t400 = -t321 / 0.2e1 + (-t325 * t441 - t565) * t601 + (t337 * t638 + t339 * t635 + t373 * t418) * t599 + t420;
t19 = 0.2e1 * mrSges(5,2) * t643 + t397 + t400;
t416 = -t19 * qJD(1) + qJD(3) * t466;
t414 = t473 * t569;
t126 = m(6) * t441 + m(7) * t475 + t425;
t393 = (t229 * t337 - t339 * t415 + t562) * t599 + t399;
t15 = t393 + t688;
t402 = (t440 + t441) * t602 + (t475 + t476) * t600 - t425;
t63 = t435 + t402;
t413 = qJD(1) * t15 - qJD(3) * t63 + qJD(4) * t126;
t271 = t310 * t531;
t293 = t363 * t330;
t46 = t271 - t293 - t436 + (t444 - t538) * t361;
t403 = -(t116 / 0.4e1 + t88 / 0.4e1) * t426 + (-t85 / 0.4e1 + t117 / 0.4e1) * t361 + t143 * t330 / 0.2e1 + t325 * t474 / 0.4e1;
t391 = -(-t309 * mrSges(7,3) / 0.2e1 + t446) * t415 - (-t538 / 0.2e1 + t445) * t229 + t149 * t589 - t310 * t152 / 0.2e1 + t112 * t573 + t403;
t405 = -t553 / 0.2e1 - t550 / 0.2e1 + t641 + t44 * t596 + t45 * t595;
t6 = t391 + t405;
t412 = -t6 * qJD(1) + t46 * qJD(3);
t389 = -t684 + (t515 / 0.2e1 - t516 / 0.2e1) * mrSges(6,3) + t52 * t457 + t53 * t458 + t607 * t467 + t612 * t566 / 0.2e1 + t496 * t683 + Ifges(5,6) * t625 + t374 * t669 + t372 * t652 * t572 + t674 * t573 + (t367 + t332) * t418 / 0.4e1 + t344 * t152 / 0.2e1 + (t143 * t566 + t309 * t52 + t310 * t53 + t344 * t40 + t345 * t41 + t679) * t599 + t149 * t578 + Ifges(5,5) * t585 + t651 * t588 + t650 * t589 + (t374 * t609 + t427 * t372 + (t430 * t569 - t663) * pkin(3)) * t601;
t390 = t555 / 0.2e1 + t44 * t456 + t45 * t459 + t552 / 0.2e1 + (t494 / 0.2e1 - t488 / 0.2e1 + t429 * t602) * qJ(5) + t469 / 0.2e1 - t468 / 0.2e1 - t689 / 0.2e1 - t499 / 0.4e1 - t501 / 0.4e1 + (t337 * t44 + t339 * t45 + t678) * t600 + t651 * t579 + t650 * t580 + (-t613 + t669) * pkin(4) + t684;
t2 = t389 + t390;
t395 = -t344 * t531 - t345 * t534 + (-mrSges(5,1) + t658) * t566 + (-mrSges(5,2) + t642) * t467;
t62 = m(7) * (t309 * t344 + t310 * t345 + t363 * t566) + (t372 * t414 + t374 * t387) * t597 + t395;
t411 = -t2 * qJD(1) - t62 * qJD(3) - t260 * t561 / 0.2e1;
t290 = t339 * t531;
t308 = t373 * t330;
t404 = t335 * t574 + t334 * t575 + (t588 + t339 / 0.2e1) * t531 - t271 / 0.2e1 - t290 / 0.2e1 + t293 / 0.2e1 + t308 / 0.2e1 + (t333 + t336) * t576;
t33 = t404 + t422;
t54 = t290 - t308 - t436 + (t444 - t535) * t361;
t392 = -(-t535 / 0.2e1 + t445) * t229 - (mrSges(7,3) * t580 + t446) * t415 + t337 * t149 / 0.2e1 + t152 * t579 + t373 * t112 / 0.2e1 + t403;
t406 = t52 * t596 + t53 * t595 - t423 + t641;
t8 = t392 + t406;
t409 = -t8 * qJD(1) - t33 * qJD(3) + t54 * qJD(4);
t329 = t330 * qJD(5);
t64 = t435 - t402;
t58 = t407 + (m(6) + m(7)) * t624;
t32 = t404 - t422;
t20 = -t397 + mrSges(5,2) * t585 + t537 / 0.2e1 + t400;
t17 = t401 + t611;
t16 = t396 + t398 + t420;
t14 = -mrSges(6,1) * t637 / 0.2e1 + t393 + t687;
t12 = -t325 * t455 + t394 + t687;
t7 = t392 - t406;
t5 = t391 - t405;
t1 = t389 - t390;
t21 = [qJD(2) * t10 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t22 + qJD(6) * t9, t517 + (t361 * t635 - t426 * t638) * t561 + t16 * qJD(3) + t20 * qJD(4) + t58 * qJD(5) + t17 * qJD(6), t539 + t16 * qJD(2) + (t685 - t555 + t495 / 0.2e1 + t489 / 0.2e1 - t552 - t44 * t531 + (t325 * t467 - t418 * t566) * mrSges(5,3) - t469 + t363 * t674 + t468 - Ifges(4,6) * t362 + Ifges(4,5) * t360 - t338 * mrSges(4,2) - t340 * mrSges(4,1) + (m(6) * t609 - t634) * t374 + t490 * t643 + t369 * t448 + t499 / 0.2e1 + t501 / 0.2e1 + m(7) * (t309 * t44 + t310 * t45 + t679) - t45 * t534 + t310 * t651 + t309 * t650 + (m(6) * t429 + t488 - t494) * t372 + (-t569 * t609 + t663) * t598) * qJD(3) + t1 * qJD(4) + t12 * qJD(5) + t5 * qJD(6), t20 * qJD(2) + t1 * qJD(3) + t14 * qJD(5) + t7 * qJD(6) + t524 + (t689 + t337 * t650 + t339 * t651 + pkin(4) * t634 - t52 * t531 - t53 * t534 + (t630 / 0.2e1 + qJ(5) * t653 + t109 * mrSges(6,3)) * t384 + (t629 / 0.2e1 - t108 * mrSges(6,3) - qJ(5) * t652) * t382 + 0.2e1 * (-pkin(4) * t609 + qJ(5) * t427) * t601 + 0.2e1 * (t337 * t52 + t339 * t53 + t678) * t599 + (t367 / 0.2e1 + t332 / 0.2e1 - Ifges(5,6)) * t418 + (-Ifges(5,5) - t484 / 0.2e1 + t490 / 0.2e1) * t325 + t685) * qJD(4), qJD(2) * t58 + qJD(3) * t12 + qJD(4) * t14 + t519, t523 + t17 * qJD(2) + t5 * qJD(3) + t7 * qJD(4) + (-mrSges(7,1) * t41 - mrSges(7,2) * t40 + t478) * qJD(6); qJD(3) * t11 - qJD(4) * t19 + qJD(5) * t59 + qJD(6) * t18 - t517, 0, t464 + t520, t416, t480, -t328 + t604; -qJD(2) * t11 + qJD(4) * t2 + qJD(5) * t13 + qJD(6) * t6 - t539, t464 - t520, qJD(4) * t62 + qJD(5) * t98 - qJD(6) * t46 (m(7) * (t337 * t344 + t339 * t345 + t373 * t566) + (-pkin(4) * t387 + qJ(5) * t414) * t597 + t395) * qJD(4) + t64 * qJD(5) + t32 * qJD(6) - t411, qJD(4) * t64 + t428, t32 * qJD(4) + (-mrSges(7,1) * t310 - mrSges(7,2) * t309 + t474) * qJD(6) - t412; qJD(2) * t19 - qJD(3) * t2 + qJD(5) * t15 + qJD(6) * t8 - t524, -t416, -t63 * qJD(5) + t33 * qJD(6) + t411, qJD(5) * t126 - qJD(6) * t54, t413 (-mrSges(7,1) * t339 - mrSges(7,2) * t337 + t474) * qJD(6) - t409; -qJD(2) * t59 - qJD(3) * t13 - qJD(4) * t15 + qJD(6) * t60 - t519, -t480, qJD(4) * t63 + t328 - t428, -t413 + t328, 0, t603; -qJD(2) * t18 - qJD(3) * t6 - qJD(4) * t8 - qJD(5) * t60 - t523, -t604, -t33 * qJD(4) - t329 + t412, -t329 + t409, -t603, 0;];
Cq  = t21;
