% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:41
% EndTime: 2019-03-09 16:44:05
% DurationCPUTime: 12.67s
% Computational Cost: add. (17725->591), mult. (34156->763), div. (0->0), fcn. (35173->6), ass. (0->346)
t370 = sin(qJ(3));
t371 = sin(qJ(2));
t624 = cos(qJ(3));
t625 = cos(qJ(2));
t329 = t370 * t371 - t624 * t625;
t330 = t370 * t625 + t371 * t624;
t372 = cos(qJ(5));
t368 = t372 ^ 2;
t369 = sin(qJ(5));
t613 = Ifges(6,4) - Ifges(7,5);
t495 = t613 * t372;
t614 = -Ifges(4,4) - Ifges(5,6);
t683 = Ifges(6,3) + Ifges(7,2);
t612 = Ifges(6,6) - Ifges(7,6);
t684 = Ifges(7,4) + Ifges(6,5);
t719 = t369 * t684;
t714 = -t612 * t372 - t719;
t731 = (t614 - t714) * t330 + (Ifges(5,3) + Ifges(4,2) - Ifges(5,2) - Ifges(4,1) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t368 + ((Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t369 + t495) * t369 - t683) * t329;
t729 = Ifges(6,1) + Ifges(7,1);
t728 = Ifges(6,2) + Ifges(7,3);
t629 = -t369 / 0.2e1;
t602 = Ifges(7,5) * t372;
t470 = Ifges(7,1) * t369 - t602;
t606 = Ifges(6,4) * t372;
t472 = Ifges(6,1) * t369 + t606;
t713 = t470 + t472;
t603 = Ifges(7,5) * t369;
t607 = Ifges(6,4) * t369;
t726 = t372 * t729 + t603 - t607;
t725 = t606 / 0.2e1 - t602 / 0.2e1 + t728 * t629;
t628 = t369 / 0.2e1;
t367 = t369 ^ 2;
t545 = t367 + t368;
t685 = mrSges(7,2) + mrSges(6,3);
t712 = t545 * t685;
t722 = mrSges(4,1) + t712;
t361 = t369 * mrSges(6,1);
t363 = t372 * mrSges(6,2);
t669 = t363 + t361;
t360 = t369 * mrSges(7,1);
t362 = t372 * mrSges(7,3);
t670 = -t362 + t360;
t695 = t669 + t670;
t721 = mrSges(5,3) + t695;
t720 = t329 * t712 / 0.2e1;
t462 = t369 * pkin(5) - qJ(6) * t372;
t447 = m(7) * t462;
t386 = qJD(5) * (-t447 - t695);
t344 = (-pkin(8) - pkin(7)) * t371;
t538 = t625 * pkin(7);
t345 = pkin(8) * t625 + t538;
t249 = -t624 * t344 + t345 * t370;
t671 = t370 * t344 + t624 * t345;
t718 = pkin(2) * (t249 * t370 + t624 * t671);
t694 = -t329 * pkin(4) + t671;
t717 = t369 * t694;
t716 = t372 * t694;
t354 = -pkin(2) * t625 - pkin(1);
t563 = t330 * qJ(4);
t430 = t354 - t563;
t207 = t329 * pkin(3) + t430;
t663 = m(5) * t207 - mrSges(5,2) * t329 - mrSges(5,3) * t330;
t715 = t329 * t713 + t330 * t684;
t464 = -Ifges(7,3) * t372 + t603;
t128 = -Ifges(7,6) * t329 + t330 * t464;
t468 = Ifges(6,2) * t372 + t607;
t129 = -Ifges(6,6) * t329 + t330 * t468;
t130 = -Ifges(7,4) * t329 + t330 * t470;
t131 = -Ifges(6,5) * t329 + t330 * t472;
t601 = Ifges(6,6) * t369;
t626 = t372 / 0.2e1;
t632 = -t329 / 0.2e1;
t693 = t372 * t725 + t628 * t726;
t392 = t128 * t628 + t129 * t629 + (Ifges(7,6) * t369 + t372 * t684 - t601) * t632 + (t130 + t131) * t626 + (-Ifges(4,5) + Ifges(5,4)) * t329 + (-Ifges(4,6) + Ifges(5,5) + t693) * t330;
t463 = pkin(5) * t372 + qJ(6) * t369;
t450 = -pkin(4) - t463;
t393 = t330 * t450 - t249;
t678 = t393 * t670;
t170 = pkin(4) * t330 + t249;
t698 = t170 * t669;
t701 = t671 * mrSges(5,2);
t702 = t671 * mrSges(4,1);
t703 = t249 * mrSges(5,3);
t704 = t249 * mrSges(4,2);
t711 = t392 + t678 + t701 + t704 - t698 - t702 - t703;
t710 = -t678 / 0.2e1 + t698 / 0.2e1 - t701 / 0.2e1 + t702 / 0.2e1 + t703 / 0.2e1 - t704 / 0.2e1;
t638 = pkin(3) + pkin(9);
t132 = t329 * t638 + t430;
t80 = -t132 * t369 + t170 * t372;
t59 = -pkin(5) * t330 - t80;
t682 = t59 + t80;
t621 = m(7) * qJ(6);
t359 = mrSges(7,3) + t621;
t700 = qJ(4) * t170;
t699 = t170 * t694;
t615 = t370 * pkin(2);
t351 = qJ(4) + t615;
t697 = t351 * t170;
t560 = t330 * t372;
t215 = mrSges(6,2) * t329 + mrSges(6,3) * t560;
t216 = mrSges(7,2) * t560 - t329 * mrSges(7,3);
t674 = t215 + t216;
t692 = -pkin(3) * t671 - qJ(4) * t249;
t537 = t624 * pkin(2);
t353 = -t537 - pkin(3);
t691 = -t351 * t249 + t353 * t671;
t689 = t369 * (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) - t372 * (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1);
t608 = mrSges(7,3) * t369;
t610 = mrSges(7,1) * t372;
t475 = t608 + t610;
t203 = t475 * t329;
t609 = mrSges(6,2) * t369;
t611 = mrSges(6,1) * t372;
t477 = -t609 + t611;
t204 = t477 * t329;
t205 = t475 * t330;
t206 = t477 * t330;
t561 = t330 * t369;
t212 = -t329 * mrSges(6,1) - mrSges(6,3) * t561;
t298 = mrSges(7,2) * t561;
t589 = t329 * mrSges(7,1);
t213 = t298 + t589;
t562 = t330 * qJ(6);
t553 = t372 * t132;
t557 = t369 * t170;
t81 = t553 + t557;
t58 = t81 + t562;
t94 = t329 * t450 + t671;
t688 = t170 * t204 + ((t129 / 0.2e1 - t128 / 0.2e1) * t372 + (t131 / 0.2e1 + t130 / 0.2e1) * t369 + (-t689 - t614) * t329) * t329 - t393 * t203 - t694 * t206 - t94 * t205 + t207 * (-mrSges(5,2) * t330 + mrSges(5,3) * t329) + t80 * t212 + t59 * t213 + t81 * t215 + t58 * t216 + t354 * (mrSges(4,1) * t330 - mrSges(4,2) * t329);
t686 = mrSges(6,1) + mrSges(7,1);
t681 = t393 * t94;
t335 = qJ(4) + t462;
t308 = t335 + t615;
t680 = t308 * t393;
t679 = t335 * t393;
t565 = t329 * t369;
t535 = mrSges(6,3) * t565;
t210 = mrSges(6,1) * t330 - t535;
t536 = mrSges(7,2) * t565;
t211 = -mrSges(7,1) * t330 + t536;
t677 = t210 - t211;
t676 = t212 - t213;
t564 = t329 * t372;
t534 = mrSges(6,3) * t564;
t214 = -t330 * mrSges(6,2) + t534;
t299 = mrSges(7,2) * t564;
t587 = t330 * mrSges(7,3);
t217 = t299 + t587;
t675 = t214 + t217;
t673 = Ifges(7,6) * t565 + t564 * t684;
t502 = t670 * t629;
t672 = t589 / 0.2e1 + t329 * t502;
t644 = -mrSges(7,1) / 0.2e1;
t668 = (t502 + t644) * t329 - t298;
t221 = pkin(3) * t330 + qJ(4) * t329;
t325 = t330 * pkin(9);
t148 = t221 + t325;
t86 = -t148 * t369 + t716;
t87 = t372 * t148 + t717;
t458 = t369 * t87 + t372 * t86;
t309 = t329 * qJ(6);
t366 = t371 * pkin(2);
t208 = t221 + t366;
t133 = t208 + t325;
t83 = t372 * t133 + t717;
t61 = -t309 + t83;
t617 = pkin(5) * t329;
t82 = -t133 * t369 + t716;
t62 = -t82 + t617;
t664 = -t369 * t61 + t372 * t62;
t650 = m(6) / 0.4e1;
t662 = t650 + m(5) / 0.4e1;
t269 = t335 * t475;
t334 = qJ(4) * t477;
t661 = -t269 / 0.2e1 - t334 / 0.2e1;
t640 = mrSges(7,3) / 0.2e1;
t643 = -mrSges(6,2) / 0.2e1;
t645 = mrSges(6,1) / 0.2e1;
t648 = m(7) / 0.2e1;
t65 = -t309 + t87;
t66 = -t86 + t617;
t659 = (-pkin(5) * t66 + qJ(6) * t65) * t648 + t87 * t643 + t86 * t645 + t66 * t644 + t65 * t640;
t658 = (-pkin(5) * t62 + qJ(6) * t61) * t648 + t83 * t643 + t82 * t645 + t62 * t644 + t61 * t640;
t400 = Ifges(7,6) * t330 + t329 * t464;
t401 = Ifges(6,6) * t330 + t329 * t468;
t431 = t462 * t329;
t657 = -t670 * t431 / 0.2e1 + t463 * t203 / 0.2e1 - t694 * t477 / 0.2e1 - t94 * t475 / 0.2e1 - t714 * t330 / 0.4e1 + t715 * t369 / 0.4e1 + (-t400 / 0.4e1 + t401 / 0.4e1) * t372 + (t713 / 0.4e1 + t725) * t565 + (-t726 / 0.2e1 - t464 / 0.4e1 + t468 / 0.4e1) * t564;
t655 = 2 * qJD(3);
t654 = m(5) / 0.2e1;
t652 = -m(6) / 0.2e1;
t651 = m(6) / 0.2e1;
t649 = -m(7) / 0.2e1;
t647 = m(7) / 0.4e1;
t646 = -mrSges(6,1) / 0.2e1;
t642 = mrSges(6,2) / 0.2e1;
t639 = m(6) + m(5);
t637 = m(7) * t62;
t636 = m(7) * t66;
t635 = t213 / 0.2e1;
t634 = -t308 / 0.2e1;
t633 = t308 / 0.2e1;
t630 = t330 / 0.2e1;
t627 = -t372 / 0.2e1;
t622 = m(5) * t671;
t620 = m(7) * t308;
t619 = m(7) * t335;
t618 = m(7) * t369;
t2 = (mrSges(4,2) * t366 + t731) * t330 + m(6) * (t80 * t82 + t81 * t83 - t699) + m(7) * (t58 * t61 + t59 * t62 + t681) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t371 + (m(4) * t354 + mrSges(4,1) * t329) * pkin(2)) * t371 + t82 * t210 + t62 * t211 + t83 * t214 + t61 * t217 + ((Ifges(3,1) - Ifges(3,2)) * t371 - pkin(1) * mrSges(3,2) + Ifges(3,4) * t625) * t625 + t688 + t663 * t208;
t600 = t2 * qJD(1);
t3 = t731 * t330 + m(6) * (t80 * t86 + t81 * t87 - t699) + m(7) * (t58 * t65 + t59 * t66 + t681) + t86 * t210 + t66 * t211 + t87 * t214 + t65 * t217 + t663 * t221 + t688;
t591 = t3 * qJD(1);
t590 = t329 * mrSges(5,1);
t588 = t330 * mrSges(5,1);
t586 = t369 * mrSges(7,2);
t583 = t372 * mrSges(7,2);
t580 = t372 * t94;
t579 = t81 * t372;
t578 = t82 * t372;
t577 = t83 * t369;
t433 = t670 * t329;
t434 = t669 * t329;
t9 = t694 * t434 + t59 * t299 + t94 * t433 - t58 * t536 + t673 * t630 + (-Ifges(6,6) * t630 + t400 / 0.2e1 - t401 / 0.2e1) * t565 + t715 * t564 / 0.2e1 + (m(7) * t94 - t203) * t431 + (m(7) * t59 - t535 - t677) * t81 + (m(7) * t58 - t534 + t675) * t80 + t693 * t329 ^ 2;
t576 = t9 * qJD(1);
t573 = qJ(4) * t206;
t15 = (-t675 * t372 + t677 * t369 + m(7) * (-t369 * t59 - t58 * t372) + m(6) * (t369 * t80 - t579) - t663) * t330;
t571 = qJD(1) * t15;
t27 = m(7) * (t330 * t58 - t565 * t94) + t203 * t565 + t330 * t217;
t570 = qJD(1) * t27;
t411 = (t58 - t81) * t372 + t682 * t369;
t380 = t210 * t629 + t211 * t628 + t411 * t648 + t675 * t626 + t685 * t329 * (-t368 / 0.2e1 - t367 / 0.2e1);
t11 = (-t447 / 0.2e1 - t363 / 0.2e1 - t360 / 0.2e1 + t362 / 0.2e1 - t361 / 0.2e1) * t330 + t380;
t569 = t11 * qJD(1);
t566 = t308 * t205;
t559 = t335 * t205;
t558 = t351 * t206;
t349 = -pkin(9) + t353;
t556 = t369 * t349;
t555 = t369 * t638;
t552 = t372 * t349;
t551 = t372 * t638;
t547 = t308 + t335;
t546 = t545 * t638 * t615;
t200 = m(7) * t561;
t544 = t200 * qJD(1);
t540 = mrSges(6,3) * t578;
t539 = mrSges(6,3) * t577;
t533 = t351 * t588;
t532 = -t615 / 0.2e1;
t531 = t615 / 0.2e1;
t522 = t646 + t644;
t521 = t642 - mrSges(7,3) / 0.2e1;
t518 = t624 * mrSges(4,2);
t517 = -t590 / 0.2e1;
t516 = t588 / 0.2e1;
t515 = -t586 / 0.2e1;
t514 = t586 / 0.2e1;
t513 = -t583 / 0.2e1;
t512 = t583 / 0.2e1;
t511 = -t579 / 0.2e1;
t510 = t579 / 0.2e1;
t509 = t624 * t351;
t505 = t561 / 0.2e1;
t504 = -t560 / 0.2e1;
t503 = t560 / 0.2e1;
t501 = t556 / 0.2e1;
t500 = t555 / 0.2e1;
t499 = t203 * t627;
t195 = t212 * t626;
t498 = -t552 / 0.2e1;
t497 = t551 / 0.2e1;
t493 = t545 * t370;
t491 = t329 * t537;
t488 = mrSges(5,2) * t615 + t537 * t721;
t481 = t547 * t649;
t310 = pkin(2) * t493;
t480 = t310 / 0.2e1 + t532;
t479 = t349 * t493;
t460 = t369 * t65 - t372 * t66;
t459 = t577 + t578;
t375 = (t372 * t210 + t369 * t675) * t531 + t353 * t517 + t66 * t512 + t65 * t515 + t213 * t498 + t349 * t195 - t710 + (-t204 - t203) * t537 / 0.2e1 - t558 / 0.2e1 + (-t697 + t458 * t349 + (t624 * t694 + (t369 * t81 + t372 * t80) * t370) * pkin(2)) * t651 - t566 / 0.2e1 + (t680 + t460 * t349 + (t624 * t94 + (t369 * t58 - t372 * t59) * t370) * pkin(2)) * t648 + t674 * t501 - t458 * mrSges(6,3) / 0.2e1 + t516 * t615 - t533 / 0.2e1 - mrSges(5,1) * t491 / 0.2e1 + t372 * t211 * t532 + (t691 + t718) * t654;
t376 = pkin(3) * t517 + t62 * t513 + t61 * t514 + qJ(4) * t516 + t212 * t497 + t710 + t559 / 0.2e1 + (t638 * t664 + t679) * t649 - t551 * t635 + (-t459 * t638 - t700) * t652 + t674 * t500 + t540 / 0.2e1 + t539 / 0.2e1 + t573 / 0.2e1 - m(5) * t692 / 0.2e1;
t4 = t375 + t376;
t41 = -t488 + (t518 - m(7) * (t308 * t624 + t479) - m(6) * (t509 + t479) - m(5) * (t353 * t370 + t509)) * pkin(2) + t722 * t615;
t457 = t4 * qJD(1) - t41 * qJD(2);
t270 = t463 * t670;
t404 = t464 * t629 + t468 * t628 + t627 * t713 + t270 - t693;
t439 = t351 * t477;
t448 = m(7) * t463;
t42 = t439 + (t448 + t475) * t308 + t404;
t449 = -pkin(5) * t213 / 0.2e1 + qJ(6) * t216 / 0.2e1;
t374 = mrSges(7,2) * t511 + Ifges(6,6) * t503 + Ifges(7,6) * t504 + t505 * t684 + t58 * t512 + t514 * t682 + t632 * t683 + t449 + t657;
t413 = -t433 / 0.2e1;
t414 = -t434 / 0.2e1;
t444 = t463 * t94;
t378 = (t308 * t431 + t349 * t411 + t444) * t649 + t308 * t413 + t351 * t414 + t210 * t501 - t211 * t556 / 0.2e1 + t675 * t498 + t349 * t720;
t6 = t374 + t378 + t658;
t456 = -t6 * qJD(1) + t42 * qJD(2);
t109 = 0.4e1 * t351 * t662 + t620 + t721;
t398 = t369 * t521 + t372 * t522;
t451 = t648 * t94 + t651 * t694;
t381 = (-t212 / 0.2e1 + t635) * t372 + (-t215 / 0.2e1 - t216 / 0.2e1) * t369 + t398 * t329 + t451;
t391 = t459 * t652 - t649 * t664;
t13 = t381 + t391;
t454 = -qJD(1) * t13 - qJD(2) * t109;
t209 = (t670 + t620) * t372;
t389 = (-t580 + (-t308 * t329 + t330 * t349) * t369) * t648 - t499;
t22 = -t637 / 0.2e1 + t389 + t668;
t453 = -qJD(1) * t22 + qJD(2) * t209;
t33 = t587 + 0.2e1 * (t562 / 0.2e1 + t553 / 0.4e1 + t557 / 0.4e1 - t81 / 0.4e1) * m(7);
t452 = qJD(1) * t33 + qJD(5) * t359;
t440 = t335 * t463;
t399 = -t613 * t369 + (-t728 + t729) * t372;
t28 = -t270 + (t351 * t646 + mrSges(7,1) * t634 + t495 + pkin(5) * t481 + (pkin(5) * t648 - t522) * t615) * t372 + (t351 * t642 + mrSges(7,3) * t634 + qJ(6) * t481 + (t621 / 0.2e1 - t521) * t615 + t399) * t369 + t661;
t44 = -m(7) * t440 + t368 * t613 + t369 * t399 - t269 - t270 - t334;
t377 = (t335 * t431 - t411 * t638 + t444) * t649 + qJ(4) * t414 + t335 * t413 - t210 * t555 / 0.2e1 + t211 * t500 + t675 * t497 - t638 * t720;
t8 = t374 + t377 + t659;
t412 = t8 * qJD(1) + t28 * qJD(2) + t44 * qJD(3);
t390 = t458 * t652 + t460 * t649;
t16 = t381 + t390;
t168 = t639 * qJ(4) + t619 + t721;
t78 = t480 * m(6) + (t480 - t462) * m(7) + (-m(7) - t639) * qJ(4) - t721;
t410 = qJD(1) * t16 - qJD(2) * t78 + qJD(3) * t168;
t110 = (t670 + (t532 + t633 + t335 / 0.2e1) * m(7)) * t372;
t388 = (-t580 + (-t329 * t335 - t330 * t638) * t369) * t648 - t499;
t23 = -t636 / 0.2e1 + t388 + t668;
t232 = (t670 + t619) * t372;
t409 = -qJD(1) * t23 + qJD(2) * t110 + qJD(3) * t232;
t387 = (pkin(5) * t586 + (-mrSges(7,2) * qJ(6) - t612) * t372 - t719) * qJD(5);
t383 = t195 + t213 * t627 + (-mrSges(5,1) + t398) * t329 + t451 + t674 * t628;
t379 = -t657 + mrSges(7,2) * t510 + t58 * t513 + t449 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t329 + t689 * t330 + t682 * t515 + (t510 + t511) * mrSges(6,3);
t355 = qJ(4) * t537;
t333 = (-m(7) * t638 - mrSges(7,2)) * t369;
t304 = (m(7) * t349 - mrSges(7,2)) * t369;
t111 = m(7) * t547 * t627 + (m(7) * t532 - t670) * t372;
t79 = 0.4e1 * t462 * t647 + m(5) * t531 + 0.2e1 * (t647 + t650) * t310 + t721 + (t647 + t662) * (0.4e1 * qJ(4) + 0.2e1 * t615);
t31 = 0.2e1 * t58 * t648 + t217;
t29 = (t611 / 0.2e1 - t609 / 0.2e1 + t610 / 0.2e1 + t608 / 0.2e1 + t448 / 0.2e1) * t615 + (t308 * t463 + t440) * t648 + t404 + t439 / 0.2e1 + t475 * t633 - t661;
t26 = t636 / 0.2e1 + t388 + t672;
t24 = t637 / 0.2e1 + t389 + t672;
t14 = t383 - t390 + t622;
t12 = t671 * t654 + t622 / 0.2e1 + t383 - t391;
t10 = mrSges(6,2) * t503 + mrSges(7,3) * t504 + t447 * t630 + t505 * t686 + t380;
t7 = -t377 + t379 + t659;
t5 = -t378 + t379 + t658;
t1 = t392 + t375 - t376;
t17 = [qJD(2) * t2 + qJD(3) * t3 + qJD(4) * t15 + qJD(5) * t9 + qJD(6) * t27, t600 + (t711 - t558 + m(7) * (-t349 * t664 + t680) - t566 - m(4) * t718 + m(6) * (t349 * t459 - t697) + t674 * t556 + t676 * t552 + t664 * mrSges(7,2) + (mrSges(3,2) * pkin(7) - Ifges(3,6)) * t371 + (-t330 * t615 + t491) * mrSges(4,3) + Ifges(3,5) * t625 - t353 * t590 - t533 - mrSges(3,1) * t538 - t540 - t539 + m(5) * t691) * qJD(2) + t1 * qJD(3) + t12 * qJD(4) + t5 * qJD(5) + t24 * qJD(6), t591 + t1 * qJD(2) + t14 * qJD(4) + t7 * qJD(5) + t26 * qJD(6) + (t692 * t654 + (-t460 * t638 + t679) * t648 + (-t458 * t638 - t700) * t651) * t655 + (-mrSges(5,1) * t563 + pkin(3) * t590 - t559 - t573 + (t66 * mrSges(7,2) - t86 * mrSges(6,3) - t638 * t676) * t372 + (-t65 * mrSges(7,2) - t87 * mrSges(6,3) - t638 * t674) * t369 + t711) * qJD(3), qJD(2) * t12 + qJD(3) * t14 + qJD(5) * t10 + t571, t5 * qJD(2) + t7 * qJD(3) + t10 * qJD(4) + t31 * qJD(6) + t576 + ((-m(7) * pkin(5) - t686) * t81 + (-mrSges(6,2) + t359) * t80 + t673 + (-mrSges(7,2) * t463 - t601) * t329) * qJD(5), qJD(2) * t24 + qJD(3) * t26 + qJD(5) * t31 + t570; qJD(3) * t4 + qJD(4) * t13 - qJD(5) * t6 + qJD(6) * t22 - t600, -qJD(3) * t41 + qJD(4) * t109 + qJD(5) * t42 - qJD(6) * t209, t79 * qJD(4) + t29 * qJD(5) + t111 * qJD(6) + ((t335 * t537 - t546) * t648 + (t355 - t546) * t651 + (-pkin(3) * t615 + t355) * t654) * t655 + t457 + (t488 + (-t370 * t722 - t518) * pkin(2)) * qJD(3), qJD(3) * t79 - t454, t29 * qJD(3) + t304 * qJD(6) + t349 * t386 + t387 + t456, qJD(3) * t111 + qJD(5) * t304 - t453; -qJD(2) * t4 + qJD(4) * t16 - qJD(5) * t8 + qJD(6) * t23 - t591, -qJD(4) * t78 - qJD(5) * t28 - qJD(6) * t110 - t457, qJD(4) * t168 - qJD(5) * t44 - qJD(6) * t232, t410, t333 * qJD(6) - t386 * t638 + t387 - t412, qJD(5) * t333 - t409; -qJD(2) * t13 - qJD(3) * t16 + qJD(5) * t11 + qJD(6) * t200 - t571, qJD(3) * t78 + t454, -t410, 0, qJD(6) * t618 + t386 + t569, qJD(5) * t618 + t544; qJD(2) * t6 + qJD(3) * t8 - qJD(4) * t11 + qJD(6) * t33 - t576, qJD(3) * t28 - t456, t412, -t569, t359 * qJD(6), t452; -qJD(2) * t22 - qJD(3) * t23 - qJD(4) * t200 - qJD(5) * t33 - t570, qJD(3) * t110 + t453, t409, -t544, -t452, 0;];
Cq  = t17;
