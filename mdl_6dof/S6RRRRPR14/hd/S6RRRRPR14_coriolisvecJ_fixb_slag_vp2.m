% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:46
% EndTime: 2019-03-10 00:10:45
% DurationCPUTime: 65.64s
% Computational Cost: add. (47563->1148), mult. (141803->1634), div. (0->0), fcn. (118662->14), ass. (0->483)
t459 = cos(qJ(2));
t451 = cos(pkin(6));
t591 = pkin(1) * t451;
t443 = t459 * t591;
t434 = qJD(1) * t443;
t455 = sin(qJ(2));
t448 = sin(pkin(6));
t450 = cos(pkin(7));
t494 = t448 * (-pkin(10) * t450 - pkin(9));
t478 = t455 * t494;
t349 = qJD(1) * t478 + t434;
t442 = t455 * t591;
t465 = t459 * t494 - t442;
t350 = t465 * qJD(1);
t447 = sin(pkin(7));
t589 = pkin(10) * t447;
t468 = (pkin(2) * t455 - t459 * t589) * t448;
t384 = qJD(1) * t468;
t454 = sin(qJ(3));
t551 = t447 * t454;
t436 = pkin(10) * t551;
t458 = cos(qJ(3));
t545 = t450 * t458;
t411 = pkin(2) * t545 - t436;
t546 = t450 * t454;
t645 = t411 * qJD(3) - t458 * t349 - t350 * t546 - t384 * t551;
t286 = -t350 * t447 + t450 * t384;
t542 = t455 * t458;
t543 = t454 * t459;
t475 = t450 * t542 + t543;
t528 = qJD(1) * t448;
t369 = t475 * t528;
t541 = t458 * t459;
t544 = t454 * t455;
t473 = -t450 * t544 + t541;
t370 = t473 * t528;
t524 = qJD(3) * t447;
t700 = -pkin(3) * t369 + pkin(11) * t370 - t286 + (pkin(3) * t454 - pkin(11) * t458) * t524;
t514 = t455 * t528;
t498 = t447 * t514;
t699 = pkin(11) * t498 - t645;
t476 = t450 * t541 - t544;
t501 = qJD(1) * t451 + qJD(2);
t479 = t447 * t501;
t316 = t458 * t479 + t476 * t528;
t313 = qJD(4) - t316;
t453 = sin(qJ(4));
t457 = cos(qJ(4));
t308 = t370 * t453 - t457 * t498;
t410 = t450 * t453 + t457 * t551;
t522 = qJD(3) * t458;
t508 = t447 * t522;
t354 = qJD(4) * t410 + t453 * t508;
t691 = t308 - t354;
t309 = t370 * t457 + t453 * t498;
t409 = -t457 * t450 + t453 * t551;
t353 = -qJD(4) * t409 + t457 * t508;
t698 = t309 - t353;
t550 = t447 * t458;
t413 = pkin(2) * t546 + pkin(10) * t550;
t697 = -t413 * qJD(3) + t454 * t349;
t392 = pkin(11) * t450 + t413;
t393 = (-pkin(3) * t458 - pkin(11) * t454 - pkin(2)) * t447;
t520 = qJD(4) * t457;
t521 = qJD(4) * t453;
t646 = -t392 * t521 + t393 * t520 + t453 * t700 - t699 * t457;
t549 = t448 * t459;
t529 = pkin(9) * t549 + t442;
t400 = t529 * qJD(1);
t513 = t459 * t528;
t310 = t400 + (t450 * t513 + t479) * pkin(10);
t466 = t451 * pkin(2) + t478;
t315 = qJD(2) * pkin(2) + qJD(1) * t466 + t434;
t378 = (-pkin(2) * t459 - t455 * t589 - pkin(1)) * t448;
t365 = qJD(1) * t378;
t483 = t315 * t450 + t365 * t447;
t220 = t458 * t310 + t454 * t483;
t696 = -qJD(5) * t453 - t220 + t313 * (pkin(4) * t453 - qJ(5) * t457);
t643 = t350 * t545 - (-pkin(3) * t514 - t384 * t458) * t447 - t697;
t474 = t450 * t543 + t542;
t467 = t474 * t448;
t317 = qJD(1) * t467 + t454 * t479;
t464 = t447 * t513 - t450 * t501 - qJD(3);
t272 = t453 * t317 + t457 * t464;
t463 = (qJD(2) * t473 + qJD(3) * t476) * t448;
t471 = qJD(3) * t479;
t274 = qJD(1) * t463 + t458 * t471;
t526 = qJD(2) * t448;
t505 = qJD(1) * t526;
t496 = t455 * t505;
t480 = t447 * t496;
t185 = -qJD(4) * t272 + t457 * t274 + t453 * t480;
t621 = t185 / 0.2e1;
t695 = Ifges(5,4) * t621;
t523 = qJD(3) * t454;
t694 = qJ(5) * t369 - (qJ(5) * t523 - qJD(5) * t458) * t447 - t646;
t693 = -pkin(4) * t691 + t698 * qJ(5) - qJD(5) * t410 + t643;
t219 = -t454 * t310 + t458 * t483;
t261 = pkin(3) * t317 - pkin(11) * t316;
t155 = t457 * t219 + t453 * t261;
t124 = qJ(5) * t317 + t155;
t446 = sin(pkin(13));
t449 = cos(pkin(13));
t518 = pkin(11) * t521;
t681 = t696 * t449 + (t124 + t518) * t446;
t692 = -t449 * t124 + t446 * t696;
t262 = -t309 * t446 + t369 * t449;
t511 = t447 * t523;
t303 = -t353 * t446 + t449 * t511;
t533 = t262 - t303;
t263 = t309 * t449 + t369 * t446;
t304 = t353 * t449 + t446 * t511;
t532 = t263 - t304;
t462 = (qJD(2) * t475 + qJD(3) * t474) * t448;
t275 = qJD(1) * t462 + t454 * t471;
t137 = -t185 * t446 + t275 * t449;
t625 = t137 / 0.2e1;
t138 = t185 * t449 + t275 * t446;
t624 = t138 / 0.2e1;
t273 = t457 * t317 - t453 * t464;
t186 = qJD(4) * t273 + t453 * t274 - t457 * t480;
t619 = t186 / 0.2e1;
t603 = -t275 / 0.2e1;
t547 = t449 * t457;
t255 = t316 * t547 + t317 * t446;
t554 = t316 * t453;
t690 = -pkin(5) * t554 + pkin(12) * t255 + (pkin(5) * t453 - pkin(12) * t547) * qJD(4) + t681;
t552 = t446 * t457;
t254 = -t316 * t552 + t317 * t449;
t548 = t449 * t453;
t689 = -pkin(12) * t254 + (-pkin(11) * t548 - pkin(12) * t552) * qJD(4) + t692;
t653 = t446 * t694 + t449 * t693;
t652 = t446 * t693 - t449 * t694;
t688 = t369 - t511;
t452 = sin(qJ(6));
t456 = cos(qJ(6));
t228 = t273 * t449 + t313 * t446;
t502 = -t273 * t446 + t449 * t313;
t670 = -t228 * t452 + t456 * t502;
t45 = qJD(6) * t670 + t137 * t452 + t138 * t456;
t132 = t228 * t456 + t452 * t502;
t46 = -qJD(6) * t132 + t137 * t456 - t138 * t452;
t10 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t186;
t352 = t465 * qJD(2);
t331 = qJD(1) * t352;
t385 = qJD(2) * t468;
t374 = qJD(1) * t385;
t430 = qJD(2) * t434;
t469 = qJD(2) * t478;
t330 = qJD(1) * t469 + t430;
t510 = t450 * t523;
t500 = -t310 * t522 - t315 * t510 - t454 * t330 - t365 * t511;
t127 = -t331 * t545 + (-pkin(3) * t496 - t374 * t458) * t447 - t500;
t509 = t450 * t522;
t139 = -t310 * t523 + t315 * t509 + t458 * t330 + t331 * t546 + t365 * t508 + t374 * t551;
t126 = pkin(11) * t480 + t139;
t280 = -t331 * t447 + t450 * t374;
t158 = pkin(3) * t275 - pkin(11) * t274 + t280;
t264 = -t315 * t447 + t450 * t365;
t195 = -pkin(3) * t316 - pkin(11) * t317 + t264;
t201 = -pkin(11) * t464 + t220;
t37 = t457 * t126 + t453 * t158 + t195 * t520 - t201 * t521;
t34 = qJ(5) * t275 + qJD(5) * t313 + t37;
t55 = pkin(4) * t186 - qJ(5) * t185 - qJD(5) * t273 + t127;
t13 = -t34 * t446 + t449 * t55;
t14 = t449 * t34 + t446 * t55;
t602 = t275 / 0.2e1;
t620 = -t186 / 0.2e1;
t635 = t46 / 0.2e1;
t636 = t45 / 0.2e1;
t5 = pkin(5) * t186 - pkin(12) * t138 + t13;
t6 = pkin(12) * t137 + t14;
t200 = pkin(3) * t464 - t219;
t111 = t272 * pkin(4) - t273 * qJ(5) + t200;
t104 = t195 * t453 + t201 * t457;
t94 = qJ(5) * t313 + t104;
t51 = t449 * t111 - t446 * t94;
t36 = pkin(5) * t272 - pkin(12) * t228 + t51;
t52 = t446 * t111 + t449 * t94;
t39 = pkin(12) * t502 + t52;
t8 = t36 * t456 - t39 * t452;
t1 = qJD(6) * t8 + t452 * t5 + t456 * t6;
t9 = t36 * t452 + t39 * t456;
t2 = -qJD(6) * t9 - t452 * t6 + t456 * t5;
t669 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t687 = -t669 - 0.2e1 * Ifges(6,5) * t624 - Ifges(7,5) * t636 - 0.2e1 * Ifges(6,6) * t625 - Ifges(7,6) * t635 - t127 * mrSges(5,1) - t13 * mrSges(6,1) - t10 / 0.2e1 + t695 + Ifges(5,2) * t620 + t14 * mrSges(6,2) - (t603 - t602) * Ifges(5,6) - ((2 * Ifges(6,3)) + Ifges(7,3) + Ifges(5,2)) * t619;
t686 = mrSges(5,2) * t127;
t685 = t200 * mrSges(5,1);
t684 = t200 * mrSges(5,2);
t683 = -pkin(5) * t691 + pkin(12) * t532 + t653;
t682 = -pkin(12) * t533 + t652;
t680 = -t449 * t518 + t692;
t647 = -t392 * t520 - t393 * t521 + t699 * t453 + t457 * t700;
t679 = Ifges(5,1) * t621 + Ifges(5,5) * t602;
t678 = -t51 * mrSges(6,1) - t8 * mrSges(7,1) + t52 * mrSges(6,2) + t9 * mrSges(7,2);
t282 = t451 * t511 + t462;
t676 = -t282 / 0.2e1;
t283 = t451 * t508 + t463;
t675 = t283 / 0.2e1;
t566 = t228 * Ifges(6,5);
t567 = t502 * Ifges(6,6);
t105 = t272 * Ifges(6,3) + t566 + t567;
t266 = qJD(6) + t272;
t565 = t266 * Ifges(7,3);
t570 = t132 * Ifges(7,5);
t571 = t670 * Ifges(7,6);
t65 = t565 + t570 + t571;
t656 = t105 + t65;
t426 = -pkin(4) * t457 - qJ(5) * t453 - pkin(3);
t418 = t449 * t426;
t332 = -pkin(12) * t548 + t418 + (-pkin(11) * t446 - pkin(5)) * t457;
t380 = pkin(11) * t547 + t446 * t426;
t553 = t446 * t453;
t355 = -pkin(12) * t553 + t380;
t271 = t332 * t452 + t355 * t456;
t674 = -qJD(6) * t271 - t452 * t689 + t456 * t690;
t270 = t332 * t456 - t355 * t452;
t673 = qJD(6) * t270 + t452 * t690 + t456 * t689;
t154 = -t453 * t219 + t261 * t457;
t125 = -pkin(4) * t317 - t154;
t515 = pkin(5) * t446 + pkin(11);
t672 = pkin(5) * t254 + t515 * t520 - t125;
t648 = pkin(4) * t688 - t647;
t103 = t195 * t457 - t453 * t201;
t671 = mrSges(5,1) * t103 - mrSges(5,2) * t104;
t668 = t103 * mrSges(5,3) - t684;
t631 = Ifges(5,4) * t620 + t679;
t667 = -t104 * mrSges(5,3) - t678 + t685;
t666 = -t37 * mrSges(5,3) - t687 - t695;
t604 = t274 / 0.2e1;
t665 = t464 / 0.2e1;
t664 = t480 / 0.2e1;
t663 = Ifges(4,5) * t675 + Ifges(4,6) * t676;
t391 = t436 + (-pkin(2) * t458 - pkin(3)) * t450;
t289 = pkin(4) * t409 - qJ(5) * t410 + t391;
t295 = t457 * t392 + t453 * t393;
t292 = -qJ(5) * t550 + t295;
t223 = t449 * t289 - t292 * t446;
t348 = t410 * t449 - t446 * t550;
t177 = pkin(5) * t409 - pkin(12) * t348 + t223;
t224 = t446 * t289 + t449 * t292;
t347 = -t410 * t446 - t449 * t550;
t199 = pkin(12) * t347 + t224;
t98 = t177 * t452 + t199 * t456;
t661 = -qJD(6) * t98 - t452 * t682 + t456 * t683;
t97 = t177 * t456 - t199 * t452;
t660 = qJD(6) * t97 + t452 * t683 + t456 * t682;
t588 = pkin(12) + qJ(5);
t427 = t588 * t446;
t428 = t588 * t449;
t357 = -t427 * t456 - t428 * t452;
t481 = t446 * t452 - t449 * t456;
t205 = pkin(4) * t273 + qJ(5) * t272;
t81 = -t103 * t446 + t449 * t205;
t62 = pkin(12) * t272 * t449 + pkin(5) * t273 + t81;
t555 = t272 * t446;
t82 = t449 * t103 + t446 * t205;
t72 = pkin(12) * t555 + t82;
t655 = -qJD(5) * t481 + qJD(6) * t357 - t452 * t62 - t456 * t72;
t358 = -t427 * t452 + t428 * t456;
t420 = t446 * t456 + t449 * t452;
t654 = -qJD(5) * t420 - qJD(6) * t358 + t452 * t72 - t456 * t62;
t651 = pkin(5) * t533 + t648;
t644 = -(t350 * t450 + t384 * t447) * t458 + t697;
t38 = -t453 * t126 + t158 * t457 - t195 * t521 - t201 * t520;
t642 = t37 * t457 - t38 * t453;
t407 = t481 * qJD(6);
t333 = (t447 * t451 + t450 * t549) * pkin(10) + t529;
t346 = t443 + t466;
t242 = -t454 * t333 + t458 * (t346 * t450 + t378 * t447);
t641 = -t38 * mrSges(5,1) + t37 * mrSges(5,2);
t107 = t228 * Ifges(6,1) + Ifges(6,4) * t502 + Ifges(6,5) * t272;
t265 = Ifges(5,4) * t272;
t579 = Ifges(6,4) * t449;
t489 = -Ifges(6,2) * t446 + t579;
t580 = Ifges(6,4) * t446;
t490 = Ifges(6,1) * t449 - t580;
t491 = mrSges(6,1) * t446 + mrSges(6,2) * t449;
t563 = t313 * Ifges(5,5);
t593 = t449 / 0.2e1;
t612 = t228 / 0.2e1;
t614 = t502 / 0.2e1;
t564 = t273 * Ifges(5,1);
t172 = -t265 + t563 + t564;
t622 = t172 / 0.2e1;
t106 = t228 * Ifges(6,4) + Ifges(6,2) * t502 + Ifges(6,6) * t272;
t630 = -t106 / 0.2e1;
t93 = -pkin(4) * t313 + qJD(5) - t103;
t640 = t489 * t614 + t490 * t612 + t93 * t491 + (-t446 * t52 - t449 * t51) * mrSges(6,3) + t622 + t563 / 0.2e1 + t446 * t630 + t107 * t593 - t265 / 0.2e1 - t668;
t639 = -0.2e1 * pkin(1);
t638 = Ifges(7,4) * t636 + Ifges(7,2) * t635 + Ifges(7,6) * t619;
t637 = Ifges(7,1) * t636 + Ifges(7,4) * t635 + Ifges(7,5) * t619;
t57 = t138 * Ifges(6,4) + t137 * Ifges(6,2) + t186 * Ifges(6,6);
t634 = t57 / 0.2e1;
t633 = Ifges(6,1) * t624 + Ifges(6,4) * t625 + Ifges(6,5) * t619;
t629 = -t670 / 0.2e1;
t628 = t670 / 0.2e1;
t627 = -t132 / 0.2e1;
t626 = t132 / 0.2e1;
t562 = t313 * Ifges(5,6);
t583 = Ifges(5,4) * t273;
t171 = -t272 * Ifges(5,2) + t562 + t583;
t623 = t171 / 0.2e1;
t191 = Ifges(4,5) * t274 - Ifges(4,6) * t275 + Ifges(4,3) * t480;
t618 = t191 / 0.2e1;
t617 = Ifges(4,1) * t604 + Ifges(4,4) * t603 + Ifges(4,5) * t664;
t615 = -t502 / 0.2e1;
t613 = -t228 / 0.2e1;
t610 = -t266 / 0.2e1;
t609 = t266 / 0.2e1;
t608 = -t272 / 0.2e1;
t607 = t272 / 0.2e1;
t606 = -t273 / 0.2e1;
t605 = t273 / 0.2e1;
t600 = -t313 / 0.2e1;
t599 = t313 / 0.2e1;
t598 = -t316 / 0.2e1;
t597 = -t317 / 0.2e1;
t596 = t317 / 0.2e1;
t594 = t447 / 0.2e1;
t590 = pkin(2) * t447;
t338 = -t448 * t476 - t451 * t550;
t435 = qJD(2) * t443;
t351 = t435 + t469;
t159 = -t333 * t523 + t346 * t509 + t458 * t351 + t352 * t546 + t378 * t508 + t385 * t551;
t512 = t455 * t526;
t497 = t447 * t512;
t150 = pkin(11) * t497 + t159;
t287 = -t352 * t447 + t450 * t385;
t176 = pkin(3) * t282 - pkin(11) * t283 + t287;
t281 = -t346 * t447 + t450 * t378;
t339 = t451 * t551 + t467;
t215 = pkin(3) * t338 - pkin(11) * t339 + t281;
t243 = t458 * t333 + t346 * t546 + t378 * t551;
t406 = -t447 * t549 + t451 * t450;
t222 = pkin(11) * t406 + t243;
t49 = t457 * t150 + t453 * t176 + t215 * t520 - t222 * t521;
t41 = qJ(5) * t282 + qJD(5) * t338 + t49;
t499 = -t333 * t522 - t346 * t510 - t454 * t351 - t378 * t511;
t151 = -t352 * t545 + (-pkin(3) * t512 - t385 * t458) * t447 - t499;
t291 = t339 * t457 + t406 * t453;
t208 = qJD(4) * t291 + t283 * t453 - t457 * t497;
t290 = t339 * t453 - t406 * t457;
t209 = -qJD(4) * t290 + t283 * t457 + t453 * t497;
t71 = pkin(4) * t208 - qJ(5) * t209 - qJD(5) * t291 + t151;
t20 = t449 * t41 + t446 * t71;
t587 = Ifges(3,4) * t455;
t586 = Ifges(4,4) * t317;
t585 = Ifges(4,4) * t454;
t584 = Ifges(4,4) * t458;
t582 = Ifges(5,4) * t453;
t581 = Ifges(5,4) * t457;
t578 = Ifges(7,4) * t132;
t577 = Ifges(3,5) * t459;
t576 = Ifges(4,5) * t317;
t575 = Ifges(6,5) * t449;
t574 = Ifges(4,6) * t316;
t573 = Ifges(6,6) * t446;
t572 = Ifges(4,3) * t450;
t569 = t139 * mrSges(4,2);
t140 = (t331 * t450 + t374 * t447) * t458 + t500;
t568 = t140 * mrSges(4,1);
t561 = t316 * mrSges(4,3);
t560 = t317 * mrSges(4,3);
t35 = -pkin(4) * t275 - t38;
t559 = t35 * t453;
t556 = Ifges(3,6) * qJD(2);
t119 = t453 * t215 + t457 * t222;
t113 = qJ(5) * t338 + t119;
t221 = -pkin(3) * t406 - t242;
t141 = pkin(4) * t290 - qJ(5) * t291 + t221;
t75 = t449 * t113 + t446 * t141;
t135 = -mrSges(6,1) * t502 + mrSges(6,2) * t228;
t231 = mrSges(5,1) * t313 - mrSges(5,3) * t273;
t540 = t135 - t231;
t174 = t254 * t456 - t255 * t452;
t307 = t407 * t453 - t420 * t520;
t539 = t174 - t307;
t175 = t254 * t452 + t255 * t456;
t408 = t420 * qJD(6);
t306 = -t408 * t453 - t481 * t520;
t538 = t175 - t306;
t276 = t347 * t456 - t348 * t452;
t178 = qJD(6) * t276 + t303 * t452 + t304 * t456;
t190 = t262 * t452 + t263 * t456;
t537 = t178 - t190;
t277 = t347 * t452 + t348 * t456;
t179 = -qJD(6) * t277 + t303 * t456 - t304 * t452;
t189 = t262 * t456 - t263 * t452;
t536 = t179 - t189;
t197 = t420 * t272;
t535 = t197 + t408;
t198 = t481 * t272;
t534 = t198 + t407;
t525 = qJD(2) * t450;
t517 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t87 = Ifges(5,5) * t185 - Ifges(5,6) * t186 + Ifges(5,3) * t275;
t16 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t503 = t524 / 0.2e1;
t19 = -t41 * t446 + t449 * t71;
t77 = -t137 * mrSges(6,1) + t138 * mrSges(6,2);
t74 = -t113 * t446 + t449 * t141;
t118 = t215 * t457 - t453 * t222;
t294 = -t453 * t392 + t393 * t457;
t293 = pkin(4) * t550 - t294;
t386 = -pkin(9) * t496 + t430;
t404 = t529 * qJD(2);
t387 = qJD(1) * t404;
t492 = -t387 * mrSges(3,1) - t386 * mrSges(3,2);
t487 = -t573 + t575;
t485 = -t13 * t446 + t14 * t449;
t251 = t291 * t449 + t338 * t446;
t48 = pkin(5) * t290 - pkin(12) * t251 + t74;
t250 = -t291 * t446 + t338 * t449;
t59 = pkin(12) * t250 + t75;
t17 = -t452 * t59 + t456 * t48;
t18 = t452 * t48 + t456 * t59;
t484 = t103 * t457 + t104 * t453;
t166 = t250 * t456 - t251 * t452;
t167 = t250 * t452 + t251 * t456;
t50 = -t453 * t150 + t176 * t457 - t215 * t521 - t222 * t520;
t114 = -pkin(4) * t338 - t118;
t47 = -pkin(4) * t282 - t50;
t460 = -t566 / 0.2e1 - t567 / 0.2e1 + t562 / 0.2e1 - t65 / 0.2e1 - t105 / 0.2e1 - t571 / 0.2e1 - t570 / 0.2e1 - t565 / 0.2e1 + t623 + t583 / 0.2e1 - t667;
t445 = -pkin(5) * t449 - pkin(4);
t432 = Ifges(3,4) * t513;
t429 = t505 * t577;
t424 = t515 * t453;
t412 = -pkin(9) * t448 * t455 + t443;
t403 = -pkin(9) * t512 + t435;
t398 = -pkin(9) * t514 + t434;
t397 = -mrSges(3,2) * t501 + mrSges(3,3) * t513;
t396 = mrSges(3,1) * t501 - mrSges(3,3) * t514;
t395 = t481 * t453;
t394 = t420 * t453;
t379 = -pkin(11) * t552 + t418;
t361 = Ifges(3,1) * t514 + Ifges(3,5) * t501 + t432;
t360 = t556 + (Ifges(3,6) * t451 + (Ifges(3,2) * t459 + t587) * t448) * qJD(1);
t312 = Ifges(4,4) * t316;
t279 = -mrSges(4,1) * t464 - t560;
t278 = mrSges(4,2) * t464 + t561;
t260 = -mrSges(4,1) * t316 + mrSges(4,2) * t317;
t258 = -mrSges(4,2) * t480 - mrSges(4,3) * t275;
t257 = mrSges(4,1) * t480 - mrSges(4,3) * t274;
t256 = -pkin(5) * t347 + t293;
t241 = Ifges(4,1) * t317 - Ifges(4,5) * t464 + t312;
t240 = Ifges(4,2) * t316 - t464 * Ifges(4,6) + t586;
t239 = -t464 * Ifges(4,3) + t574 + t576;
t230 = -mrSges(5,2) * t313 - mrSges(5,3) * t272;
t207 = mrSges(4,1) * t275 + mrSges(4,2) * t274;
t206 = mrSges(5,1) * t272 + mrSges(5,2) * t273;
t192 = Ifges(4,4) * t274 - Ifges(4,2) * t275 + Ifges(4,6) * t480;
t170 = t273 * Ifges(5,5) - t272 * Ifges(5,6) + t313 * Ifges(5,3);
t169 = mrSges(6,1) * t272 - mrSges(6,3) * t228;
t168 = -mrSges(6,2) * t272 + mrSges(6,3) * t502;
t162 = t209 * t449 + t282 * t446;
t161 = -t209 * t446 + t282 * t449;
t160 = (t352 * t450 + t385 * t447) * t458 + t499;
t143 = -mrSges(5,2) * t275 - mrSges(5,3) * t186;
t142 = mrSges(5,1) * t275 - mrSges(5,3) * t185;
t128 = Ifges(7,4) * t670;
t101 = mrSges(7,1) * t266 - mrSges(7,3) * t132;
t100 = -mrSges(7,2) * t266 + mrSges(7,3) * t670;
t99 = mrSges(5,1) * t186 + mrSges(5,2) * t185;
t91 = -pkin(5) * t250 + t114;
t90 = -pkin(5) * t555 + t104;
t86 = mrSges(6,1) * t186 - mrSges(6,3) * t138;
t85 = -mrSges(6,2) * t186 + mrSges(6,3) * t137;
t78 = -pkin(5) * t502 + t93;
t76 = -mrSges(7,1) * t670 + mrSges(7,2) * t132;
t67 = Ifges(7,1) * t132 + Ifges(7,5) * t266 + t128;
t66 = Ifges(7,2) * t670 + Ifges(7,6) * t266 + t578;
t61 = -qJD(6) * t167 + t161 * t456 - t162 * t452;
t60 = qJD(6) * t166 + t161 * t452 + t162 * t456;
t33 = -mrSges(7,2) * t186 + mrSges(7,3) * t46;
t32 = mrSges(7,1) * t186 - mrSges(7,3) * t45;
t28 = -pkin(5) * t161 + t47;
t27 = -pkin(5) * t137 + t35;
t15 = pkin(12) * t161 + t20;
t7 = pkin(5) * t208 - pkin(12) * t162 + t19;
t4 = -qJD(6) * t18 - t15 * t452 + t456 * t7;
t3 = qJD(6) * t17 + t15 * t456 + t452 * t7;
t11 = [(Ifges(7,1) * t167 + Ifges(7,4) * t166) * t636 + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t626 + ((t386 * t459 + t387 * t455) * mrSges(3,3) + ((t361 / 0.2e1 - t398 * mrSges(3,3) + Ifges(3,5) * qJD(2) / 0.2e1) * t459 + (-t360 / 0.2e1 - t400 * mrSges(3,3) - t556 / 0.2e1 + (-t220 * mrSges(4,2) + t219 * mrSges(4,1) + t576 / 0.2e1 + t574 / 0.2e1 + t239 / 0.2e1 + (qJD(3) / 0.2e1 + t525 / 0.2e1) * Ifges(4,3)) * t447) * t455) * qJD(2)) * t448 + m(3) * (t386 * t529 - t387 * t412 - t398 * t404 + t400 * t403) + (Ifges(6,4) * t251 + Ifges(6,2) * t250) * t625 + (Ifges(6,4) * t162 + Ifges(6,2) * t161) * t614 + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t628 + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t635 + (Ifges(6,5) * t162 + Ifges(6,6) * t161) * t607 + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t609 + (Ifges(5,4) * t291 + Ifges(5,6) * t338) * t620 + (Ifges(5,4) * t209 + Ifges(5,6) * t282) * t608 + (t1 * t166 - t167 * t2 - t60 * t8 + t61 * t9) * mrSges(7,3) + (t429 / 0.2e1 + t492) * t451 + (-t104 * t282 + t127 * t291 + t200 * t209 - t338 * t37) * mrSges(5,2) + m(7) * (t1 * t18 + t17 * t2 + t27 * t91 + t28 * t78 + t3 * t9 + t4 * t8) + m(6) * (t114 * t35 + t13 * t74 + t14 * t75 + t19 * t51 + t20 * t52 + t47 * t93) + m(5) * (t103 * t50 + t104 * t49 + t118 * t38 + t119 * t37 + t127 * t221 + t151 * t200) + m(4) * (t139 * t243 + t140 * t242 + t159 * t220 + t160 * t219 + t264 * t287 + t280 * t281) + (Ifges(6,5) * t251 + Ifges(7,5) * t167 + Ifges(6,6) * t250 + Ifges(7,6) * t166) * t619 - t406 * t569 + (-t139 * t338 - t140 * t339 - t219 * t283 - t220 * t282) * mrSges(4,3) + (-t13 * t251 + t14 * t250 + t161 * t52 - t162 * t51) * mrSges(6,3) + (Ifges(5,1) * t291 + Ifges(5,5) * t338) * t621 + (Ifges(5,1) * t209 + Ifges(5,5) * t282) * t605 + t316 * (Ifges(4,4) * t283 - Ifges(4,2) * t282) / 0.2e1 + (Ifges(5,5) * t209 + Ifges(5,3) * t282) * t599 + (Ifges(5,5) * t291 + Ifges(5,3) * t338) * t602 + (Ifges(6,1) * t251 + Ifges(6,4) * t250) * t624 + (Ifges(6,1) * t162 + Ifges(6,4) * t161) * t612 + t403 * t397 - t404 * t396 + (qJD(3) + t525) * t663 + t167 * t637 + t166 * t638 + t291 * t631 + t251 * t633 + t250 * t634 + t339 * t617 + t406 * t618 + t209 * t622 + (Ifges(4,4) * t339 - Ifges(4,2) * t338 + Ifges(4,6) * t406) * t603 + (Ifges(4,1) * t339 - Ifges(4,4) * t338 + Ifges(4,5) * t406) * t604 + (Ifges(4,1) * t283 - Ifges(4,4) * t282) * t596 + t406 * t568 + t280 * (mrSges(4,1) * t338 + mrSges(4,2) * t339) + t38 * (mrSges(5,1) * t338 - mrSges(5,3) * t291) + t338 * t87 / 0.2e1 - t338 * t192 / 0.2e1 + t17 * t32 + t18 * t33 + t666 * t290 + (Ifges(7,6) * t628 + Ifges(7,5) * t626 - Ifges(5,4) * t605 + Ifges(6,3) * t607 - Ifges(5,2) * t608 + Ifges(7,3) * t609 + Ifges(6,5) * t612 + Ifges(6,6) * t614 - Ifges(5,6) * t599 + t656 / 0.2e1 - t171 / 0.2e1 + t667) * t208 + (t406 * t663 + ((-t412 * mrSges(3,3) + Ifges(3,5) * t451 + (mrSges(3,2) * t639 + 0.3e1 / 0.2e1 * Ifges(3,4) * t459) * t448) * t459 + ((Ifges(4,5) * t339 - Ifges(4,6) * t338 + Ifges(4,3) * t406) * t594 - t529 * mrSges(3,3) + (t572 * t594 - 0.3e1 / 0.2e1 * Ifges(3,6)) * t451 + (mrSges(3,1) * t639 - 0.3e1 / 0.2e1 * t587 + (-t447 ^ 2 * Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t459) * t448) * t455) * t526) * qJD(1) + t61 * t66 / 0.2e1 + t60 * t67 / 0.2e1 + t28 * t76 + t78 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + t75 * t85 + t74 * t86 + t91 * t16 + t3 * t100 + t4 * t101 + t114 * t77 + t47 * t135 + t118 * t142 + t119 * t143 + t161 * t106 / 0.2e1 + t162 * t107 / 0.2e1 + t93 * (-mrSges(6,1) * t161 + mrSges(6,2) * t162) + t27 * (-mrSges(7,1) * t166 + mrSges(7,2) * t167) + t20 * t168 + t19 * t169 + t241 * t675 + t240 * t676 + t151 * t206 + t221 * t99 + t49 * t230 + t50 * t231 + t35 * (-mrSges(6,1) * t250 + mrSges(6,2) * t251) + t242 * t257 + t243 * t258 + t159 * t278 + t160 * t279 + t281 * t207 + t282 * t170 / 0.2e1 + t103 * (mrSges(5,1) * t282 - mrSges(5,3) * t209) + t264 * (mrSges(4,1) * t282 + mrSges(4,2) * t283) + t287 * t260; (-t87 / 0.2e1 + t192 / 0.2e1 - Ifges(5,6) * t620 - Ifges(5,5) * t621 - Ifges(5,3) * t602 + t641) * t550 + (Ifges(5,4) * t353 + Ifges(6,5) * t263 - Ifges(5,2) * t354 + Ifges(5,6) * t511 + Ifges(6,6) * t262 + Ifges(6,3) * t308) * t608 + t668 * t698 - ((-Ifges(3,2) * t514 + t361 + t432) * t459 + t501 * (-Ifges(3,6) * t455 + t577)) * t528 / 0.2e1 + (t1 * t276 - t2 * t277 + t536 * t9 - t537 * t8) * mrSges(7,3) + (-t455 * (Ifges(3,1) * t459 - t587) / 0.2e1 + pkin(1) * (mrSges(3,1) * t455 + mrSges(3,2) * t459)) * qJD(1) ^ 2 * t448 ^ 2 + t492 + (t170 - t240) * (t454 * t503 - t369 / 0.2e1) + (-t13 * t348 + t14 * t347 + t51 * t532 - t52 * t533) * mrSges(6,3) + (-t139 * t450 - t264 * t370 + (t220 * t514 + t264 * t522 + t280 * t454) * t447) * mrSges(4,2) + (t140 * t450 - t264 * t369 + (-t219 * t514 + t264 * t523 - t280 * t458) * t447) * mrSges(4,1) + (Ifges(6,4) * t348 + Ifges(6,2) * t347) * t625 + (t304 / 0.2e1 - t263 / 0.2e1) * t107 + t429 + (Ifges(7,4) * t277 + Ifges(7,2) * t276) * t635 - t207 * t590 + (Ifges(5,4) * t309 + Ifges(6,5) * t304 - Ifges(5,2) * t308 + Ifges(5,6) * t369 + Ifges(6,6) * t303 + Ifges(6,3) * t354) * t607 + t643 * t206 + t644 * t279 + (t139 * t413 + t140 * t411 + t219 * t644 + t220 * t645 - t264 * t286 - t280 * t590) * m(4) + t645 * t278 + t646 * t230 + t647 * t231 + (t103 * t647 + t104 * t646 + t127 * t391 + t200 * t643 + t294 * t38 + t295 * t37) * m(5) - t239 * t498 / 0.2e1 + (t353 / 0.2e1 - t309 / 0.2e1) * t172 - t464 * (Ifges(4,5) * t458 - Ifges(4,6) * t454) * t524 / 0.2e1 + (t316 * (-Ifges(4,2) * t454 + t584) + t317 * (Ifges(4,1) * t458 - t585) + t458 * t241) * t503 - Ifges(3,6) * t496 + t360 * t514 / 0.2e1 + (mrSges(6,1) * t533 - mrSges(6,2) * t532) * t93 + (-mrSges(7,1) * t536 + mrSges(7,2) * t537) * t78 + t648 * t135 + (t219 * t370 + t220 * t369 + (t139 * t458 - t140 * t454 + (-t219 * t458 - t220 * t454) * qJD(3)) * t447) * mrSges(4,3) + t411 * t257 + t413 * t258 - t398 * t397 + t400 * t396 + t391 * t99 - t370 * t241 / 0.2e1 + (t572 + (Ifges(4,5) * t454 + Ifges(4,6) * t458) * t447) * t664 + (Ifges(4,5) * t370 - Ifges(4,6) * t369 + Ifges(4,3) * t498) * t665 + t277 * t637 + t276 * t638 + (Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t308) * t627 + (Ifges(7,4) * t178 + Ifges(7,2) * t179 + Ifges(7,6) * t354) * t628 + (Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t308) * t629 + t348 * t633 + t347 * t634 + (Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t308) * t615 + t551 * t617 + t450 * t618 + (Ifges(7,1) * t178 + Ifges(7,4) * t179 + Ifges(7,5) * t354) * t626 + (Ifges(4,6) * t450 + (Ifges(4,2) * t458 + t585) * t447) * t603 + (Ifges(4,5) * t450 + (Ifges(4,1) * t454 + t584) * t447) * t604 + (Ifges(5,1) * t353 - Ifges(5,4) * t354 + Ifges(5,5) * t511) * t605 + (Ifges(5,1) * t309 - Ifges(5,4) * t308 + Ifges(5,5) * t369) * t606 + (Ifges(7,5) * t178 + Ifges(7,6) * t179 + Ifges(7,3) * t354) * t609 + (Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t308) * t610 + (Ifges(6,1) * t304 + Ifges(6,4) * t303 + Ifges(6,5) * t354) * t612 + (Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t308) * t613 + (Ifges(6,4) * t304 + Ifges(6,2) * t303 + Ifges(6,6) * t354) * t614 + (Ifges(4,1) * t370 - Ifges(4,4) * t369 + Ifges(4,5) * t498) * t597 + (Ifges(4,4) * t370 - Ifges(4,2) * t369 + Ifges(4,6) * t498) * t598 + (Ifges(5,5) * t353 - Ifges(5,6) * t354 + Ifges(5,3) * t511) * t599 + (Ifges(5,5) * t309 - Ifges(5,6) * t308 + Ifges(5,3) * t369) * t600 + t35 * (-mrSges(6,1) * t347 + mrSges(6,2) * t348) + t293 * t77 + t294 * t142 + t295 * t143 + (t303 / 0.2e1 - t262 / 0.2e1) * t106 + (t179 / 0.2e1 - t189 / 0.2e1) * t66 - t671 * t688 + t651 * t76 + t652 * t168 + t653 * t169 + (t13 * t223 + t14 * t224 + t293 * t35 + t51 * t653 + t52 * t652 + t648 * t93) * m(6) + (-t171 + t656) * (t354 / 0.2e1 - t308 / 0.2e1) + (t398 * t459 + t400 * t455) * mrSges(3,3) * t528 + t666 * t409 + t660 * t100 + t661 * t101 + (t1 * t98 + t2 * t97 + t256 * t27 + t651 * t78 + t660 * t9 + t661 * t8) * m(7) + (-mrSges(5,3) * t38 + 0.2e1 * t631 + t686) * t410 + (t178 / 0.2e1 - t190 / 0.2e1) * t67 - t667 * t691 + t97 * t32 + t98 * t33 + (Ifges(6,1) * t348 + Ifges(6,4) * t347) * t624 + (Ifges(6,5) * t348 + Ifges(7,5) * t277 + Ifges(6,6) * t347 + Ifges(7,6) * t276) * t619 + (Ifges(7,1) * t277 + Ifges(7,4) * t276) * t636 + t223 * t86 + t224 * t85 + t256 * t16 + t27 * (-mrSges(7,1) * t276 + mrSges(7,2) * t277) - t286 * t260; (Ifges(6,1) * t255 + Ifges(6,4) * t254) * t613 + (-Ifges(7,1) * t395 - Ifges(7,4) * t394) * t636 + (Ifges(7,1) * t175 + Ifges(7,4) * t174) * t627 + (-t13 * t548 - t14 * t553 - t254 * t52 + t255 * t51) * mrSges(6,3) + t642 * mrSges(5,3) + (Ifges(7,5) * t627 + Ifges(7,6) * t629 + Ifges(6,6) * t615 + t623 + Ifges(6,3) * t608 + Ifges(7,3) * t610 + Ifges(6,5) * t613 - t656 / 0.2e1 + t678) * t554 + (-t27 * t395 - t538 * t78) * mrSges(7,2) + (-pkin(3) * t127 - t103 * t154 - t104 * t155 - t200 * t220 + (-qJD(4) * t484 + t642) * pkin(11)) * m(5) + t582 * t620 + (-t264 * mrSges(4,2) + t484 * mrSges(5,3) + Ifges(4,1) * t597 + Ifges(4,5) * t665 + t581 * t607 - t582 * t606) * t316 + (-Ifges(7,5) * t395 - Ifges(7,6) * t394) * t619 + (-t1 * t394 + t2 * t395 + t538 * t8 - t539 * t9) * mrSges(7,3) + t568 - t569 - t57 * t553 / 0.2e1 + t581 * t621 + (-t586 + t170) * t597 + (Ifges(6,4) * t255 + Ifges(6,2) * t254) * t615 + (-t206 + t279 + t560) * t220 + (-t278 + t561) * t219 + (t27 * t394 + t539 * t78) * mrSges(7,1) + t191 + (t312 + t241) * t598 + t424 * t16 + t379 * t86 + t380 * t85 - t395 * t637 - t394 * t638 + (Ifges(7,1) * t306 + Ifges(7,4) * t307) * t626 + (Ifges(7,4) * t306 + Ifges(7,2) * t307) * t628 + t254 * t630 + t548 * t633 + (Ifges(7,5) * t306 + Ifges(7,6) * t307) * t609 + t240 * t596 + t491 * t559 + ((-t684 + Ifges(5,1) * t606 + Ifges(5,5) * t600 - t172 / 0.2e1) * t316 + t143 * pkin(11) + (t564 / 0.2e1 + t487 * t607 + (m(6) * t93 + t540) * pkin(11) + t640) * qJD(4) + t687) * t457 + (-Ifges(7,4) * t395 - Ifges(7,2) * t394) * t635 + (Ifges(7,4) * t175 + Ifges(7,2) * t174) * t629 + t680 * t168 + t681 * t169 + (pkin(11) * t559 - t125 * t93 + t13 * t379 + t14 * t380 + t51 * t681 + t52 * t680) * m(6) + ((-t142 + t77) * pkin(11) + (-Ifges(5,2) * t607 - Ifges(5,6) * t600 - t685) * t316 + (-pkin(11) * t230 + t272 * t517 - t460) * qJD(4) + t686 + t631 + t490 * t624 + t489 * t625 + t487 * t619 + t679) * t453 + (t307 / 0.2e1 - t174 / 0.2e1) * t66 - pkin(3) * t99 + (-t264 * mrSges(4,1) + Ifges(5,5) * t606 - Ifges(4,2) * t598 - Ifges(4,6) * t665 + Ifges(5,6) * t607 + Ifges(5,3) * t600 - t671) * t317 - t125 * t135 + t672 * t76 + t673 * t100 + t674 * t101 + (t1 * t271 + t2 * t270 + t27 * t424 + t672 * t78 + t673 * t9 + t674 * t8) * m(7) + (t306 / 0.2e1 - t175 / 0.2e1) * t67 - t155 * t230 - t154 * t231 + (Ifges(6,5) * t255 + Ifges(6,6) * t254) * t608 + (Ifges(7,5) * t175 + Ifges(7,6) * t174) * t610 - t255 * t107 / 0.2e1 - t93 * (-mrSges(6,1) * t254 + mrSges(6,2) * t255) + t270 * t32 + t271 * t33; -t641 + (Ifges(6,5) * t446 + Ifges(7,5) * t420 + Ifges(6,6) * t449 - Ifges(7,6) * t481) * t619 + (-t1 * t481 - t2 * t420 + t534 * t8 - t535 * t9) * mrSges(7,3) + t27 * (mrSges(7,1) * t481 + mrSges(7,2) * t420) + (Ifges(7,1) * t420 - Ifges(7,4) * t481) * t636 + (Ifges(7,4) * t420 - Ifges(7,2) * t481) * t635 - t481 * t638 + (-t446 * t86 + t449 * t85) * qJ(5) + (-t407 / 0.2e1 - t198 / 0.2e1) * t67 + (-Ifges(7,4) * t407 - Ifges(7,2) * t408) * t628 + (-Ifges(7,1) * t407 - Ifges(7,4) * t408) * t626 + (-Ifges(7,5) * t407 - Ifges(7,6) * t408) * t609 + (-t408 / 0.2e1 - t197 / 0.2e1) * t66 - t540 * t104 + t460 * t273 + ((t575 / 0.2e1 - t573 / 0.2e1) * t272 + (Ifges(5,1) / 0.2e1 - t517) * t273 + t640) * t272 + (t168 * t449 - t169 * t446) * qJD(5) + (mrSges(7,1) * t535 - mrSges(7,2) * t534) * t78 + t35 * (-mrSges(6,1) * t449 + mrSges(6,2) * t446) + t445 * t16 + t420 * t637 + (Ifges(7,1) * t198 + Ifges(7,4) * t197) * t627 + (Ifges(7,4) * t198 + Ifges(7,2) * t197) * t629 + t446 * t633 + (Ifges(6,1) * t446 + t579) * t624 + (Ifges(6,2) * t449 + t580) * t625 + (Ifges(7,5) * t198 + Ifges(7,6) * t197) * t610 + t57 * t593 + t357 * t32 + t358 * t33 + t485 * mrSges(6,3) + t654 * t101 + t655 * t100 + (t1 * t358 + t2 * t357 + t27 * t445 + t654 * t8 + t655 * t9 - t78 * t90) * m(7) + t87 - pkin(4) * t77 + (-pkin(4) * t35 + (-t446 * t51 + t449 * t52) * qJD(5) + t485 * qJ(5) - t104 * t93 - t51 * t81 - t52 * t82) * m(6) - t90 * t76 - t82 * t168 - t81 * t169 - t103 * t230; -t670 * t100 + t132 * t101 - t502 * t168 + t228 * t169 + t16 + t77 + (t132 * t8 - t670 * t9 + t27) * m(7) + (t228 * t51 - t502 * t52 + t35) * m(6); -t78 * (mrSges(7,1) * t132 + mrSges(7,2) * t670) + (Ifges(7,5) * t670 - Ifges(7,6) * t132) * t610 - t8 * t100 + t9 * t101 + (Ifges(7,1) * t670 - t578) * t627 + t66 * t626 + (t132 * t9 + t670 * t8) * mrSges(7,3) + t10 + (-Ifges(7,2) * t132 + t128 + t67) * t629 + t669;];
tauc  = t11(:);
