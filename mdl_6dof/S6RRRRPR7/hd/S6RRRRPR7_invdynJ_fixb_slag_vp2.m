% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:35
% EndTime: 2019-03-09 22:28:23
% DurationCPUTime: 61.13s
% Computational Cost: add. (41019->1092), mult. (98949->1469), div. (0->0), fcn. (80748->18), ass. (0->497)
t472 = cos(qJ(2));
t467 = sin(qJ(2));
t461 = sin(pkin(6));
t574 = qJD(1) * t461;
t543 = t467 * t574;
t463 = cos(pkin(6));
t573 = qJD(1) * t463;
t559 = pkin(1) * t573;
t369 = -pkin(8) * t543 + t472 * t559;
t495 = (pkin(2) * t467 - pkin(9) * t472) * t461;
t370 = qJD(1) * t495;
t466 = sin(qJ(3));
t471 = cos(qJ(3));
t263 = -t369 * t466 + t471 * t370;
t474 = -pkin(10) - pkin(9);
t544 = qJD(3) * t474;
t587 = t471 * t472;
t794 = -(pkin(3) * t467 - pkin(10) * t587) * t574 - t263 + t471 * t544;
t264 = t471 * t369 + t466 * t370;
t542 = t472 * t574;
t520 = t466 * t542;
t793 = -pkin(10) * t520 - t466 * t544 + t264;
t465 = sin(qJ(4));
t470 = cos(qJ(4));
t399 = -t465 * t466 + t470 * t471;
t733 = qJD(3) + qJD(4);
t301 = t733 * t399;
t322 = t399 * t542;
t792 = t301 - t322;
t400 = t465 * t471 + t466 * t470;
t302 = t733 * t400;
t321 = t400 * t542;
t789 = t302 - t321;
t422 = t474 * t466;
t423 = t474 * t471;
t566 = qJD(4) * t470;
t567 = qJD(4) * t465;
t742 = t422 * t566 + t423 * t567 + t465 * t794 - t793 * t470;
t324 = t465 * t422 - t470 * t423;
t741 = -qJD(4) * t324 + t793 * t465 + t470 * t794;
t791 = -pkin(4) * t543 - qJ(5) * t792 - qJD(5) * t400 + t741;
t790 = qJ(5) * t789 - qJD(5) * t399 - t742;
t783 = -mrSges(7,3) + mrSges(6,2);
t464 = sin(qJ(6));
t469 = cos(qJ(6));
t512 = -mrSges(7,1) * t469 + mrSges(7,2) * t464;
t781 = t512 - mrSges(6,1);
t434 = pkin(8) * t542;
t373 = t467 * t559 + t434;
t569 = qJD(3) * t466;
t738 = -t373 + (-t520 + t569) * pkin(3);
t419 = qJD(3) - t542;
t411 = qJD(4) + t419;
t664 = -t411 / 0.2e1;
t439 = qJD(2) + t573;
t340 = t439 * t471 - t466 * t543;
t341 = t439 * t466 + t471 * t543;
t251 = t340 * t465 + t341 * t470;
t460 = sin(pkin(12));
t462 = cos(pkin(12));
t526 = t470 * t340 - t341 * t465;
t762 = t462 * t251 + t460 * t526;
t678 = -t762 / 0.2e1;
t763 = -t251 * t460 + t462 * t526;
t680 = -t763 / 0.2e1;
t788 = -Ifges(6,4) * t678 - Ifges(6,2) * t680 - Ifges(6,6) * t664;
t175 = qJD(6) - t763;
t681 = t175 / 0.2e1;
t162 = t411 * t464 + t469 * t762;
t687 = t162 / 0.2e1;
t161 = t411 * t469 - t464 * t762;
t689 = t161 / 0.2e1;
t787 = Ifges(7,5) * t687 + Ifges(7,6) * t689 + Ifges(7,3) * t681;
t773 = Ifges(6,3) + Ifges(5,3);
t786 = m(7) * pkin(5) - t781;
t525 = -m(7) * pkin(11) + t783;
t752 = t460 * t791 - t790 * t462;
t739 = pkin(4) * t789 + t738;
t114 = Ifges(6,1) * t762 + Ifges(6,4) * t763 + t411 * Ifges(6,5);
t663 = t411 / 0.2e1;
t677 = t762 / 0.2e1;
t679 = t763 / 0.2e1;
t319 = -pkin(2) * t439 - t369;
t252 = -pkin(3) * t340 + t319;
t188 = -pkin(4) * t526 + qJD(5) + t252;
t320 = pkin(9) * t439 + t373;
t359 = (-pkin(2) * t472 - pkin(9) * t467 - pkin(1)) * t461;
t331 = qJD(1) * t359;
t236 = -t320 * t466 + t471 * t331;
t204 = -pkin(10) * t341 + t236;
t193 = pkin(3) * t419 + t204;
t237 = t320 * t471 + t331 * t466;
t205 = pkin(10) * t340 + t237;
t201 = t465 * t205;
t135 = t470 * t193 - t201;
t769 = qJ(5) * t251;
t115 = t135 - t769;
t104 = pkin(4) * t411 + t115;
t203 = t470 * t205;
t136 = t193 * t465 + t203;
t749 = qJ(5) * t526;
t116 = t136 + t749;
t606 = t116 * t460;
t61 = t104 * t462 - t606;
t646 = t61 * mrSges(6,3);
t728 = mrSges(6,2) * t188 - t646;
t785 = Ifges(6,1) * t677 + Ifges(6,4) * t679 + Ifges(6,5) * t663 + t114 / 0.2e1 + t728;
t110 = t462 * t116;
t62 = t460 * t104 + t110;
t60 = pkin(11) * t411 + t62;
t95 = -pkin(5) * t763 - pkin(11) * t762 + t188;
t29 = -t464 * t60 + t469 * t95;
t30 = t464 * t95 + t469 * t60;
t784 = mrSges(6,1) * t188 + mrSges(7,1) * t29 - mrSges(7,2) * t30 - mrSges(6,3) * t62 + t787 - t788;
t782 = -pkin(11) * t543 + t752;
t216 = t301 * t460 + t462 * t302;
t217 = t301 * t462 - t302 * t460;
t225 = t462 * t321 + t322 * t460;
t226 = -t321 * t460 + t322 * t462;
t780 = t739 + (-t217 + t226) * pkin(11) + (t216 - t225) * pkin(5);
t570 = qJD(2) * t472;
t377 = (qJD(1) * t570 + qJDD(1) * t467) * t461;
t562 = qJDD(1) * t463;
t438 = qJDD(2) + t562;
t240 = qJD(3) * t340 + t377 * t471 + t438 * t466;
t241 = -qJD(3) * t341 - t377 * t466 + t438 * t471;
t140 = qJD(4) * t526 + t240 * t470 + t241 * t465;
t571 = qJD(2) * t461;
t541 = t467 * t571;
t563 = qJDD(1) * t461;
t376 = -qJD(1) * t541 + t472 * t563;
t362 = qJDD(3) - t376;
t351 = qJDD(4) + t362;
t659 = pkin(1) * t463;
t558 = qJD(2) * t659;
t521 = qJD(1) * t558;
t555 = pkin(1) * t562;
t271 = pkin(8) * t376 + t467 * t555 + t472 * t521;
t256 = pkin(9) * t438 + t271;
t261 = -pkin(1) * t563 - pkin(2) * t376 - pkin(9) * t377;
t147 = -qJD(3) * t237 - t256 * t466 + t471 * t261;
t109 = pkin(3) * t362 - pkin(10) * t240 + t147;
t568 = qJD(3) * t471;
t146 = t471 * t256 + t466 * t261 - t320 * t569 + t331 * t568;
t125 = pkin(10) * t241 + t146;
t42 = -qJD(4) * t136 + t470 * t109 - t125 * t465;
t24 = pkin(4) * t351 - qJ(5) * t140 - qJD(5) * t251 + t42;
t141 = -qJD(4) * t251 - t240 * t465 + t241 * t470;
t41 = t465 * t109 + t470 * t125 + t193 * t566 - t205 * t567;
t26 = qJ(5) * t141 + qJD(5) * t526 + t41;
t11 = t460 * t24 + t462 * t26;
t89 = t140 * t462 + t141 * t460;
t57 = qJD(6) * t161 + t351 * t464 + t469 * t89;
t58 = -qJD(6) * t162 + t351 * t469 - t464 * t89;
t88 = -t140 * t460 + t141 * t462;
t87 = qJDD(6) - t88;
t14 = Ifges(7,5) * t57 + Ifges(7,6) * t58 + Ifges(7,3) * t87;
t666 = t351 / 0.2e1;
t701 = t89 / 0.2e1;
t702 = t88 / 0.2e1;
t703 = t87 / 0.2e1;
t709 = t58 / 0.2e1;
t710 = t57 / 0.2e1;
t597 = t461 * t467;
t440 = pkin(8) * t597;
t272 = -qJD(2) * t434 - qJDD(1) * t440 - t467 * t521 + t472 * t555;
t257 = -pkin(2) * t438 - t272;
t186 = -pkin(3) * t241 + t257;
t99 = -pkin(4) * t141 + qJDD(5) + t186;
t22 = -pkin(5) * t88 - pkin(11) * t89 + t99;
t8 = pkin(11) * t351 + t11;
t2 = qJD(6) * t29 + t22 * t464 + t469 * t8;
t3 = -qJD(6) * t30 + t22 * t469 - t464 * t8;
t727 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t778 = t727 + mrSges(6,1) * t99 - mrSges(6,3) * t11 + Ifges(7,5) * t710 + Ifges(7,6) * t709 + Ifges(7,3) * t703 + t14 / 0.2e1 + (-t666 - t351 / 0.2e1) * Ifges(6,6) + (-t702 - t88 / 0.2e1) * Ifges(6,2) + (-t701 - t89 / 0.2e1) * Ifges(6,4);
t776 = -Ifges(6,4) * t677 - Ifges(6,2) * t679 - Ifges(6,6) * t663 + t784 + t787;
t775 = pkin(5) * t762 - pkin(11) * t763;
t16 = Ifges(7,1) * t57 + Ifges(7,4) * t58 + Ifges(7,5) * t87;
t511 = mrSges(7,1) * t464 + mrSges(7,2) * t469;
t59 = -pkin(5) * t411 - t61;
t489 = t59 * t511;
t504 = Ifges(7,5) * t469 - Ifges(7,6) * t464;
t627 = Ifges(7,4) * t469;
t506 = -Ifges(7,2) * t464 + t627;
t628 = Ifges(7,4) * t464;
t509 = Ifges(7,1) * t469 - t628;
t565 = qJD(6) * t464;
t532 = -t565 / 0.2e1;
t160 = Ifges(7,4) * t161;
t81 = t162 * Ifges(7,1) + t175 * Ifges(7,5) + t160;
t609 = t469 * t81;
t548 = t609 / 0.2e1;
t564 = qJD(6) * t469;
t635 = mrSges(7,3) * t469;
t636 = mrSges(7,3) * t464;
t661 = t464 / 0.2e1;
t682 = -t175 / 0.2e1;
t688 = -t162 / 0.2e1;
t690 = -t161 / 0.2e1;
t10 = t24 * t462 - t26 * t460;
t7 = -pkin(5) * t351 - t10;
t15 = Ifges(7,4) * t57 + Ifges(7,2) * t58 + Ifges(7,6) * t87;
t714 = t15 / 0.2e1;
t719 = Ifges(7,5) * t688 + Ifges(7,6) * t690 + Ifges(7,3) * t682 - t784 + t788;
t722 = Ifges(6,1) * t678 + Ifges(6,4) * t680 + Ifges(6,5) * t664 - t114 / 0.2e1;
t764 = -t42 * mrSges(5,1) - t10 * mrSges(6,1) + t41 * mrSges(5,2) + t11 * mrSges(6,2);
t767 = Ifges(5,5) * t140 + Ifges(6,5) * t89 + Ifges(5,6) * t141 + Ifges(6,6) * t88 + t351 * t773;
t629 = Ifges(7,4) * t162;
t80 = Ifges(7,2) * t161 + Ifges(7,6) * t175 + t629;
t774 = (t29 * t635 + t30 * t636 + t504 * t682 + t506 * t690 + t509 * t688 - t489 - t609 / 0.2e1 + t80 * t661 + t722 - t728) * t763 + t719 * t762 - t3 * t636 + t2 * t635 + t16 * t661 + t469 * t714 + (Ifges(7,1) * t464 + t627) * t710 + (Ifges(7,2) * t469 + t628) * t709 + t7 * t512 + t80 * t532 + (Ifges(7,5) * t464 + Ifges(7,6) * t469) * t703 + (-t29 * t564 - t30 * t565) * mrSges(7,3) + (t489 + t548) * qJD(6) + (t161 * t506 + t162 * t509 + t175 * t504) * qJD(6) / 0.2e1 - t764 + t767;
t652 = pkin(4) * t251;
t771 = t271 * mrSges(3,2);
t754 = t790 * t460 + t462 * t791;
t770 = -m(4) * pkin(9) + m(5) * t474 - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t101 = -mrSges(7,1) * t161 + mrSges(7,2) * t162;
t168 = mrSges(6,1) * t411 - mrSges(6,3) * t762;
t746 = -t101 + t168;
t459 = qJ(3) + qJ(4);
t454 = pkin(12) + t459;
t446 = sin(t454);
t447 = cos(t454);
t768 = -t446 * t525 + t447 * t786;
t473 = cos(qJ(1));
t586 = t472 * t473;
t468 = sin(qJ(1));
t590 = t467 * t468;
t390 = -t463 * t590 + t586;
t455 = sin(t459);
t456 = cos(t459);
t596 = t461 * t468;
t298 = -t390 * t455 + t456 * t596;
t766 = -t455 * t597 + t456 * t463;
t765 = -t147 * mrSges(4,1) + t146 * mrSges(4,2);
t760 = t99 * mrSges(6,2) - mrSges(6,3) * t10 + 0.2e1 * Ifges(6,1) * t701 + 0.2e1 * Ifges(6,4) * t702 + 0.2e1 * Ifges(6,5) * t666;
t757 = -m(7) - m(6);
t694 = t140 / 0.2e1;
t693 = t141 / 0.2e1;
t674 = t240 / 0.2e1;
t673 = t241 / 0.2e1;
t665 = t362 / 0.2e1;
t296 = -t462 * t399 + t400 * t460;
t297 = t399 * t460 + t400 * t462;
t653 = pkin(3) * t471;
t452 = pkin(2) + t653;
t352 = -pkin(4) * t399 - t452;
t196 = pkin(5) * t296 - pkin(11) * t297 + t352;
t266 = qJ(5) * t399 + t324;
t323 = t470 * t422 + t423 * t465;
t490 = -qJ(5) * t400 + t323;
t198 = t462 * t266 + t460 * t490;
t131 = t196 * t469 - t198 * t464;
t756 = qJD(6) * t131 + t464 * t780 + t469 * t782;
t132 = t196 * t464 + t198 * t469;
t755 = -qJD(6) * t132 - t464 * t782 + t469 * t780;
t753 = pkin(5) * t543 - t754;
t716 = m(5) * pkin(3);
t750 = -t716 - mrSges(4,1);
t167 = -mrSges(6,2) * t411 + mrSges(6,3) * t763;
t105 = -mrSges(7,2) * t175 + mrSges(7,3) * t161;
t106 = mrSges(7,1) * t175 - mrSges(7,3) * t162;
t501 = t105 * t469 - t106 * t464;
t745 = -t167 - t501;
t211 = -t226 * t464 + t469 * t543;
t488 = t464 * t217 + t297 * t564;
t744 = t211 + t488;
t212 = t226 * t469 + t464 * t543;
t487 = -t469 * t217 + t297 * t565;
t743 = t212 + t487;
t594 = t461 * t472;
t396 = pkin(8) * t594 + t467 * t659;
t358 = pkin(9) * t463 + t396;
t259 = -t358 * t466 + t471 * t359;
t595 = t461 * t471;
t386 = t463 * t466 + t467 * t595;
t215 = -pkin(3) * t594 - pkin(10) * t386 + t259;
t260 = t471 * t358 + t466 * t359;
t385 = t463 * t471 - t466 * t597;
t223 = pkin(10) * t385 + t260;
t152 = t465 * t215 + t470 * t223;
t523 = mrSges(3,3) * t543;
t740 = -mrSges(3,1) * t439 - mrSges(4,1) * t340 + mrSges(4,2) * t341 + t523;
t658 = pkin(1) * t472;
t395 = t463 * t658 - t440;
t480 = mrSges(3,2) + t770;
t737 = -t480 + t511;
t514 = -mrSges(4,1) * t471 + mrSges(4,2) * t466;
t736 = -m(4) * pkin(2) - m(5) * t452 - t456 * mrSges(5,1) + t455 * mrSges(5,2) + t514;
t734 = t146 * t471 - t147 * t466;
t732 = mrSges(3,1) - t736;
t333 = -t446 * t597 + t447 * t463;
t334 = t446 * t463 + t447 * t597;
t731 = -t766 * mrSges(5,1) - (-t455 * t463 - t456 * t597) * mrSges(5,2) + t783 * t334 + t781 * t333;
t289 = t390 * t446 - t447 * t596;
t290 = t390 * t447 + t446 * t596;
t299 = t390 * t456 + t455 * t596;
t730 = -t298 * mrSges(5,1) + t299 * mrSges(5,2) - t781 * t289 + t290 * t783;
t588 = t468 * t472;
t589 = t467 * t473;
t388 = t463 * t589 + t588;
t593 = t461 * t473;
t285 = -t388 * t446 - t447 * t593;
t286 = t388 * t447 - t446 * t593;
t415 = t455 * t593;
t493 = -t388 * t455 - t456 * t593;
t729 = -t493 * mrSges(5,1) - (-t388 * t456 + t415) * mrSges(5,2) + t783 * t286 + t781 * t285;
t725 = t732 + t768;
t27 = mrSges(7,1) * t87 - mrSges(7,3) * t57;
t28 = -mrSges(7,2) * t87 + mrSges(7,3) * t58;
t724 = m(7) * (t2 * t469 - t3 * t464 + (-t29 * t469 - t30 * t464) * qJD(6)) - t106 * t564 - t105 * t565 + t28 * t469 - t27 * t464;
t634 = Ifges(3,4) * t467;
t508 = Ifges(3,2) * t472 + t634;
t610 = t439 * Ifges(3,6);
t721 = t136 * mrSges(5,2) + t508 * t574 / 0.2e1 + t610 / 0.2e1 - t135 * mrSges(5,1) - t61 * mrSges(6,1);
t717 = t461 ^ 2;
t713 = t16 / 0.2e1;
t708 = Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t666;
t707 = Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t666;
t700 = pkin(1) * mrSges(3,1);
t699 = pkin(1) * mrSges(3,2);
t692 = Ifges(4,4) * t674 + Ifges(4,2) * t673 + Ifges(4,6) * t665;
t691 = Ifges(4,1) * t674 + Ifges(4,4) * t673 + Ifges(4,5) * t665;
t630 = Ifges(5,4) * t251;
t171 = Ifges(5,2) * t526 + t411 * Ifges(5,6) + t630;
t686 = -t171 / 0.2e1;
t685 = t171 / 0.2e1;
t246 = Ifges(5,4) * t526;
t172 = t251 * Ifges(5,1) + t411 * Ifges(5,5) + t246;
t684 = -t172 / 0.2e1;
t683 = t172 / 0.2e1;
t614 = t341 * Ifges(4,4);
t230 = t340 * Ifges(4,2) + t419 * Ifges(4,6) + t614;
t676 = t230 / 0.2e1;
t332 = Ifges(4,4) * t340;
t231 = t341 * Ifges(4,1) + t419 * Ifges(4,5) + t332;
t675 = t231 / 0.2e1;
t672 = -t526 / 0.2e1;
t671 = t526 / 0.2e1;
t670 = -t251 / 0.2e1;
t669 = t251 / 0.2e1;
t667 = t341 / 0.2e1;
t662 = t463 / 0.2e1;
t657 = pkin(3) * t341;
t656 = pkin(3) * t385;
t655 = pkin(3) * t466;
t654 = pkin(3) * t470;
t651 = pkin(4) * t460;
t650 = pkin(4) * t462;
t269 = t385 * t470 - t386 * t465;
t540 = t461 * t570;
t303 = -qJD(3) * t386 - t466 * t540;
t304 = qJD(3) * t385 + t471 * t540;
t183 = qJD(4) * t269 + t303 * t465 + t304 * t470;
t270 = t385 * t465 + t386 * t470;
t371 = qJD(2) * t495;
t374 = t395 * qJD(2);
t195 = -qJD(3) * t260 + t471 * t371 - t374 * t466;
t159 = pkin(3) * t541 - pkin(10) * t304 + t195;
t194 = -t358 * t569 + t359 * t568 + t466 * t371 + t471 * t374;
t169 = pkin(10) * t303 + t194;
t68 = -qJD(4) * t152 + t470 * t159 - t169 * t465;
t47 = pkin(4) * t541 - qJ(5) * t183 - qJD(5) * t270 + t68;
t184 = -qJD(4) * t270 + t303 * t470 - t304 * t465;
t67 = t465 * t159 + t470 * t169 + t215 * t566 - t223 * t567;
t49 = qJ(5) * t184 + qJD(5) * t269 + t67;
t20 = t460 * t47 + t462 * t49;
t642 = mrSges(4,2) * t471;
t640 = mrSges(5,3) * t135;
t639 = mrSges(5,3) * t136;
t638 = mrSges(5,3) * t526;
t637 = mrSges(5,3) * t251;
t633 = Ifges(3,4) * t472;
t632 = Ifges(4,4) * t466;
t631 = Ifges(4,4) * t471;
t626 = pkin(3) * qJD(4);
t623 = t763 * Ifges(6,6);
t622 = t762 * Ifges(6,5);
t621 = t236 * mrSges(4,3);
t620 = t237 * mrSges(4,3);
t619 = t526 * Ifges(5,6);
t618 = t251 * Ifges(5,5);
t615 = t340 * Ifges(4,6);
t613 = t341 * Ifges(4,5);
t612 = t419 * Ifges(4,3);
t611 = t439 * Ifges(3,5);
t603 = t297 * t464;
t602 = t297 * t469;
t598 = t460 * t465;
t592 = t462 * t465;
t151 = t470 * t215 - t223 * t465;
t129 = -pkin(4) * t594 - qJ(5) * t270 + t151;
t133 = qJ(5) * t269 + t152;
t72 = t460 * t129 + t462 * t133;
t143 = t470 * t204 - t201;
t413 = pkin(4) * t455 + t655;
t414 = pkin(4) * t456 + t653;
t579 = -t390 * t413 + t414 * t596;
t576 = -t413 * t597 + t463 * t414;
t375 = pkin(8) * t540 + t467 * t558;
t451 = pkin(4) + t654;
t379 = pkin(3) * t592 + t460 * t451;
t575 = t473 * pkin(1) + pkin(8) * t596;
t552 = t464 * t594;
t551 = t469 * t594;
t546 = Ifges(4,5) * t240 + Ifges(4,6) * t241 + Ifges(4,3) * t362;
t545 = Ifges(3,5) * t377 + Ifges(3,6) * t376 + Ifges(3,3) * t438;
t37 = -t88 * mrSges(6,1) + t89 * mrSges(6,2);
t531 = t468 * pkin(1) - pkin(8) * t593;
t530 = t285 * pkin(5) + pkin(11) * t286;
t529 = -t289 * pkin(5) + pkin(11) * t290;
t528 = t333 * pkin(5) + pkin(11) * t334;
t142 = -t204 * t465 - t203;
t524 = -m(5) * t655 - mrSges(3,3);
t522 = mrSges(3,3) * t542;
t517 = t298 * pkin(4);
t255 = -pkin(3) * t303 + t375;
t210 = t657 + t652;
t516 = -t388 * t413 - t414 * t593;
t515 = mrSges(4,1) * t385 - mrSges(4,2) * t386;
t510 = Ifges(4,1) * t471 - t632;
t507 = -Ifges(4,2) * t466 + t631;
t505 = Ifges(4,5) * t471 - Ifges(4,6) * t466;
t503 = -t29 * t464 + t30 * t469;
t19 = -t460 * t49 + t462 * t47;
t389 = t463 * t588 + t589;
t407 = pkin(2) + t414;
t458 = -qJ(5) + t474;
t502 = -t389 * t458 + t390 * t407 + t413 * t596 + t575;
t199 = -t462 * t269 + t270 * t460;
t200 = t269 * t460 + t270 * t462;
t357 = t440 + (-pkin(2) - t658) * t463;
t280 = t357 - t656;
t208 = -pkin(4) * t269 + t280;
t102 = pkin(5) * t199 - pkin(11) * t200 + t208;
t70 = -pkin(11) * t594 + t72;
t38 = t102 * t469 - t464 * t70;
t39 = t102 * t464 + t469 * t70;
t71 = t129 * t462 - t133 * t460;
t499 = t766 * pkin(4);
t378 = -pkin(3) * t598 + t451 * t462;
t189 = -t200 * t464 - t551;
t494 = -t200 * t469 + t552;
t306 = -t390 * t466 + t468 * t595;
t491 = t142 - t749;
t148 = -pkin(4) * t184 + t255;
t481 = t493 * pkin(4);
t387 = -t463 * t586 + t590;
t479 = -g(1) * t389 - g(2) * t387 + g(3) * t594;
t449 = -pkin(5) - t650;
t432 = Ifges(3,4) * t542;
t427 = t466 * t593;
t391 = (-mrSges(3,1) * t472 + mrSges(3,2) * t467) * t461;
t365 = -pkin(5) - t378;
t364 = -mrSges(3,2) * t439 + t522;
t317 = Ifges(3,1) * t543 + t432 + t611;
t307 = t390 * t471 + t466 * t596;
t284 = mrSges(4,1) * t419 - mrSges(4,3) * t341;
t283 = -mrSges(4,2) * t419 + mrSges(4,3) * t340;
t233 = t290 * t469 + t389 * t464;
t232 = -t290 * t464 + t389 * t469;
t229 = t612 + t613 + t615;
t222 = mrSges(5,1) * t411 - t637;
t221 = -mrSges(5,2) * t411 + t638;
t207 = -mrSges(4,2) * t362 + mrSges(4,3) * t241;
t206 = mrSges(4,1) * t362 - mrSges(4,3) * t240;
t197 = t266 * t460 - t462 * t490;
t187 = -mrSges(5,1) * t526 + mrSges(5,2) * t251;
t173 = -mrSges(4,1) * t241 + mrSges(4,2) * t240;
t170 = t411 * Ifges(5,3) + t618 + t619;
t127 = -mrSges(5,2) * t351 + mrSges(5,3) * t141;
t126 = mrSges(5,1) * t351 - mrSges(5,3) * t140;
t121 = t143 - t769;
t120 = t183 * t462 + t184 * t460;
t119 = t183 * t460 - t462 * t184;
t118 = -mrSges(6,1) * t763 + mrSges(6,2) * t762;
t112 = t411 * Ifges(6,3) + t622 + t623;
t100 = t652 + t775;
t98 = t210 + t775;
t92 = qJD(6) * t494 - t120 * t464 + t469 * t541;
t91 = qJD(6) * t189 + t120 * t469 + t464 * t541;
t90 = -mrSges(5,1) * t141 + mrSges(5,2) * t140;
t74 = mrSges(6,1) * t351 - mrSges(6,3) * t89;
t73 = -mrSges(6,2) * t351 + mrSges(6,3) * t88;
t69 = pkin(5) * t594 - t71;
t66 = t462 * t121 + t460 * t491;
t64 = t115 * t462 - t606;
t63 = t115 * t460 + t110;
t43 = pkin(5) * t119 - pkin(11) * t120 + t148;
t34 = t100 * t464 + t469 * t64;
t33 = t100 * t469 - t464 * t64;
t32 = t464 * t98 + t469 * t66;
t31 = -t464 * t66 + t469 * t98;
t21 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t18 = pkin(11) * t541 + t20;
t17 = -pkin(5) * t541 - t19;
t5 = -qJD(6) * t39 - t18 * t464 + t43 * t469;
t4 = qJD(6) * t38 + t18 * t469 + t43 * t464;
t1 = [(-Ifges(7,5) * t494 + Ifges(7,6) * t189) * t703 + (-Ifges(7,4) * t494 + Ifges(7,2) * t189) * t709 + (-Ifges(7,1) * t494 + Ifges(7,4) * t189) * t710 + (t189 * t2 - t29 * t91 + t3 * t494 + t30 * t92) * mrSges(7,3) - t494 * t713 + t7 * (-mrSges(7,1) * t189 - mrSges(7,2) * t494) + t785 * t120 + (Ifges(5,1) * t270 + Ifges(5,4) * t269) * t694 + (Ifges(4,1) * t386 + Ifges(4,4) * t385) * t674 + (Ifges(4,5) * t386 + Ifges(4,6) * t385) * t665 + t740 * t375 + m(4) * (t146 * t260 + t147 * t259 + t194 * t237 + t195 * t236 + t257 * t357 + t319 * t375) + m(5) * (t135 * t68 + t136 * t67 + t151 * t42 + t152 * t41 + t186 * t280 + t252 * t255) + m(7) * (t17 * t59 + t2 * t39 + t29 * t5 + t3 * t38 + t30 * t4 + t69 * t7) + m(6) * (t10 * t71 + t11 * t72 + t148 * t188 + t19 * t61 + t20 * t62 + t208 * t99) + (t146 * t385 - t147 * t386 - t236 * t304 + t237 * t303) * mrSges(4,3) + t189 * t714 + (Ifges(7,4) * t91 + Ifges(7,2) * t92) * t689 + (Ifges(5,5) * t270 + Ifges(5,6) * t269) * t666 + (Ifges(5,4) * t270 + Ifges(5,2) * t269) * t693 + (-t135 * t183 + t136 * t184 + t269 * t41 - t270 * t42) * mrSges(5,3) + (Ifges(7,1) * t91 + Ifges(7,4) * t92) * t687 + (-pkin(1) * t391 * t461 + Ifges(2,3)) * qJDD(1) + (-m(4) * (pkin(2) * t390 + t575) - t307 * mrSges(4,1) - t306 * mrSges(4,2) - m(6) * t502 - t290 * mrSges(6,1) - m(7) * (pkin(5) * t290 + t502) - t233 * mrSges(7,1) - t232 * mrSges(7,2) - m(3) * t575 - t390 * mrSges(3,1) - mrSges(2,1) * t473 - m(5) * (t390 * t452 + t575) - t299 * mrSges(5,1) - t298 * mrSges(5,2) + (t461 * t524 + mrSges(2,2)) * t468 + t525 * t289 + t480 * t389) * g(2) + t272 * (mrSges(3,1) * t463 - mrSges(3,3) * t597) + (Ifges(7,5) * t91 + Ifges(7,6) * t92) * t681 + t340 * (Ifges(4,4) * t304 + Ifges(4,2) * t303) / 0.2e1 + (t395 * mrSges(3,1) - t396 * mrSges(3,2) + Ifges(3,3) * t662 + (Ifges(3,5) * t467 + Ifges(3,6) * t472) * t461) * t438 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t717 + t271 * t396 + t272 * t395 - t369 * t375 + t373 * t374) + ((t611 / 0.2e1 + t317 / 0.2e1 - t369 * mrSges(3,3) + (-t699 + t633 / 0.2e1) * t574) * t472 + (-t721 + t613 / 0.2e1 + t622 / 0.2e1 + t615 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t411 + t612 / 0.2e1 + t618 / 0.2e1 + t619 / 0.2e1 + t112 / 0.2e1 + t623 / 0.2e1 + t229 / 0.2e1 - t610 / 0.2e1 + t170 / 0.2e1 - t373 * mrSges(3,3) - t237 * mrSges(4,2) + t236 * mrSges(4,1) - t62 * mrSges(6,2) + (-t700 - t634 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t472) * t574) * t467) * t571 + t419 * (Ifges(4,5) * t304 + Ifges(4,6) * t303) / 0.2e1 + (Ifges(4,4) * t386 + Ifges(4,2) * t385) * t673 + t270 * t707 + t269 * t708 + t386 * t691 + t385 * t692 + t304 * t675 + t303 * t676 + t183 * t683 + t184 * t685 + t545 * t662 + (Ifges(5,5) * t183 + Ifges(5,6) * t184) * t663 + (Ifges(4,1) * t304 + Ifges(4,4) * t303) * t667 + (Ifges(5,1) * t183 + Ifges(5,4) * t184) * t669 + (Ifges(5,4) * t183 + Ifges(5,2) * t184) * t671 - t257 * t515 + t374 * t364 + t357 * t173 + t319 * (-mrSges(4,1) * t303 + mrSges(4,2) * t304) + t280 * t90 + t194 * t283 + t195 * t284 + t186 * (-mrSges(5,1) * t269 + mrSges(5,2) * t270) + t252 * (-mrSges(5,1) * t184 + mrSges(5,2) * t183) + t255 * t187 + t259 * t206 + t260 * t207 + t67 * t221 + t68 * t222 + t208 * t37 + t20 * t167 + t19 * t168 + t151 * t126 + t152 * t127 + t148 * t118 + t760 * t200 + (mrSges(2,1) * t468 - t427 * mrSges(4,1) - t415 * mrSges(5,1) + t732 * t388 + (mrSges(2,2) + (-mrSges(5,2) * t456 + t524 - t642) * t461) * t473 + t525 * t285 + t786 * t286 + t737 * t387 - t757 * (-t387 * t458 + t388 * t407 - t413 * t593 + t531) + (m(3) + m(4) + m(5)) * t531) * g(1) + t776 * t119 + (-t395 * mrSges(3,3) + Ifges(3,5) * t662 + (t467 * Ifges(3,1) + t633 - t699) * t461) * t377 + (t396 * mrSges(3,3) + Ifges(3,6) * t662 + (t508 + t700) * t461) * t376 + t778 * t199 + t38 * t27 + t39 * t28 - (t546 + t767) * t594 / 0.2e1 + t69 * t21 + t72 * t73 + t71 * t74 + t91 * t81 / 0.2e1 - t463 * t771 + t92 * t80 / 0.2e1 + t59 * (-mrSges(7,1) * t92 + mrSges(7,2) * t91) + t17 * t101 + t4 * t105 + t5 * t106 + (mrSges(3,3) * t271 - Ifges(4,5) * t674 - Ifges(5,5) * t694 - Ifges(6,5) * t701 - Ifges(4,6) * t673 - Ifges(5,6) * t693 - Ifges(6,6) * t702 - Ifges(4,3) * t665 - t666 * t773 + t764 + t765) * t594; -((t229 + t170 + t112) * t467 + (-Ifges(3,2) * t543 + t471 * t231 + t317 + t432) * t472 + t419 * (Ifges(4,3) * t467 + t472 * t505) + t341 * (Ifges(4,5) * t467 + t472 * t510) + t340 * (Ifges(4,6) * t467 + t472 * t507) + t439 * (Ifges(3,5) * t472 - Ifges(3,6) * t467)) * t574 / 0.2e1 + (t548 + t785) * t217 + (t757 * (-t389 * t407 - t390 * t458) - t737 * t390 + t725 * t389) * g(1) + (t757 * (-t387 * t407 - t388 * t458) - t737 * t388 + t725 * t387) * g(2) + t419 * t319 * (mrSges(4,1) * t466 + t642) - (Ifges(5,4) * t669 + Ifges(5,2) * t671 + Ifges(5,6) * t663 + t639 + t685) * t302 + (-t621 + t675) * t568 + (-Ifges(7,5) * t487 - Ifges(7,6) * t488) * t681 + (Ifges(7,5) * t212 + Ifges(7,6) * t211) * t682 + t752 * t167 + t753 * t101 + (-t10 * t197 + t11 * t198 + t188 * t739 + t352 * t99 + t61 * t754 + t62 * t752) * m(6) + t754 * t168 + (-t284 * t568 - t283 * t569 + t471 * t207 + m(4) * ((-t236 * t471 - t237 * t466) * qJD(3) + t734) - t466 * t206) * pkin(9) + t734 * mrSges(4,3) + (-pkin(2) * t257 - t236 * t263 - t237 * t264) * m(4) + (-t620 - t230 / 0.2e1) * t569 + (-t188 * t226 + t543 * t62) * mrSges(6,2) + t738 * t187 + t739 * t118 + (-m(4) * t319 + t523 - t740) * t373 + t741 * t222 + (t135 * t741 + t136 * t742 - t186 * t452 + t252 * t738 + t323 * t42 + t324 * t41) * m(5) + t742 * t221 - t744 * t80 / 0.2e1 + (mrSges(7,1) * t744 - mrSges(7,2) * t743) * t59 + (-t2 * t603 + t29 * t743 - t3 * t602 - t30 * t744) * mrSges(7,3) + (Ifges(5,1) * t669 + Ifges(5,4) * t671 + Ifges(5,5) * t663 - t640 + t683) * t301 + (-Ifges(7,4) * t487 - Ifges(7,2) * t488) * t689 + (Ifges(7,4) * t212 + Ifges(7,2) * t211) * t690 + (t21 - t74) * t197 + t602 * t713 + (t722 + t646) * t226 + (-t364 + t522) * t369 + (Ifges(7,1) * t212 + Ifges(7,4) * t211) * t688 + (-Ifges(7,1) * t487 - Ifges(7,4) * t488) * t687 + (-t237 * (-mrSges(4,3) * t466 * t472 - mrSges(4,2) * t467) - t236 * (mrSges(4,1) * t467 - mrSges(4,3) * t587)) * t574 - t15 * t603 / 0.2e1 + (t340 * t507 + t341 * t510 + t419 * t505) * qJD(3) / 0.2e1 + (t789 * mrSges(5,1) + mrSges(5,2) * t792) * t252 + (Ifges(5,1) * t322 - Ifges(5,4) * t321) * t670 + t755 * t106 + (t131 * t3 + t132 * t2 + t197 * t7 + t29 * t755 + t30 * t756 + t59 * t753) * m(7) + t756 * t105 + (-t467 * (Ifges(3,1) * t472 - t634) / 0.2e1 + pkin(1) * (mrSges(3,1) * t467 + mrSges(3,2) * t472)) * qJD(1) ^ 2 * t717 - t771 + (Ifges(5,5) * t322 - Ifges(5,6) * t321) * t664 + (Ifges(5,4) * t322 - Ifges(5,2) * t321) * t672 + (t135 * t322 + t136 * t321 + t399 * t41 - t400 * t42) * mrSges(5,3) + t400 * t707 + t399 * t708 + t466 * t691 + t471 * t692 + (Ifges(5,4) * t400 + Ifges(5,2) * t399) * t693 + (Ifges(5,1) * t400 + Ifges(5,4) * t399) * t694 + (Ifges(4,2) * t471 + t632) * t673 + (Ifges(4,1) * t466 + t631) * t674 + t520 * t676 + t322 * t684 - t321 * t686 + (Ifges(4,5) * t466 + Ifges(4,6) * t471) * t665 + (Ifges(5,5) * t400 + Ifges(5,6) * t399) * t666 + t257 * t514 + t186 * (-mrSges(5,1) * t399 + mrSges(5,2) * t400) + t352 * t37 + t323 * t126 + t324 * t127 - t264 * t283 - t263 * t284 + t272 * mrSges(3,1) - t212 * t81 / 0.2e1 + t198 * t73 - pkin(2) * t173 + t131 * t27 + t132 * t28 - t452 * t90 + (t504 * t703 + t506 * t709 + t509 * t710 + t511 * t7 + t532 * t81 + t760) * t297 + t776 * t216 + t778 * t296 + t545 + t719 * t225 + (t757 * t407 * t594 + t391 + ((t736 - t768) * t472 + (-t458 * t757 - t511 + t770) * t467) * t461) * g(3) + (Ifges(5,5) * t670 + Ifges(6,5) * t678 + Ifges(5,6) * t672 + Ifges(6,6) * t680 + t773 * t664 + t721) * t543; -(mrSges(5,1) * t252 + Ifges(5,4) * t670 + Ifges(5,2) * t672 + Ifges(5,6) * t664 - t639 + t686) * t251 + (m(6) * t62 + m(7) * t503 - t745) * (t462 * t470 - t598) * t626 + (-m(7) * (t516 + t530) - m(6) * t516 - (-t388 * t471 + t427) * mrSges(4,2) + t750 * (-t388 * t466 - t471 * t593) + t729) * g(2) + (-m(7) * (t529 + t579) - m(6) * t579 + mrSges(4,2) * t307 + t750 * t306 + t730) * g(1) + (-t29 * t31 - t30 * t32 + t365 * t7) * m(7) + (t10 * t378 + t11 * t379 - t188 * t210 - t62 * t66) * m(6) + (t127 * t465 + t221 * t566 - t222 * t567) * pkin(3) + (t41 * t465 + t42 * t470 + (-t135 * t465 + t136 * t470) * qJD(4)) * t716 - t419 * (Ifges(4,5) * t340 - Ifges(4,6) * t341) / 0.2e1 - t187 * t657 - m(5) * (t135 * t142 + t136 * t143 + t252 * t657) - (-Ifges(4,2) * t341 + t231 + t332) * t340 / 0.2e1 + t724 * (pkin(11) + t379) + (-mrSges(5,2) * t252 + Ifges(5,1) * t670 + Ifges(5,4) * t672 + Ifges(5,5) * t664 + t640 + t684) * t526 - t765 + (-m(5) * t656 - m(7) * (t528 + t576) - m(6) * t576 - t515 + t731) * g(3) + t230 * t667 + t774 + t378 * t74 + t379 * t73 + t365 * t21 - t319 * (mrSges(4,1) * t341 + mrSges(4,2) * t340) - t236 * t283 + t237 * t284 - t143 * t221 - t142 * t222 - t210 * t118 - t66 * t167 - t341 * (Ifges(4,1) * t340 - t614) / 0.2e1 + t341 * t620 + t340 * t621 + t126 * t654 + t546 + (-m(6) * t61 + m(7) * t59 - t746) * (-t121 * t460 + t462 * t491 + (t460 * t470 + t592) * t626) - t32 * t105 - t31 * t106; t746 * t63 + (-t221 + t638) * t135 + t724 * (pkin(11) + t651) + (Ifges(5,5) * t526 - Ifges(5,6) * t251) * t664 + (Ifges(5,1) * t526 - t630) * t670 - t252 * (mrSges(5,1) * t251 + mrSges(5,2) * t526) + ((t10 * t462 + t11 * t460) * pkin(4) - t188 * t652 + t61 * t63 - t62 * t64) * m(6) - t118 * t652 + (-t29 * t33 - t30 * t34 + t449 * t7 - t59 * t63) * m(7) + (-m(7) * (t481 + t530) - m(6) * t481 + t729) * g(2) + (-m(7) * (t517 + t529) - m(6) * t517 + t730) * g(1) + (-m(7) * (t499 + t528) - m(6) * t499 + t731) * g(3) + t171 * t669 + (-Ifges(5,2) * t251 + t172 + t246) * t672 + t774 - t64 * t167 + t449 * t21 + (t222 + t637) * t136 + t74 * t650 + t73 * t651 - t34 * t105 - t33 * t106; t469 * t27 + t464 * t28 + t746 * t762 + t501 * qJD(6) + t745 * t763 + t37 + (t175 * t503 + t2 * t464 + t3 * t469 - t59 * t762 + t479) * m(7) + (t61 * t762 - t62 * t763 + t479 + t99) * m(6); -t59 * (mrSges(7,1) * t162 + mrSges(7,2) * t161) + (Ifges(7,1) * t161 - t629) * t688 + t80 * t687 + (Ifges(7,5) * t161 - Ifges(7,6) * t162) * t682 - t29 * t105 + t30 * t106 - g(1) * (mrSges(7,1) * t232 - mrSges(7,2) * t233) - g(2) * ((-t286 * t464 + t387 * t469) * mrSges(7,1) + (-t286 * t469 - t387 * t464) * mrSges(7,2)) - g(3) * ((-t334 * t464 - t551) * mrSges(7,1) + (-t334 * t469 + t552) * mrSges(7,2)) + (t161 * t29 + t162 * t30) * mrSges(7,3) + t14 + (-Ifges(7,2) * t162 + t160 + t81) * t690 + t727;];
tau  = t1;
