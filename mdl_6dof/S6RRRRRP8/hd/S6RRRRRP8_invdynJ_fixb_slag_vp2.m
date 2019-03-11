% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:19
% EndTime: 2019-03-10 01:49:08
% DurationCPUTime: 64.58s
% Computational Cost: add. (30773->1117), mult. (72996->1466), div. (0->0), fcn. (58676->14), ass. (0->491)
t472 = cos(qJ(2));
t467 = sin(qJ(2));
t462 = sin(pkin(6));
t613 = qJD(1) * t462;
t574 = t467 * t613;
t463 = cos(pkin(6));
t612 = qJD(1) * t463;
t599 = pkin(1) * t612;
t384 = -pkin(8) * t574 + t472 * t599;
t497 = (pkin(2) * t467 - pkin(9) * t472) * t462;
t385 = qJD(1) * t497;
t466 = sin(qJ(3));
t471 = cos(qJ(3));
t260 = -t384 * t466 + t471 * t385;
t473 = -pkin(10) - pkin(9);
t578 = qJD(3) * t473;
t622 = t471 * t472;
t846 = -(pkin(3) * t467 - pkin(10) * t622) * t613 - t260 + t471 * t578;
t261 = t471 * t384 + t466 * t385;
t573 = t472 * t613;
t540 = t466 * t573;
t845 = -pkin(10) * t540 - t466 * t578 + t261;
t465 = sin(qJ(4));
t470 = cos(qJ(4));
t406 = t465 * t471 + t466 * t470;
t761 = qJD(3) + qJD(4);
t306 = t761 * t406;
t319 = t406 * t573;
t835 = t306 - t319;
t801 = Ifges(6,1) + Ifges(7,1);
t799 = Ifges(6,5) + Ifges(7,4);
t429 = t473 * t466;
t430 = t473 * t471;
t770 = t470 * t429 + t430 * t465;
t778 = qJD(4) * t770 + t465 * t846 - t845 * t470;
t442 = pkin(8) * t573;
t387 = t467 * t599 + t442;
t608 = qJD(3) * t466;
t763 = -t387 + (-t540 + t608) * pkin(3);
t445 = qJD(2) + t612;
t349 = t445 * t471 - t466 * t574;
t609 = qJD(2) * t472;
t391 = (qJD(1) * t609 + qJDD(1) * t467) * t462;
t601 = qJDD(1) * t463;
t444 = qJDD(2) + t601;
t228 = qJD(3) * t349 + t391 * t471 + t444 * t466;
t350 = t445 * t466 + t471 * t574;
t229 = -qJD(3) * t350 - t391 * t466 + t444 * t471;
t506 = t349 * t465 + t470 * t350;
t127 = -qJD(4) * t506 - t228 * t465 + t229 * t470;
t124 = qJDD(5) - t127;
t735 = t124 / 0.2e1;
t551 = t470 * t349 - t350 * t465;
t126 = qJD(4) * t551 + t228 * t470 + t229 * t465;
t464 = sin(qJ(5));
t469 = cos(qJ(5));
t424 = qJD(3) - t573;
t498 = -qJD(4) - t424;
t201 = t464 * t506 + t469 * t498;
t610 = qJD(2) * t462;
t572 = t467 * t610;
t602 = qJDD(1) * t462;
t390 = -qJD(1) * t572 + t472 * t602;
t371 = qJDD(3) - t390;
t363 = qJDD(4) + t371;
t78 = -qJD(5) * t201 + t469 * t126 + t464 * t363;
t742 = t78 / 0.2e1;
t844 = t799 * t735 + t801 * t742;
t202 = -t464 * t498 + t469 * t506;
t79 = qJD(5) * t202 + t464 * t126 - t469 * t363;
t740 = t79 / 0.2e1;
t800 = -Ifges(6,4) + Ifges(7,5);
t827 = Ifges(6,6) - Ifges(7,6);
t797 = Ifges(6,3) + Ifges(7,2);
t802 = mrSges(6,3) + mrSges(7,2);
t831 = mrSges(5,2) - t802;
t843 = pkin(11) * t574 - t778;
t405 = t465 * t466 - t470 * t471;
t305 = t761 * t405;
t320 = t405 * t573;
t842 = t763 + (t305 - t320) * pkin(11) + t835 * pkin(4);
t632 = t462 * t467;
t446 = pkin(8) * t632;
t544 = qJD(2) * t599;
t594 = pkin(1) * t601;
t268 = -qJD(2) * t442 - qJDD(1) * t446 - t467 * t544 + t472 * t594;
t251 = -pkin(2) * t444 - t268;
t170 = -pkin(3) * t229 + t251;
t318 = pkin(9) * t445 + t387;
t368 = (-pkin(2) * t472 - pkin(9) * t467 - pkin(1)) * t462;
t332 = qJD(1) * t368;
t224 = -t318 * t466 + t471 * t332;
t189 = -pkin(10) * t350 + t224;
t175 = pkin(3) * t424 + t189;
t225 = t318 * t471 + t332 * t466;
t190 = pkin(10) * t349 + t225;
t605 = qJD(4) * t470;
t606 = qJD(4) * t465;
t267 = pkin(8) * t390 + t467 * t594 + t472 * t544;
t250 = pkin(9) * t444 + t267;
t255 = -pkin(1) * t602 - pkin(2) * t390 - pkin(9) * t391;
t131 = -qJD(3) * t225 - t250 * t466 + t471 * t255;
t87 = pkin(3) * t371 - pkin(10) * t228 + t131;
t607 = qJD(3) * t471;
t130 = t471 * t250 + t466 * t255 - t318 * t608 + t332 * t607;
t94 = pkin(10) * t229 + t130;
t30 = t175 * t605 - t190 * t606 + t465 * t87 + t470 * t94;
t703 = t363 / 0.2e1;
t733 = t127 / 0.2e1;
t734 = t126 / 0.2e1;
t741 = -t79 / 0.2e1;
t240 = qJD(5) - t551;
t182 = t470 * t190;
t117 = t465 * t175 + t182;
t111 = -pkin(11) * t498 + t117;
t317 = -pkin(2) * t445 - t384;
t246 = -pkin(3) * t349 + t317;
t133 = -pkin(4) * t551 - pkin(11) * t506 + t246;
t27 = pkin(11) * t363 + t30;
t35 = -pkin(4) * t127 - pkin(11) * t126 + t170;
t603 = qJD(5) * t469;
t604 = qJD(5) * t464;
t6 = -t111 * t604 + t133 * t603 + t469 * t27 + t464 * t35;
t2 = qJ(6) * t124 + qJD(6) * t240 + t6;
t53 = t111 * t469 + t133 * t464;
t7 = -qJD(5) * t53 - t27 * t464 + t35 * t469;
t4 = -pkin(5) * t124 + qJDD(6) - t7;
t751 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t796 = t124 * t797 + t78 * t799 - t79 * t827;
t841 = t751 + Ifges(6,6) * t741 + Ifges(7,6) * t740 + t170 * mrSges(5,1) + t735 * t797 + t742 * t799 + t796 / 0.2e1 - t30 * mrSges(5,3) + (-t363 / 0.2e1 - t703) * Ifges(5,6) + (-t127 / 0.2e1 - t733) * Ifges(5,2) + (-t126 / 0.2e1 - t734) * Ifges(5,4);
t31 = -t175 * t606 - t190 * t605 - t465 * t94 + t470 * t87;
t840 = mrSges(5,2) * t170 - mrSges(5,3) * t31 + 0.2e1 * Ifges(5,1) * t734 + 0.2e1 * Ifges(5,4) * t733 + 0.2e1 * Ifges(5,5) * t703;
t676 = Ifges(6,4) * t202;
t107 = -t201 * Ifges(6,2) + t240 * Ifges(6,6) + t676;
t737 = -t107 / 0.2e1;
t838 = t800 * t740 + t844;
t837 = mrSges(6,1) + mrSges(7,1);
t836 = mrSges(6,2) - mrSges(7,3);
t105 = t202 * Ifges(6,5) - t201 * Ifges(6,6) + t240 * Ifges(6,3);
t106 = t202 * Ifges(7,4) + t240 * Ifges(7,2) + t201 * Ifges(7,6);
t823 = t106 + t105;
t199 = Ifges(6,4) * t201;
t672 = Ifges(7,5) * t201;
t822 = t202 * t801 + t799 * t240 - t199 + t672;
t322 = t429 * t465 - t430 * t470;
t779 = -qJD(4) * t322 + t845 * t465 + t470 * t846;
t52 = -t111 * t464 + t133 * t469;
t783 = qJD(6) - t52;
t44 = -pkin(5) * t240 + t783;
t535 = t2 * t469 + t4 * t464;
t834 = t44 * t603 + t535;
t45 = qJ(6) * t240 + t53;
t479 = Ifges(5,6) * t498;
t833 = t117 * mrSges(5,3) - t246 * mrSges(5,1) - t52 * mrSges(6,1) + t44 * mrSges(7,1) + t53 * mrSges(6,2) - t45 * mrSges(7,3) - t479 / 0.2e1;
t181 = t465 * t190;
t116 = t470 * t175 - t181;
t480 = Ifges(5,5) * t498;
t832 = t480 / 0.2e1 - t246 * mrSges(5,2) + t116 * mrSges(5,3) - t464 * t737;
t677 = Ifges(5,4) * t506;
t164 = Ifges(5,2) * t551 - t479 + t677;
t830 = t164 / 0.2e1 + t833;
t393 = t463 * t471 - t466 * t632;
t504 = pkin(3) * t393;
t456 = pkin(3) * t471 + pkin(2);
t283 = pkin(4) * t405 - pkin(11) * t406 - t456;
t793 = t283 * t603 - t322 * t604 + t464 * t842 - t469 * t843;
t772 = t464 * t283 + t469 * t322;
t792 = -qJD(5) * t772 + t464 * t843 + t469 * t842;
t825 = t267 * mrSges(3,2);
t824 = t31 * mrSges(5,1) - t30 * mrSges(5,2);
t128 = t189 * t465 + t182;
t821 = pkin(3) * t606 - t128;
t780 = pkin(4) * t574 - t779;
t516 = pkin(5) * t464 - qJ(6) * t469;
t820 = pkin(5) * t604 - qJ(6) * t603 - qJD(6) * t464 - t516 * t551;
t110 = pkin(4) * t498 - t116;
t528 = mrSges(7,1) * t464 - mrSges(7,3) * t469;
t529 = mrSges(6,1) * t464 + mrSges(6,2) * t469;
t58 = t201 * pkin(5) - t202 * qJ(6) + t110;
t819 = t110 * t529 + t58 * t528;
t697 = cos(qJ(1));
t575 = t697 * t472;
t468 = sin(qJ(1));
t625 = t467 * t468;
t398 = -t463 * t625 + t575;
t630 = t462 * t471;
t311 = -t398 * t466 + t468 * t630;
t818 = t131 * mrSges(4,1) - t130 * mrSges(4,2);
t534 = -t464 * t7 + t469 * t6;
t711 = -t506 / 0.2e1;
t713 = -t551 / 0.2e1;
t716 = -t240 / 0.2e1;
t723 = -t202 / 0.2e1;
t724 = t201 / 0.2e1;
t725 = -t201 / 0.2e1;
t817 = -Ifges(5,4) * t711 - Ifges(5,2) * t713 + Ifges(6,6) * t724 + Ifges(7,6) * t725 + t716 * t797 + t723 * t799 + t830;
t198 = Ifges(7,5) * t202;
t104 = t240 * Ifges(7,6) + t201 * Ifges(7,3) + t198;
t627 = t464 * t104;
t698 = -t469 / 0.2e1;
t239 = Ifges(5,4) * t551;
t165 = Ifges(5,1) * t506 + t239 - t480;
t728 = -t165 / 0.2e1;
t816 = -t627 / 0.2e1 + t728 - t819 + t822 * t698 + t832;
t36 = -mrSges(7,2) * t79 + mrSges(7,3) * t124;
t37 = mrSges(6,1) * t124 - mrSges(6,3) * t78;
t38 = -t124 * mrSges(7,1) + t78 * mrSges(7,2);
t39 = -mrSges(6,2) * t124 - mrSges(6,3) * t79;
t815 = (t36 + t39) * t469 + (-t37 + t38) * t464;
t576 = t697 * t467;
t624 = t468 * t472;
t396 = t463 * t576 + t624;
t461 = qJ(3) + qJ(4);
t458 = sin(t461);
t459 = cos(t461);
t577 = t462 * t697;
t299 = t396 * t459 - t458 * t577;
t395 = -t463 * t575 + t625;
t230 = t299 * t464 - t395 * t469;
t813 = t299 * t469 + t395 * t464;
t812 = -t104 / 0.2e1 - t58 * mrSges(7,1);
t19 = t78 * Ifges(7,5) + t124 * Ifges(7,6) + t79 * Ifges(7,3);
t22 = t78 * Ifges(6,4) - t79 * Ifges(6,2) + t124 * Ifges(6,6);
t811 = -mrSges(7,2) * t2 - mrSges(6,3) * t6 - t22 / 0.2e1 + t19 / 0.2e1;
t810 = m(7) + m(6);
t719 = t228 / 0.2e1;
t718 = t229 / 0.2e1;
t702 = t371 / 0.2e1;
t807 = pkin(5) * t506;
t795 = qJ(6) * t835 + qJD(6) * t405 + t793;
t794 = -pkin(5) * t835 - t792;
t512 = t464 * t53 + t469 * t52;
t790 = t512 * mrSges(6,3);
t789 = m(4) * pkin(9) + mrSges(4,3) + mrSges(5,3);
t788 = t820 + t821;
t787 = -t117 + t820;
t262 = -t320 * t464 - t469 * t574;
t263 = -t320 * t469 + t464 * t574;
t517 = pkin(5) * t469 + qJ(6) * t464;
t786 = -pkin(5) * t262 + qJ(6) * t263 - t516 * t305 + (qJD(5) * t517 - qJD(6) * t469) * t406 + t780;
t785 = qJ(6) * t506;
t629 = t462 * t472;
t402 = t463 * t467 * pkin(1) + pkin(8) * t629;
t367 = pkin(9) * t463 + t402;
t253 = -t367 * t466 + t471 * t368;
t394 = t463 * t466 + t467 * t630;
t197 = -pkin(3) * t629 - pkin(10) * t394 + t253;
t254 = t471 * t367 + t466 * t368;
t213 = pkin(10) * t393 + t254;
t138 = t465 * t197 + t470 * t213;
t136 = -pkin(11) * t629 + t138;
t265 = t393 * t465 + t394 * t470;
t696 = pkin(1) * t472;
t366 = t446 + (-pkin(2) - t696) * t463;
t270 = t366 - t504;
t505 = t470 * t393 - t394 * t465;
t162 = -pkin(4) * t505 - pkin(11) * t265 + t270;
t782 = t469 * t136 + t464 * t162;
t141 = mrSges(6,1) * t201 + mrSges(6,2) * t202;
t212 = -mrSges(5,1) * t498 - mrSges(5,3) * t506;
t781 = t141 - t212;
t546 = mrSges(3,3) * t574;
t777 = -mrSges(3,1) * t445 - mrSges(4,1) * t349 + mrSges(4,2) * t350 + t546;
t570 = t406 * t603;
t495 = -t305 * t464 + t570;
t776 = t262 - t495;
t648 = t305 * t469;
t494 = t406 * t604 + t648;
t775 = t263 + t494;
t298 = -t396 * t458 - t459 * t577;
t651 = t298 * t469;
t652 = t298 * t464;
t774 = -pkin(5) * t651 - qJ(6) * t652;
t631 = t462 * t468;
t302 = t398 * t458 - t459 * t631;
t649 = t302 * t469;
t650 = t302 * t464;
t773 = -pkin(5) * t649 - qJ(6) * t650;
t364 = -t458 * t632 + t459 * t463;
t645 = t364 * t469;
t646 = t364 * t464;
t771 = -pkin(5) * t645 - qJ(6) * t646;
t401 = t463 * t696 - t446;
t519 = Ifges(6,5) * t469 - Ifges(6,6) * t464;
t521 = Ifges(7,4) * t469 + Ifges(7,6) * t464;
t769 = t519 + t521;
t671 = Ifges(7,5) * t464;
t525 = Ifges(7,1) * t469 + t671;
t675 = Ifges(6,4) * t464;
t526 = Ifges(6,1) * t469 - t675;
t768 = t526 + t525;
t531 = -mrSges(4,1) * t471 + mrSges(4,2) * t466;
t767 = -m(4) * pkin(2) - t459 * mrSges(5,1) + mrSges(5,2) * t458 + t531;
t694 = pkin(3) * t465;
t454 = pkin(11) + t694;
t596 = pkin(3) * t605;
t766 = -t454 * t604 + t469 * t596;
t765 = -t454 * t603 - t464 * t596;
t764 = -t52 * t603 - t53 * t604;
t760 = mrSges(3,2) - t789;
t759 = mrSges(3,1) - t767;
t365 = t458 * t463 + t459 * t632;
t758 = -t364 * mrSges(5,1) + t365 * t831 - t645 * t837 + t646 * t836;
t303 = t398 * t459 + t458 * t631;
t757 = t302 * mrSges(5,1) + t303 * t831 + t649 * t837 - t650 * t836;
t756 = -t298 * mrSges(5,1) + t299 * t831 - t651 * t837 + t652 * t836;
t386 = qJD(2) * t497;
t388 = t401 * qJD(2);
t177 = -qJD(3) * t254 + t471 * t386 - t388 * t466;
t571 = t462 * t609;
t308 = qJD(3) * t393 + t471 * t571;
t149 = pkin(3) * t572 - pkin(10) * t308 + t177;
t176 = -t367 * t608 + t368 * t607 + t466 * t386 + t471 * t388;
t307 = -qJD(3) * t394 - t466 * t571;
t160 = pkin(10) * t307 + t176;
t46 = t465 * t149 + t470 * t160 + t197 * t605 - t213 * t606;
t42 = pkin(11) * t572 + t46;
t168 = qJD(4) * t505 + t307 * t465 + t308 * t470;
t169 = qJD(4) * t265 - t470 * t307 + t308 * t465;
t389 = t402 * qJD(2);
t249 = -pkin(3) * t307 + t389;
t81 = pkin(4) * t169 - pkin(11) * t168 + t249;
t13 = -qJD(5) * t782 - t42 * t464 + t469 * t81;
t549 = m(7) * pkin(5) + t837;
t541 = -m(7) * qJ(6) + t836;
t172 = pkin(4) * t506 - pkin(11) * t551;
t747 = t462 ^ 2;
t739 = pkin(1) * mrSges(3,1);
t738 = pkin(1) * mrSges(3,2);
t732 = Ifges(4,4) * t719 + Ifges(4,2) * t718 + Ifges(4,6) * t702;
t731 = Ifges(4,1) * t719 + Ifges(4,4) * t718 + Ifges(4,5) * t702;
t730 = -t164 / 0.2e1;
t727 = t165 / 0.2e1;
t722 = t202 / 0.2e1;
t662 = t350 * Ifges(4,4);
t220 = t349 * Ifges(4,2) + t424 * Ifges(4,6) + t662;
t721 = t220 / 0.2e1;
t333 = Ifges(4,4) * t349;
t221 = t350 * Ifges(4,1) + t424 * Ifges(4,5) + t333;
t720 = t221 / 0.2e1;
t715 = t240 / 0.2e1;
t712 = t551 / 0.2e1;
t710 = t506 / 0.2e1;
t704 = t350 / 0.2e1;
t700 = t463 / 0.2e1;
t695 = pkin(3) * t350;
t693 = pkin(3) * t470;
t692 = pkin(4) * t459;
t684 = mrSges(4,3) * t349;
t683 = mrSges(6,3) * t201;
t682 = mrSges(6,3) * t202;
t681 = Ifges(3,4) * t467;
t680 = Ifges(3,4) * t472;
t679 = Ifges(4,4) * t466;
t678 = Ifges(4,4) * t471;
t674 = Ifges(6,4) * t469;
t673 = Ifges(5,5) * t506;
t670 = Ifges(7,5) * t469;
t669 = Ifges(5,6) * t551;
t664 = t225 * mrSges(4,3);
t663 = t349 * Ifges(4,6);
t661 = t350 * Ifges(4,5);
t660 = t424 * Ifges(4,3);
t659 = t445 * Ifges(3,5);
t658 = t445 * Ifges(3,6);
t657 = t45 * t464;
t656 = t130 * t471;
t655 = t131 * t466;
t654 = t551 * t464;
t653 = t551 * t469;
t644 = t395 * t458;
t397 = t463 * t624 + t576;
t641 = t397 * t458;
t638 = t406 * t469;
t634 = t459 * t464;
t633 = t459 * t469;
t623 = t469 * t472;
t64 = t469 * t116 + t464 * t172;
t129 = t189 * t470 - t181;
t150 = t172 + t695;
t62 = t469 * t129 + t464 * t150;
t153 = -mrSges(7,2) * t201 + mrSges(7,3) * t240;
t154 = -mrSges(6,2) * t240 - t683;
t621 = -t153 - t154;
t155 = mrSges(6,1) * t240 - t682;
t156 = -mrSges(7,1) * t240 + mrSges(7,2) * t202;
t620 = -t155 + t156;
t617 = -t395 * t456 - t396 * t473;
t616 = -t397 * t456 - t398 * t473;
t614 = t697 * pkin(1) + pkin(8) * t631;
t448 = pkin(4) * t629;
t592 = t45 * t604;
t589 = t458 * t629;
t587 = t466 * t631;
t585 = t462 * t623;
t434 = t464 * t629;
t584 = Ifges(5,5) * t126 + Ifges(5,6) * t127 + Ifges(5,3) * t363;
t583 = Ifges(4,5) * t228 + Ifges(4,6) * t229 + Ifges(4,3) * t371;
t579 = Ifges(3,5) * t391 + Ifges(3,6) * t390 + Ifges(3,3) * t444;
t563 = t627 / 0.2e1;
t557 = -t604 / 0.2e1;
t556 = t603 / 0.2e1;
t555 = -pkin(1) * t468 + pkin(8) * t577;
t554 = t298 * pkin(4) + t299 * pkin(11);
t553 = -t302 * pkin(4) + pkin(11) * t303;
t552 = t364 * pkin(4) + pkin(11) * t365;
t137 = t197 * t470 - t465 * t213;
t435 = t466 * t577;
t550 = -t396 * t471 + t435;
t545 = mrSges(3,3) * t573;
t536 = t311 * pkin(3);
t135 = -t137 + t448;
t532 = mrSges(4,1) * t393 - mrSges(4,2) * t394;
t527 = Ifges(4,1) * t471 - t679;
t524 = Ifges(3,2) * t472 + t681;
t523 = -Ifges(4,2) * t466 + t678;
t522 = -Ifges(6,2) * t464 + t674;
t520 = Ifges(4,5) * t471 - Ifges(4,6) * t466;
t518 = Ifges(7,3) * t464 + t670;
t513 = t44 * t469 - t657;
t511 = pkin(3) * t587 - t397 * t473 + t398 * t456 + t614;
t63 = -t116 * t464 + t172 * t469;
t61 = -t129 * t464 + t150 * t469;
t71 = -t136 * t464 + t162 * t469;
t206 = t283 * t469 - t322 * t464;
t417 = -pkin(4) - t517;
t47 = t149 * t470 - t465 * t160 - t197 * t606 - t213 * t605;
t500 = -pkin(11) * t810 + t831;
t237 = t265 * t464 + t585;
t493 = pkin(3) * t435 + t395 * t473 - t396 * t456 + t555;
t491 = t201 * t522;
t490 = t201 * t518;
t489 = t202 * t526;
t488 = t202 * t525;
t487 = t240 * t521;
t486 = t240 * t519;
t12 = -t136 * t604 + t162 * t603 + t469 * t42 + t464 * t81;
t482 = t396 * t466 + t471 * t577;
t478 = t536 + t553;
t28 = -pkin(4) * t363 - t31;
t477 = t482 * pkin(3);
t476 = t504 + t552;
t43 = -pkin(4) * t572 - t47;
t475 = -t477 + t554;
t9 = pkin(5) * t79 - qJ(6) * t78 - qJD(6) * t202 + t28;
t474 = t584 + (-Ifges(7,3) * t469 + t671) * t740 + (Ifges(6,2) * t469 + t675) * t741 + t107 * t557 + t28 * (-mrSges(6,1) * t469 + mrSges(6,2) * t464) + t469 * t22 / 0.2e1 + t9 * (-mrSges(7,1) * t469 - mrSges(7,3) * t464) + t19 * t698 + (t464 * t801 - t670 + t674) * t742 + (t464 * t799 + t469 * t827) * t735 + t464 * t838 + t822 * t556 + t534 * mrSges(6,3) + t834 * mrSges(7,2) + (t563 - t491 / 0.2e1 + t819) * qJD(5) + (t490 + t489 + t488 + t487 + t486) * qJD(5) / 0.2e1 + t824;
t455 = -pkin(4) - t693;
t440 = Ifges(3,4) * t573;
t408 = t456 * t629;
t403 = t417 - t693;
t399 = (-mrSges(3,1) * t472 + mrSges(3,2) * t467) * t462;
t383 = -mrSges(3,2) * t445 + t545;
t315 = Ifges(3,1) * t574 + t440 + t659;
t314 = t524 * t613 + t658;
t312 = t398 * t471 + t587;
t296 = t365 * t464 + t585;
t282 = mrSges(4,1) * t424 - mrSges(4,3) * t350;
t281 = -mrSges(4,2) * t424 + t684;
t238 = t265 * t469 - t434;
t235 = t303 * t469 + t397 * t464;
t234 = t303 * t464 - t397 * t469;
t219 = t660 + t661 + t663;
t216 = t406 * t516 - t770;
t211 = mrSges(5,2) * t498 + mrSges(5,3) * t551;
t192 = -mrSges(4,2) * t371 + mrSges(4,3) * t229;
t191 = mrSges(4,1) * t371 - mrSges(4,3) * t228;
t188 = -pkin(5) * t405 - t206;
t187 = qJ(6) * t405 + t772;
t171 = -mrSges(5,1) * t551 + mrSges(5,2) * t506;
t166 = -mrSges(4,1) * t229 + mrSges(4,2) * t228;
t163 = -t498 * Ifges(5,3) + t669 + t673;
t140 = mrSges(7,1) * t201 - mrSges(7,3) * t202;
t139 = pkin(5) * t202 + qJ(6) * t201;
t121 = -qJD(5) * t434 + t168 * t464 + t265 * t603 - t469 * t572;
t120 = -qJD(5) * t237 + t168 * t469 + t464 * t572;
t103 = -mrSges(5,2) * t363 + mrSges(5,3) * t127;
t102 = mrSges(5,1) * t363 - mrSges(5,3) * t126;
t82 = pkin(5) * t237 - qJ(6) * t238 + t135;
t60 = pkin(5) * t505 - t71;
t59 = -qJ(6) * t505 + t782;
t57 = -t63 - t807;
t56 = t64 + t785;
t54 = -mrSges(5,1) * t127 + mrSges(5,2) * t126;
t51 = -t61 - t807;
t50 = t62 + t785;
t33 = mrSges(6,1) * t79 + mrSges(6,2) * t78;
t32 = mrSges(7,1) * t79 - mrSges(7,3) * t78;
t14 = pkin(5) * t121 - qJ(6) * t120 - qJD(6) * t238 + t43;
t11 = -pkin(5) * t169 - t13;
t10 = qJ(6) * t169 - qJD(6) * t505 + t12;
t1 = [-(t583 + t584) * t629 / 0.2e1 + (t110 * t120 - t169 * t53) * mrSges(6,2) + (-t120 * t58 + t169 * t45) * mrSges(7,3) + (-t116 * t168 - t117 * t169) * mrSges(5,3) + (t28 * mrSges(6,2) + t4 * mrSges(7,2) - t7 * mrSges(6,3) - t9 * mrSges(7,3) + Ifges(6,4) * t741 + Ifges(7,5) * t740 + t838 + t844) * t238 + (Ifges(4,5) * t394 + Ifges(4,6) * t393) * t702 + t840 * t265 - t841 * t505 + t424 * (Ifges(4,5) * t308 + Ifges(4,6) * t307) / 0.2e1 + t308 * t720 + t307 * t721 + (Ifges(5,1) * t168 - Ifges(5,4) * t169) * t710 + (Ifges(5,4) * t168 - Ifges(5,2) * t169) * t712 + (Ifges(3,3) * t700 + t401 * mrSges(3,1) - t402 * mrSges(3,2) + (Ifges(3,5) * t467 + Ifges(3,6) * t472) * t462) * t444 + t268 * (mrSges(3,1) * t463 - mrSges(3,3) * t632) + (t130 * t393 - t131 * t394 - t224 * t308 + t225 * t307) * mrSges(4,3) + t777 * t389 + (Ifges(6,4) * t120 + Ifges(6,6) * t169) * t725 + (Ifges(4,4) * t394 + Ifges(4,2) * t393) * t718 + (Ifges(7,5) * t120 + Ifges(7,6) * t169) * t724 + t168 * t727 + t169 * t730 + t394 * t731 + t393 * t732 + (-pkin(1) * t399 * t462 + Ifges(2,3)) * qJDD(1) + m(4) * (t130 * t254 + t131 * t253 + t176 * t225 + t177 * t224 + t251 * t366 + t317 * t389) + m(5) * (t116 * t47 + t117 * t46 + t137 * t31 + t138 * t30 + t170 * t270 + t246 * t249) + m(7) * (t10 * t45 + t11 * t44 + t14 * t58 + t2 * t59 + t4 * t60 + t82 * t9) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t747 + t267 * t402 + t268 * t401 - t384 * t389 + t387 * t388) - t251 * t532 + t349 * (Ifges(4,4) * t308 + Ifges(4,2) * t307) / 0.2e1 + t82 * t32 + t71 * t37 + t59 * t36 + t60 * t38 + (Ifges(4,1) * t394 + Ifges(4,4) * t393) * t719 + t388 * t383 + (-m(4) * (-pkin(2) * t396 + t555) - t550 * mrSges(4,1) - t482 * mrSges(4,2) - m(3) * t555 + t396 * mrSges(3,1) - mrSges(3,3) * t577 + t468 * mrSges(2,1) + t697 * mrSges(2,2) - m(5) * t493 + t299 * mrSges(5,1) + t549 * t813 - t541 * t230 - t760 * t395 + t500 * t298 + t810 * (pkin(4) * t299 - t493)) * g(1) + (-t697 * mrSges(2,1) - m(3) * t614 - t398 * mrSges(3,1) - m(5) * t511 - t303 * mrSges(5,1) - m(4) * (pkin(2) * t398 + t614) - t312 * mrSges(4,1) - t311 * mrSges(4,2) + (-mrSges(3,3) * t462 + mrSges(2,2)) * t468 - t549 * t235 + t541 * t234 + t760 * t397 + t500 * t302 - t810 * (t303 * pkin(4) + t511)) * g(2) + t366 * t166 + t317 * (-mrSges(4,1) * t307 + mrSges(4,2) * t308) + t176 * t281 + t177 * t282 + t270 * t54 + t246 * (mrSges(5,1) * t169 + mrSges(5,2) * t168) + t249 * t171 + t253 * t191 + t254 * t192 + t822 * t120 / 0.2e1 + t823 * t169 / 0.2e1 + m(6) * (t110 * t43 + t12 * t53 + t13 * t52 + t135 * t28 + t6 * t782 + t7 * t71) + t782 * t39 + t135 * t33 + t137 * t102 + t138 * t103 + t14 * t140 + t43 * t141 + (mrSges(3,3) * t267 - Ifges(4,5) * t719 - Ifges(5,5) * t734 - Ifges(4,6) * t718 - Ifges(5,6) * t733 - Ifges(4,3) * t702 - Ifges(5,3) * t703 - t818 - t824) * t629 + t10 * t153 + t12 * t154 + t13 * t155 + t11 * t156 - t463 * t825 + t52 * (mrSges(6,1) * t169 - mrSges(6,3) * t120) + t44 * (-mrSges(7,1) * t169 + mrSges(7,2) * t120) - t498 * (Ifges(5,5) * t168 - Ifges(5,6) * t169) / 0.2e1 + (Ifges(3,5) * t700 - t401 * mrSges(3,3) + (t467 * Ifges(3,1) + t680 - t738) * t462) * t391 + (Ifges(3,6) * t700 + t402 * mrSges(3,3) + (t524 + t739) * t462) * t390 + (t120 * t799 + t169 * t797) * t715 + (t110 * mrSges(6,1) - t45 * mrSges(7,2) - t53 * mrSges(6,3) - Ifges(6,2) * t725 + Ifges(7,3) * t724 - t715 * t827 + t722 * t800 + t737 - t812) * t121 + (mrSges(6,1) * t28 + mrSges(7,1) * t9 - Ifges(6,2) * t741 + Ifges(7,3) * t740 - t735 * t827 + t742 * t800 + t811) * t237 + (t120 * t801 + t169 * t799) * t722 + t46 * t211 + t47 * t212 + ((t659 / 0.2e1 + t315 / 0.2e1 - t384 * mrSges(3,3) + (-t738 + t680 / 0.2e1) * t613) * t472 + (-t658 / 0.2e1 + t761 * Ifges(5,3) / 0.2e1 + t660 / 0.2e1 + t663 / 0.2e1 + t661 / 0.2e1 - t225 * mrSges(4,2) + t224 * mrSges(4,1) + t673 / 0.2e1 + t669 / 0.2e1 + t116 * mrSges(5,1) - t117 * mrSges(5,2) + t163 / 0.2e1 + t219 / 0.2e1 - t314 / 0.2e1 - t387 * mrSges(3,3) + (-t739 - t681 / 0.2e1 + (-Ifges(5,3) / 0.2e1 + Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t472) * t613) * t467) * t610 + t579 * t700 + (Ifges(4,1) * t308 + Ifges(4,4) * t307) * t704; (t545 - t383) * t384 + t638 * t838 + (-t467 * (Ifges(3,1) * t472 - t681) / 0.2e1 + pkin(1) * (mrSges(3,1) * t467 + mrSges(3,2) * t472)) * qJD(1) ^ 2 * t747 + (-pkin(2) * t251 - t224 * t260 - t225 * t261) * m(4) + (t349 * t523 + t350 * t527 + t424 * t520) * qJD(3) / 0.2e1 + (t104 * t556 + t28 * t529 + t464 * t811 + t518 * t740 + t522 * t741 + t528 * t9 + t735 * t769 + t742 * t768 + t840) * t406 + t841 * t405 - mrSges(4,3) * t655 + (Ifges(7,3) * t725 - Ifges(6,2) * t724 + t107 / 0.2e1 + t800 * t723 - t827 * t716 + t812) * t262 + (-t494 * t799 - t495 * t827) * t715 + (Ifges(4,2) * t471 + t679) * t718 + (Ifges(4,1) * t466 + t678) * t719 + t607 * t720 + t540 * t721 + (-(mrSges(4,1) * t467 - mrSges(4,3) * t622) * t613 - mrSges(4,3) * t607) * t224 + t579 + (t52 * t775 + t53 * t776 - t638 * t7) * mrSges(6,3) + (t4 * t638 - t44 * t775 + t45 * t776) * mrSges(7,2) + (-mrSges(6,1) * t776 - mrSges(6,2) * t775) * t110 + (-m(4) * t317 + t546 - t777) * t387 + t778 * t211 + t779 * t212 + t780 * t141 - t825 + t314 * t574 / 0.2e1 - t225 * (-mrSges(4,3) * t466 * t472 - mrSges(4,2) * t467) * t613 + (-t466 * t191 + m(4) * (t656 - t655 + (-t224 * t471 - t225 * t466) * qJD(3)) + t471 * t192 - t282 * t607 - t281 * t608) * pkin(9) + (-Ifges(7,5) * t494 + Ifges(7,3) * t495) * t724 - t320 * t728 + t466 * t731 + t471 * t732 + t570 * t737 + t817 * t319 + t251 * t531 + (-Ifges(6,4) * t494 - Ifges(6,2) * t495) * t725 + t424 * t317 * (mrSges(4,1) * t466 + mrSges(4,2) * t471) + t58 * (mrSges(7,1) * t495 + mrSges(7,3) * t494) + t763 * t171 + (-m(5) * t408 + t399 - t802 * t589 - t810 * (pkin(11) * t589 + t459 * t448 - t473 * t632 + t408) + t541 * (t434 * t459 - t469 * t632) + (t767 * t472 + (m(5) * t473 - t789) * t467 - t549 * (t459 * t623 + t464 * t467)) * t462) * g(3) - t456 * t54 - (t33 - t102) * t770 + (t116 * t779 + t117 * t778 - t170 * t456 + t246 * t763 + t30 * t322 + t31 * t770) * m(5) + (t110 * t780 + t206 * t7 - t28 * t770 + t52 * t792 + t53 * t793 + t6 * t772) * m(6) - ((t219 + t163) * t467 + t424 * (Ifges(4,3) * t467 + t472 * t520) + t350 * (Ifges(4,5) * t467 + t472 * t527) + t349 * (Ifges(4,6) * t467 + t472 * t523) + t445 * (Ifges(3,5) * t472 - Ifges(3,6) * t467) + (-Ifges(3,2) * t574 + t471 * t221 + t315 + t440) * t472) * t613 / 0.2e1 + (-t664 - t220 / 0.2e1) * t608 + t322 * t103 - t261 * t281 - t260 * t282 + t268 * mrSges(3,1) - (Ifges(5,1) * t710 + Ifges(5,4) * t712 + t563 + t727 - t832) * t305 + (-Ifges(5,4) * t710 - Ifges(5,2) * t712 + Ifges(6,6) * t725 + Ifges(7,6) * t724 + t715 * t797 + t722 * t799 + t730 - t833) * t306 + t822 * (t406 * t557 - t648 / 0.2e1 - t263 / 0.2e1) + t823 * (-t319 / 0.2e1 + t306 / 0.2e1) + t772 * t39 + t786 * t140 - pkin(2) * t166 + t792 * t155 + t793 * t154 + t794 * t156 + t187 * t36 + t188 * t38 + t795 * t153 + (t187 * t2 + t188 * t4 + t216 * t9 + t44 * t794 + t45 * t795 + t58 * t786) * m(7) + (t58 * mrSges(7,3) + Ifges(6,4) * t724 + Ifges(7,5) * t725 + t716 * t799 + t723 * t801) * t263 + (-t494 * t801 + t495 * t800) * t722 + t206 * t37 + t216 * t32 + (-m(5) * t616 + t802 * t641 - t810 * (-pkin(11) * t641 - t397 * t692 + t616) - t549 * (-t397 * t633 + t398 * t464) + t541 * (-t397 * t634 - t398 * t469) + t760 * t398 + t759 * t397) * g(1) + (-m(5) * t617 + t802 * t644 - t810 * (-pkin(11) * t644 - t395 * t692 + t617) - t549 * (-t395 * t633 + t396 * t464) + t541 * (-t395 * t634 - t396 * t469) + t760 * t396 + t759 * t395) * g(2) + (t117 * t574 + t246 * t320) * mrSges(5,2) + t498 * (-Ifges(5,5) * t320 + Ifges(5,3) * t574) / 0.2e1 - t116 * (mrSges(5,1) * t574 + mrSges(5,3) * t320) + (-Ifges(5,4) * t320 + Ifges(5,6) * t574) * t713 + (-Ifges(5,1) * t320 + Ifges(5,5) * t574) * t711 + mrSges(4,3) * t656 + (Ifges(4,5) * t466 + Ifges(4,6) * t471) * t702; (-t281 + t684) * t224 + (t596 - t129) * t211 + (-t44 * t653 + t45 * t654 - t592) * mrSges(7,2) - (-Ifges(4,2) * t350 + t221 + t333) * t349 / 0.2e1 + (-t110 * t128 - t52 * t61 - t53 * t62 + t28 * t455 + (t110 * t465 + (-t464 * t52 + t469 * t53) * t470) * qJD(4) * pkin(3)) * m(6) + (-t44 * t51 - t45 * t50 + t403 * t9 + (t44 * t464 + t45 * t469) * t596 + t788 * t58) * m(7) + (t116 * t128 - t117 * t129 - t246 * t695 + (t30 * t465 + t31 * t470 + (-t116 * t465 + t117 * t470) * qJD(4)) * pkin(3)) * m(5) - t171 * t695 + t818 + t583 + (-m(5) * t504 - m(6) * t476 - t532 - m(7) * (t476 - t771) + t758) * g(3) + (-m(5) * t536 - m(6) * t478 - m(7) * (t478 + t773) - mrSges(4,1) * t311 + mrSges(4,2) * t312 + t757) * g(1) + (m(5) * t477 - m(6) * t475 - m(7) * (t475 - t774) + mrSges(4,1) * t482 - mrSges(4,2) * t550 + t756) * g(2) + ((-qJD(5) * t512 + t534) * m(6) + (qJD(5) * t513 + t535) * m(7) + t815) * t454 + t474 + (Ifges(5,1) * t711 + Ifges(5,4) * t713 + t518 * t725 + t522 * t724 + t716 * t769 + t723 * t768 + t816) * t551 + t817 * t506 - t350 * (Ifges(4,1) * t349 - t662) / 0.2e1 + (t52 * t653 + t53 * t654 + t764) * mrSges(6,3) + (-t51 - t765) * t156 + (-t61 + t765) * t155 + (-t62 + t766) * t154 + (-t50 + t766) * t153 + t455 * t33 - t424 * (Ifges(4,5) * t349 - Ifges(4,6) * t350) / 0.2e1 + t403 * t32 - t317 * (mrSges(4,1) * t350 + mrSges(4,2) * t349) + t225 * t282 + t821 * t781 + t823 * t711 + t788 * t140 + t350 * t664 + t102 * t693 + t103 * t694 + t220 * t704; t756 * g(2) + t757 * g(1) + t758 * g(3) + t474 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t506 + t790 - t513 * mrSges(7,2) - t487 / 0.2e1 - t488 / 0.2e1 - t489 / 0.2e1 - t490 / 0.2e1 + t491 / 0.2e1 - t486 / 0.2e1 - t239 / 0.2e1 + t816) * t551 + (-mrSges(7,2) * t657 - t790) * qJD(5) - t781 * t117 - pkin(4) * t33 + (-t105 / 0.2e1 - t106 / 0.2e1 + t677 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1) * t240 + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t202 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t201 + t830) * t506 + t417 * t32 - t56 * t153 - t64 * t154 - t63 * t155 - t57 * t156 + t787 * t140 - t116 * t211 + (t417 * t9 + (-t554 + t774) * g(2) + (-t553 - t773) * g(1) + (-t552 + t771) * g(3) - t44 * t57 - t45 * t56 + t787 * t58) * m(7) + (-pkin(4) * t28 - g(1) * t553 - g(2) * t554 - g(3) * t552 - t110 * t117 - t52 * t63 - t53 * t64) * m(6) + (m(7) * (-t592 + t834) + m(6) * (t534 + t764) + (t464 * t621 + t469 * t620) * qJD(5) + t815) * pkin(11); (-t201 * t801 + t104 + t198 - t676) * t723 + (-pkin(5) * t4 + qJ(6) * t2 - t139 * t58 - t44 * t53 + t45 * t783) * m(7) + (-t201 * t799 - t202 * t827) * t716 + (-Ifges(6,2) * t202 - t199 + t822) * t724 + t751 + (t549 * t230 + t541 * t813) * g(2) + t107 * t722 + (Ifges(7,3) * t202 - t672) * t725 + (-t620 + t682) * t53 + (t621 - t683) * t52 + (t234 * t549 + t235 * t541) * g(1) + (t541 * (t365 * t469 - t434) + t549 * t296) * g(3) + (t201 * t44 + t202 * t45) * mrSges(7,2) + qJ(6) * t36 - pkin(5) * t38 - t139 * t140 + qJD(6) * t153 - t58 * (mrSges(7,1) * t202 + mrSges(7,3) * t201) - t110 * (mrSges(6,1) * t202 - mrSges(6,2) * t201) + t796; t202 * t140 - t240 * t153 + (-g(1) * t234 - g(2) * t230 - g(3) * t296 + t202 * t58 - t240 * t45 + t4) * m(7) + t38;];
tau  = t1;
