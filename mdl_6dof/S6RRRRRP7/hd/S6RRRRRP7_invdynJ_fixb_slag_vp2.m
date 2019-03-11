% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP7
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:32
% EndTime: 2019-03-10 01:36:18
% DurationCPUTime: 61.58s
% Computational Cost: add. (30454->1057), mult. (72252->1379), div. (0->0), fcn. (58101->14), ass. (0->488)
t456 = cos(pkin(6));
t585 = qJD(1) * t456;
t436 = qJD(2) + t585;
t460 = sin(qJ(3));
t465 = cos(qJ(3));
t461 = sin(qJ(2));
t455 = sin(pkin(6));
t586 = qJD(1) * t455;
t552 = t461 * t586;
t347 = t436 * t465 - t460 * t552;
t348 = t436 * t460 + t465 * t552;
t459 = sin(qJ(4));
t464 = cos(qJ(4));
t527 = t464 * t347 - t348 * t459;
t709 = t527 / 0.2e1;
t466 = cos(qJ(2));
t551 = t466 * t586;
t414 = qJD(3) - t551;
t406 = qJD(4) + t414;
t695 = t406 / 0.2e1;
t495 = t347 * t459 + t464 * t348;
t707 = t495 / 0.2e1;
t855 = Ifges(5,4) * t707 + Ifges(5,6) * t695;
t861 = -Ifges(5,2) * t709 - t855;
t572 = pkin(1) * t585;
t377 = -pkin(8) * t552 + t466 * t572;
t490 = (pkin(2) * t461 - pkin(9) * t466) * t455;
t378 = qJD(1) * t490;
t265 = -t377 * t460 + t465 * t378;
t467 = -pkin(10) - pkin(9);
t556 = qJD(3) * t467;
t603 = t465 * t466;
t860 = -(pkin(3) * t461 - pkin(10) * t603) * t586 - t265 + t465 * t556;
t266 = t465 * t377 + t460 * t378;
t520 = t460 * t551;
t859 = -pkin(10) * t520 - t460 * t556 + t266;
t808 = Ifges(6,4) + Ifges(7,4);
t400 = t459 * t465 + t460 * t464;
t763 = qJD(3) + qJD(4);
t311 = t763 * t400;
t324 = t400 * t551;
t857 = t311 - t324;
t458 = sin(qJ(5));
t463 = cos(qJ(5));
t673 = -mrSges(7,1) - mrSges(6,1);
t810 = mrSges(6,2) + mrSges(7,2);
t856 = -t458 * t810 - t463 * t673 + mrSges(5,1);
t809 = Ifges(6,1) + Ifges(7,1);
t807 = Ifges(6,5) + Ifges(7,5);
t806 = Ifges(6,2) + Ifges(7,2);
t805 = Ifges(6,6) + Ifges(7,6);
t421 = t467 * t460;
t422 = t467 * t465;
t774 = t464 * t421 + t422 * t459;
t781 = qJD(4) * t774 + t459 * t860 - t859 * t464;
t432 = pkin(8) * t551;
t380 = t461 * t572 + t432;
t581 = qJD(3) * t460;
t775 = -t380 + (-t520 + t581) * pkin(3);
t323 = pkin(9) * t436 + t380;
t366 = (-pkin(2) * t466 - pkin(9) * t461 - pkin(1)) * t455;
t334 = qJD(1) * t366;
t229 = -t323 * t460 + t465 * t334;
t194 = -pkin(10) * t348 + t229;
t180 = pkin(3) * t414 + t194;
t230 = t323 * t465 + t334 * t460;
t195 = pkin(10) * t347 + t230;
t189 = t464 * t195;
t116 = t180 * t459 + t189;
t322 = -pkin(2) * t436 - t377;
t251 = -pkin(3) * t347 + t322;
t245 = qJD(5) - t527;
t208 = t406 * t458 + t463 * t495;
t111 = pkin(11) * t406 + t116;
t134 = -pkin(4) * t527 - pkin(11) * t495 + t251;
t54 = -t111 * t458 + t463 * t134;
t41 = -qJ(6) * t208 + t54;
t33 = pkin(5) * t245 + t41;
t207 = t406 * t463 - t458 * t495;
t55 = t111 * t463 + t134 * t458;
t42 = qJ(6) * t207 + t55;
t839 = -t251 * mrSges(5,1) - t54 * mrSges(6,1) - t33 * mrSges(7,1) + t55 * mrSges(6,2) + t42 * mrSges(7,2) + t116 * mrSges(5,3) - t861;
t457 = -qJ(6) - pkin(11);
t744 = m(6) * pkin(11);
t835 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t475 = m(7) * t457 - t744 + t835;
t583 = qJD(2) * t455;
t550 = t461 * t583;
t575 = qJDD(1) * t455;
t383 = -qJD(1) * t550 + t466 * t575;
t369 = qJDD(3) - t383;
t361 = qJDD(4) + t369;
t699 = t361 / 0.2e1;
t854 = 0.2e1 * t699;
t582 = qJD(2) * t466;
t384 = (qJD(1) * t582 + qJDD(1) * t461) * t455;
t574 = qJDD(1) * t456;
t435 = qJDD(2) + t574;
t233 = qJD(3) * t347 + t384 * t465 + t435 * t460;
t234 = -qJD(3) * t348 - t384 * t460 + t435 * t465;
t127 = -qJD(4) * t495 - t233 * t459 + t234 * t464;
t731 = t127 / 0.2e1;
t853 = 0.2e1 * t731;
t126 = qJD(4) * t527 + t233 * t464 + t234 * t459;
t732 = t126 / 0.2e1;
t852 = 0.2e1 * t732;
t849 = pkin(11) * t552 - t781;
t399 = t459 * t460 - t464 * t465;
t310 = t763 * t399;
t325 = t399 * t551;
t848 = t775 + (t310 - t325) * pkin(11) + t857 * pkin(4);
t847 = t808 * t207;
t846 = t808 * t208;
t712 = t245 / 0.2e1;
t720 = t208 / 0.2e1;
t722 = t207 / 0.2e1;
t804 = Ifges(6,3) + Ifges(7,3);
t845 = t712 * t804 + t720 * t807 + t722 * t805 - t839 + t861;
t690 = cos(qJ(1));
t554 = t690 * t461;
t462 = sin(qJ(1));
t605 = t462 * t466;
t390 = t456 * t554 + t605;
t454 = qJ(3) + qJ(4);
t450 = sin(t454);
t451 = cos(t454);
t555 = t455 * t690;
t301 = t390 * t450 + t451 * t555;
t302 = t390 * t451 - t450 * t555;
t843 = t301 * t856 + t302 * t835;
t553 = t690 * t466;
t606 = t461 * t462;
t392 = -t456 * t606 + t553;
t611 = t455 * t462;
t305 = t392 * t450 - t451 * t611;
t306 = t392 * t451 + t450 * t611;
t842 = t305 * t856 + t306 * t835;
t612 = t455 * t461;
t362 = t450 * t612 - t456 * t451;
t363 = t450 * t456 + t451 * t612;
t841 = t362 * t856 + t363 * t835;
t104 = t208 * Ifges(7,5) + t207 * Ifges(7,6) + t245 * Ifges(7,3);
t105 = t208 * Ifges(6,5) + t207 * Ifges(6,6) + t245 * Ifges(6,3);
t786 = t105 + t104;
t785 = t207 * t806 + t245 * t805 + t846;
t784 = t208 * t809 + t245 * t807 + t847;
t328 = t421 * t459 - t422 * t464;
t780 = -qJD(4) * t328 + t859 * t459 + t464 * t860;
t268 = t325 * t458 + t463 * t552;
t576 = qJD(5) * t463;
t548 = t400 * t576;
t628 = t310 * t458;
t487 = t548 - t628;
t777 = t268 + t487;
t578 = qJD(4) * t464;
t579 = qJD(4) * t459;
t521 = qJD(2) * t572;
t566 = pkin(1) * t574;
t275 = pkin(8) * t383 + t461 * t566 + t466 * t521;
t255 = pkin(9) * t435 + t275;
t260 = -pkin(1) * t575 - pkin(2) * t383 - pkin(9) * t384;
t131 = -qJD(3) * t230 - t255 * t460 + t465 * t260;
t86 = pkin(3) * t369 - pkin(10) * t233 + t131;
t580 = qJD(3) * t465;
t130 = t465 * t255 + t460 * t260 - t323 * t581 + t334 * t580;
t93 = pkin(10) * t234 + t130;
t28 = t180 * t578 - t195 * t579 + t459 * t86 + t464 * t93;
t25 = pkin(11) * t361 + t28;
t437 = pkin(8) * t612;
t276 = -qJD(2) * t432 - qJDD(1) * t437 - t461 * t521 + t466 * t566;
t256 = -pkin(2) * t435 - t276;
t174 = -pkin(3) * t234 + t256;
t36 = -pkin(4) * t127 - pkin(11) * t126 + t174;
t6 = -qJD(5) * t55 - t25 * t458 + t463 * t36;
t124 = qJDD(5) - t127;
t78 = qJD(5) * t207 + t126 * t463 + t361 * t458;
t1 = pkin(5) * t124 - qJ(6) * t78 - qJD(6) * t208 + t6;
t736 = t1 * mrSges(7,3);
t838 = t6 * mrSges(6,3) + t736;
t244 = Ifges(5,4) * t527;
t645 = t406 * Ifges(5,5);
t167 = Ifges(5,1) * t495 + t244 + t645;
t188 = t459 * t195;
t115 = t180 * t464 - t188;
t754 = -t251 * mrSges(5,2) + t115 * mrSges(5,3);
t837 = -t754 + Ifges(5,1) * t707 + Ifges(5,4) * t709 + Ifges(5,5) * t695 + t167 / 0.2e1;
t733 = t124 / 0.2e1;
t79 = -qJD(5) * t208 - t126 * t458 + t361 * t463;
t739 = t79 / 0.2e1;
t740 = t78 / 0.2e1;
t577 = qJD(5) * t458;
t5 = -t111 * t577 + t134 * t576 + t463 * t25 + t458 * t36;
t3 = qJ(6) * t79 + qJD(6) * t207 + t5;
t750 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t803 = t124 * t804 + t78 * t807 + t79 * t805;
t833 = -t803 / 0.2e1 - t174 * mrSges(5,1) - t1 * mrSges(7,1) + t28 * mrSges(5,3) + Ifges(5,4) * t852 + Ifges(5,2) * t853 + Ifges(5,6) * t854 - t804 * t733 - t805 * t739 - t807 * t740 - t750;
t29 = -t180 * t579 - t195 * t578 - t459 * t93 + t464 * t86;
t26 = -pkin(4) * t361 - t29;
t12 = -pkin(5) * t79 + qJDD(6) + t26;
t832 = t26 * mrSges(6,2) + t12 * mrSges(7,2) + t807 * t733 + t809 * t740;
t802 = t124 * t805 + t78 * t808 + t79 * t806;
t831 = -t26 * mrSges(6,1) - t12 * mrSges(7,1) + t805 * t733 + t806 * t739 + t802 / 0.2e1 + t3 * mrSges(7,3);
t743 = m(7) * pkin(5);
t828 = t784 / 0.2e1;
t387 = t456 * t465 - t460 * t612;
t491 = pkin(3) * t387;
t787 = t527 * t458;
t826 = pkin(5) * t787;
t801 = t124 * t807 + t78 * t809 + t79 * t808;
t569 = pkin(3) * t578;
t129 = t194 * t464 - t188;
t176 = pkin(4) * t495 - pkin(11) * t527;
t688 = pkin(3) * t348;
t151 = t176 + t688;
t61 = t463 * t129 + t458 * t151;
t825 = t463 * t569 - t61;
t824 = t458 * t849 + t463 * t848;
t448 = pkin(3) * t465 + pkin(2);
t288 = pkin(4) * t399 - pkin(11) * t400 - t448;
t823 = t288 * t576 + t458 * t848 - t463 * t849;
t779 = pkin(4) * t552 - t780;
t110 = -pkin(4) * t406 - t115;
t509 = mrSges(7,1) * t458 + mrSges(7,2) * t463;
t510 = mrSges(6,1) * t458 + mrSges(6,2) * t463;
t84 = -pkin(5) * t207 + qJD(6) + t110;
t822 = t110 * t510 + t84 * t509;
t610 = t455 * t465;
t316 = -t392 * t460 + t462 * t610;
t128 = t194 * t459 + t189;
t821 = pkin(3) * t579 - t128;
t696 = -t406 / 0.2e1;
t708 = -t495 / 0.2e1;
t710 = -t527 / 0.2e1;
t713 = -t245 / 0.2e1;
t721 = -t208 / 0.2e1;
t723 = -t207 / 0.2e1;
t820 = -Ifges(5,4) * t708 - Ifges(5,2) * t710 - Ifges(5,6) * t696 + t713 * t804 + t721 * t807 + t723 * t805 + t839;
t446 = pkin(5) * t463 + pkin(4);
t512 = -mrSges(4,1) * t465 + mrSges(4,2) * t460;
t819 = -m(4) * pkin(2) + t512 + (-m(6) * pkin(4) - m(7) * t446 - mrSges(5,1)) * t451 + t475 * t450;
t818 = qJ(6) * t787 + t463 * qJD(6);
t389 = -t456 * t553 + t606;
t817 = t302 * t458 - t389 * t463;
t816 = -t302 * t463 - t389 * t458;
t766 = -t54 * mrSges(6,3) - t33 * mrSges(7,3);
t767 = -t55 * mrSges(6,3) - t42 * mrSges(7,3);
t813 = t174 * mrSges(5,2) - t29 * mrSges(5,3) + Ifges(5,1) * t852 + Ifges(5,4) * t853 + Ifges(5,5) * t854;
t717 = t233 / 0.2e1;
t716 = t234 / 0.2e1;
t698 = t369 / 0.2e1;
t269 = -t325 * t463 + t458 * t552;
t315 = t463 * t328;
t492 = qJ(6) * t310 - qJD(6) * t400;
t800 = t492 * t463 + (-t315 + (qJ(6) * t400 - t288) * t458) * qJD(5) + qJ(6) * t269 + t824 + t857 * pkin(5);
t799 = (-qJD(5) * t328 + t492) * t458 + t823 + (-t548 - t268) * qJ(6);
t213 = t458 * t288 + t315;
t798 = -qJD(5) * t213 + t824;
t797 = -t328 * t577 + t823;
t796 = t275 * mrSges(3,2);
t795 = -m(4) * pkin(9) - mrSges(4,3) - mrSges(5,3);
t794 = t29 * mrSges(5,1) - t28 * mrSges(5,2);
t687 = pkin(3) * t459;
t445 = pkin(11) + t687;
t602 = -qJ(6) - t445;
t525 = qJD(5) * t602;
t793 = t458 * t525 + t818 + t825;
t60 = -t129 * t458 + t463 * t151;
t452 = t463 * qJ(6);
t751 = pkin(5) * t495 - t452 * t527;
t792 = (-qJD(6) - t569) * t458 + t463 * t525 - t60 - t751;
t529 = qJD(5) * t457;
t64 = t463 * t115 + t458 * t176;
t791 = t458 * t529 - t64 + t818;
t63 = -t115 * t458 + t463 * t176;
t790 = -qJD(6) * t458 + t463 * t529 - t63 - t751;
t568 = pkin(5) * t577;
t789 = t568 + t821 - t826;
t788 = t743 + mrSges(7,1);
t783 = pkin(5) * t777 + t779;
t609 = t455 * t466;
t396 = t456 * t461 * pkin(1) + pkin(8) * t609;
t365 = pkin(9) * t456 + t396;
t258 = -t365 * t460 + t465 * t366;
t388 = t456 * t460 + t461 * t610;
t203 = -pkin(3) * t609 - pkin(10) * t388 + t258;
t259 = t465 * t365 + t460 * t366;
t219 = pkin(10) * t387 + t259;
t139 = t459 * t203 + t464 * t219;
t137 = -pkin(11) * t609 + t139;
t271 = t387 * t459 + t388 * t464;
t689 = pkin(1) * t466;
t364 = t437 + (-pkin(2) - t689) * t456;
t279 = t364 - t491;
t494 = t464 * t387 - t388 * t459;
t164 = -pkin(4) * t494 - pkin(11) * t271 + t279;
t71 = t463 * t137 + t458 * t164;
t141 = -mrSges(6,1) * t207 + mrSges(6,2) * t208;
t218 = mrSges(5,1) * t406 - mrSges(5,3) * t495;
t782 = t141 - t218;
t523 = mrSges(3,3) * t552;
t778 = -mrSges(3,1) * t436 - mrSges(4,1) * t347 + mrSges(4,2) * t348 + t523;
t627 = t310 * t463;
t486 = t400 * t577 + t627;
t776 = t269 + t486;
t395 = t456 * t689 - t437;
t499 = Ifges(7,5) * t463 - Ifges(7,6) * t458;
t500 = Ifges(6,5) * t463 - Ifges(6,6) * t458;
t772 = t499 + t500;
t659 = Ifges(7,4) * t463;
t502 = -Ifges(7,2) * t458 + t659;
t661 = Ifges(6,4) * t463;
t503 = -Ifges(6,2) * t458 + t661;
t771 = t502 + t503;
t660 = Ifges(7,4) * t458;
t506 = Ifges(7,1) * t463 - t660;
t662 = Ifges(6,4) * t458;
t507 = Ifges(6,1) * t463 - t662;
t770 = t506 + t507;
t769 = t130 * t465 - t131 * t460;
t768 = t131 * mrSges(4,1) - t130 * mrSges(4,2);
t764 = m(6) + m(7) + m(5);
t762 = -mrSges(6,1) - t788;
t761 = mrSges(3,2) + t795;
t753 = -t458 * t743 + t795;
t752 = mrSges(3,1) - t819;
t693 = t458 / 0.2e1;
t726 = -t167 / 0.2e1;
t749 = t726 + t754 - t822 + t785 * t693 - t784 * t463 / 0.2e1;
t745 = t455 ^ 2;
t738 = pkin(1) * mrSges(3,1);
t737 = pkin(1) * mrSges(3,2);
t730 = Ifges(4,4) * t717 + Ifges(4,2) * t716 + Ifges(4,6) * t698;
t729 = Ifges(4,1) * t717 + Ifges(4,4) * t716 + Ifges(4,5) * t698;
t647 = t348 * Ifges(4,4);
t225 = t347 * Ifges(4,2) + t414 * Ifges(4,6) + t647;
t719 = t225 / 0.2e1;
t336 = Ifges(4,4) * t347;
t226 = t348 * Ifges(4,1) + t414 * Ifges(4,5) + t336;
t718 = t226 / 0.2e1;
t700 = t348 / 0.2e1;
t694 = t456 / 0.2e1;
t686 = pkin(3) * t464;
t683 = pkin(11) * t463;
t679 = t463 * t5;
t676 = t6 * t458;
t671 = mrSges(6,3) * t207;
t670 = mrSges(6,3) * t208;
t669 = mrSges(7,3) * t207;
t668 = mrSges(7,3) * t208;
t667 = Ifges(3,4) * t461;
t666 = Ifges(3,4) * t466;
t665 = Ifges(4,4) * t460;
t664 = Ifges(4,4) * t465;
t652 = t229 * mrSges(4,3);
t651 = t230 * mrSges(4,3);
t650 = t527 * Ifges(5,6);
t649 = t495 * Ifges(5,5);
t648 = t347 * Ifges(4,6);
t646 = t348 * Ifges(4,5);
t643 = t406 * Ifges(5,3);
t642 = t414 * Ifges(4,3);
t641 = t436 * Ifges(3,5);
t640 = t436 * Ifges(3,6);
t633 = t527 * t463;
t621 = t390 * t458;
t620 = t392 * t458;
t618 = t400 * t458;
t617 = t400 * t463;
t615 = t445 * t463;
t614 = t451 * t458;
t613 = t451 * t463;
t607 = t458 * t466;
t604 = t463 * t466;
t601 = -t301 * t446 - t302 * t457;
t600 = -t305 * t446 - t306 * t457;
t593 = -t362 * t446 - t363 * t457;
t587 = t690 * pkin(1) + pkin(8) * t611;
t565 = t455 * t607;
t563 = t460 * t611;
t561 = t455 * t604;
t560 = Ifges(5,5) * t126 + Ifges(5,6) * t127 + Ifges(5,3) * t361;
t558 = Ifges(4,5) * t233 + Ifges(4,6) * t234 + Ifges(4,3) * t369;
t557 = Ifges(3,5) * t384 + Ifges(3,6) * t383 + Ifges(3,3) * t435;
t549 = t455 * t582;
t31 = -t79 * mrSges(7,1) + t78 * mrSges(7,2);
t535 = -t577 / 0.2e1;
t533 = -pkin(1) * t462 + pkin(8) * t555;
t532 = -t301 * pkin(4) + t302 * pkin(11);
t531 = -t305 * pkin(4) + pkin(11) * t306;
t530 = -t362 * pkin(4) + pkin(11) * t363;
t70 = -t137 * t458 + t463 * t164;
t138 = t203 * t464 - t459 * t219;
t212 = t463 * t288 - t328 * t458;
t425 = t460 * t555;
t526 = -t390 * t465 + t425;
t522 = mrSges(3,3) * t551;
t515 = t316 * pkin(3);
t136 = pkin(4) * t609 - t138;
t513 = mrSges(4,1) * t387 - mrSges(4,2) * t388;
t508 = Ifges(4,1) * t465 - t665;
t505 = Ifges(3,2) * t466 + t667;
t504 = -Ifges(4,2) * t460 + t664;
t501 = Ifges(4,5) * t465 - Ifges(4,6) * t460;
t498 = t458 * t55 + t463 * t54;
t391 = t456 * t605 + t554;
t497 = pkin(3) * t563 - t391 * t467 + t392 * t448 + t587;
t239 = -t306 * t458 + t391 * t463;
t379 = qJD(2) * t490;
t381 = t395 * qJD(2);
t182 = -qJD(3) * t259 + t465 * t379 - t381 * t460;
t313 = qJD(3) * t387 + t465 * t549;
t150 = pkin(3) * t550 - pkin(10) * t313 + t182;
t181 = -t365 * t581 + t366 * t580 + t460 * t379 + t465 * t381;
t312 = -qJD(3) * t388 - t460 * t549;
t161 = pkin(10) * t312 + t181;
t48 = t150 * t464 - t459 * t161 - t203 * t579 - t219 * t578;
t242 = -t271 * t458 - t561;
t489 = -t271 * t463 + t565;
t485 = pkin(3) * t425 + t389 * t467 - t390 * t448 + t533;
t483 = t207 * t503;
t482 = t207 * t502;
t481 = t208 * t507;
t480 = t208 * t506;
t479 = t245 * t500;
t478 = t245 * t499;
t47 = t459 * t150 + t464 * t161 + t203 * t578 - t219 * t579;
t44 = pkin(11) * t550 + t47;
t172 = qJD(4) * t494 + t312 * t459 + t313 * t464;
t173 = qJD(4) * t271 - t464 * t312 + t313 * t459;
t382 = t396 * qJD(2);
t254 = -pkin(3) * t312 + t382;
t83 = pkin(4) * t173 - pkin(11) * t172 + t254;
t9 = -t137 * t577 + t164 * t576 + t463 * t44 + t458 * t83;
t474 = t390 * t460 + t465 * t555;
t473 = -mrSges(3,2) - t753;
t472 = t474 * pkin(3);
t45 = -pkin(4) * t550 - t48;
t10 = -qJD(5) * t71 - t44 * t458 + t463 * t83;
t469 = -qJD(5) * t498 - t676 + t679;
t468 = t794 + t560 + mrSges(6,3) * t679 + (t659 + t661) * t740 + (t660 + t662) * t739 + t801 * t693 + t785 * t535 + t576 * t828 + t822 * qJD(5) + (t483 + t482 + t481 + t480 + t479 + t478) * qJD(5) / 0.2e1 + t832 * t458 + t831 * t463;
t447 = -pkin(4) - t686;
t430 = Ifges(3,4) * t551;
t419 = t452 + t683;
t418 = t457 * t458;
t410 = -t446 - t686;
t398 = t452 + t615;
t397 = t602 * t458;
t393 = (-mrSges(3,1) * t466 + mrSges(3,2) * t461) * t455;
t376 = -mrSges(3,2) * t436 + t522;
t320 = Ifges(3,1) * t552 + t430 + t641;
t319 = t505 * t586 + t640;
t317 = t392 * t465 + t563;
t287 = mrSges(4,1) * t414 - mrSges(4,3) * t348;
t286 = -mrSges(4,2) * t414 + mrSges(4,3) * t347;
t267 = pkin(5) * t618 - t774;
t240 = t306 * t463 + t391 * t458;
t224 = t642 + t646 + t648;
t217 = -mrSges(5,2) * t406 + mrSges(5,3) * t527;
t197 = -mrSges(4,2) * t369 + mrSges(4,3) * t234;
t196 = mrSges(4,1) * t369 - mrSges(4,3) * t233;
t187 = -qJ(6) * t618 + t213;
t177 = pkin(5) * t399 - t400 * t452 + t212;
t175 = -mrSges(5,1) * t527 + mrSges(5,2) * t495;
t168 = -mrSges(4,1) * t234 + mrSges(4,2) * t233;
t165 = t643 + t649 + t650;
t157 = mrSges(6,1) * t245 - t670;
t156 = mrSges(7,1) * t245 - t668;
t155 = -mrSges(6,2) * t245 + t671;
t154 = -mrSges(7,2) * t245 + t669;
t140 = -mrSges(7,1) * t207 + mrSges(7,2) * t208;
t121 = qJD(5) * t489 - t172 * t458 + t463 * t550;
t120 = qJD(5) * t242 + t172 * t463 + t458 * t550;
t103 = -mrSges(5,2) * t361 + mrSges(5,3) * t127;
t102 = mrSges(5,1) * t361 - mrSges(5,3) * t126;
t97 = -pkin(5) * t242 + t136;
t87 = t116 + t826;
t58 = qJ(6) * t242 + t71;
t56 = -mrSges(5,1) * t127 + mrSges(5,2) * t126;
t52 = -pkin(5) * t494 + qJ(6) * t489 + t70;
t40 = -mrSges(6,2) * t124 + mrSges(6,3) * t79;
t39 = -mrSges(7,2) * t124 + mrSges(7,3) * t79;
t38 = mrSges(6,1) * t124 - mrSges(6,3) * t78;
t37 = mrSges(7,1) * t124 - mrSges(7,3) * t78;
t32 = -mrSges(6,1) * t79 + mrSges(6,2) * t78;
t30 = -pkin(5) * t121 + t45;
t8 = qJ(6) * t121 + qJD(6) * t242 + t9;
t7 = pkin(5) * t173 - qJ(6) * t120 + qJD(6) * t489 + t10;
t2 = [-(t558 + t560) * t609 / 0.2e1 + (-m(7) * (-t302 * t446 + t485) - m(5) * t485 + t302 * mrSges(5,1) - m(6) * (-pkin(4) * t302 + t485) - m(4) * (-pkin(2) * t390 + t533) - t526 * mrSges(4,1) - t474 * mrSges(4,2) - m(3) * t533 + t390 * mrSges(3,1) - mrSges(3,3) * t555 + t462 * mrSges(2,1) + t690 * mrSges(2,2) + t673 * t816 - t810 * t817 + t473 * t389 - t475 * t301) * g(1) + (-m(3) * t587 - t392 * mrSges(3,1) - m(7) * (t306 * t446 + t497) - m(5) * t497 - t306 * mrSges(5,1) - m(6) * (pkin(4) * t306 + t497) - m(4) * (pkin(2) * t392 + t587) - t317 * mrSges(4,1) - t316 * mrSges(4,2) - t690 * mrSges(2,1) + (-mrSges(3,3) * t455 + mrSges(2,2)) * t462 + t673 * t240 - t810 * t239 - t473 * t391 + t475 * t305) * g(2) + (Ifges(4,5) * t388 + Ifges(4,6) * t387) * t698 + (t805 * t712 + t806 * t722 + t785 / 0.2e1 + t808 * t720 - mrSges(7,1) * t84 - mrSges(6,1) * t110 - t767) * t121 + (mrSges(6,2) * t110 + mrSges(7,2) * t84 + t807 * t712 + t809 * t720 + t808 * t722 + t766 + t828) * t120 + (t130 * t387 - t131 * t388 - t229 * t313 + t230 * t312) * mrSges(4,3) + m(7) * (t1 * t52 + t12 * t97 + t3 * t58 + t30 * t84 + t33 * t7 + t42 * t8) + m(6) * (t10 * t54 + t110 * t45 + t136 * t26 + t5 * t71 + t55 * t9 + t6 * t70) + m(5) * (t115 * t48 + t116 * t47 + t138 * t29 + t139 * t28 + t174 * t279 + t251 * t254) + m(4) * (t130 * t259 + t131 * t258 + t181 * t230 + t182 * t229 + t256 * t364 + t322 * t382) + (mrSges(3,3) * t275 - Ifges(4,5) * t717 - Ifges(5,5) * t732 - Ifges(4,6) * t716 - Ifges(5,6) * t731 - Ifges(4,3) * t698 - Ifges(5,3) * t699 - t768 - t794) * t609 + t778 * t382 + (Ifges(4,1) * t388 + Ifges(4,4) * t387) * t717 + t813 * t271 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t745 + t275 * t396 + t276 * t395 - t377 * t382 + t380 * t381) + (Ifges(4,4) * t388 + Ifges(4,2) * t387) * t716 + t70 * t38 + t71 * t40 + (Ifges(3,3) * t694 - t396 * mrSges(3,2) + t395 * mrSges(3,1) + (Ifges(3,5) * t461 + Ifges(3,6) * t466) * t455) * t435 + t276 * (mrSges(3,1) * t456 - mrSges(3,3) * t612) + t414 * (Ifges(4,5) * t313 + Ifges(4,6) * t312) / 0.2e1 - t456 * t796 + (-pkin(1) * t393 * t455 + Ifges(2,3)) * qJDD(1) + (t786 / 0.2e1 + t845) * t173 + (Ifges(3,5) * t694 - t395 * mrSges(3,3) + (t461 * Ifges(3,1) + t666 - t737) * t455) * t384 + ((-t377 * mrSges(3,3) + t641 / 0.2e1 + t320 / 0.2e1 + (-t737 + t666 / 0.2e1) * t586) * t466 + (-t380 * mrSges(3,3) - t640 / 0.2e1 + t643 / 0.2e1 + t642 / 0.2e1 + t646 / 0.2e1 + t648 / 0.2e1 - t230 * mrSges(4,2) + t229 * mrSges(4,1) - t116 * mrSges(5,2) + t649 / 0.2e1 + t650 / 0.2e1 + t115 * mrSges(5,1) + t224 / 0.2e1 - t319 / 0.2e1 + t165 / 0.2e1 + (-t738 - t667 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t466) * t586) * t461) * t583 + (Ifges(3,6) * t694 + t396 * mrSges(3,3) + (t505 + t738) * t455) * t383 + t58 * t39 + t52 * t37 + t381 * t376 + t364 * t168 + t322 * (-mrSges(4,1) * t312 + mrSges(4,2) * t313) + t181 * t286 + t182 * t287 + t279 * t56 + t254 * t175 + t258 * t196 + t259 * t197 + (t5 * mrSges(6,3) + t808 * t740 + t831) * t242 + t97 * t31 - t256 * t513 + t136 * t32 + t138 * t102 + t139 * t103 + t30 * t140 + t45 * t141 + t8 * t154 + t9 * t155 + t7 * t156 + t10 * t157 + t833 * t494 + t837 * t172 + (-t808 * t739 - t801 / 0.2e1 - t832 + t838) * t489 + t47 * t217 + t48 * t218 + t347 * (Ifges(4,4) * t313 + Ifges(4,2) * t312) / 0.2e1 + t557 * t694 + (Ifges(4,1) * t313 + Ifges(4,4) * t312) * t700 + t313 * t718 + t312 * t719 + t388 * t729 + t387 * t730; (t347 * t504 + t348 * t508 + t414 * t501) * qJD(3) / 0.2e1 - ((-Ifges(3,2) * t552 + t465 * t226 + t320 + t430) * t466 + t414 * (Ifges(4,3) * t461 + t466 * t501) + t348 * (Ifges(4,5) * t461 + t466 * t508) + t347 * (Ifges(4,6) * t461 + t466 * t504) + t436 * (Ifges(3,5) * t466 - Ifges(3,6) * t461) + (t224 + t165) * t461) * t586 / 0.2e1 - (t32 - t102) * t774 + (t110 * t779 + t212 * t6 + t213 * t5 - t26 * t774 + t54 * t798 + t55 * t797) * m(6) + (t115 * t780 + t116 * t781 - t174 * t448 + t251 * t775 + t28 * t328 + t29 * t774) * m(5) + t820 * t324 + (-t764 * t448 * t609 + t393 + (t673 * (t451 * t604 + t458 * t461) - t810 * (-t451 * t607 + t461 * t463) + t819 * t466 + (t467 * t764 + t753) * t461) * t455) * g(3) + (pkin(1) * (mrSges(3,1) * t461 + mrSges(3,2) * t466) - t461 * (Ifges(3,1) * t466 - t667) / 0.2e1) * qJD(1) ^ 2 * t745 + (-t652 + t718) * t580 + (-Ifges(5,1) * t325 + Ifges(5,5) * t552) * t708 + (-Ifges(5,4) * t325 + Ifges(5,6) * t552) * t710 + (t116 * t552 + t251 * t325) * mrSges(5,2) - t115 * (mrSges(5,1) * t552 + mrSges(5,3) * t325) + (-Ifges(5,5) * t325 + Ifges(5,3) * t552) * t696 + t797 * t155 + (t268 * t805 + t269 * t807) * t713 + (-t486 * t807 - t487 * t805) * t712 + (t268 * t806 + t269 * t808) * t723 + (-t486 * t808 - t487 * t806) * t722 + t798 * t157 + t799 * t154 + (t1 * t177 + t12 * t267 + t187 * t3 + t33 * t800 + t42 * t799 + t783 * t84) * m(7) + t800 * t156 + t801 * t617 / 0.2e1 - t802 * t618 / 0.2e1 + t775 * t175 + (mrSges(7,1) * t777 - mrSges(7,2) * t776) * t84 + (mrSges(6,1) * t777 - mrSges(6,2) * t776) * t110 + (-t5 * t618 + t54 * t776 - t55 * t777 - t6 * t617) * mrSges(6,3) + (-t1 * t617 - t3 * t618 + t33 * t776 - t42 * t777) * mrSges(7,3) + (-m(4) * t322 + t523 - t778) * t380 + t779 * t141 + (-t287 * t580 + m(4) * ((-t229 * t465 - t230 * t460) * qJD(3) + t769) - t286 * t581 + t465 * t197 - t460 * t196) * pkin(9) + t769 * mrSges(4,3) + (-t229 * (mrSges(4,1) * t461 - mrSges(4,3) * t603) - t230 * (-mrSges(4,3) * t460 * t466 - mrSges(4,2) * t461)) * t586 + (t12 * t509 + t26 * t510 + t733 * t772 + t739 * t771 + t740 * t770 + t813) * t400 - t796 + t780 * t218 + t781 * t217 + t783 * t140 + t414 * t322 * (mrSges(4,1) * t460 + mrSges(4,2) * t465) + (-t376 + t522) * t377 + t319 * t552 / 0.2e1 + (-t621 * t743 - t764 * (-t389 * t448 - t390 * t467) + t673 * (-t389 * t613 + t621) - t810 * (t389 * t614 + t390 * t463) + t761 * t390 + t752 * t389) * g(2) + (-t620 * t743 - t764 * (-t391 * t448 - t392 * t467) + t673 * (-t391 * t613 + t620) - t810 * (t391 * t614 + t392 * t463) + t761 * t392 + t752 * t391) * g(1) + (-pkin(2) * t256 - t229 * t265 - t230 * t266) * m(4) + (t268 * t808 + t269 * t809) * t721 + (-t486 * t809 - t487 * t808) * t720 + t845 * t311 + t557 - t448 * t56 + t328 * t103 - t266 * t286 - t265 * t287 + t276 * mrSges(3,1) + t267 * t31 + (-t651 - t225 / 0.2e1) * t581 + t256 * t512 - t833 * t399 - pkin(2) * t168 + t177 * t37 + t187 * t39 - t837 * t310 + t784 * (t400 * t535 - t627 / 0.2e1 - t269 / 0.2e1) + t785 * (t628 / 0.2e1 - t548 / 0.2e1 - t268 / 0.2e1) + t786 * (-t324 / 0.2e1 + t311 / 0.2e1) + t212 * t38 + t213 * t40 + (Ifges(4,5) * t460 + Ifges(4,6) * t465) * t698 + (Ifges(4,2) * t465 + t665) * t716 + (Ifges(4,1) * t460 + t664) * t717 + t520 * t719 - t325 * t726 + t460 * t729 + t465 * t730; -(-Ifges(4,2) * t348 + t226 + t336) * t347 / 0.2e1 + (t115 * t128 - t116 * t129 - t251 * t688 + (t28 * t459 + t29 * t464 + (-t115 * t459 + t116 * t464) * qJD(4)) * pkin(3)) * m(5) + (-t445 * t577 + t825) * t155 + t820 * t495 + t821 * t782 + (t569 - t129) * t217 + t40 * t615 + (-t445 * t576 - t458 * t569 - t60) * t157 + (-t110 * t128 - t54 * t60 - t55 * t61 + t26 * t447 + (t110 * t459 + (-t458 * t54 + t463 * t55) * t464) * qJD(4) * pkin(3) + t469 * t445) * m(6) + t789 * t140 + t792 * t156 + (t1 * t397 + t12 * t410 + t3 * t398 + t33 * t792 + t42 * t793 + t789 * t84) * m(7) + t793 * t154 + (Ifges(5,1) * t708 + Ifges(5,4) * t710 + Ifges(5,5) * t696 + t713 * t772 + t721 * t770 + t723 * t771 + t749) * t527 + (-m(5) * t491 - m(7) * (t491 + t593) - m(6) * (t491 + t530) - t513 + t841) * g(3) + (-mrSges(4,1) * t316 + mrSges(4,2) * t317 - m(5) * t515 - m(7) * (t515 + t600) - m(6) * (t515 + t531) + t842) * g(1) + (mrSges(4,1) * t474 - mrSges(4,2) * t526 + m(5) * t472 - m(7) * (-t472 + t601) - m(6) * (-t472 + t532) + t843) * g(2) + (t54 * t633 + t55 * t787 - t676) * mrSges(6,3) + (t33 * t633 + t42 * t787) * mrSges(7,3) + t768 - t175 * t688 + t766 * t576 + t767 * t577 + (-t38 * t445 - t736) * t458 + t468 + t558 + t447 * t32 + t410 * t31 - t414 * (Ifges(4,5) * t347 - Ifges(4,6) * t348) / 0.2e1 + t397 * t37 + t398 * t39 - t322 * (mrSges(4,1) * t348 + mrSges(4,2) * t347) - t229 * t286 + t230 * t287 + t348 * t651 + t347 * t652 - t348 * (Ifges(4,1) * t347 - t647) / 0.2e1 + t786 * t708 + t102 * t686 + t103 * t687 + t225 * t700; t469 * t744 + t790 * t156 + t842 * g(1) + t843 * g(2) + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t495 - t645 / 0.2e1 - t483 / 0.2e1 - t481 / 0.2e1 - t482 / 0.2e1 - t478 / 0.2e1 - t479 / 0.2e1 - t480 / 0.2e1 + (t33 * t463 + t42 * t458) * mrSges(7,3) + t749 + t498 * mrSges(6,3) - t244 / 0.2e1) * t527 + (-t104 / 0.2e1 - t105 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t245 + (-Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1) * t208 + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t207 + t839 + t855) * t495 + ((-pkin(11) * t157 + t766) * t463 + (pkin(5) * t140 - pkin(11) * t155 + t767) * t458) * qJD(5) + t468 - t782 * t116 + t791 * t154 + (-pkin(11) * t38 - t838) * t458 - pkin(4) * t32 - t446 * t31 + t418 * t37 + t419 * t39 + t841 * g(3) - t87 * t140 - t64 * t155 - t63 * t157 - t115 * t217 + t40 * t683 + (-t600 * g(1) - t601 * g(2) - t593 * g(3) + t1 * t418 - t12 * t446 + t3 * t419 + (t568 - t87) * t84 + t791 * t42 + t790 * t33) * m(7) + (-pkin(4) * t26 - g(1) * t531 - g(2) * t532 - g(3) * t530 - t110 * t116 - t54 * t63 - t55 * t64) * m(6); (-t762 * t817 - t810 * t816) * g(2) + (-t810 * (-t363 * t463 + t565) + t762 * (-t363 * t458 - t561)) * g(3) + (t239 * t762 + t240 * t810) * g(1) + t750 + (-t208 * t806 + t784 + t847) * t723 + t785 * t720 + (t207 * t809 - t846) * t721 + (t157 + t670) * t55 + (-t155 + t671) * t54 + (-m(7) * (-t33 + t41) + t156 + t668) * t42 + (t207 * t807 - t208 * t805) * t713 + t788 * t1 - t41 * t154 - t84 * (mrSges(7,1) * t208 + mrSges(7,2) * t207) - t110 * (mrSges(6,1) * t208 + mrSges(6,2) * t207) + t33 * t669 + t803 + ((-m(7) * t84 - t140) * t208 + t37) * pkin(5); -t207 * t154 + t208 * t156 + (-g(1) * t305 - g(2) * t301 - g(3) * t362 - t207 * t42 + t208 * t33 + t12) * m(7) + t31;];
tau  = t2;
