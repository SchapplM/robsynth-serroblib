% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:33
% EndTime: 2019-03-09 11:05:07
% DurationCPUTime: 18.91s
% Computational Cost: add. (37015->928), mult. (88163->1273), div. (0->0), fcn. (98882->10), ass. (0->471)
t544 = sin(qJ(6));
t538 = t544 ^ 2;
t547 = cos(qJ(6));
t539 = t547 ^ 2;
t672 = t538 + t539;
t855 = Ifges(6,4) - Ifges(5,5);
t854 = Ifges(6,5) - Ifges(5,6);
t810 = -mrSges(7,3) / 0.2e1;
t853 = t672 * t810;
t540 = sin(pkin(11));
t542 = cos(pkin(11));
t543 = cos(pkin(6));
t541 = sin(pkin(6));
t546 = sin(qJ(2));
t698 = t541 * t546;
t492 = t540 * t543 + t542 * t698;
t545 = sin(qJ(4));
t664 = t540 * t698;
t696 = t542 * t543;
t591 = -t664 + t696;
t767 = cos(qJ(4));
t569 = t492 * t767 + t545 * t591;
t373 = t492 * t545 - t767 * t591;
t848 = qJ(5) * t373;
t234 = pkin(4) * t569 + t848;
t548 = cos(qJ(2));
t495 = (pkin(2) * t546 - qJ(3) * t548) * t541;
t529 = pkin(8) * t698;
t765 = pkin(1) * t548;
t496 = t543 * t765 - t529;
t383 = t542 * t495 - t496 * t540;
t695 = t542 * t548;
t307 = (pkin(3) * t546 - pkin(9) * t695) * t541 + t383;
t384 = t540 * t495 + t542 * t496;
t697 = t541 * t548;
t663 = t540 * t697;
t335 = -pkin(9) * t663 + t384;
t170 = t545 * t307 + t767 * t335;
t154 = -qJ(5) * t698 - t170;
t509 = t540 * t767 + t545 * t542;
t456 = t509 * t697;
t136 = -pkin(5) * t456 - t154;
t624 = -t307 * t767 + t545 * t335;
t155 = -pkin(4) * t698 + t624;
t397 = t456 * t547 - t544 * t698;
t398 = t456 * t544 + t547 * t698;
t656 = t767 * t542;
t457 = -t545 * t663 + t656 * t697;
t724 = t457 * mrSges(6,1);
t411 = mrSges(6,2) * t698 + t724;
t517 = Ifges(7,5) * t547 - Ifges(7,6) * t544;
t752 = Ifges(7,4) * t544;
t521 = Ifges(7,1) * t547 - t752;
t804 = pkin(4) + pkin(10);
t128 = t457 * pkin(5) - t698 * t804 + t624;
t497 = t543 * t546 * pkin(1) + pkin(8) * t697;
t453 = pkin(3) * t663 + t497;
t598 = -qJ(5) * t457 + t453;
t193 = t456 * t804 + t598;
t78 = t128 * t547 - t193 * t544;
t79 = t128 * t544 + t193 * t547;
t615 = t79 * t544 + t78 * t547;
t666 = Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t751 = Ifges(7,4) * t547;
t519 = -Ifges(7,2) * t544 + t751;
t779 = -t519 / 0.4e1;
t713 = t547 * mrSges(7,2);
t717 = t544 * mrSges(7,1);
t516 = t713 + t717;
t781 = -t516 / 0.2e1;
t801 = -t155 / 0.2e1;
t811 = mrSges(6,3) / 0.2e1;
t817 = -m(7) / 0.2e1;
t820 = -m(6) / 0.2e1;
t846 = mrSges(5,1) / 0.2e1;
t852 = -(t517 / 0.4e1 + t666) * t457 + (-pkin(4) * t155 - qJ(5) * t154) * t820 + (qJ(5) * t136 - t615 * t804) * t817 + pkin(4) * t411 / 0.2e1 + t136 * t781 + t154 * t811 + mrSges(6,2) * t801 + t624 * t846 + t170 * mrSges(5,2) / 0.2e1 + t397 * t779 - t398 * t521 / 0.4e1;
t639 = m(7) * t672;
t597 = t639 / 0.4e1 + m(6) / 0.4e1;
t851 = 0.2e1 * t597 * t569;
t790 = t397 / 0.2e1;
t789 = t398 / 0.2e1;
t850 = t457 / 0.2e1;
t761 = mrSges(5,3) + mrSges(6,1);
t849 = Ifges(5,4) + Ifges(6,6);
t786 = -t456 / 0.2e1;
t785 = t456 / 0.2e1;
t792 = t569 / 0.2e1;
t650 = -t697 / 0.2e1;
t842 = -t698 / 0.2e1;
t651 = t698 / 0.2e1;
t840 = -mrSges(5,1) + mrSges(6,2);
t839 = mrSges(5,2) - mrSges(6,3);
t838 = Ifges(6,1) + Ifges(5,3);
t640 = t538 / 0.2e1 + t539 / 0.2e1;
t837 = t640 * mrSges(7,3);
t836 = mrSges(6,3) + t516;
t835 = t540 * t492 + t542 * t591;
t507 = t540 * t545 - t656;
t834 = t853 * t507;
t832 = t507 * t855 + t854 * t509;
t831 = t373 * t855 + t854 * t569;
t701 = t507 * t544;
t414 = mrSges(7,1) * t509 - mrSges(7,3) * t701;
t684 = t547 * t414;
t700 = t507 * t547;
t416 = -mrSges(7,2) * t509 + mrSges(7,3) * t700;
t689 = t544 * t416;
t830 = t684 + t689;
t731 = t569 * mrSges(7,1);
t315 = -t544 * t373 + t547 * t697;
t733 = t315 * mrSges(7,3);
t208 = t731 + t733;
t686 = t547 * t208;
t730 = t569 * mrSges(7,2);
t314 = t373 * t547 + t544 * t697;
t735 = t314 * mrSges(7,3);
t207 = -t730 + t735;
t693 = t544 * t207;
t829 = t686 + t693;
t828 = mrSges(4,1) * t540 + mrSges(4,2) * t542;
t743 = Ifges(7,6) * t547;
t749 = Ifges(7,5) * t544;
t619 = t743 + t749;
t827 = -t619 / 0.4e1 + t749 / 0.2e1 + t743 / 0.2e1;
t826 = Ifges(7,5) * t789 + Ifges(7,6) * t790 + Ifges(7,3) * t850;
t687 = t547 * t207;
t692 = t544 * t208;
t588 = t692 / 0.2e1 - t687 / 0.2e1;
t773 = -t547 / 0.2e1;
t774 = t544 / 0.2e1;
t825 = (t314 * t773 + t315 * t774) * mrSges(7,3) - t588;
t824 = 0.2e1 * t509;
t823 = 2 * qJD(4);
t822 = m(4) / 0.2e1;
t821 = m(5) / 0.2e1;
t819 = m(6) / 0.2e1;
t816 = m(7) / 0.2e1;
t815 = -mrSges(5,2) / 0.2e1;
t814 = -mrSges(6,2) / 0.2e1;
t813 = -mrSges(7,2) / 0.2e1;
t808 = Ifges(4,5) / 0.2e1;
t806 = t78 / 0.2e1;
t600 = -pkin(2) * t548 - qJ(3) * t546 - pkin(1);
t572 = -t548 * pkin(3) + t542 * t600;
t570 = t572 * t541;
t625 = qJ(3) * t543 + t497;
t593 = t540 * t625;
t574 = -t492 * pkin(9) - t593;
t559 = t570 + t574;
t583 = t600 * t541;
t337 = t540 * t583 + t542 * t625;
t276 = pkin(9) * t591 + t337;
t657 = t767 * t276;
t764 = t373 * pkin(5);
t557 = t545 * t559 + t657 - t764;
t682 = t548 * qJ(5);
t662 = t541 * t682;
t95 = t557 - t662;
t805 = t95 / 0.2e1;
t803 = -qJ(5) / 0.2e1;
t802 = qJ(5) / 0.2e1;
t800 = t207 / 0.2e1;
t706 = qJ(5) * t507;
t347 = t509 * t804 + t706;
t760 = pkin(9) + qJ(3);
t513 = t760 * t542;
t638 = t760 * t540;
t436 = t513 * t767 - t545 * t638;
t363 = -t507 * pkin(5) + t436;
t211 = -t347 * t544 + t363 * t547;
t799 = t211 / 0.2e1;
t798 = t234 / 0.2e1;
t797 = t314 / 0.2e1;
t796 = -t315 / 0.2e1;
t526 = mrSges(6,3) * t697;
t732 = t373 * mrSges(6,1);
t329 = t526 + t732;
t795 = t329 / 0.2e1;
t794 = t363 / 0.2e1;
t793 = -t569 / 0.2e1;
t714 = t547 * mrSges(7,1);
t716 = t544 * mrSges(7,2);
t515 = -t714 + t716;
t401 = t515 * t507;
t788 = t401 / 0.2e1;
t787 = t416 / 0.2e1;
t784 = -t507 / 0.2e1;
t783 = -t509 / 0.2e1;
t782 = t509 / 0.2e1;
t780 = -t517 / 0.2e1;
t778 = -t540 / 0.2e1;
t777 = t542 / 0.2e1;
t776 = -t544 / 0.2e1;
t775 = -t544 / 0.4e1;
t772 = t547 / 0.2e1;
t771 = t547 / 0.4e1;
t768 = t804 / 0.2e1;
t421 = pkin(4) * t509 + t706;
t766 = m(6) * t421;
t534 = -pkin(3) * t542 - pkin(2);
t395 = pkin(3) * t664 + t529 + (t534 - t765) * t543;
t566 = -qJ(5) * t569 + t395;
t109 = t373 * t804 + t566;
t134 = t276 * t545 - t767 * t559;
t100 = -pkin(5) * t569 - t134;
t530 = pkin(4) * t697;
t89 = pkin(10) * t697 - t100 + t530;
t58 = -t109 * t544 + t547 * t89;
t763 = t58 * mrSges(7,3);
t59 = t109 * t547 + t544 * t89;
t762 = t59 * mrSges(7,3);
t757 = Ifges(4,4) * t540;
t756 = Ifges(4,4) * t542;
t755 = Ifges(5,4) * t569;
t754 = Ifges(5,4) * t509;
t753 = Ifges(7,4) * t315;
t748 = Ifges(4,2) * t540;
t747 = Ifges(4,6) * t540;
t746 = Ifges(6,6) * t569;
t745 = Ifges(6,6) * t509;
t742 = Ifges(7,3) * t373;
t740 = Ifges(7,3) * t507;
t599 = -qJ(5) * t509 + t534;
t322 = t507 * t804 + t599;
t435 = t513 * t545 + t767 * t638;
t361 = pkin(5) * t509 + t435;
t198 = -t322 * t544 + t361 * t547;
t739 = t198 * mrSges(7,3);
t199 = t322 * t547 + t361 * t544;
t738 = t199 * mrSges(7,3);
t567 = -t545 * t574 - t657;
t123 = (-t545 * t572 + t682) * t541 + t567;
t124 = t134 + t530;
t728 = t569 * Ifges(7,6);
t126 = Ifges(7,2) * t314 + t728 - t753;
t312 = Ifges(7,4) * t314;
t729 = t569 * Ifges(7,5);
t127 = -Ifges(7,1) * t315 + t312 + t729;
t135 = t545 * t570 - t567;
t151 = t373 * pkin(4) + t566;
t734 = t315 * mrSges(7,2);
t736 = t314 * mrSges(7,1);
t172 = -t734 - t736;
t191 = Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t457;
t192 = Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t457;
t236 = -mrSges(6,2) * t373 - mrSges(6,3) * t569;
t247 = -mrSges(7,1) * t397 + mrSges(7,2) * t398;
t254 = pkin(4) * t456 + t598;
t277 = -mrSges(7,2) * t457 + mrSges(7,3) * t397;
t278 = mrSges(7,1) * t457 - mrSges(7,3) * t398;
t723 = t457 * mrSges(5,2);
t726 = t456 * mrSges(5,1);
t316 = t723 + t726;
t722 = t457 * mrSges(6,3);
t725 = t456 * mrSges(6,2);
t317 = -t722 - t725;
t330 = mrSges(6,1) * t569 - mrSges(6,2) * t697;
t331 = mrSges(5,2) * t697 - t373 * mrSges(5,3);
t332 = -mrSges(5,1) * t697 - mrSges(5,3) * t569;
t336 = t542 * t583 - t593;
t410 = mrSges(6,1) * t456 - mrSges(6,3) * t698;
t412 = -mrSges(5,2) * t698 - mrSges(5,3) * t456;
t413 = mrSges(5,1) * t698 - mrSges(5,3) * t457;
t433 = (Ifges(4,6) * t546 + (-t748 + t756) * t548) * t541;
t451 = mrSges(4,2) * t697 + mrSges(4,3) * t591;
t452 = -mrSges(4,1) * t697 - mrSges(4,3) * t492;
t465 = t828 * t697;
t469 = t529 + (-pkin(2) - t765) * t543;
t485 = (-mrSges(4,3) * t540 * t548 - mrSges(4,2) * t546) * t541;
t486 = (mrSges(4,1) * t546 - mrSges(4,3) * t695) * t541;
t527 = Ifges(3,5) * t697;
t582 = t453 * mrSges(5,2) + Ifges(6,4) * t842 + Ifges(5,5) * t651 + t826 + (Ifges(5,1) + Ifges(6,2)) * t850 + t849 * t786;
t601 = t453 * mrSges(5,1) + Ifges(6,5) * t651 + Ifges(5,6) * t842 - t849 * t457 / 0.2e1 + (Ifges(5,2) + Ifges(6,3)) * t785;
t125 = -Ifges(7,5) * t315 + Ifges(7,6) * t314 + Ifges(7,3) * t569;
t364 = Ifges(6,6) * t373;
t219 = -Ifges(6,4) * t697 - Ifges(6,2) * t569 + t364;
t369 = Ifges(5,4) * t373;
t221 = Ifges(5,1) * t569 - Ifges(5,5) * t697 - t369;
t623 = t221 / 0.2e1 + t125 / 0.2e1 - t219 / 0.2e1;
t665 = Ifges(5,6) / 0.2e1 - Ifges(6,5) / 0.2e1;
t631 = t665 * t456;
t632 = t497 * mrSges(4,2) + (Ifges(4,5) * t546 + (t542 * Ifges(4,1) - t757) * t548) * t541 / 0.2e1;
t633 = t497 * mrSges(4,1) - t433 / 0.2e1;
t636 = Ifges(4,6) * t777 - Ifges(3,6);
t218 = -Ifges(6,5) * t697 + Ifges(6,3) * t373 - t746;
t220 = -Ifges(5,2) * t373 - Ifges(5,6) * t697 + t755;
t645 = t218 / 0.2e1 - t220 / 0.2e1;
t721 = t496 * mrSges(3,2);
t3 = m(5) * (t134 * t624 + t135 * t170 + t395 * t453) - t624 * t332 + t127 * t789 + t126 * t790 + t192 * t796 + t191 * t797 + (t433 * t777 - t721 + t527 / 0.2e1 + (-t542 * mrSges(4,1) - mrSges(3,1)) * t497) * t543 + t645 * t456 + t632 * t492 + m(6) * (t123 * t154 + t124 * t155 + t151 * t254) + m(7) * (t136 * t95 + t58 * t78 + t59 * t79) + m(4) * (t336 * t383 + t337 * t384 + t469 * t497) + t623 * t457 + t601 * t373 + t337 * t485 + t336 * t486 + t469 * t465 + t384 * t451 + t383 * t452 + t123 * t410 + t124 * t411 + t135 * t412 - t134 * t413 + t395 * t316 + t154 * t329 + t155 * t330 + t170 * t331 + t151 * t317 + t582 * t569 + t59 * t277 + t58 * t278 + t254 * t236 + t95 * t247 + t79 * t207 + t78 * t208 + t136 * t172 + ((t492 * t808 + t636 * t543 + t633 * t540 + (-pkin(1) * mrSges(3,1) + (-Ifges(3,4) - t747 / 0.2e1) * t546) * t541 + t666 * t569 - t665 * t373) * t546 + ((Ifges(4,1) * t492 + Ifges(4,4) * t696) * t777 + (Ifges(4,4) * t492 + Ifges(4,2) * t696) * t778 + Ifges(3,4) * t697 + Ifges(3,5) * t543 / 0.2e1 + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t542 + t747) * t548) * t541 - t666 * t457 + t631 + (Ifges(3,1) - Ifges(3,2) - Ifges(4,3) + (-t756 / 0.2e1 + t748 / 0.2e1) * t540 - t838) * t698) * t548) * t541;
t737 = t3 * qJD(1);
t140 = t569 * t619 - t742;
t620 = Ifges(7,2) * t547 + t752;
t141 = -Ifges(7,6) * t373 + t569 * t620;
t621 = Ifges(7,1) * t544 + t751;
t142 = -Ifges(7,5) * t373 + t569 * t621;
t228 = t515 * t569;
t715 = t544 * mrSges(7,3);
t231 = -mrSges(7,1) * t373 - t569 * t715;
t712 = t547 * mrSges(7,3);
t232 = mrSges(7,2) * t373 + t569 * t712;
t233 = -mrSges(6,2) * t569 + t373 * mrSges(6,3);
t235 = mrSges(5,1) * t569 - t373 * mrSges(5,2);
t237 = Ifges(6,3) * t569 + t364;
t238 = Ifges(6,2) * t373 + t746;
t239 = -Ifges(5,2) * t569 - t369;
t240 = -Ifges(5,1) * t373 - t755;
t678 = t330 - t332;
t679 = t329 - t331;
t688 = t547 * t126;
t694 = t544 * t127;
t101 = t135 - t764;
t145 = t569 * t804 + t848;
t70 = t101 * t547 - t145 * t544;
t71 = t101 * t544 + t145 * t547;
t4 = t395 * t235 + t142 * t796 + t141 * t797 + t58 * t231 + t59 * t232 + t151 * t233 + t234 * t236 + t95 * t228 + t71 * t207 + t70 * t208 + t100 * t172 + t678 * t135 + t679 * t134 + m(6) * (t123 * t134 + t124 * t135 + t151 * t234) + m(7) * (t100 * t95 + t58 * t70 + t59 * t71) + (-t134 * mrSges(5,3) - t124 * mrSges(6,1) + t237 / 0.2e1 - t239 / 0.2e1 - t623) * t373 + (-t135 * mrSges(5,3) + t123 * mrSges(6,1) - t238 / 0.2e1 + t240 / 0.2e1 + t140 / 0.2e1 + t694 / 0.2e1 + t688 / 0.2e1 + t645) * t569 + t831 * t650;
t727 = t4 * qJD(1);
t720 = t507 * mrSges(6,1);
t719 = t509 * Ifges(7,5);
t718 = t509 * Ifges(7,6);
t171 = -mrSges(7,1) * t315 + mrSges(7,2) * t314;
t173 = Ifges(7,5) * t314 + Ifges(7,6) * t315;
t174 = Ifges(7,2) * t315 + t312;
t175 = Ifges(7,1) * t314 + t753;
t9 = -t59 * t208 + t95 * t171 + t58 * t207 + t173 * t792 + (t762 - t175 / 0.2e1 + t126 / 0.2e1) * t315 + (-t763 + t127 / 0.2e1 + t174 / 0.2e1) * t314;
t709 = t9 * qJD(1);
t618 = t544 * t59 + t547 * t58;
t12 = t591 * t451 - t492 * t452 - (t172 - t679) * t373 + (t678 + t829) * t569 + m(7) * (-t373 * t95 + t569 * t618) + m(6) * (t123 * t373 + t124 * t569) + m(5) * (t134 * t569 - t135 * t373) + m(4) * (-t336 * t492 + t337 * t591);
t705 = qJD(1) * t12;
t617 = t544 * t58 - t547 * t59;
t17 = (-m(6) * t123 + m(7) * t95 + t172 - t329) * t697 + (m(6) * t151 - m(7) * t617 + t236 + t687 - t692) * t569;
t704 = qJD(1) * t17;
t703 = t211 * t547;
t212 = t347 * t547 + t363 * t544;
t702 = t212 * t544;
t490 = Ifges(7,4) * t700;
t327 = Ifges(7,1) * t701 + t490 + t719;
t691 = t544 * t327;
t690 = t544 * t414;
t325 = t507 * t620 + t718;
t685 = t547 * t325;
t683 = t547 * t416;
t681 = t804 * t277;
t680 = t804 * t414;
t671 = qJD(4) * t544;
t670 = qJD(4) * t547;
t669 = -mrSges(6,1) / 0.2e1 - mrSges(5,3) / 0.2e1;
t668 = t814 + t846;
t667 = t811 + t815;
t648 = t804 * t772;
t647 = -t126 / 0.4e1 + t175 / 0.4e1;
t646 = t127 / 0.4e1 + t174 / 0.4e1;
t406 = t507 * t521;
t644 = t325 / 0.4e1 - t406 / 0.4e1;
t405 = -Ifges(7,2) * t701 + t490;
t643 = t327 / 0.4e1 + t405 / 0.4e1;
t642 = t621 / 0.4e1 + t519 / 0.4e1;
t641 = t521 / 0.4e1 - t620 / 0.4e1;
t637 = t672 * t804;
t630 = -t639 / 0.2e1;
t404 = Ifges(7,5) * t700 - Ifges(7,6) * t701;
t627 = t208 / 0.2e1 - t733 / 0.2e1;
t626 = -t735 / 0.2e1 + t800;
t323 = Ifges(7,3) * t509 + t507 * t619;
t500 = Ifges(6,6) * t507;
t427 = -Ifges(6,2) * t509 + t500;
t505 = Ifges(5,4) * t507;
t431 = Ifges(5,1) * t509 - t505;
t622 = t431 / 0.2e1 + t323 / 0.2e1 - t427 / 0.2e1;
t616 = t544 * t71 + t547 * t70;
t324 = t619 * t509 - t740;
t326 = -Ifges(7,6) * t507 + t509 * t620;
t328 = -Ifges(7,5) * t507 + t509 * t621;
t400 = pkin(4) * t507 + t599;
t402 = t515 * t509;
t415 = -mrSges(7,1) * t507 - t509 * t715;
t417 = mrSges(7,2) * t507 + t509 * t712;
t498 = t507 * mrSges(6,3);
t420 = -t509 * mrSges(6,2) + t498;
t499 = t507 * mrSges(5,2);
t422 = t509 * mrSges(5,1) - t499;
t423 = -mrSges(6,2) * t507 - mrSges(6,3) * t509;
t424 = Ifges(6,3) * t509 + t500;
t425 = Ifges(6,3) * t507 - t745;
t426 = Ifges(6,2) * t507 + t745;
t428 = -Ifges(5,2) * t509 - t505;
t429 = -Ifges(5,2) * t507 + t754;
t430 = -Ifges(5,1) * t507 - t754;
t15 = t198 * t415 + t199 * t417 + t211 * t414 + t212 * t416 - t361 * t401 + t363 * t402 + t421 * t423 + t534 * t422 + m(7) * (t198 * t211 + t199 * t212 - t361 * t363) + (-t429 / 0.2e1 - t426 / 0.2e1 + t430 / 0.2e1 + t425 / 0.2e1 + t324 / 0.2e1 + t691 / 0.2e1 + t685 / 0.2e1) * t509 + (-t428 / 0.2e1 + t424 / 0.2e1 + t328 * t774 + t326 * t772 - t622) * t507 + (t420 + t766) * t400;
t587 = t691 / 0.4e1 + t685 / 0.4e1;
t589 = t694 / 0.4e1 + t688 / 0.4e1;
t550 = (t100 * t363 + t198 * t70 + t199 * t71 + t211 * t58 + t212 * t59 - t361 * t95) * t816 - t361 * t172 / 0.2e1 + t71 * t787 + t100 * t788 + t228 * t794 + t423 * t798 + t208 * t799 + t212 * t800 + (t795 - t331 / 0.2e1) * t435 + (-t239 / 0.4e1 - t221 / 0.4e1 - t125 / 0.4e1 + t237 / 0.4e1 + t219 / 0.4e1 + t544 * t142 / 0.4e1 + t141 * t771 + (t134 / 0.2e1 - t124 / 0.2e1) * mrSges(6,1)) * t507 + (t427 / 0.4e1 - t428 / 0.4e1 - t431 / 0.4e1 + t424 / 0.4e1 - t323 / 0.4e1 + t669 * t435) * t373 + (-t220 / 0.4e1 + t240 / 0.4e1 + t218 / 0.4e1 + t140 / 0.4e1 - t238 / 0.4e1 + (t123 / 0.2e1 + t135 / 0.2e1) * mrSges(6,1) + t589) * t509 + t534 * t235 / 0.2e1 + t58 * t415 / 0.2e1 + t59 * t417 / 0.2e1 + t151 * t420 / 0.2e1 + t421 * t236 / 0.2e1 + t395 * t422 / 0.2e1 + t70 * t414 / 0.2e1 + t400 * t233 / 0.2e1 + (t330 / 0.2e1 - t332 / 0.2e1) * t436 + t314 * t326 / 0.4e1 - t315 * t328 / 0.4e1 + (-t429 / 0.4e1 + t430 / 0.4e1 + t425 / 0.4e1 - t426 / 0.4e1 + t324 / 0.4e1 + t669 * t436 + t587) * t569 + t198 * t231 / 0.2e1 + t199 * t232 / 0.2e1 + t402 * t805 + (t151 * t421 + t234 * t400 + (t124 - t134) * t436 + (t123 + t135) * t435) * t819;
t573 = -t832 * t548 / 0.4e1;
t2 = ((-Ifges(6,1) / 0.2e1 - Ifges(5,3) / 0.2e1) * t546 + t573) * t541 + (-t192 / 0.4e1 + mrSges(7,3) * t806 + t278 * t768) * t547 + t631 + (t191 / 0.4e1 + t79 * mrSges(7,3) / 0.2e1 + t681 / 0.2e1) * t544 + (t410 / 0.2e1 - t247 / 0.2e1) * qJ(5) + t550 + t852;
t614 = t2 * qJD(1) + t15 * qJD(2);
t403 = t507 * t516;
t24 = t363 * t403 + t404 * t782 + t198 * t416 - t199 * t414 + ((-t739 + t405 / 0.2e1 + t327 / 0.2e1) * t547 + (-t738 - t325 / 0.2e1 + t406 / 0.2e1) * t544) * t507;
t553 = ((-t763 / 0.2e1 + t646) * t547 + (-t762 / 0.2e1 + t647) * t544) * t507 + (-t739 / 0.2e1 + t643) * t314 + (t738 / 0.2e1 + t644) * t315 + t198 * t800 - t199 * t208 / 0.2e1 + t171 * t794 + t569 * t404 / 0.4e1 + t509 * t173 / 0.4e1 + t58 * t787 - t59 * t414 / 0.2e1 + t403 * t805;
t568 = mrSges(7,1) * t806 + t79 * t813 + t826;
t5 = t553 - t568;
t613 = t5 * qJD(1) + t24 * qJD(2);
t607 = t198 * t547 + t199 * t544;
t37 = (t507 * t761 - t401) * t507 + (t509 * t761 + t830) * t509 + m(7) * (-t363 * t507 + t509 * t607) + (m(4) * qJ(3) + mrSges(4,3)) * (t540 ^ 2 + t542 ^ 2) + (m(6) + m(5)) * (t435 * t509 - t436 * t507);
t604 = -t373 * t436 + t435 * t569;
t551 = (qJ(3) * t835 - t540 * t336 + t542 * t337) * t822 + (t134 * t509 - t135 * t507 + t604) * t821 + (t123 * t507 + t124 * t509 + t604) * t819 + (-t363 * t373 - t507 * t95 + t509 * t618 + t569 * t607) * t816 - t373 * t788 + t507 * t795 + t332 * t783 + t452 * t778 + t451 * t777 + t835 * mrSges(4,3) / 0.2e1 + t830 * t792 + (t172 + t331) * t784 + (t330 + t829) * t782 + t761 * (-t373 * t784 + t569 * t782);
t555 = -m(4) * t497 / 0.2e1 - m(5) * t453 / 0.2e1 + t254 * t820 + (-t544 * t78 + t547 * t79) * t817 - t726 / 0.2e1 + t725 / 0.2e1 - t723 / 0.2e1 + t722 / 0.2e1 + t278 * t774 + t277 * t773 + t828 * t650;
t7 = t551 + t555;
t612 = -qJD(1) * t7 - qJD(2) * t37;
t586 = t690 / 0.2e1 - t683 / 0.2e1;
t606 = t198 * t544 - t199 * t547;
t554 = (-t236 / 0.2e1 + t588) * t509 + (-t423 / 0.2e1 + t586) * t569 + (-t151 * t509 - t400 * t569 - t436 * t697) * t819 + (-t363 * t697 + t509 * t617 + t569 * t606) * t816;
t561 = m(6) * t801 + t615 * t817 - t724 / 0.2e1 + t277 * t776 + t278 * t773;
t592 = (-t401 / 0.2e1 + t720 / 0.2e1) * t548;
t14 = (t546 * t814 + t592) * t541 + t554 + t561;
t48 = (-m(6) * t400 + m(7) * t606 - t423 - t683 + t690) * t509;
t611 = -qJD(1) * t14 - qJD(2) * t48;
t558 = t668 * t569 + t667 * t373 + m(6) * t798 + (-t544 * t70 + t547 * t71) * t816 + t231 * t776 + t232 * t772;
t562 = -t234 * t820 + (-t569 * t637 - t848) * t817;
t20 = -(t781 - t667) * t373 + (t837 + t668) * t569 + t558 + t562;
t560 = -t766 / 0.2e1 + (-t509 * t637 - t706) * t816 + t507 * t781;
t564 = t766 / 0.2e1 + (-t211 * t544 + t212 * t547) * t816 + t415 * t776 + t417 * t772;
t44 = t499 - t498 + (-t837 + t840) * t509 + t560 - t564;
t610 = qJD(1) * t20 - qJD(2) * t44;
t27 = (mrSges(7,1) * t792 + t627) * t547 + (mrSges(7,2) * t793 + t626) * t544;
t595 = t716 / 0.2e1 - t714 / 0.2e1;
t580 = t595 * t509;
t585 = -t689 / 0.2e1 - t684 / 0.2e1;
t68 = -t580 - t585;
t609 = qJD(1) * t27 + qJD(2) * t68;
t28 = (t730 / 0.2e1 - t626) * t547 + (t731 / 0.2e1 + t627) * t544;
t594 = t713 / 0.2e1 + t717 / 0.2e1;
t579 = t594 * t509;
t60 = t507 * t837 + t579 + t586;
t608 = qJD(1) * t28 + qJD(2) * t60;
t605 = t702 + t703;
t167 = (m(7) * t640 + t819) * t824;
t72 = 0.2e1 * t851;
t602 = qJD(1) * t72 + qJD(2) * t167;
t596 = mrSges(7,1) * t799 + t212 * t813;
t590 = t403 * t803 + t515 * t794;
t584 = t415 * t773 + t417 * t776;
t563 = t642 * t315 + t641 * t314 + t171 * t802 - t95 * t515 / 0.2e1;
t575 = t742 / 0.2e1 - t70 * mrSges(7,1) / 0.2e1 + t71 * mrSges(7,2) / 0.2e1;
t10 = (-0.3e1 / 0.4e1 * t728 - t626 * t804 + t647) * t547 + (-0.3e1 / 0.4e1 * t729 + t627 * t804 - t646) * t544 + t563 + t575;
t25 = (-Ifges(7,3) / 0.2e1 - t804 * t837) * t507 + (0.3e1 / 0.4e1 * t718 + t416 * t768 - t641 * t507 + t644) * t547 + (0.3e1 / 0.4e1 * t719 - t680 / 0.2e1 + t642 * t507 + t643) * t544 + t590 + t596;
t279 = -qJ(5) * t515 + (-t519 / 0.2e1 - t621 / 0.2e1) * t547 + (t620 / 0.2e1 - t521 / 0.2e1) * t544;
t577 = t10 * qJD(1) - t25 * qJD(2) + t279 * qJD(4);
t525 = -0.2e1 * t662;
t552 = -t526 + (t525 + t135) * t819 + (t525 + t557) * t816 - t736 / 0.2e1 - t734 / 0.2e1 + t516 * t650;
t565 = t135 * t820 + t231 * t773 + t232 * t776 + t616 * t817;
t22 = t552 + t565;
t494 = (m(6) + m(7)) * qJ(5) + t836;
t62 = t595 * t507 + 0.2e1 * (t363 / 0.4e1 - t703 / 0.4e1 - t702 / 0.4e1) * m(7) + t584;
t576 = qJD(1) * t22 + qJD(2) * t62 + qJD(4) * t494;
t571 = qJD(4) * (-qJ(5) * mrSges(6,1) + t519 * t772 + t521 * t774);
t168 = m(6) * t783 + t509 * t630 + t597 * t824;
t73 = t851 + (t630 + t820) * t569;
t69 = -t580 + t585;
t61 = t579 - t586 + t834;
t53 = t605 * t816 + m(7) * t794 + m(6) * t436 + (-mrSges(6,1) + t595) * t507 - t584;
t46 = -t509 * t837 + t560 + t564;
t30 = -t686 / 0.2e1 + t715 * t797 - t693 / 0.2e1 + t315 * t712 / 0.2e1 - t595 * t569;
t29 = t594 * t569 + t825;
t26 = t406 * t771 + t405 * t775 - t416 * t648 - t680 * t776 - t740 / 0.2e1 - t587 - t590 + t596 + (t779 - t621 / 0.4e1) * t701 + (-t620 + t521) * t700 / 0.4e1 - t834 * t804 + t827 * t509;
t21 = t552 - t565 - t732;
t19 = mrSges(5,1) * t793 + mrSges(6,2) * t792 + t558 - t562 + t853 * t569 + (-t815 - t836 / 0.2e1) * t373;
t13 = mrSges(6,2) * t651 + t541 * t592 + t554 - t561;
t11 = t174 * t775 + t175 * t771 + t569 * t827 - t804 * t825 + t563 - t575 - t589;
t8 = t551 - t555;
t6 = t553 + t568;
t1 = Ifges(6,5) * t785 + Ifges(5,6) * t786 + t191 * t775 + t192 * t771 + t247 * t802 - t278 * t648 + t410 * t803 + t573 * t541 + t615 * t810 + t838 * t651 - t681 * t774 + t550 - t852;
t16 = [qJD(2) * t3 + qJD(3) * t12 + qJD(4) * t4 - qJD(5) * t17 + qJD(6) * t9, t8 * qJD(3) + t1 * qJD(4) + t13 * qJD(5) + t6 * qJD(6) + t737 + ((t155 * mrSges(6,1) + mrSges(5,3) * t624 + t582) * t509 + 0.2e1 * (t170 * t436 + t435 * t624 + t453 * t534) * t821 + t527 + t327 * t789 + t325 * t790 + t425 * t785 + t429 * t786 + t534 * t316 - t497 * mrSges(3,1) - pkin(2) * t465 - t435 * t413 + t436 * t412 - t436 * t410 + t435 * t411 + t79 * t416 + t254 * t423 + t78 * t414 + t400 * t317 + t136 * t401 + t363 * t247 + t199 * t277 + t198 * t278 + (t384 * mrSges(4,3) + qJ(3) * t485 - t633) * t542 + (-t383 * mrSges(4,3) - qJ(3) * t486 + t632) * t540 + t622 * t457 + (t154 * mrSges(6,1) - t170 * mrSges(5,3) + t191 * t772 + t192 * t774 + t601) * t507 + ((-t507 * t665 + t509 * t666 + t540 * t808 + t636) * t546 + ((Ifges(4,1) * t540 + t756) * t777 + (Ifges(4,2) * t542 + t757) * t778) * t548) * t541 + 0.2e1 * (t136 * t363 + t198 * t78 + t199 * t79) * t816 + 0.2e1 * (-t154 * t436 + t155 * t435 + t254 * t400) * t819 + 0.2e1 * (-pkin(2) * t497 + (-t383 * t540 + t384 * t542) * qJ(3)) * t822 - t721) * qJD(2), qJD(2) * t8 + qJD(4) * t19 + qJD(5) * t73 + qJD(6) * t30 + t705, t727 + t1 * qJD(2) + t19 * qJD(3) + (pkin(4) * t732 + qJ(5) * t228 + t100 * t516 + t134 * t839 + t135 * t840 + t373 * t780 + t831) * qJD(4) + t21 * qJD(5) + t11 * qJD(6) + (t142 / 0.2e1 - t804 * t231 - t70 * mrSges(7,3)) * t670 + (-t141 / 0.2e1 - t71 * mrSges(7,3) - t804 * t232) * t671 + ((qJ(5) * t100 - t616 * t804) * t816 + (-pkin(4) * t135 - qJ(5) * t134) * t819) * t823 + t569 * t571, qJD(2) * t13 + qJD(3) * t73 + qJD(4) * t21 + qJD(6) * t29 - t704, t709 + t6 * qJD(2) + t30 * qJD(3) + t11 * qJD(4) + t29 * qJD(5) + (-mrSges(7,1) * t59 - mrSges(7,2) * t58 + t173) * qJD(6); qJD(3) * t7 + qJD(4) * t2 + qJD(5) * t14 + qJD(6) * t5 - t737, qJD(3) * t37 + qJD(4) * t15 + qJD(5) * t48 + qJD(6) * t24, qJD(4) * t46 + qJD(5) * t168 + qJD(6) * t69 - t612, t46 * qJD(3) + (pkin(4) * t720 + qJ(5) * t402 - t361 * t516 + t435 * t839 + t436 * t840 + t507 * t780 + t832) * qJD(4) + t53 * qJD(5) + t26 * qJD(6) + (t328 / 0.2e1 - t804 * t415 - t211 * mrSges(7,3)) * t670 + (-t326 / 0.2e1 - t804 * t417 - t212 * mrSges(7,3)) * t671 + ((-qJ(5) * t361 - t605 * t804) * t816 + (-pkin(4) * t436 - qJ(5) * t435) * t819) * t823 + t509 * t571 + t614, qJD(3) * t168 + qJD(4) * t53 + qJD(6) * t61 - t611, t69 * qJD(3) + t26 * qJD(4) + t61 * qJD(5) + (-mrSges(7,1) * t199 - mrSges(7,2) * t198 + t404) * qJD(6) + t613; -qJD(2) * t7 + qJD(4) * t20 - qJD(5) * t72 - qJD(6) * t27 - t705, -qJD(4) * t44 - qJD(5) * t167 - qJD(6) * t68 + t612, 0, t610, -t602, qJD(6) * t515 - t609; -qJD(2) * t2 - qJD(3) * t20 + qJD(5) * t22 + qJD(6) * t10 - t727, qJD(3) * t44 + qJD(5) * t62 - qJD(6) * t25 - t614, -t610, qJD(5) * t494 + qJD(6) * t279, t576 ((mrSges(7,2) * t804 - Ifges(7,6)) * t547 + (mrSges(7,1) * t804 - Ifges(7,5)) * t544) * qJD(6) + t577; -qJD(2) * t14 + qJD(3) * t72 - qJD(4) * t22 - qJD(6) * t28 + t704, qJD(3) * t167 - qJD(4) * t62 - qJD(6) * t60 + t611, t602, -t576, 0, -qJD(6) * t516 - t608; -qJD(2) * t5 + qJD(3) * t27 - qJD(4) * t10 + qJD(5) * t28 - t709, qJD(3) * t68 + qJD(4) * t25 + qJD(5) * t60 - t613, t609, -t577, t608, 0;];
Cq  = t16;
