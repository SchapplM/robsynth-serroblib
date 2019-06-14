% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 19:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:58:11
% EndTime: 2019-05-07 18:58:32
% DurationCPUTime: 10.70s
% Computational Cost: add. (164296->354), mult. (348698->434), div. (0->0), fcn. (269477->10), ass. (0->144)
t805 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t782 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t781 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t804 = Ifges(5,2) + Ifges(7,2) + Ifges(6,3);
t780 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t803 = -Ifges(5,3) - Ifges(6,2) - Ifges(7,3);
t752 = sin(pkin(6));
t756 = sin(qJ(2));
t758 = cos(qJ(2));
t783 = qJD(1) * qJD(2);
t739 = (-qJDD(1) * t758 + t756 * t783) * t752;
t753 = cos(pkin(6));
t749 = t753 * qJD(1) + qJD(2);
t755 = sin(qJ(3));
t785 = qJD(1) * t752;
t775 = t756 * t785;
t799 = cos(qJ(3));
t727 = t799 * t749 - t755 * t775;
t738 = (qJDD(1) * t756 + t758 * t783) * t752;
t748 = t753 * qJDD(1) + qJDD(2);
t707 = t727 * qJD(3) + t799 * t738 + t755 * t748;
t728 = t755 * t749 + t799 * t775;
t784 = qJD(1) * t758;
t774 = t752 * t784;
t744 = qJD(3) - t774;
t754 = sin(qJ(4));
t798 = cos(qJ(4));
t713 = t754 * t728 - t798 * t744;
t731 = qJDD(3) + t739;
t658 = -t713 * qJD(4) + t798 * t707 + t754 * t731;
t737 = (-pkin(2) * t758 - pkin(9) * t756) * t785;
t747 = t749 ^ 2;
t757 = sin(qJ(1));
t759 = cos(qJ(1));
t745 = t757 * g(1) - t759 * g(2);
t760 = qJD(1) ^ 2;
t797 = pkin(8) * t752;
t734 = qJDD(1) * pkin(1) + t760 * t797 + t745;
t746 = -t759 * g(1) - t757 * g(2);
t735 = -t760 * pkin(1) + qJDD(1) * t797 + t746;
t790 = t753 * t756;
t786 = t734 * t790 + t758 * t735;
t681 = -t747 * pkin(2) + t748 * pkin(9) + (-g(3) * t756 + t737 * t784) * t752 + t786;
t796 = t753 * g(3);
t682 = t739 * pkin(2) - t738 * pkin(9) - t796 + (-t734 + (pkin(2) * t756 - pkin(9) * t758) * t749 * qJD(1)) * t752;
t649 = -t755 * t681 + t799 * t682;
t711 = -t727 * pkin(3) - t728 * pkin(10);
t743 = t744 ^ 2;
t764 = t731 * pkin(3) + t743 * pkin(10) - t728 * t711 + t649;
t725 = qJD(4) - t727;
t794 = t713 * t725;
t802 = (-t658 + t794) * qJ(5) - t764;
t714 = t798 * t728 + t754 * t744;
t801 = -0.2e1 * t714;
t800 = 2 * qJD(5);
t795 = -mrSges(5,3) - mrSges(6,2);
t691 = t725 * mrSges(7,2) + t713 * mrSges(7,3);
t793 = t725 * t691;
t792 = t752 * t756;
t791 = t752 * t758;
t789 = t753 * t758;
t709 = -g(3) * t792 + t786;
t732 = t749 * mrSges(3,1) - mrSges(3,3) * t775;
t736 = (-mrSges(3,1) * t758 + mrSges(3,2) * t756) * t785;
t650 = t799 * t681 + t755 * t682;
t706 = -t728 * qJD(3) - t755 * t738 + t799 * t748;
t710 = -t727 * mrSges(4,1) + t728 * mrSges(4,2);
t716 = t744 * mrSges(4,1) - t728 * mrSges(4,3);
t646 = -t743 * pkin(3) + t731 * pkin(10) + t727 * t711 + t650;
t708 = -g(3) * t791 + t734 * t789 - t756 * t735;
t680 = -t748 * pkin(2) - t747 * pkin(9) + t737 * t775 - t708;
t648 = (-t727 * t744 - t707) * pkin(10) + (t728 * t744 - t706) * pkin(3) + t680;
t641 = -t754 * t646 + t798 * t648;
t692 = -t725 * mrSges(5,2) - t713 * mrSges(5,3);
t704 = qJDD(4) - t706;
t685 = t713 * pkin(4) - t714 * qJ(5);
t724 = t725 ^ 2;
t639 = -t704 * pkin(4) - t724 * qJ(5) + t714 * t685 + qJDD(5) - t641;
t690 = -t713 * mrSges(6,2) + t725 * mrSges(6,3);
t632 = qJD(6) * t801 + (-t658 - t794) * qJ(6) + (t713 * t714 - t704) * pkin(5) + t639;
t687 = -t713 * mrSges(7,1) + t714 * mrSges(7,2);
t769 = -m(7) * t632 + t658 * mrSges(7,3) + t714 * t687;
t763 = -m(6) * t639 + t704 * mrSges(6,1) + t725 * t690 + t769;
t686 = t713 * mrSges(6,1) - t714 * mrSges(6,3);
t787 = -t713 * mrSges(5,1) - t714 * mrSges(5,2) - t686;
t627 = m(5) * t641 + (t691 + t692) * t725 + t787 * t714 + (mrSges(5,1) + mrSges(7,1)) * t704 + t795 * t658 + t763;
t642 = t798 * t646 + t754 * t648;
t657 = t714 * qJD(4) + t754 * t707 - t798 * t731;
t694 = -t725 * mrSges(7,1) - t714 * mrSges(7,3);
t695 = t725 * mrSges(5,1) - t714 * mrSges(5,3);
t638 = -t724 * pkin(4) + t704 * qJ(5) - t713 * t685 + t725 * t800 + t642;
t696 = -t725 * mrSges(6,1) + t714 * mrSges(6,2);
t693 = -t725 * pkin(5) - t714 * qJ(6);
t712 = t713 ^ 2;
t634 = -t712 * pkin(5) + t657 * qJ(6) + 0.2e1 * qJD(6) * t713 + t725 * t693 + t638;
t779 = m(7) * t634 + t657 * mrSges(7,3) + t713 * t687;
t767 = m(6) * t638 + t704 * mrSges(6,3) + t725 * t696 + t779;
t629 = m(5) * t642 + (t694 - t695) * t725 + t787 * t713 + (-mrSges(5,2) + mrSges(7,2)) * t704 + t795 * t657 + t767;
t771 = -t754 * t627 + t798 * t629;
t623 = m(4) * t650 - t731 * mrSges(4,2) + t706 * mrSges(4,3) + t727 * t710 - t744 * t716 + t771;
t715 = -t744 * mrSges(4,2) + t727 * mrSges(4,3);
t640 = qJD(5) * t801 + (t714 * t725 + t657) * pkin(4) + t802;
t636 = -t712 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t657 + (-pkin(4) * t725 + t693 + t800) * t714 - t802;
t768 = -m(7) * t636 + t657 * mrSges(7,1) - t658 * mrSges(7,2) + t713 * t691 - t714 * t694;
t630 = m(6) * t640 + t657 * mrSges(6,1) - t658 * mrSges(6,3) + t713 * t690 - t714 * t696 + t768;
t761 = m(5) * t764 - t657 * mrSges(5,1) - t658 * mrSges(5,2) - t713 * t692 - t714 * t695 - t630;
t626 = m(4) * t649 + t731 * mrSges(4,1) - t707 * mrSges(4,3) - t728 * t710 + t744 * t715 + t761;
t772 = t799 * t623 - t755 * t626;
t613 = m(3) * t709 - t748 * mrSges(3,2) - t739 * mrSges(3,3) - t749 * t732 + t736 * t774 + t772;
t616 = t755 * t623 + t799 * t626;
t720 = -t752 * t734 - t796;
t733 = -t749 * mrSges(3,2) + mrSges(3,3) * t774;
t615 = m(3) * t720 + t739 * mrSges(3,1) + t738 * mrSges(3,2) + (t732 * t756 - t733 * t758) * t785 + t616;
t624 = t798 * t627 + t754 * t629;
t762 = -m(4) * t680 + t706 * mrSges(4,1) - t707 * mrSges(4,2) + t727 * t715 - t728 * t716 - t624;
t620 = m(3) * t708 + t748 * mrSges(3,1) - t738 * mrSges(3,3) + t749 * t733 - t736 * t775 + t762;
t602 = t613 * t790 - t752 * t615 + t620 * t789;
t600 = m(2) * t745 + qJDD(1) * mrSges(2,1) - t760 * mrSges(2,2) + t602;
t608 = t758 * t613 - t756 * t620;
t607 = m(2) * t746 - t760 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t608;
t788 = t759 * t600 + t757 * t607;
t601 = t613 * t792 + t753 * t615 + t620 * t791;
t778 = t713 * t780 - t714 * t781 + t725 * t803;
t777 = t713 * t804 - t714 * t782 - t725 * t780;
t776 = t782 * t713 - t714 * t805 - t781 * t725;
t773 = -t757 * t600 + t759 * t607;
t609 = mrSges(5,1) * t764 + mrSges(5,3) * t642 - mrSges(6,1) * t640 + mrSges(6,2) * t638 + mrSges(7,1) * t636 - mrSges(7,3) * t634 - pkin(5) * t768 - qJ(6) * t779 - pkin(4) * t630 + (-qJ(6) * t694 - t776) * t725 + t778 * t714 + (-qJ(6) * mrSges(7,2) + t780) * t704 + t782 * t658 - t804 * t657;
t631 = -t704 * mrSges(7,1) - t769 - t793;
t617 = -mrSges(5,2) * t764 + mrSges(6,2) * t639 + mrSges(7,2) * t636 - mrSges(5,3) * t641 - mrSges(6,3) * t640 - mrSges(7,3) * t632 - qJ(5) * t630 - qJ(6) * t631 - t782 * t657 + t658 * t805 + t781 * t704 + t778 * t713 + t777 * t725;
t700 = Ifges(4,5) * t728 + Ifges(4,6) * t727 + Ifges(4,3) * t744;
t701 = Ifges(4,4) * t728 + Ifges(4,2) * t727 + Ifges(4,6) * t744;
t603 = mrSges(4,2) * t680 - mrSges(4,3) * t649 + Ifges(4,1) * t707 + Ifges(4,4) * t706 + Ifges(4,5) * t731 - pkin(10) * t624 - t754 * t609 + t798 * t617 + t727 * t700 - t744 * t701;
t702 = Ifges(4,1) * t728 + Ifges(4,4) * t727 + Ifges(4,5) * t744;
t604 = -qJ(5) * (t725 * t694 + t767) - pkin(4) * (t763 + t793) + t744 * t702 - t728 * t700 + Ifges(4,6) * t731 + Ifges(4,2) * t706 + Ifges(4,4) * t707 - mrSges(4,1) * t680 + mrSges(4,3) * t650 + mrSges(5,2) * t642 + mrSges(6,1) * t639 - mrSges(5,1) * t641 - mrSges(6,3) * t638 - mrSges(7,2) * t634 + mrSges(7,1) * t632 + pkin(5) * t631 - pkin(3) * t624 + (pkin(4) * t686 + t777) * t714 + (qJ(5) * t686 + t776) * t713 + (-pkin(4) * mrSges(7,1) - qJ(5) * mrSges(7,2) + t803) * t704 + (pkin(4) * mrSges(6,2) - t781) * t658 + (qJ(5) * mrSges(6,2) + t780) * t657;
t717 = Ifges(3,3) * t749 + (Ifges(3,5) * t756 + Ifges(3,6) * t758) * t785;
t718 = Ifges(3,6) * t749 + (Ifges(3,4) * t756 + Ifges(3,2) * t758) * t785;
t597 = mrSges(3,2) * t720 - mrSges(3,3) * t708 + Ifges(3,1) * t738 - Ifges(3,4) * t739 + Ifges(3,5) * t748 - pkin(9) * t616 + t799 * t603 - t755 * t604 + t717 * t774 - t749 * t718;
t719 = Ifges(3,5) * t749 + (Ifges(3,1) * t756 + Ifges(3,4) * t758) * t785;
t598 = Ifges(3,4) * t738 - Ifges(3,2) * t739 + Ifges(3,6) * t748 - t717 * t775 + t749 * t719 - mrSges(3,1) * t720 + mrSges(3,3) * t709 - Ifges(4,5) * t707 - Ifges(4,6) * t706 - Ifges(4,3) * t731 - t728 * t701 + t727 * t702 - mrSges(4,1) * t649 + mrSges(4,2) * t650 - t754 * t617 - t798 * t609 - pkin(3) * t761 - pkin(10) * t771 - pkin(2) * t616;
t765 = pkin(8) * t608 + t597 * t756 + t598 * t758;
t596 = Ifges(3,5) * t738 - Ifges(3,6) * t739 + Ifges(3,3) * t748 + mrSges(3,1) * t708 - mrSges(3,2) * t709 + t755 * t603 + t799 * t604 + pkin(2) * t762 + pkin(9) * t772 + (t718 * t756 - t719 * t758) * t785;
t595 = -mrSges(2,2) * g(3) - mrSges(2,3) * t745 + Ifges(2,5) * qJDD(1) - t760 * Ifges(2,6) + t758 * t597 - t756 * t598 + (-t601 * t752 - t602 * t753) * pkin(8);
t594 = mrSges(2,1) * g(3) + mrSges(2,3) * t746 + t760 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t601 - t752 * t596 + t765 * t753;
t1 = [-m(1) * g(1) + t773; -m(1) * g(2) + t788; (-m(1) - m(2)) * g(3) + t601; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t788 - t757 * t594 + t759 * t595; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t773 + t759 * t594 + t757 * t595; -mrSges(1,1) * g(2) + mrSges(2,1) * t745 + mrSges(1,2) * g(1) - mrSges(2,2) * t746 + Ifges(2,3) * qJDD(1) + pkin(1) * t602 + t753 * t596 + t765 * t752;];
tauB  = t1;
