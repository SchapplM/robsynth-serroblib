% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 13:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:59:53
% EndTime: 2019-05-07 13:03:01
% DurationCPUTime: 136.39s
% Computational Cost: add. (2099257->416), mult. (5506761->555), div. (0->0), fcn. (4676163->16), ass. (0->180)
t790 = cos(pkin(6));
t783 = qJD(1) * t790 + qJD(2);
t786 = sin(pkin(7));
t789 = cos(pkin(7));
t787 = sin(pkin(6));
t799 = cos(qJ(2));
t818 = qJD(1) * t799;
t814 = t787 * t818;
t768 = (t783 * t786 + t789 * t814) * pkin(10);
t794 = sin(qJ(2));
t820 = qJD(1) * t787;
t834 = pkin(10) * t786;
t772 = (-pkin(2) * t799 - t794 * t834) * t820;
t816 = qJD(1) * qJD(2);
t778 = (qJDD(1) * t794 + t799 * t816) * t787;
t782 = qJDD(1) * t790 + qJDD(2);
t795 = sin(qJ(1));
t800 = cos(qJ(1));
t780 = t795 * g(1) - g(2) * t800;
t801 = qJD(1) ^ 2;
t835 = pkin(9) * t787;
t775 = qJDD(1) * pkin(1) + t801 * t835 + t780;
t781 = -g(1) * t800 - g(2) * t795;
t776 = -pkin(1) * t801 + qJDD(1) * t835 + t781;
t823 = t790 * t799;
t808 = t775 * t823 - t794 * t776;
t819 = qJD(1) * t794;
t833 = pkin(10) * t789;
t728 = -t778 * t833 + t782 * pkin(2) + t783 * t768 + (-g(3) * t799 - t772 * t819) * t787 + t808;
t815 = t787 * t819;
t771 = pkin(2) * t783 - t815 * t833;
t779 = (qJDD(1) * t799 - t794 * t816) * t787;
t807 = t779 * t789 + t782 * t786;
t824 = t790 * t794;
t821 = t775 * t824 + t799 * t776;
t729 = -t783 * t771 + (-g(3) * t794 + t772 * t818) * t787 + t807 * pkin(10) + t821;
t832 = t790 * g(3);
t733 = -t778 * t834 - t779 * pkin(2) - t832 + (-t775 + (-t768 * t799 + t771 * t794) * qJD(1)) * t787;
t793 = sin(qJ(3));
t798 = cos(qJ(3));
t826 = t789 * t798;
t830 = t786 * t798;
t694 = t728 * t826 - t793 * t729 + t733 * t830;
t825 = t789 * t799;
t758 = t783 * t830 + (-t793 * t794 + t798 * t825) * t820;
t744 = t758 * qJD(3) + t798 * t778 + t793 * t807;
t831 = t786 * t793;
t759 = t783 * t831 + (t793 * t825 + t794 * t798) * t820;
t760 = -t779 * t786 + t782 * t789 + qJDD(3);
t769 = t783 * t789 - t786 * t814 + qJD(3);
t684 = (t758 * t769 - t744) * qJ(4) + (t758 * t759 + t760) * pkin(3) + t694;
t827 = t789 * t793;
t695 = t728 * t827 + t798 * t729 + t733 * t831;
t743 = -t759 * qJD(3) - t793 * t778 + t798 * t807;
t753 = pkin(3) * t769 - qJ(4) * t759;
t757 = t758 ^ 2;
t687 = -pkin(3) * t757 + qJ(4) * t743 - t753 * t769 + t695;
t785 = sin(pkin(13));
t788 = cos(pkin(13));
t750 = t758 * t785 + t759 * t788;
t676 = -0.2e1 * qJD(4) * t750 + t788 * t684 - t785 * t687;
t829 = t787 * t794;
t828 = t787 * t799;
t749 = t758 * t788 - t759 * t785;
t677 = 0.2e1 * qJD(4) * t749 + t785 * t684 + t788 * t687;
t715 = t743 * t788 - t744 * t785;
t723 = -mrSges(5,1) * t749 + mrSges(5,2) * t750;
t738 = mrSges(5,1) * t769 - mrSges(5,3) * t750;
t724 = -pkin(4) * t749 - pkin(11) * t750;
t767 = t769 ^ 2;
t675 = -pkin(4) * t767 + pkin(11) * t760 + t724 * t749 + t677;
t707 = -t786 * t728 + t789 * t733;
t693 = -t743 * pkin(3) - t757 * qJ(4) + t759 * t753 + qJDD(4) + t707;
t716 = t743 * t785 + t744 * t788;
t679 = (-t749 * t769 - t716) * pkin(11) + (t750 * t769 - t715) * pkin(4) + t693;
t792 = sin(qJ(5));
t797 = cos(qJ(5));
t671 = t797 * t675 + t792 * t679;
t736 = t750 * t797 + t769 * t792;
t699 = -qJD(5) * t736 - t716 * t792 + t760 * t797;
t735 = -t750 * t792 + t769 * t797;
t709 = -mrSges(6,1) * t735 + mrSges(6,2) * t736;
t712 = qJDD(5) - t715;
t748 = qJD(5) - t749;
t718 = mrSges(6,1) * t748 - mrSges(6,3) * t736;
t710 = -pkin(5) * t735 - pkin(12) * t736;
t747 = t748 ^ 2;
t669 = -pkin(5) * t747 + pkin(12) * t712 + t710 * t735 + t671;
t674 = -t760 * pkin(4) - t767 * pkin(11) + t750 * t724 - t676;
t700 = qJD(5) * t735 + t716 * t797 + t760 * t792;
t672 = (-t735 * t748 - t700) * pkin(12) + (t736 * t748 - t699) * pkin(5) + t674;
t791 = sin(qJ(6));
t796 = cos(qJ(6));
t666 = -t669 * t791 + t672 * t796;
t713 = -t736 * t791 + t748 * t796;
t682 = qJD(6) * t713 + t700 * t796 + t712 * t791;
t714 = t736 * t796 + t748 * t791;
t696 = -mrSges(7,1) * t713 + mrSges(7,2) * t714;
t698 = qJDD(6) - t699;
t734 = qJD(6) - t735;
t701 = -mrSges(7,2) * t734 + mrSges(7,3) * t713;
t664 = m(7) * t666 + mrSges(7,1) * t698 - mrSges(7,3) * t682 - t696 * t714 + t701 * t734;
t667 = t669 * t796 + t672 * t791;
t681 = -qJD(6) * t714 - t700 * t791 + t712 * t796;
t702 = mrSges(7,1) * t734 - mrSges(7,3) * t714;
t665 = m(7) * t667 - mrSges(7,2) * t698 + mrSges(7,3) * t681 + t696 * t713 - t702 * t734;
t810 = -t664 * t791 + t796 * t665;
t657 = m(6) * t671 - mrSges(6,2) * t712 + mrSges(6,3) * t699 + t709 * t735 - t718 * t748 + t810;
t670 = -t675 * t792 + t679 * t797;
t717 = -mrSges(6,2) * t748 + mrSges(6,3) * t735;
t668 = -pkin(5) * t712 - pkin(12) * t747 + t710 * t736 - t670;
t804 = -m(7) * t668 + t681 * mrSges(7,1) - mrSges(7,2) * t682 + t713 * t701 - t702 * t714;
t662 = m(6) * t670 + mrSges(6,1) * t712 - mrSges(6,3) * t700 - t709 * t736 + t717 * t748 + t804;
t811 = t797 * t657 - t662 * t792;
t649 = m(5) * t677 - mrSges(5,2) * t760 + mrSges(5,3) * t715 + t723 * t749 - t738 * t769 + t811;
t737 = -mrSges(5,2) * t769 + mrSges(5,3) * t749;
t658 = t664 * t796 + t665 * t791;
t802 = -m(6) * t674 + t699 * mrSges(6,1) - mrSges(6,2) * t700 + t735 * t717 - t718 * t736 - t658;
t654 = m(5) * t676 + mrSges(5,1) * t760 - mrSges(5,3) * t716 - t723 * t750 + t737 * t769 + t802;
t644 = t785 * t649 + t788 * t654;
t751 = -mrSges(4,1) * t758 + mrSges(4,2) * t759;
t752 = -mrSges(4,2) * t769 + mrSges(4,3) * t758;
t642 = m(4) * t694 + mrSges(4,1) * t760 - mrSges(4,3) * t744 - t751 * t759 + t752 * t769 + t644;
t754 = mrSges(4,1) * t769 - mrSges(4,3) * t759;
t812 = t788 * t649 - t654 * t785;
t643 = m(4) * t695 - mrSges(4,2) * t760 + mrSges(4,3) * t743 + t751 * t758 - t754 * t769 + t812;
t652 = t792 * t657 + t797 * t662;
t803 = m(5) * t693 - mrSges(5,1) * t715 + t716 * mrSges(5,2) - t737 * t749 + t750 * t738 + t652;
t651 = m(4) * t707 - mrSges(4,1) * t743 + mrSges(4,2) * t744 - t752 * t758 + t754 * t759 + t803;
t629 = t642 * t826 + t643 * t827 - t651 * t786;
t755 = -g(3) * t828 + t808;
t774 = -mrSges(3,2) * t783 + mrSges(3,3) * t814;
t777 = (-mrSges(3,1) * t799 + mrSges(3,2) * t794) * t820;
t625 = m(3) * t755 + mrSges(3,1) * t782 - mrSges(3,3) * t778 + t774 * t783 - t777 * t815 + t629;
t628 = t642 * t830 + t643 * t831 + t789 * t651;
t764 = -t787 * t775 - t832;
t773 = mrSges(3,1) * t783 - mrSges(3,3) * t815;
t627 = m(3) * t764 - t779 * mrSges(3,1) + t778 * mrSges(3,2) + (t773 * t794 - t774 * t799) * t820 + t628;
t635 = -t642 * t793 + t798 * t643;
t756 = -g(3) * t829 + t821;
t634 = m(3) * t756 - mrSges(3,2) * t782 + mrSges(3,3) * t779 - t773 * t783 + t777 * t814 + t635;
t615 = t625 * t823 - t627 * t787 + t634 * t824;
t613 = m(2) * t780 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t801 + t615;
t621 = -t625 * t794 + t799 * t634;
t620 = m(2) * t781 - mrSges(2,1) * t801 - qJDD(1) * mrSges(2,2) + t621;
t822 = t800 * t613 + t795 * t620;
t614 = t625 * t828 + t790 * t627 + t634 * t829;
t813 = -t613 * t795 + t800 * t620;
t688 = Ifges(7,5) * t714 + Ifges(7,6) * t713 + Ifges(7,3) * t734;
t690 = Ifges(7,1) * t714 + Ifges(7,4) * t713 + Ifges(7,5) * t734;
t659 = -mrSges(7,1) * t668 + mrSges(7,3) * t667 + Ifges(7,4) * t682 + Ifges(7,2) * t681 + Ifges(7,6) * t698 - t688 * t714 + t690 * t734;
t689 = Ifges(7,4) * t714 + Ifges(7,2) * t713 + Ifges(7,6) * t734;
t660 = mrSges(7,2) * t668 - mrSges(7,3) * t666 + Ifges(7,1) * t682 + Ifges(7,4) * t681 + Ifges(7,5) * t698 + t688 * t713 - t689 * t734;
t703 = Ifges(6,5) * t736 + Ifges(6,6) * t735 + Ifges(6,3) * t748;
t704 = Ifges(6,4) * t736 + Ifges(6,2) * t735 + Ifges(6,6) * t748;
t645 = mrSges(6,2) * t674 - mrSges(6,3) * t670 + Ifges(6,1) * t700 + Ifges(6,4) * t699 + Ifges(6,5) * t712 - pkin(12) * t658 - t659 * t791 + t660 * t796 + t703 * t735 - t704 * t748;
t705 = Ifges(6,1) * t736 + Ifges(6,4) * t735 + Ifges(6,5) * t748;
t646 = -mrSges(6,1) * t674 - mrSges(7,1) * t666 + mrSges(7,2) * t667 + mrSges(6,3) * t671 + Ifges(6,4) * t700 - Ifges(7,5) * t682 + Ifges(6,2) * t699 + Ifges(6,6) * t712 - Ifges(7,6) * t681 - Ifges(7,3) * t698 - pkin(5) * t658 - t689 * t714 + t690 * t713 - t703 * t736 + t705 * t748;
t720 = Ifges(5,4) * t750 + Ifges(5,2) * t749 + Ifges(5,6) * t769;
t721 = Ifges(5,1) * t750 + Ifges(5,4) * t749 + Ifges(5,5) * t769;
t740 = Ifges(4,4) * t759 + Ifges(4,2) * t758 + Ifges(4,6) * t769;
t741 = Ifges(4,1) * t759 + Ifges(4,4) * t758 + Ifges(4,5) * t769;
t622 = Ifges(4,5) * t744 + Ifges(4,6) * t743 + t759 * t740 - t758 * t741 + mrSges(4,1) * t694 - mrSges(4,2) * t695 + Ifges(5,5) * t716 + Ifges(5,6) * t715 + t750 * t720 - t749 * t721 + mrSges(5,1) * t676 - mrSges(5,2) * t677 + t792 * t645 + t797 * t646 + pkin(4) * t802 + pkin(11) * t811 + pkin(3) * t644 + (Ifges(4,3) + Ifges(5,3)) * t760;
t761 = Ifges(3,3) * t783 + (Ifges(3,5) * t794 + Ifges(3,6) * t799) * t820;
t763 = Ifges(3,5) * t783 + (Ifges(3,1) * t794 + Ifges(3,4) * t799) * t820;
t719 = Ifges(5,5) * t750 + Ifges(5,6) * t749 + Ifges(5,3) * t769;
t630 = mrSges(5,2) * t693 - mrSges(5,3) * t676 + Ifges(5,1) * t716 + Ifges(5,4) * t715 + Ifges(5,5) * t760 - pkin(11) * t652 + t645 * t797 - t646 * t792 + t719 * t749 - t720 * t769;
t636 = Ifges(5,4) * t716 + Ifges(5,2) * t715 + Ifges(5,6) * t760 - t750 * t719 + t769 * t721 - mrSges(5,1) * t693 + mrSges(5,3) * t677 - Ifges(6,5) * t700 - Ifges(6,6) * t699 - Ifges(6,3) * t712 - t736 * t704 + t735 * t705 - mrSges(6,1) * t670 + mrSges(6,2) * t671 - t791 * t660 - t796 * t659 - pkin(5) * t804 - pkin(12) * t810 - pkin(4) * t652;
t739 = Ifges(4,5) * t759 + Ifges(4,6) * t758 + Ifges(4,3) * t769;
t616 = -mrSges(4,1) * t707 + mrSges(4,3) * t695 + Ifges(4,4) * t744 + Ifges(4,2) * t743 + Ifges(4,6) * t760 - pkin(3) * t803 + qJ(4) * t812 + t785 * t630 + t788 * t636 - t759 * t739 + t769 * t741;
t617 = mrSges(4,2) * t707 - mrSges(4,3) * t694 + Ifges(4,1) * t744 + Ifges(4,4) * t743 + Ifges(4,5) * t760 - qJ(4) * t644 + t630 * t788 - t636 * t785 + t739 * t758 - t740 * t769;
t805 = pkin(10) * t635 + t616 * t798 + t617 * t793;
t610 = -mrSges(3,1) * t764 + mrSges(3,3) * t756 + Ifges(3,4) * t778 + Ifges(3,2) * t779 + Ifges(3,6) * t782 - pkin(2) * t628 - t786 * t622 - t761 * t815 + t783 * t763 + t789 * t805;
t762 = Ifges(3,6) * t783 + (Ifges(3,4) * t794 + Ifges(3,2) * t799) * t820;
t611 = t761 * t814 + mrSges(3,2) * t764 - mrSges(3,3) * t755 + Ifges(3,1) * t778 + Ifges(3,4) * t779 + Ifges(3,5) * t782 - t793 * t616 + t798 * t617 - t783 * t762 + (-t628 * t786 - t629 * t789) * pkin(10);
t806 = pkin(9) * t621 + t610 * t799 + t611 * t794;
t609 = mrSges(3,1) * t755 - mrSges(3,2) * t756 + Ifges(3,5) * t778 + Ifges(3,6) * t779 + Ifges(3,3) * t782 + pkin(2) * t629 + t789 * t622 + (t762 * t794 - t763 * t799) * t820 + t805 * t786;
t608 = -mrSges(2,2) * g(3) - mrSges(2,3) * t780 + Ifges(2,5) * qJDD(1) - t801 * Ifges(2,6) - t794 * t610 + t799 * t611 + (-t614 * t787 - t615 * t790) * pkin(9);
t607 = mrSges(2,1) * g(3) + mrSges(2,3) * t781 + t801 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t614 - t787 * t609 + t790 * t806;
t1 = [-m(1) * g(1) + t813; -m(1) * g(2) + t822; (-m(1) - m(2)) * g(3) + t614; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t822 - t795 * t607 + t800 * t608; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t813 + t800 * t607 + t795 * t608; -mrSges(1,1) * g(2) + mrSges(2,1) * t780 + mrSges(1,2) * g(1) - mrSges(2,2) * t781 + Ifges(2,3) * qJDD(1) + pkin(1) * t615 + t790 * t609 + t787 * t806;];
tauB  = t1;
