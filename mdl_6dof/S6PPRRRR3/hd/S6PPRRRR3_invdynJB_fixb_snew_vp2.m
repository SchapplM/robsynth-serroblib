% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 21:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:14:32
% EndTime: 2019-05-04 21:15:16
% DurationCPUTime: 45.24s
% Computational Cost: add. (898715->306), mult. (1709194->416), div. (0->0), fcn. (1448434->18), ass. (0->150)
t805 = sin(pkin(13));
t810 = cos(pkin(13));
t798 = -t810 * g(1) - t805 * g(2);
t804 = sin(pkin(14));
t809 = cos(pkin(14));
t797 = t805 * g(1) - t810 * g(2);
t803 = -g(3) + qJDD(1);
t808 = sin(pkin(6));
t813 = cos(pkin(6));
t830 = t797 * t813 + t803 * t808;
t772 = -t804 * t798 + t830 * t809;
t773 = t809 * t798 + t830 * t804;
t785 = -t808 * t797 + t813 * t803 + qJDD(2);
t817 = sin(qJ(3));
t812 = cos(pkin(7));
t821 = cos(qJ(3));
t840 = t812 * t821;
t807 = sin(pkin(7));
t843 = t807 * t821;
t747 = t772 * t840 - t817 * t773 + t785 * t843;
t822 = qJD(3) ^ 2;
t806 = sin(pkin(8));
t849 = pkin(10) * t806;
t742 = qJDD(3) * pkin(3) + t822 * t849 + t747;
t841 = t812 * t817;
t844 = t807 * t817;
t748 = t772 * t841 + t821 * t773 + t785 * t844;
t743 = -t822 * pkin(3) + qJDD(3) * t849 + t748;
t756 = -t807 * t772 + t812 * t785;
t820 = cos(qJ(4));
t811 = cos(pkin(8));
t816 = sin(qJ(4));
t842 = t811 * t816;
t845 = t806 * t816;
t735 = t742 * t842 + t820 * t743 + t756 * t845;
t801 = t811 * qJD(3) + qJD(4);
t838 = qJD(3) * t806;
t836 = t816 * t838;
t789 = t801 * mrSges(5,1) - mrSges(5,3) * t836;
t791 = (-mrSges(5,1) * t820 + mrSges(5,2) * t816) * t838;
t837 = qJD(3) * qJD(4);
t794 = (-qJDD(3) * t820 + t816 * t837) * t806;
t800 = t811 * qJDD(3) + qJDD(4);
t792 = (-pkin(4) * t820 - pkin(11) * t816) * t838;
t799 = t801 ^ 2;
t835 = t820 * t838;
t733 = -t799 * pkin(4) + t800 * pkin(11) + t792 * t835 + t735;
t755 = t811 * t756;
t793 = (qJDD(3) * t816 + t820 * t837) * t806;
t737 = t794 * pkin(4) - t793 * pkin(11) + t755 + (-t742 + (pkin(4) * t816 - pkin(11) * t820) * t801 * qJD(3)) * t806;
t815 = sin(qJ(5));
t819 = cos(qJ(5));
t729 = t819 * t733 + t815 * t737;
t786 = t819 * t801 - t815 * t836;
t787 = t815 * t801 + t819 * t836;
t768 = -t786 * pkin(5) - t787 * pkin(12);
t788 = qJDD(5) + t794;
t796 = qJD(5) - t835;
t795 = t796 ^ 2;
t727 = -t795 * pkin(5) + t788 * pkin(12) + t786 * t768 + t729;
t734 = -t816 * t743 + (t742 * t811 + t756 * t806) * t820;
t732 = -t800 * pkin(4) - t799 * pkin(11) + t792 * t836 - t734;
t765 = -t787 * qJD(5) - t815 * t793 + t819 * t800;
t766 = t786 * qJD(5) + t819 * t793 + t815 * t800;
t730 = (-t786 * t796 - t766) * pkin(12) + (t787 * t796 - t765) * pkin(5) + t732;
t814 = sin(qJ(6));
t818 = cos(qJ(6));
t723 = -t814 * t727 + t818 * t730;
t774 = -t814 * t787 + t818 * t796;
t746 = t774 * qJD(6) + t818 * t766 + t814 * t788;
t775 = t818 * t787 + t814 * t796;
t753 = -t774 * mrSges(7,1) + t775 * mrSges(7,2);
t784 = qJD(6) - t786;
t757 = -t784 * mrSges(7,2) + t774 * mrSges(7,3);
t763 = qJDD(6) - t765;
t721 = m(7) * t723 + t763 * mrSges(7,1) - t746 * mrSges(7,3) - t775 * t753 + t784 * t757;
t724 = t818 * t727 + t814 * t730;
t745 = -t775 * qJD(6) - t814 * t766 + t818 * t788;
t758 = t784 * mrSges(7,1) - t775 * mrSges(7,3);
t722 = m(7) * t724 - t763 * mrSges(7,2) + t745 * mrSges(7,3) + t774 * t753 - t784 * t758;
t715 = -t814 * t721 + t818 * t722;
t767 = -t786 * mrSges(6,1) + t787 * mrSges(6,2);
t777 = t796 * mrSges(6,1) - t787 * mrSges(6,3);
t713 = m(6) * t729 - t788 * mrSges(6,2) + t765 * mrSges(6,3) + t786 * t767 - t796 * t777 + t715;
t728 = -t815 * t733 + t819 * t737;
t726 = -t788 * pkin(5) - t795 * pkin(12) + t787 * t768 - t728;
t725 = -m(7) * t726 + t745 * mrSges(7,1) - t746 * mrSges(7,2) + t774 * t757 - t775 * t758;
t776 = -t796 * mrSges(6,2) + t786 * mrSges(6,3);
t719 = m(6) * t728 + t788 * mrSges(6,1) - t766 * mrSges(6,3) - t787 * t767 + t796 * t776 + t725;
t833 = t819 * t713 - t815 * t719;
t704 = m(5) * t735 - t800 * mrSges(5,2) - t794 * mrSges(5,3) - t801 * t789 + t791 * t835 + t833;
t707 = t815 * t713 + t819 * t719;
t738 = -t806 * t742 + t755;
t790 = -t801 * mrSges(5,2) + mrSges(5,3) * t835;
t706 = m(5) * t738 + t794 * mrSges(5,1) + t793 * mrSges(5,2) + (t789 * t816 - t790 * t820) * t838 + t707;
t714 = t818 * t721 + t814 * t722;
t825 = -m(6) * t732 + t765 * mrSges(6,1) - t766 * mrSges(6,2) + t786 * t776 - t787 * t777 - t714;
t710 = m(5) * t734 + t800 * mrSges(5,1) - t793 * mrSges(5,3) + t801 * t790 - t791 * t836 + t825;
t846 = t710 * t820;
t693 = t704 * t842 - t806 * t706 + t811 * t846;
t689 = m(4) * t747 + qJDD(3) * mrSges(4,1) - t822 * mrSges(4,2) + t693;
t692 = t704 * t845 + t811 * t706 + t806 * t846;
t691 = m(4) * t756 + t692;
t698 = t820 * t704 - t816 * t710;
t697 = m(4) * t748 - t822 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t698;
t678 = t689 * t840 - t807 * t691 + t697 * t841;
t674 = m(3) * t772 + t678;
t683 = -t817 * t689 + t821 * t697;
t682 = m(3) * t773 + t683;
t850 = t674 * t809 + t682 * t804;
t677 = t689 * t843 + t812 * t691 + t697 * t844;
t676 = m(3) * t785 + t677;
t664 = -t808 * t676 + t850 * t813;
t662 = m(2) * t797 + t664;
t671 = -t804 * t674 + t809 * t682;
t670 = m(2) * t798 + t671;
t839 = t810 * t662 + t805 * t670;
t663 = t813 * t676 + t850 * t808;
t834 = -t805 * t662 + t810 * t670;
t832 = m(2) * t803 + t663;
t749 = Ifges(7,5) * t775 + Ifges(7,6) * t774 + Ifges(7,3) * t784;
t751 = Ifges(7,1) * t775 + Ifges(7,4) * t774 + Ifges(7,5) * t784;
t716 = -mrSges(7,1) * t726 + mrSges(7,3) * t724 + Ifges(7,4) * t746 + Ifges(7,2) * t745 + Ifges(7,6) * t763 - t775 * t749 + t784 * t751;
t750 = Ifges(7,4) * t775 + Ifges(7,2) * t774 + Ifges(7,6) * t784;
t717 = mrSges(7,2) * t726 - mrSges(7,3) * t723 + Ifges(7,1) * t746 + Ifges(7,4) * t745 + Ifges(7,5) * t763 + t774 * t749 - t784 * t750;
t759 = Ifges(6,5) * t787 + Ifges(6,6) * t786 + Ifges(6,3) * t796;
t760 = Ifges(6,4) * t787 + Ifges(6,2) * t786 + Ifges(6,6) * t796;
t699 = mrSges(6,2) * t732 - mrSges(6,3) * t728 + Ifges(6,1) * t766 + Ifges(6,4) * t765 + Ifges(6,5) * t788 - pkin(12) * t714 - t814 * t716 + t818 * t717 + t786 * t759 - t796 * t760;
t761 = Ifges(6,1) * t787 + Ifges(6,4) * t786 + Ifges(6,5) * t796;
t824 = mrSges(7,1) * t723 - mrSges(7,2) * t724 + Ifges(7,5) * t746 + Ifges(7,6) * t745 + Ifges(7,3) * t763 + t775 * t750 - t774 * t751;
t700 = -mrSges(6,1) * t732 + mrSges(6,3) * t729 + Ifges(6,4) * t766 + Ifges(6,2) * t765 + Ifges(6,6) * t788 - pkin(5) * t714 - t787 * t759 + t796 * t761 - t824;
t781 = Ifges(5,6) * t801 + (Ifges(5,4) * t816 + Ifges(5,2) * t820) * t838;
t782 = Ifges(5,5) * t801 + (Ifges(5,1) * t816 + Ifges(5,4) * t820) * t838;
t684 = Ifges(5,5) * t793 - Ifges(5,6) * t794 + Ifges(5,3) * t800 + mrSges(5,1) * t734 - mrSges(5,2) * t735 + t815 * t699 + t819 * t700 + pkin(4) * t825 + pkin(11) * t833 + (t781 * t816 - t782 * t820) * t838;
t780 = Ifges(5,3) * t801 + (Ifges(5,5) * t816 + Ifges(5,6) * t820) * t838;
t685 = mrSges(5,2) * t738 - mrSges(5,3) * t734 + Ifges(5,1) * t793 - Ifges(5,4) * t794 + Ifges(5,5) * t800 - pkin(11) * t707 + t819 * t699 - t815 * t700 + t780 * t835 - t801 * t781;
t823 = mrSges(6,1) * t728 - mrSges(6,2) * t729 + Ifges(6,5) * t766 + Ifges(6,6) * t765 + Ifges(6,3) * t788 + pkin(5) * t725 + pkin(12) * t715 + t818 * t716 + t814 * t717 + t787 * t760 - t786 * t761;
t686 = -mrSges(5,1) * t738 + mrSges(5,3) * t735 + Ifges(5,4) * t793 - Ifges(5,2) * t794 + Ifges(5,6) * t800 - pkin(4) * t707 - t780 * t836 + t801 * t782 - t823;
t827 = pkin(10) * t698 + t685 * t816 + t686 * t820;
t666 = -mrSges(4,1) * t756 + mrSges(4,3) * t748 + t822 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t692 - t806 * t684 + t827 * t811;
t667 = mrSges(4,2) * t756 - mrSges(4,3) * t747 + Ifges(4,5) * qJDD(3) - t822 * Ifges(4,6) + t820 * t685 - t816 * t686 + (-t692 * t806 - t693 * t811) * pkin(10);
t828 = pkin(9) * t683 + t666 * t821 + t667 * t817;
t665 = mrSges(4,1) * t747 - mrSges(4,2) * t748 + Ifges(4,3) * qJDD(3) + pkin(3) * t693 + t811 * t684 + t827 * t806;
t659 = -mrSges(3,1) * t785 + mrSges(3,3) * t773 - pkin(2) * t677 - t807 * t665 + t828 * t812;
t660 = mrSges(3,2) * t785 - mrSges(3,3) * t772 - t817 * t666 + t821 * t667 + (-t677 * t807 - t678 * t812) * pkin(9);
t826 = qJ(2) * t671 + t659 * t809 + t660 * t804;
t658 = mrSges(3,1) * t772 - mrSges(3,2) * t773 + pkin(2) * t678 + t812 * t665 + t828 * t807;
t657 = mrSges(2,2) * t803 - mrSges(2,3) * t797 - t804 * t659 + t809 * t660 + (-t663 * t808 - t664 * t813) * qJ(2);
t656 = -mrSges(2,1) * t803 + mrSges(2,3) * t798 - pkin(1) * t663 - t808 * t658 + t826 * t813;
t1 = [-m(1) * g(1) + t834; -m(1) * g(2) + t839; -m(1) * g(3) + t832; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t839 - t805 * t656 + t810 * t657; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t834 + t810 * t656 + t805 * t657; -mrSges(1,1) * g(2) + mrSges(2,1) * t797 + mrSges(1,2) * g(1) - mrSges(2,2) * t798 + pkin(1) * t664 + t813 * t658 + t826 * t808; t832; t676; t665; t684; t823; t824;];
tauJB  = t1;
