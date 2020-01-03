% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:49:10
% EndTime: 2019-12-31 22:49:53
% DurationCPUTime: 42.38s
% Computational Cost: add. (658545->345), mult. (1636442->467), div. (0->0), fcn. (1352335->14), ass. (0->155)
t839 = cos(pkin(5));
t833 = t839 * qJD(1) + qJD(2);
t836 = sin(pkin(6));
t838 = cos(pkin(6));
t837 = sin(pkin(5));
t848 = cos(qJ(2));
t867 = qJD(1) * t848;
t863 = t837 * t867;
t817 = (t833 * t836 + t838 * t863) * pkin(9);
t843 = sin(qJ(2));
t869 = qJD(1) * t837;
t882 = pkin(9) * t836;
t821 = (-pkin(2) * t848 - t843 * t882) * t869;
t866 = qJD(1) * qJD(2);
t827 = (qJDD(1) * t843 + t848 * t866) * t837;
t832 = t839 * qJDD(1) + qJDD(2);
t844 = sin(qJ(1));
t849 = cos(qJ(1));
t830 = t844 * g(1) - t849 * g(2);
t850 = qJD(1) ^ 2;
t883 = pkin(8) * t837;
t824 = qJDD(1) * pkin(1) + t850 * t883 + t830;
t831 = -t849 * g(1) - t844 * g(2);
t825 = -t850 * pkin(1) + qJDD(1) * t883 + t831;
t872 = t839 * t848;
t860 = t824 * t872 - t843 * t825;
t868 = qJD(1) * t843;
t881 = pkin(9) * t838;
t778 = -t827 * t881 + t832 * pkin(2) + t833 * t817 + (-g(3) * t848 - t821 * t868) * t837 + t860;
t864 = t837 * t868;
t820 = t833 * pkin(2) - t864 * t881;
t828 = (qJDD(1) * t848 - t843 * t866) * t837;
t858 = t828 * t838 + t832 * t836;
t873 = t839 * t843;
t870 = t824 * t873 + t848 * t825;
t779 = -t833 * t820 + (-g(3) * t843 + t821 * t867) * t837 + t858 * pkin(9) + t870;
t880 = t839 * g(3);
t784 = -t827 * t882 - t828 * pkin(2) - t880 + (-t824 + (-t817 * t848 + t820 * t843) * qJD(1)) * t837;
t842 = sin(qJ(3));
t847 = cos(qJ(3));
t757 = -t842 * t779 + (t778 * t838 + t784 * t836) * t847;
t874 = t838 * t848;
t879 = t836 * t842;
t808 = t833 * t879 + (t842 * t874 + t843 * t847) * t869;
t794 = -t808 * qJD(3) - t842 * t827 + t858 * t847;
t878 = t836 * t847;
t807 = (-t842 * t843 + t847 * t874) * t869 + t833 * t878;
t877 = t837 * t843;
t876 = t837 * t848;
t875 = t838 * t842;
t758 = t778 * t875 + t847 * t779 + t784 * t879;
t796 = -t807 * mrSges(4,1) + t808 * mrSges(4,2);
t818 = t838 * t833 - t836 * t863 + qJD(3);
t802 = t818 * mrSges(4,1) - t808 * mrSges(4,3);
t809 = -t836 * t828 + t838 * t832 + qJDD(3);
t797 = -t807 * pkin(3) - t808 * pkin(10);
t816 = t818 ^ 2;
t751 = -t816 * pkin(3) + t809 * pkin(10) + t807 * t797 + t758;
t763 = -t836 * t778 + t838 * t784;
t795 = t807 * qJD(3) + t847 * t827 + t858 * t842;
t753 = (-t807 * t818 - t795) * pkin(10) + (t808 * t818 - t794) * pkin(3) + t763;
t841 = sin(qJ(4));
t846 = cos(qJ(4));
t748 = t846 * t751 + t841 * t753;
t799 = -t841 * t808 + t846 * t818;
t800 = t846 * t808 + t841 * t818;
t781 = -t799 * pkin(4) - t800 * pkin(11);
t793 = qJDD(4) - t794;
t806 = qJD(4) - t807;
t805 = t806 ^ 2;
t745 = -t805 * pkin(4) + t793 * pkin(11) + t799 * t781 + t748;
t750 = -t809 * pkin(3) - t816 * pkin(10) + t808 * t797 - t757;
t767 = -t800 * qJD(4) - t841 * t795 + t846 * t809;
t768 = t799 * qJD(4) + t846 * t795 + t841 * t809;
t746 = (-t799 * t806 - t768) * pkin(11) + (t800 * t806 - t767) * pkin(4) + t750;
t840 = sin(qJ(5));
t845 = cos(qJ(5));
t742 = -t840 * t745 + t845 * t746;
t786 = -t840 * t800 + t845 * t806;
t756 = t786 * qJD(5) + t845 * t768 + t840 * t793;
t787 = t845 * t800 + t840 * t806;
t764 = -t786 * mrSges(6,1) + t787 * mrSges(6,2);
t766 = qJDD(5) - t767;
t798 = qJD(5) - t799;
t769 = -t798 * mrSges(6,2) + t786 * mrSges(6,3);
t739 = m(6) * t742 + t766 * mrSges(6,1) - t756 * mrSges(6,3) - t787 * t764 + t798 * t769;
t743 = t845 * t745 + t840 * t746;
t755 = -t787 * qJD(5) - t840 * t768 + t845 * t793;
t770 = t798 * mrSges(6,1) - t787 * mrSges(6,3);
t740 = m(6) * t743 - t766 * mrSges(6,2) + t755 * mrSges(6,3) + t786 * t764 - t798 * t770;
t733 = -t840 * t739 + t845 * t740;
t780 = -t799 * mrSges(5,1) + t800 * mrSges(5,2);
t789 = t806 * mrSges(5,1) - t800 * mrSges(5,3);
t731 = m(5) * t748 - t793 * mrSges(5,2) + t767 * mrSges(5,3) + t799 * t780 - t806 * t789 + t733;
t747 = -t841 * t751 + t846 * t753;
t744 = -t793 * pkin(4) - t805 * pkin(11) + t800 * t781 - t747;
t741 = -m(6) * t744 + t755 * mrSges(6,1) - t756 * mrSges(6,2) + t786 * t769 - t787 * t770;
t788 = -t806 * mrSges(5,2) + t799 * mrSges(5,3);
t737 = m(5) * t747 + t793 * mrSges(5,1) - t768 * mrSges(5,3) - t800 * t780 + t806 * t788 + t741;
t861 = t846 * t731 - t841 * t737;
t722 = m(4) * t758 - t809 * mrSges(4,2) + t794 * mrSges(4,3) + t807 * t796 - t818 * t802 + t861;
t725 = t841 * t731 + t846 * t737;
t801 = -t818 * mrSges(4,2) + t807 * mrSges(4,3);
t724 = m(4) * t763 - t794 * mrSges(4,1) + t795 * mrSges(4,2) - t807 * t801 + t808 * t802 + t725;
t732 = t845 * t739 + t840 * t740;
t853 = -m(5) * t750 + t767 * mrSges(5,1) - t768 * mrSges(5,2) + t799 * t788 - t800 * t789 - t732;
t728 = m(4) * t757 + t809 * mrSges(4,1) - t795 * mrSges(4,3) - t808 * t796 + t818 * t801 + t853;
t711 = t838 * t847 * t728 + t722 * t875 - t836 * t724;
t803 = -g(3) * t876 + t860;
t823 = -t833 * mrSges(3,2) + mrSges(3,3) * t863;
t826 = (-mrSges(3,1) * t848 + mrSges(3,2) * t843) * t869;
t707 = m(3) * t803 + t832 * mrSges(3,1) - t827 * mrSges(3,3) + t833 * t823 - t826 * t864 + t711;
t710 = t722 * t879 + t838 * t724 + t728 * t878;
t813 = -t837 * t824 - t880;
t822 = t833 * mrSges(3,1) - mrSges(3,3) * t864;
t709 = m(3) * t813 - t828 * mrSges(3,1) + t827 * mrSges(3,2) + (t822 * t843 - t823 * t848) * t869 + t710;
t716 = t847 * t722 - t842 * t728;
t804 = -g(3) * t877 + t870;
t715 = m(3) * t804 - t832 * mrSges(3,2) + t828 * mrSges(3,3) - t833 * t822 + t826 * t863 + t716;
t696 = t707 * t872 - t837 * t709 + t715 * t873;
t693 = m(2) * t830 + qJDD(1) * mrSges(2,1) - t850 * mrSges(2,2) + t696;
t701 = -t843 * t707 + t848 * t715;
t699 = m(2) * t831 - t850 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t701;
t871 = t849 * t693 + t844 * t699;
t695 = t707 * t876 + t839 * t709 + t715 * t877;
t862 = -t844 * t693 + t849 * t699;
t759 = Ifges(6,5) * t787 + Ifges(6,6) * t786 + Ifges(6,3) * t798;
t761 = Ifges(6,1) * t787 + Ifges(6,4) * t786 + Ifges(6,5) * t798;
t734 = -mrSges(6,1) * t744 + mrSges(6,3) * t743 + Ifges(6,4) * t756 + Ifges(6,2) * t755 + Ifges(6,6) * t766 - t787 * t759 + t798 * t761;
t760 = Ifges(6,4) * t787 + Ifges(6,2) * t786 + Ifges(6,6) * t798;
t735 = mrSges(6,2) * t744 - mrSges(6,3) * t742 + Ifges(6,1) * t756 + Ifges(6,4) * t755 + Ifges(6,5) * t766 + t786 * t759 - t798 * t760;
t771 = Ifges(5,5) * t800 + Ifges(5,6) * t799 + Ifges(5,3) * t806;
t772 = Ifges(5,4) * t800 + Ifges(5,2) * t799 + Ifges(5,6) * t806;
t717 = mrSges(5,2) * t750 - mrSges(5,3) * t747 + Ifges(5,1) * t768 + Ifges(5,4) * t767 + Ifges(5,5) * t793 - pkin(11) * t732 - t840 * t734 + t845 * t735 + t799 * t771 - t806 * t772;
t773 = Ifges(5,1) * t800 + Ifges(5,4) * t799 + Ifges(5,5) * t806;
t852 = mrSges(6,1) * t742 - mrSges(6,2) * t743 + Ifges(6,5) * t756 + Ifges(6,6) * t755 + Ifges(6,3) * t766 + t787 * t760 - t786 * t761;
t718 = -mrSges(5,1) * t750 + mrSges(5,3) * t748 + Ifges(5,4) * t768 + Ifges(5,2) * t767 + Ifges(5,6) * t793 - pkin(4) * t732 - t800 * t771 + t806 * t773 - t852;
t790 = Ifges(4,5) * t808 + Ifges(4,6) * t807 + Ifges(4,3) * t818;
t791 = Ifges(4,4) * t808 + Ifges(4,2) * t807 + Ifges(4,6) * t818;
t703 = mrSges(4,2) * t763 - mrSges(4,3) * t757 + Ifges(4,1) * t795 + Ifges(4,4) * t794 + Ifges(4,5) * t809 - pkin(10) * t725 + t846 * t717 - t841 * t718 + t807 * t790 - t818 * t791;
t792 = Ifges(4,1) * t808 + Ifges(4,4) * t807 + Ifges(4,5) * t818;
t851 = mrSges(5,1) * t747 - mrSges(5,2) * t748 + Ifges(5,5) * t768 + Ifges(5,6) * t767 + Ifges(5,3) * t793 + pkin(4) * t741 + pkin(11) * t733 + t845 * t734 + t840 * t735 + t800 * t772 - t799 * t773;
t704 = -mrSges(4,1) * t763 + mrSges(4,3) * t758 + Ifges(4,4) * t795 + Ifges(4,2) * t794 + Ifges(4,6) * t809 - pkin(3) * t725 - t808 * t790 + t818 * t792 - t851;
t855 = pkin(9) * t716 + t703 * t842 + t704 * t847;
t702 = mrSges(4,1) * t757 - mrSges(4,2) * t758 + Ifges(4,5) * t795 + Ifges(4,6) * t794 + Ifges(4,3) * t809 + pkin(3) * t853 + pkin(10) * t861 + t841 * t717 + t846 * t718 + t808 * t791 - t807 * t792;
t811 = Ifges(3,6) * t833 + (Ifges(3,4) * t843 + Ifges(3,2) * t848) * t869;
t812 = Ifges(3,5) * t833 + (Ifges(3,1) * t843 + Ifges(3,4) * t848) * t869;
t687 = mrSges(3,1) * t803 - mrSges(3,2) * t804 + Ifges(3,5) * t827 + Ifges(3,6) * t828 + Ifges(3,3) * t832 + pkin(2) * t711 + t838 * t702 + (t811 * t843 - t812 * t848) * t869 + t855 * t836;
t810 = Ifges(3,3) * t833 + (Ifges(3,5) * t843 + Ifges(3,6) * t848) * t869;
t689 = -mrSges(3,1) * t813 + mrSges(3,3) * t804 + Ifges(3,4) * t827 + Ifges(3,2) * t828 + Ifges(3,6) * t832 - pkin(2) * t710 - t836 * t702 - t810 * t864 + t833 * t812 + t855 * t838;
t691 = t810 * t863 + mrSges(3,2) * t813 - mrSges(3,3) * t803 + Ifges(3,1) * t827 + Ifges(3,4) * t828 + Ifges(3,5) * t832 + t847 * t703 - t842 * t704 - t833 * t811 + (-t710 * t836 - t711 * t838) * pkin(9);
t854 = mrSges(2,1) * t830 - mrSges(2,2) * t831 + Ifges(2,3) * qJDD(1) + pkin(1) * t696 + t839 * t687 + t689 * t876 + t691 * t877 + t701 * t883;
t685 = -mrSges(2,2) * g(3) - mrSges(2,3) * t830 + Ifges(2,5) * qJDD(1) - t850 * Ifges(2,6) - t843 * t689 + t848 * t691 + (-t695 * t837 - t696 * t839) * pkin(8);
t684 = mrSges(2,1) * g(3) + mrSges(2,3) * t831 + t850 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t695 - t837 * t687 + (pkin(8) * t701 + t689 * t848 + t691 * t843) * t839;
t1 = [-m(1) * g(1) + t862; -m(1) * g(2) + t871; (-m(1) - m(2)) * g(3) + t695; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t871 - t844 * t684 + t849 * t685; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t862 + t849 * t684 + t844 * t685; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t854; t854; t687; t702; t851; t852;];
tauJB = t1;
