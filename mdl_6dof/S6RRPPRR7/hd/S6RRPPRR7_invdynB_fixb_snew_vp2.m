% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:03:22
% EndTime: 2019-05-06 11:03:35
% DurationCPUTime: 8.02s
% Computational Cost: add. (84736->367), mult. (195261->448), div. (0->0), fcn. (132493->10), ass. (0->155)
t879 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t856 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t855 = Ifges(5,5) + Ifges(3,6) - Ifges(4,6);
t878 = Ifges(4,2) + Ifges(3,3) + Ifges(5,3);
t854 = Ifges(5,6) + Ifges(3,5) + Ifges(4,4);
t877 = -Ifges(4,3) - Ifges(5,1) - Ifges(3,2);
t817 = cos(pkin(6));
t809 = qJD(1) * t817 + qJD(2);
t876 = t809 ^ 2;
t875 = pkin(3) + pkin(9);
t874 = pkin(2) * t809;
t820 = sin(qJ(2));
t873 = g(3) * t820;
t808 = qJDD(1) * t817 + qJDD(2);
t872 = t808 * pkin(2);
t871 = mrSges(3,3) + mrSges(4,2);
t816 = sin(pkin(6));
t826 = qJD(1) ^ 2;
t870 = t816 ^ 2 * t826;
t869 = t816 * t820;
t824 = cos(qJ(2));
t868 = t816 * t824;
t867 = t817 * t820;
t866 = t817 * t824;
t821 = sin(qJ(1));
t825 = cos(qJ(1));
t803 = t821 * g(1) - g(2) * t825;
t782 = pkin(8) * t816 * t826 + qJDD(1) * pkin(1) + t803;
t750 = -t817 * g(3) - t816 * t782;
t861 = qJD(1) * t816;
t847 = t820 * t861;
t778 = -mrSges(4,1) * t809 + mrSges(4,2) * t847;
t860 = qJD(1) * t824;
t846 = t816 * t860;
t780 = -mrSges(3,2) * t809 + mrSges(3,3) * t846;
t781 = mrSges(4,2) * t846 + mrSges(4,3) * t809;
t789 = (qJD(2) * t860 + qJDD(1) * t820) * t816;
t857 = qJDD(1) * t816;
t790 = -qJD(2) * t847 + t824 * t857;
t842 = t809 * t846;
t841 = -t790 * pkin(2) + t750 + (-t789 - t842) * qJ(3);
t716 = (-(2 * qJD(3)) + t874) * t847 + t841;
t775 = -pkin(3) * t809 - qJ(4) * t847;
t852 = t824 ^ 2 * t870;
t830 = -qJ(4) * t852 + qJDD(4) - t841 + ((2 * qJD(3)) + t775) * t847;
t702 = t789 * pkin(4) + t875 * t790 + (pkin(4) * t824 + (-pkin(2) - pkin(9)) * t820) * t809 * t861 + t830;
t788 = (pkin(4) * t820 + pkin(9) * t824) * t861;
t804 = -g(1) * t825 - g(2) * t821;
t783 = -pkin(1) * t826 + pkin(8) * t857 + t804;
t733 = -g(3) * t868 + t782 * t866 - t820 * t783;
t784 = (-pkin(2) * t824 - qJ(3) * t820) * t861;
t834 = -qJ(3) * t876 + t784 * t847 + qJDD(3) - t733;
t858 = qJD(1) * qJD(4);
t843 = -0.2e1 * t816 * t858;
t828 = t820 * t843 + t834 + (-t789 + t842) * qJ(4);
t851 = t824 * t870;
t706 = -t876 * pkin(4) + (-pkin(3) * t851 - t788 * t861) * t820 + (-pkin(2) - t875) * t808 + t828;
t819 = sin(qJ(5));
t823 = cos(qJ(5));
t699 = t819 * t702 + t823 * t706;
t765 = -t809 * t819 - t823 * t846;
t731 = -qJD(5) * t765 + t790 * t819 - t808 * t823;
t764 = -t809 * t823 + t819 * t846;
t735 = -mrSges(6,1) * t764 + mrSges(6,2) * t765;
t795 = qJD(5) + t847;
t740 = mrSges(6,1) * t795 - mrSges(6,3) * t765;
t773 = qJDD(5) + t789;
t736 = -pkin(5) * t764 - pkin(10) * t765;
t793 = t795 ^ 2;
t697 = -pkin(5) * t793 + pkin(10) * t773 + t736 * t764 + t699;
t859 = qJD(3) * t809;
t794 = 0.2e1 * t859;
t863 = t782 * t867 + t824 * t783;
t840 = pkin(2) * t876 - t808 * qJ(3) - t784 * t846 - t863;
t831 = pkin(3) * t852 + t790 * qJ(4) - t809 * t775 + t840;
t704 = t808 * pkin(4) - t876 * pkin(9) + t794 + t824 * t843 + (-t788 * t860 - t873) * t816 - t831;
t732 = qJD(5) * t764 - t790 * t823 - t808 * t819;
t700 = t704 + (-t764 * t795 - t732) * pkin(10) + (t765 * t795 - t731) * pkin(5);
t818 = sin(qJ(6));
t822 = cos(qJ(6));
t694 = -t697 * t818 + t700 * t822;
t737 = -t765 * t818 + t795 * t822;
t712 = qJD(6) * t737 + t732 * t822 + t773 * t818;
t738 = t765 * t822 + t795 * t818;
t722 = -mrSges(7,1) * t737 + mrSges(7,2) * t738;
t762 = qJD(6) - t764;
t723 = -mrSges(7,2) * t762 + mrSges(7,3) * t737;
t729 = qJDD(6) - t731;
t692 = m(7) * t694 + mrSges(7,1) * t729 - mrSges(7,3) * t712 - t722 * t738 + t723 * t762;
t695 = t697 * t822 + t700 * t818;
t711 = -qJD(6) * t738 - t732 * t818 + t773 * t822;
t724 = mrSges(7,1) * t762 - mrSges(7,3) * t738;
t693 = m(7) * t695 - mrSges(7,2) * t729 + mrSges(7,3) * t711 + t722 * t737 - t724 * t762;
t844 = -t692 * t818 + t822 * t693;
t683 = m(6) * t699 - mrSges(6,2) * t773 + mrSges(6,3) * t731 + t735 * t764 - t740 * t795 + t844;
t698 = t702 * t823 - t706 * t819;
t739 = -mrSges(6,2) * t795 + mrSges(6,3) * t764;
t696 = -pkin(5) * t773 - pkin(10) * t793 + t736 * t765 - t698;
t837 = -m(7) * t696 + t711 * mrSges(7,1) - mrSges(7,2) * t712 + t737 * t723 - t724 * t738;
t688 = m(6) * t698 + mrSges(6,1) * t773 - mrSges(6,3) * t732 - t735 * t765 + t739 * t795 + t837;
t677 = t819 * t683 + t823 * t688;
t709 = t790 * pkin(3) - t847 * t874 + t830;
t779 = -mrSges(5,1) * t809 + mrSges(5,3) * t846;
t836 = -m(5) * t709 - t789 * mrSges(5,1) + t779 * t846 - t677;
t833 = m(4) * t716 - t790 * mrSges(4,1) + t836;
t776 = mrSges(5,2) * t809 - mrSges(5,3) * t847;
t862 = -mrSges(3,1) * t809 + mrSges(3,3) * t847 + t776;
t671 = m(3) * t750 + (-mrSges(3,1) + mrSges(5,2)) * t790 + (mrSges(3,2) - mrSges(4,3)) * t789 + ((-t780 - t781) * t824 + (-t778 - t862) * t820) * t861 + t833;
t785 = (-mrSges(4,1) * t824 - mrSges(4,3) * t820) * t861;
t786 = (-mrSges(3,1) * t824 + mrSges(3,2) * t820) * t861;
t721 = t834 - t872;
t708 = -t872 + (-t820 * t851 - t808) * pkin(3) + t828;
t787 = (mrSges(5,1) * t820 - mrSges(5,2) * t824) * t861;
t864 = t823 * t683 - t819 * t688;
t839 = m(5) * t708 + t808 * mrSges(5,2) + t809 * t779 - t787 * t847 + t864;
t832 = -m(4) * t721 + t808 * mrSges(4,1) + t809 * t781 - t839;
t674 = m(3) * t733 + t808 * mrSges(3,1) + t809 * t780 + (-t785 - t786) * t847 + (mrSges(5,3) - t871) * t789 + t832;
t853 = g(3) * t869;
t734 = -t853 + t863;
t715 = t794 - t840 - t853;
t707 = -0.2e1 * t859 + (0.2e1 * t824 * t858 + t873) * t816 + t831;
t684 = t822 * t692 + t818 * t693;
t835 = -m(6) * t704 + t731 * mrSges(6,1) - t732 * mrSges(6,2) + t764 * t739 - t765 * t740 - t684;
t829 = -m(5) * t707 - t790 * mrSges(5,3) - t835;
t827 = m(4) * t715 + t808 * mrSges(4,3) + t809 * t778 + t785 * t846 + t829;
t681 = t827 + (t786 - t787) * t846 + t862 * t809 + (-mrSges(3,2) + mrSges(5,1)) * t808 + t871 * t790 + m(3) * t734;
t663 = -t671 * t816 + t674 * t866 + t681 * t867;
t661 = m(2) * t803 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t826 + t663;
t667 = -t674 * t820 + t824 * t681;
t666 = m(2) * t804 - mrSges(2,1) * t826 - qJDD(1) * mrSges(2,2) + t667;
t865 = t825 * t661 + t821 * t666;
t662 = t817 * t671 + t674 * t868 + t681 * t869;
t850 = (t820 * t854 + t824 * t855) * t861 + t878 * t809;
t849 = (-t820 * t856 + t824 * t877) * t861 - t855 * t809;
t848 = (-t820 * t879 - t856 * t824) * t861 - t854 * t809;
t845 = -t661 * t821 + t825 * t666;
t717 = Ifges(7,5) * t738 + Ifges(7,6) * t737 + Ifges(7,3) * t762;
t719 = Ifges(7,1) * t738 + Ifges(7,4) * t737 + Ifges(7,5) * t762;
t685 = -mrSges(7,1) * t696 + mrSges(7,3) * t695 + Ifges(7,4) * t712 + Ifges(7,2) * t711 + Ifges(7,6) * t729 - t717 * t738 + t719 * t762;
t718 = Ifges(7,4) * t738 + Ifges(7,2) * t737 + Ifges(7,6) * t762;
t686 = mrSges(7,2) * t696 - mrSges(7,3) * t694 + Ifges(7,1) * t712 + Ifges(7,4) * t711 + Ifges(7,5) * t729 + t717 * t737 - t718 * t762;
t725 = Ifges(6,5) * t765 + Ifges(6,6) * t764 + Ifges(6,3) * t795;
t726 = Ifges(6,4) * t765 + Ifges(6,2) * t764 + Ifges(6,6) * t795;
t668 = mrSges(6,2) * t704 - mrSges(6,3) * t698 + Ifges(6,1) * t732 + Ifges(6,4) * t731 + Ifges(6,5) * t773 - pkin(10) * t684 - t685 * t818 + t686 * t822 + t725 * t764 - t726 * t795;
t727 = Ifges(6,1) * t765 + Ifges(6,4) * t764 + Ifges(6,5) * t795;
t669 = -mrSges(6,1) * t704 - mrSges(7,1) * t694 + mrSges(7,2) * t695 + mrSges(6,3) * t699 + Ifges(6,4) * t732 - Ifges(7,5) * t712 + Ifges(6,2) * t731 + Ifges(6,6) * t773 - Ifges(7,6) * t711 - Ifges(7,3) * t729 - pkin(5) * t684 - t718 * t738 + t719 * t737 - t725 * t765 + t727 * t795;
t675 = t790 * mrSges(5,2) - t789 * mrSges(4,3) + (-t781 * t824 + (-t776 - t778) * t820) * t861 + t833;
t658 = -t823 * t668 - pkin(3) * t836 - qJ(4) * t829 + t819 * t669 - mrSges(3,1) * t750 + mrSges(3,3) * t734 + mrSges(4,2) * t715 - mrSges(4,1) * t716 + mrSges(5,3) * t707 - mrSges(5,2) * t709 - pkin(2) * t675 + pkin(9) * t677 + (-qJ(4) * t776 - t848) * t809 + (-mrSges(5,1) * qJ(4) + t855) * t808 + (-mrSges(5,2) * pkin(3) - t877) * t790 + t856 * t789 + (qJ(4) * t787 * t824 + (pkin(3) * t776 - t850) * t820) * t861;
t676 = -t789 * mrSges(5,3) + t839;
t659 = pkin(10) * t844 + t849 * t809 + t879 * t789 + pkin(5) * t837 + t822 * t685 + t818 * t686 + Ifges(6,3) * t773 - t764 * t727 + t765 * t726 + mrSges(3,2) * t750 + Ifges(6,6) * t731 + Ifges(6,5) * t732 - mrSges(3,3) * t733 + mrSges(4,2) * t721 - mrSges(4,3) * t716 - mrSges(5,3) * t708 + mrSges(5,1) * t709 - mrSges(6,2) * t699 + mrSges(6,1) * t698 - qJ(3) * t675 - qJ(4) * t676 + pkin(4) * t677 + t850 * t846 + t854 * t808 + t856 * t790;
t838 = pkin(8) * t667 + t658 * t824 + t659 * t820;
t657 = -t819 * t668 + pkin(2) * t832 - t823 * t669 + qJ(3) * (t809 * t776 + t827) - pkin(4) * t835 + mrSges(3,1) * t733 - mrSges(3,2) * t734 - mrSges(4,1) * t721 + mrSges(4,3) * t715 - mrSges(5,1) * t707 + mrSges(5,2) * t708 - pkin(3) * t676 - pkin(9) * t864 + (mrSges(5,1) * qJ(3) + t878) * t808 + (mrSges(4,2) * qJ(3) + t855) * t790 + (pkin(2) * (-mrSges(4,2) + mrSges(5,3)) + t854) * t789 + ((-qJ(3) * t787 + t848) * t824 + (-pkin(2) * t785 - t849) * t820) * t861;
t656 = -mrSges(2,2) * g(3) - mrSges(2,3) * t803 + Ifges(2,5) * qJDD(1) - t826 * Ifges(2,6) - t820 * t658 + t824 * t659 + (-t662 * t816 - t663 * t817) * pkin(8);
t655 = mrSges(2,1) * g(3) + mrSges(2,3) * t804 + t826 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t662 - t816 * t657 + t817 * t838;
t1 = [-m(1) * g(1) + t845; -m(1) * g(2) + t865; (-m(1) - m(2)) * g(3) + t662; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t865 - t821 * t655 + t825 * t656; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t845 + t825 * t655 + t821 * t656; -mrSges(1,1) * g(2) + mrSges(2,1) * t803 + mrSges(1,2) * g(1) - mrSges(2,2) * t804 + Ifges(2,3) * qJDD(1) + pkin(1) * t663 + t817 * t657 + t816 * t838;];
tauB  = t1;
