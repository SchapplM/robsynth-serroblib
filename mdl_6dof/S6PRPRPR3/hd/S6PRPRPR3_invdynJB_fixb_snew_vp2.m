% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 22:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:33:27
% EndTime: 2019-05-04 22:33:34
% DurationCPUTime: 6.05s
% Computational Cost: add. (74174->284), mult. (138345->347), div. (0->0), fcn. (87090->12), ass. (0->134)
t857 = Ifges(5,1) + Ifges(6,2);
t847 = Ifges(5,4) + Ifges(6,6);
t846 = Ifges(5,5) - Ifges(6,4);
t856 = Ifges(5,2) + Ifges(6,3);
t845 = Ifges(5,6) - Ifges(6,5);
t855 = Ifges(5,3) + Ifges(6,1);
t806 = sin(qJ(4));
t809 = cos(qJ(4));
t776 = (mrSges(6,2) * t809 - mrSges(6,3) * t806) * qJD(2);
t832 = qJD(2) * t809;
t785 = -mrSges(6,1) * t832 - qJD(4) * mrSges(6,3);
t800 = sin(pkin(10));
t803 = cos(pkin(10));
t781 = g(1) * t800 - g(2) * t803;
t798 = -g(3) + qJDD(1);
t801 = sin(pkin(6));
t804 = cos(pkin(6));
t754 = -t781 * t801 + t804 * t798;
t753 = qJDD(3) + t754;
t831 = qJD(2) * qJD(4);
t830 = t809 * t831;
t778 = qJDD(2) * t806 + t830;
t782 = -g(1) * t803 - g(2) * t800;
t807 = sin(qJ(2));
t810 = cos(qJ(2));
t840 = t804 * t810;
t842 = t801 * t810;
t737 = t781 * t840 - t782 * t807 + t798 * t842;
t735 = qJDD(2) * pkin(2) + t737;
t841 = t804 * t807;
t843 = t801 * t807;
t738 = t781 * t841 + t810 * t782 + t798 * t843;
t812 = qJD(2) ^ 2;
t736 = -pkin(2) * t812 + t738;
t799 = sin(pkin(11));
t802 = cos(pkin(11));
t730 = t799 * t735 + t802 * t736;
t728 = -pkin(3) * t812 + qJDD(2) * pkin(8) + t730;
t725 = t806 * t728;
t775 = (-pkin(4) * t809 - qJ(5) * t806) * qJD(2);
t811 = qJD(4) ^ 2;
t833 = qJD(2) * t806;
t824 = -t811 * qJ(5) + t775 * t833 + qJDD(5) + t725;
t849 = pkin(9) * t812;
t850 = -pkin(4) - pkin(9);
t718 = t778 * pkin(5) + t850 * qJDD(4) + (-pkin(5) * t831 - t806 * t849 - t753) * t809 + t824;
t829 = t806 * t831;
t779 = qJDD(2) * t809 - t829;
t789 = pkin(5) * t833 - qJD(4) * pkin(9);
t797 = t809 ^ 2;
t729 = t802 * t735 - t799 * t736;
t823 = -qJDD(2) * pkin(3) - t729;
t851 = -2 * qJD(5);
t815 = pkin(4) * t829 + t833 * t851 + (-t778 - t830) * qJ(5) + t823;
t719 = -t789 * t833 + (-pkin(5) * t797 - pkin(8)) * t812 + t850 * t779 + t815;
t805 = sin(qJ(6));
t808 = cos(qJ(6));
t714 = t718 * t808 - t719 * t805;
t773 = -qJD(4) * t805 - t808 * t832;
t747 = qJD(6) * t773 + qJDD(4) * t808 - t779 * t805;
t774 = qJD(4) * t808 - t805 * t832;
t748 = -mrSges(7,1) * t773 + mrSges(7,2) * t774;
t791 = qJD(6) + t833;
t751 = -mrSges(7,2) * t791 + mrSges(7,3) * t773;
t771 = qJDD(6) + t778;
t711 = m(7) * t714 + mrSges(7,1) * t771 - t747 * mrSges(7,3) - t748 * t774 + t751 * t791;
t715 = t718 * t805 + t719 * t808;
t746 = -qJD(6) * t774 - qJDD(4) * t805 - t779 * t808;
t752 = mrSges(7,1) * t791 - mrSges(7,3) * t774;
t712 = m(7) * t715 - mrSges(7,2) * t771 + t746 * mrSges(7,3) + t748 * t773 - t752 * t791;
t702 = t808 * t711 + t805 * t712;
t839 = t809 * t753;
t721 = -qJDD(4) * pkin(4) + t824 - t839;
t819 = -m(6) * t721 - t778 * mrSges(6,1) - t702;
t701 = qJDD(4) * mrSges(6,2) + qJD(4) * t785 + t776 * t833 - t819;
t724 = t809 * t728 + t806 * t753;
t817 = -t811 * pkin(4) + qJDD(4) * qJ(5) + t775 * t832 + t724;
t717 = -t797 * t849 + t779 * pkin(5) + ((2 * qJD(5)) + t789) * qJD(4) + t817;
t739 = Ifges(7,5) * t774 + Ifges(7,6) * t773 + Ifges(7,3) * t791;
t741 = Ifges(7,1) * t774 + Ifges(7,4) * t773 + Ifges(7,5) * t791;
t703 = -mrSges(7,1) * t717 + mrSges(7,3) * t715 + Ifges(7,4) * t747 + Ifges(7,2) * t746 + Ifges(7,6) * t771 - t739 * t774 + t741 * t791;
t740 = Ifges(7,4) * t774 + Ifges(7,2) * t773 + Ifges(7,6) * t791;
t704 = mrSges(7,2) * t717 - mrSges(7,3) * t714 + Ifges(7,1) * t747 + Ifges(7,4) * t746 + Ifges(7,5) * t771 + t739 * t773 - t740 * t791;
t720 = qJD(4) * t851 - t817;
t723 = -t725 + t839;
t786 = mrSges(6,1) * t833 + qJD(4) * mrSges(6,2);
t820 = -m(7) * t717 + t746 * mrSges(7,1) - t747 * mrSges(7,2) + t773 * t751 - t774 * t752;
t816 = -m(6) * t720 + qJDD(4) * mrSges(6,3) + qJD(4) * t786 + t776 * t832 - t820;
t834 = t846 * qJD(4) + (t806 * t857 + t847 * t809) * qJD(2);
t835 = t845 * qJD(4) + (t806 * t847 + t809 * t856) * qJD(2);
t854 = (t835 * t806 - t834 * t809) * qJD(2) + t855 * qJDD(4) + t846 * t778 + t845 * t779 + mrSges(5,1) * t723 - mrSges(5,2) * t724 + mrSges(6,2) * t721 - mrSges(6,3) * t720 - pkin(4) * t701 - pkin(9) * t702 + qJ(5) * (mrSges(6,1) * t779 + t816) - t805 * t703 + t808 * t704;
t848 = t812 * pkin(8);
t777 = (-mrSges(5,1) * t809 + mrSges(5,2) * t806) * qJD(2);
t784 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t832;
t699 = m(5) * t723 - t778 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t784 - t785) * qJD(4) + (-t776 - t777) * t833 + t819;
t783 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t833;
t707 = t777 * t832 + m(5) * t724 - qJDD(4) * mrSges(5,2) - qJD(4) * t783 + (mrSges(5,3) + mrSges(6,1)) * t779 + t816;
t826 = -t699 * t806 + t809 * t707;
t691 = m(4) * t730 - mrSges(4,1) * t812 - qJDD(2) * mrSges(4,2) + t826;
t727 = t823 - t848;
t722 = -t779 * pkin(4) + t815 - t848;
t837 = -t805 * t711 + t808 * t712;
t822 = -m(6) * t722 - t779 * mrSges(6,2) + t786 * t833 - t837;
t814 = -m(5) * t727 + t784 * t832 + t779 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t778 + (-t783 * t806 - t785 * t809) * qJD(2) + t822;
t697 = m(4) * t729 + qJDD(2) * mrSges(4,1) - t812 * mrSges(4,2) + t814;
t688 = t799 * t691 + t802 * t697;
t686 = m(3) * t737 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t812 + t688;
t827 = t802 * t691 - t697 * t799;
t687 = m(3) * t738 - mrSges(3,1) * t812 - qJDD(2) * mrSges(3,2) + t827;
t695 = t809 * t699 + t806 * t707;
t694 = m(4) * t753 + t695;
t693 = m(3) * t754 + t694;
t673 = t686 * t840 + t687 * t841 - t693 * t801;
t671 = m(2) * t781 + t673;
t677 = -t686 * t807 + t810 * t687;
t676 = m(2) * t782 + t677;
t838 = t803 * t671 + t800 * t676;
t836 = t855 * qJD(4) + (t806 * t846 + t809 * t845) * qJD(2);
t672 = t686 * t842 + t687 * t843 + t804 * t693;
t828 = -t671 * t800 + t803 * t676;
t825 = m(2) * t798 + t672;
t700 = -t778 * mrSges(6,3) + t785 * t832 - t822;
t679 = -mrSges(5,1) * t727 - mrSges(6,1) * t720 + mrSges(6,2) * t722 + mrSges(5,3) * t724 - pkin(4) * t700 - pkin(5) * t820 - pkin(9) * t837 + t834 * qJD(4) + t845 * qJDD(4) - t808 * t703 - t805 * t704 + t847 * t778 + t779 * t856 - t836 * t833;
t818 = mrSges(7,1) * t714 - mrSges(7,2) * t715 + Ifges(7,5) * t747 + Ifges(7,6) * t746 + Ifges(7,3) * t771 + t774 * t740 - t773 * t741;
t680 = mrSges(6,1) * t721 + mrSges(5,2) * t727 - mrSges(5,3) * t723 - mrSges(6,3) * t722 + pkin(5) * t702 - qJ(5) * t700 - t835 * qJD(4) + t846 * qJDD(4) + t778 * t857 + t847 * t779 + t836 * t832 + t818;
t669 = mrSges(4,2) * t753 - mrSges(4,3) * t729 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t812 - pkin(8) * t695 - t679 * t806 + t680 * t809;
t678 = -mrSges(4,1) * t753 + mrSges(4,3) * t730 + t812 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t695 - t854;
t666 = -mrSges(3,1) * t754 + mrSges(3,3) * t738 + t812 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t694 + qJ(3) * t827 + t799 * t669 + t802 * t678;
t667 = mrSges(3,2) * t754 - mrSges(3,3) * t737 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t812 - qJ(3) * t688 + t669 * t802 - t678 * t799;
t821 = pkin(7) * t677 + t666 * t810 + t667 * t807;
t668 = mrSges(3,1) * t737 - mrSges(3,2) * t738 + mrSges(4,1) * t729 - mrSges(4,2) * t730 + t806 * t680 + t809 * t679 + pkin(3) * t814 + pkin(8) * t826 + pkin(2) * t688 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t665 = mrSges(2,2) * t798 - mrSges(2,3) * t781 - t807 * t666 + t810 * t667 + (-t672 * t801 - t673 * t804) * pkin(7);
t664 = -mrSges(2,1) * t798 + mrSges(2,3) * t782 - pkin(1) * t672 - t801 * t668 + t821 * t804;
t1 = [-m(1) * g(1) + t828; -m(1) * g(2) + t838; -m(1) * g(3) + t825; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t838 - t800 * t664 + t803 * t665; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t828 + t803 * t664 + t800 * t665; -mrSges(1,1) * g(2) + mrSges(2,1) * t781 + mrSges(1,2) * g(1) - mrSges(2,2) * t782 + pkin(1) * t673 + t804 * t668 + t821 * t801; t825; t668; t694; t854; t701; t818;];
tauJB  = t1;
