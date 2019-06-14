% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:05:17
% EndTime: 2019-05-05 18:05:25
% DurationCPUTime: 5.72s
% Computational Cost: add. (64569->314), mult. (134052->370), div. (0->0), fcn. (84103->8), ass. (0->126)
t867 = Ifges(6,1) + Ifges(7,1);
t857 = Ifges(6,4) - Ifges(7,5);
t855 = -Ifges(6,5) - Ifges(7,4);
t866 = Ifges(6,2) + Ifges(7,3);
t854 = Ifges(6,6) - Ifges(7,6);
t865 = -Ifges(6,3) - Ifges(7,2);
t822 = sin(qJ(1));
t824 = cos(qJ(1));
t806 = -t824 * g(1) - t822 * g(2);
t864 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t806;
t818 = sin(pkin(9));
t819 = cos(pkin(9));
t823 = cos(qJ(3));
t846 = qJD(1) * t823;
t792 = qJD(3) * t819 - t818 * t846;
t793 = qJD(3) * t818 + t819 * t846;
t820 = sin(qJ(5));
t861 = cos(qJ(5));
t758 = -t861 * t792 + t820 * t793;
t821 = sin(qJ(3));
t845 = qJD(1) * qJD(3);
t842 = t821 * t845;
t801 = qJDD(1) * t823 - t842;
t771 = qJDD(3) * t819 - t801 * t818;
t772 = qJDD(3) * t818 + t801 * t819;
t728 = -t758 * qJD(5) + t820 * t771 + t861 * t772;
t759 = t820 * t792 + t861 * t793;
t739 = mrSges(7,1) * t758 - mrSges(7,3) * t759;
t826 = qJD(1) ^ 2;
t862 = (-pkin(1) - pkin(7));
t768 = (t862 * t826) - t864;
t841 = t823 * t845;
t800 = qJDD(1) * t821 + t841;
t743 = (-t801 + t842) * qJ(4) + (t800 + t841) * pkin(3) + t768;
t805 = t822 * g(1) - t824 * g(2);
t834 = -t826 * qJ(2) + qJDD(2) - t805;
t774 = t862 * qJDD(1) + t834;
t763 = -g(3) * t823 + t821 * t774;
t798 = (pkin(3) * t821 - qJ(4) * t823) * qJD(1);
t825 = qJD(3) ^ 2;
t847 = qJD(1) * t821;
t748 = -pkin(3) * t825 + qJDD(3) * qJ(4) - t798 * t847 + t763;
t721 = -0.2e1 * qJD(4) * t793 + t819 * t743 - t818 * t748;
t718 = (t792 * t847 - t772) * pkin(8) + (t792 * t793 + t800) * pkin(4) + t721;
t722 = 0.2e1 * qJD(4) * t792 + t818 * t743 + t819 * t748;
t773 = pkin(4) * t847 - pkin(8) * t793;
t791 = t792 ^ 2;
t720 = -pkin(4) * t791 + pkin(8) * t771 - t773 * t847 + t722;
t713 = t861 * t718 - t820 * t720;
t738 = pkin(5) * t758 - qJ(6) * t759;
t797 = qJDD(5) + t800;
t808 = qJD(5) + t847;
t807 = t808 ^ 2;
t712 = -t797 * pkin(5) - t807 * qJ(6) + t759 * t738 + qJDD(6) - t713;
t750 = -mrSges(7,2) * t758 + mrSges(7,3) * t808;
t836 = -m(7) * t712 + t797 * mrSges(7,1) + t808 * t750;
t708 = t728 * mrSges(7,2) + t759 * t739 - t836;
t714 = t820 * t718 + t861 * t720;
t711 = -pkin(5) * t807 + qJ(6) * t797 + 0.2e1 * qJD(6) * t808 - t738 * t758 + t714;
t727 = t759 * qJD(5) - t861 * t771 + t820 * t772;
t752 = -mrSges(7,1) * t808 + mrSges(7,2) * t759;
t843 = m(7) * t711 + t797 * mrSges(7,3) + t808 * t752;
t849 = -t857 * t758 + t867 * t759 - t855 * t808;
t851 = t866 * t758 - t857 * t759 - t854 * t808;
t863 = t854 * t727 + t855 * t728 - t849 * t758 + t851 * t759 + t865 * t797 - mrSges(6,1) * t713 + mrSges(7,1) * t712 + mrSges(6,2) * t714 - mrSges(7,3) * t711 + pkin(5) * t708 - qJ(6) * (-t727 * mrSges(7,2) - t758 * t739 + t843);
t860 = mrSges(2,1) - mrSges(3,2);
t859 = -mrSges(6,3) - mrSges(7,2);
t858 = -Ifges(3,4) + Ifges(2,5);
t856 = (Ifges(3,5) - Ifges(2,6));
t799 = (mrSges(4,1) * t821 + mrSges(4,2) * t823) * qJD(1);
t804 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t846;
t751 = mrSges(6,1) * t808 - mrSges(6,3) * t759;
t848 = -mrSges(6,1) * t758 - mrSges(6,2) * t759 - t739;
t701 = m(6) * t714 - t797 * mrSges(6,2) + t859 * t727 - t808 * t751 + t848 * t758 + t843;
t749 = -mrSges(6,2) * t808 - mrSges(6,3) * t758;
t703 = m(6) * t713 + t797 * mrSges(6,1) + t859 * t728 + t808 * t749 + t848 * t759 + t836;
t698 = t820 * t701 + t861 * t703;
t760 = -mrSges(5,1) * t792 + mrSges(5,2) * t793;
t769 = -mrSges(5,2) * t847 + mrSges(5,3) * t792;
t694 = m(5) * t721 + mrSges(5,1) * t800 - mrSges(5,3) * t772 - t760 * t793 + t769 * t847 + t698;
t770 = mrSges(5,1) * t847 - mrSges(5,3) * t793;
t837 = t861 * t701 - t703 * t820;
t695 = m(5) * t722 - mrSges(5,2) * t800 + mrSges(5,3) * t771 + t760 * t792 - t770 * t847 + t837;
t838 = -t694 * t818 + t819 * t695;
t688 = m(4) * t763 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t800 - qJD(3) * t804 - t799 * t847 + t838;
t762 = t821 * g(3) + t823 * t774;
t745 = -qJDD(3) * pkin(3) - t825 * qJ(4) + t798 * t846 + qJDD(4) - t762;
t723 = -t771 * pkin(4) - t791 * pkin(8) + t793 * t773 + t745;
t716 = -0.2e1 * qJD(6) * t759 + (t758 * t808 - t728) * qJ(6) + (t759 * t808 + t727) * pkin(5) + t723;
t709 = m(7) * t716 + t727 * mrSges(7,1) - t728 * mrSges(7,3) + t758 * t750 - t759 * t752;
t829 = m(6) * t723 + t727 * mrSges(6,1) + t728 * mrSges(6,2) + t758 * t749 + t759 * t751 + t709;
t706 = m(5) * t745 - t771 * mrSges(5,1) + t772 * mrSges(5,2) - t792 * t769 + t793 * t770 + t829;
t803 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t847;
t704 = m(4) * t762 + qJDD(3) * mrSges(4,1) - t801 * mrSges(4,3) + qJD(3) * t803 - t799 * t846 - t706;
t682 = t821 * t688 + t823 * t704;
t779 = -qJDD(1) * pkin(1) + t834;
t833 = -m(3) * t779 + (t826 * mrSges(3,3)) - t682;
t678 = m(2) * t805 - (t826 * mrSges(2,2)) + t860 * qJDD(1) + t833;
t777 = t826 * pkin(1) + t864;
t690 = t819 * t694 + t818 * t695;
t832 = -m(4) * t768 - mrSges(4,1) * t800 - t801 * mrSges(4,2) - t803 * t847 - t804 * t846 - t690;
t830 = -m(3) * t777 + (t826 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t832;
t685 = m(2) * t806 - (mrSges(2,1) * t826) - qJDD(1) * mrSges(2,2) + t830;
t852 = t824 * t678 + t822 * t685;
t850 = t854 * t758 + t855 * t759 + t865 * t808;
t840 = -t678 * t822 + t824 * t685;
t839 = t823 * t688 - t821 * t704;
t696 = -mrSges(6,1) * t723 - mrSges(7,1) * t716 + mrSges(7,2) * t711 + mrSges(6,3) * t714 - pkin(5) * t709 - t866 * t727 + t857 * t728 + t850 * t759 + t854 * t797 + t849 * t808;
t697 = mrSges(6,2) * t723 + mrSges(7,2) * t712 - mrSges(6,3) * t713 - mrSges(7,3) * t716 - qJ(6) * t709 - t857 * t727 + t867 * t728 + t850 * t758 - t855 * t797 + t851 * t808;
t753 = Ifges(5,5) * t793 + Ifges(5,6) * t792 + Ifges(5,3) * t847;
t755 = Ifges(5,1) * t793 + Ifges(5,4) * t792 + Ifges(5,5) * t847;
t674 = -mrSges(5,1) * t745 + mrSges(5,3) * t722 + Ifges(5,4) * t772 + Ifges(5,2) * t771 + Ifges(5,6) * t800 - pkin(4) * t829 + pkin(8) * t837 + t861 * t696 + t820 * t697 - t793 * t753 + t755 * t847;
t754 = Ifges(5,4) * t793 + Ifges(5,2) * t792 + Ifges(5,6) * t847;
t676 = mrSges(5,2) * t745 - mrSges(5,3) * t721 + Ifges(5,1) * t772 + Ifges(5,4) * t771 + Ifges(5,5) * t800 - pkin(8) * t698 - t820 * t696 + t861 * t697 + t792 * t753 - t754 * t847;
t785 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t823 - Ifges(4,2) * t821) * qJD(1);
t786 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t823 - Ifges(4,4) * t821) * qJD(1);
t831 = mrSges(4,1) * t762 - mrSges(4,2) * t763 + Ifges(4,5) * t801 - Ifges(4,6) * t800 + Ifges(4,3) * qJDD(3) - pkin(3) * t706 + qJ(4) * t838 + t819 * t674 + t818 * t676 + t785 * t846 + t786 * t847;
t784 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t823 - Ifges(4,6) * t821) * qJD(1);
t671 = mrSges(4,2) * t768 - mrSges(4,3) * t762 + Ifges(4,1) * t801 - Ifges(4,4) * t800 + Ifges(4,5) * qJDD(3) - qJ(4) * t690 - qJD(3) * t785 - t674 * t818 + t676 * t819 - t784 * t847;
t672 = t863 - pkin(4) * t698 - t784 * t846 - pkin(3) * t690 + Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(4,2)) * t800 - mrSges(5,1) * t721 + mrSges(5,2) * t722 + mrSges(4,3) * t763 - mrSges(4,1) * t768 - Ifges(5,6) * t771 - Ifges(5,5) * t772 + qJD(3) * t786 + t792 * t755 - t793 * t754 + Ifges(4,4) * t801;
t680 = qJDD(1) * mrSges(3,2) - t833;
t827 = mrSges(2,1) * t805 - mrSges(2,2) * t806 + mrSges(3,2) * t779 - mrSges(3,3) * t777 - pkin(1) * t680 - pkin(7) * t682 + qJ(2) * t830 + t823 * t671 - t672 * t821 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t681 = -m(3) * g(3) + t839;
t669 = (t856 * t826) + t831 + pkin(2) * t682 - qJ(2) * t681 + t858 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + mrSges(3,1) * t779 - mrSges(2,3) * t805;
t668 = -mrSges(3,1) * t777 + mrSges(2,3) * t806 - pkin(1) * t681 - pkin(2) * t832 - pkin(7) * t839 + t860 * g(3) - t856 * qJDD(1) - t821 * t671 - t823 * t672 + t858 * t826;
t1 = [-m(1) * g(1) + t840; -m(1) * g(2) + t852; (-m(1) - m(2) - m(3)) * g(3) + t839; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t852 - t822 * t668 + t824 * t669; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t840 + t824 * t668 + t822 * t669; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t827; t827; t680; t831; t706; -t863; t708;];
tauJB  = t1;
