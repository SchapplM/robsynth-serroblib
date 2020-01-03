% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:15:11
% EndTime: 2019-12-31 22:15:21
% DurationCPUTime: 7.69s
% Computational Cost: add. (116577->303), mult. (249130->386), div. (0->0), fcn. (189232->10), ass. (0->130)
t812 = Ifges(5,1) + Ifges(6,1);
t804 = Ifges(5,4) - Ifges(6,5);
t803 = -Ifges(5,5) - Ifges(6,4);
t811 = Ifges(5,2) + Ifges(6,3);
t802 = Ifges(5,6) - Ifges(6,6);
t810 = -Ifges(5,3) - Ifges(6,2);
t767 = sin(pkin(5));
t771 = sin(qJ(2));
t774 = cos(qJ(2));
t788 = qJD(1) * qJD(2);
t754 = (-qJDD(1) * t774 + t771 * t788) * t767;
t768 = cos(pkin(5));
t763 = t768 * qJD(1) + qJD(2);
t770 = sin(qJ(3));
t773 = cos(qJ(3));
t790 = qJD(1) * t767;
t786 = t771 * t790;
t742 = t773 * t763 - t770 * t786;
t753 = (qJDD(1) * t771 + t774 * t788) * t767;
t762 = t768 * qJDD(1) + qJDD(2);
t724 = t742 * qJD(3) + t773 * t753 + t770 * t762;
t743 = t770 * t763 + t773 * t786;
t789 = qJD(1) * t774;
t785 = t767 * t789;
t758 = qJD(3) - t785;
t769 = sin(qJ(4));
t808 = cos(qJ(4));
t729 = t769 * t743 - t808 * t758;
t746 = qJDD(3) + t754;
t690 = -t729 * qJD(4) + t808 * t724 + t769 * t746;
t730 = t808 * t743 + t769 * t758;
t708 = t729 * mrSges(6,1) - t730 * mrSges(6,3);
t752 = (-pkin(2) * t774 - pkin(8) * t771) * t790;
t761 = t763 ^ 2;
t772 = sin(qJ(1));
t775 = cos(qJ(1));
t759 = t772 * g(1) - t775 * g(2);
t776 = qJD(1) ^ 2;
t807 = pkin(7) * t767;
t749 = qJDD(1) * pkin(1) + t776 * t807 + t759;
t760 = -t775 * g(1) - t772 * g(2);
t750 = -t776 * pkin(1) + qJDD(1) * t807 + t760;
t798 = t768 * t771;
t791 = t749 * t798 + t774 * t750;
t703 = -t761 * pkin(2) + t762 * pkin(8) + (-g(3) * t771 + t752 * t789) * t767 + t791;
t806 = t768 * g(3);
t704 = t754 * pkin(2) - t753 * pkin(8) - t806 + (-t749 + (pkin(2) * t771 - pkin(8) * t774) * t763 * qJD(1)) * t767;
t686 = t773 * t703 + t770 * t704;
t728 = -t742 * pkin(3) - t743 * pkin(9);
t757 = t758 ^ 2;
t682 = -t757 * pkin(3) + t746 * pkin(9) + t742 * t728 + t686;
t797 = t768 * t774;
t799 = t767 * t774;
t725 = -g(3) * t799 + t749 * t797 - t771 * t750;
t702 = -t762 * pkin(2) - t761 * pkin(8) + t752 * t786 - t725;
t723 = -t743 * qJD(3) - t770 * t753 + t773 * t762;
t684 = (-t742 * t758 - t724) * pkin(9) + (t743 * t758 - t723) * pkin(3) + t702;
t678 = -t769 * t682 + t808 * t684;
t707 = t729 * pkin(4) - t730 * qJ(5);
t721 = qJDD(4) - t723;
t740 = qJD(4) - t742;
t739 = t740 ^ 2;
t676 = -t721 * pkin(4) - t739 * qJ(5) + t730 * t707 + qJDD(5) - t678;
t711 = -t729 * mrSges(6,2) + t740 * mrSges(6,3);
t782 = -m(6) * t676 + t721 * mrSges(6,1) + t740 * t711;
t672 = t690 * mrSges(6,2) + t730 * t708 - t782;
t679 = t808 * t682 + t769 * t684;
t675 = -t739 * pkin(4) + t721 * qJ(5) + 0.2e1 * qJD(5) * t740 - t729 * t707 + t679;
t689 = t730 * qJD(4) + t769 * t724 - t808 * t746;
t714 = -t740 * mrSges(6,1) + t730 * mrSges(6,2);
t787 = m(6) * t675 + t721 * mrSges(6,3) + t740 * t714;
t793 = t804 * t729 - t812 * t730 + t803 * t740;
t794 = t811 * t729 - t804 * t730 - t802 * t740;
t809 = -t802 * t689 - t803 * t690 - t810 * t721 - t793 * t729 - t794 * t730 + mrSges(5,1) * t678 - mrSges(6,1) * t676 - mrSges(5,2) * t679 + mrSges(6,3) * t675 - pkin(4) * t672 + qJ(5) * (-t689 * mrSges(6,2) - t729 * t708 + t787);
t805 = -mrSges(5,3) - mrSges(6,2);
t800 = t767 * t771;
t726 = -g(3) * t800 + t791;
t747 = t763 * mrSges(3,1) - mrSges(3,3) * t786;
t751 = (-mrSges(3,1) * t774 + mrSges(3,2) * t771) * t790;
t713 = t740 * mrSges(5,1) - t730 * mrSges(5,3);
t792 = -t729 * mrSges(5,1) - t730 * mrSges(5,2) - t708;
t668 = m(5) * t679 - t721 * mrSges(5,2) + t805 * t689 - t740 * t713 + t792 * t729 + t787;
t712 = -t740 * mrSges(5,2) - t729 * mrSges(5,3);
t669 = m(5) * t678 + t721 * mrSges(5,1) + t805 * t690 + t740 * t712 + t792 * t730 + t782;
t664 = t808 * t668 - t769 * t669;
t727 = -t742 * mrSges(4,1) + t743 * mrSges(4,2);
t732 = t758 * mrSges(4,1) - t743 * mrSges(4,3);
t661 = m(4) * t686 - t746 * mrSges(4,2) + t723 * mrSges(4,3) + t742 * t727 - t758 * t732 + t664;
t685 = -t770 * t703 + t773 * t704;
t681 = -t746 * pkin(3) - t757 * pkin(9) + t743 * t728 - t685;
t677 = -0.2e1 * qJD(5) * t730 + (t729 * t740 - t690) * qJ(5) + (t730 * t740 + t689) * pkin(4) + t681;
t673 = m(6) * t677 + t689 * mrSges(6,1) - t690 * mrSges(6,3) + t729 * t711 - t730 * t714;
t670 = -m(5) * t681 - t689 * mrSges(5,1) - t690 * mrSges(5,2) - t729 * t712 - t730 * t713 - t673;
t731 = -t758 * mrSges(4,2) + t742 * mrSges(4,3);
t666 = m(4) * t685 + t746 * mrSges(4,1) - t724 * mrSges(4,3) - t743 * t727 + t758 * t731 + t670;
t783 = t773 * t661 - t770 * t666;
t651 = m(3) * t726 - t762 * mrSges(3,2) - t754 * mrSges(3,3) - t763 * t747 + t751 * t785 + t783;
t654 = t770 * t661 + t773 * t666;
t736 = -t767 * t749 - t806;
t748 = -t763 * mrSges(3,2) + mrSges(3,3) * t785;
t653 = m(3) * t736 + t754 * mrSges(3,1) + t753 * mrSges(3,2) + (t747 * t771 - t748 * t774) * t790 + t654;
t663 = t769 * t668 + t808 * t669;
t778 = -m(4) * t702 + t723 * mrSges(4,1) - t724 * mrSges(4,2) + t742 * t731 - t743 * t732 - t663;
t657 = m(3) * t725 + t762 * mrSges(3,1) - t753 * mrSges(3,3) + t763 * t748 - t751 * t786 + t778;
t640 = t651 * t798 - t767 * t653 + t657 * t797;
t637 = m(2) * t759 + qJDD(1) * mrSges(2,1) - t776 * mrSges(2,2) + t640;
t646 = t774 * t651 - t771 * t657;
t644 = m(2) * t760 - t776 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t646;
t796 = t775 * t637 + t772 * t644;
t795 = t802 * t729 + t803 * t730 + t810 * t740;
t639 = t651 * t800 + t768 * t653 + t657 * t799;
t784 = -t772 * t637 + t775 * t644;
t658 = -mrSges(5,1) * t681 - mrSges(6,1) * t677 + mrSges(6,2) * t675 + mrSges(5,3) * t679 - pkin(4) * t673 - t811 * t689 + t804 * t690 + t802 * t721 + t795 * t730 - t793 * t740;
t662 = mrSges(5,2) * t681 + mrSges(6,2) * t676 - mrSges(5,3) * t678 - mrSges(6,3) * t677 - qJ(5) * t673 - t804 * t689 + t812 * t690 - t803 * t721 + t795 * t729 + t794 * t740;
t717 = Ifges(4,5) * t743 + Ifges(4,6) * t742 + Ifges(4,3) * t758;
t718 = Ifges(4,4) * t743 + Ifges(4,2) * t742 + Ifges(4,6) * t758;
t641 = mrSges(4,2) * t702 - mrSges(4,3) * t685 + Ifges(4,1) * t724 + Ifges(4,4) * t723 + Ifges(4,5) * t746 - pkin(9) * t663 - t769 * t658 + t808 * t662 + t742 * t717 - t758 * t718;
t719 = Ifges(4,1) * t743 + Ifges(4,4) * t742 + Ifges(4,5) * t758;
t647 = -mrSges(4,1) * t702 + mrSges(4,3) * t686 + Ifges(4,4) * t724 + Ifges(4,2) * t723 + Ifges(4,6) * t746 - pkin(3) * t663 - t743 * t717 + t758 * t719 - t809;
t734 = Ifges(3,6) * t763 + (Ifges(3,4) * t771 + Ifges(3,2) * t774) * t790;
t735 = Ifges(3,5) * t763 + (Ifges(3,1) * t771 + Ifges(3,4) * t774) * t790;
t631 = Ifges(3,5) * t753 - Ifges(3,6) * t754 + Ifges(3,3) * t762 + mrSges(3,1) * t725 - mrSges(3,2) * t726 + t770 * t641 + t773 * t647 + pkin(2) * t778 + pkin(8) * t783 + (t734 * t771 - t735 * t774) * t790;
t733 = Ifges(3,3) * t763 + (Ifges(3,5) * t771 + Ifges(3,6) * t774) * t790;
t633 = mrSges(3,2) * t736 - mrSges(3,3) * t725 + Ifges(3,1) * t753 - Ifges(3,4) * t754 + Ifges(3,5) * t762 - pkin(8) * t654 + t773 * t641 - t770 * t647 + t733 * t785 - t763 * t734;
t777 = mrSges(4,1) * t685 - mrSges(4,2) * t686 + Ifges(4,5) * t724 + Ifges(4,6) * t723 + Ifges(4,3) * t746 + pkin(3) * t670 + pkin(9) * t664 + t808 * t658 + t769 * t662 + t743 * t718 - t742 * t719;
t635 = -mrSges(3,1) * t736 + mrSges(3,3) * t726 + Ifges(3,4) * t753 - Ifges(3,2) * t754 + Ifges(3,6) * t762 - pkin(2) * t654 - t733 * t786 + t763 * t735 - t777;
t780 = mrSges(2,1) * t759 - mrSges(2,2) * t760 + Ifges(2,3) * qJDD(1) + pkin(1) * t640 + t768 * t631 + t633 * t800 + t635 * t799 + t646 * t807;
t629 = -mrSges(2,2) * g(3) - mrSges(2,3) * t759 + Ifges(2,5) * qJDD(1) - t776 * Ifges(2,6) + t774 * t633 - t771 * t635 + (-t639 * t767 - t640 * t768) * pkin(7);
t628 = mrSges(2,1) * g(3) + mrSges(2,3) * t760 + t776 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t639 - t767 * t631 + (pkin(7) * t646 + t633 * t771 + t635 * t774) * t768;
t1 = [-m(1) * g(1) + t784; -m(1) * g(2) + t796; (-m(1) - m(2)) * g(3) + t639; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t796 - t772 * t628 + t775 * t629; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t784 + t775 * t628 + t772 * t629; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t780; t780; t631; t777; t809; t672;];
tauJB = t1;
