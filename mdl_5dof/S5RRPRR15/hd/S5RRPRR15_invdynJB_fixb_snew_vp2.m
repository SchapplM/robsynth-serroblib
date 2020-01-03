% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR15_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:36
% EndTime: 2019-12-31 20:41:42
% DurationCPUTime: 4.04s
% Computational Cost: add. (40856->295), mult. (83913->356), div. (0->0), fcn. (48513->8), ass. (0->121)
t814 = -2 * qJD(3);
t813 = Ifges(3,1) + Ifges(4,2);
t806 = Ifges(3,4) + Ifges(4,6);
t805 = Ifges(3,5) - Ifges(4,4);
t812 = Ifges(3,2) + Ifges(4,3);
t804 = Ifges(3,6) - Ifges(4,5);
t811 = Ifges(3,3) + Ifges(4,1);
t773 = cos(qJ(2));
t769 = sin(qJ(2));
t796 = qJD(1) * qJD(2);
t793 = t769 * t796;
t743 = qJDD(1) * t773 - t793;
t759 = t769 * qJD(1);
t751 = pkin(3) * t759 - qJD(2) * pkin(7);
t766 = t773 ^ 2;
t776 = qJD(1) ^ 2;
t794 = t773 * t796;
t742 = qJDD(1) * t769 + t794;
t770 = sin(qJ(1));
t774 = cos(qJ(1));
t752 = g(1) * t770 - t774 * g(2);
t789 = -qJDD(1) * pkin(1) - t752;
t782 = pkin(2) * t793 + t759 * t814 + (-t742 - t794) * qJ(3) + t789;
t680 = -t751 * t759 + (-pkin(3) * t766 - pkin(6)) * t776 + (-pkin(2) - pkin(7)) * t743 + t782;
t753 = -g(1) * t774 - g(2) * t770;
t727 = -pkin(1) * t776 + qJDD(1) * pkin(6) + t753;
t711 = -t773 * g(3) - t769 * t727;
t739 = (-pkin(2) * t773 - qJ(3) * t769) * qJD(1);
t775 = qJD(2) ^ 2;
t694 = -qJDD(2) * pkin(2) - qJ(3) * t775 + t739 * t759 + qJDD(3) - t711;
t685 = (-t769 * t773 * t776 - qJDD(2)) * pkin(7) + (t742 - t794) * pkin(3) + t694;
t768 = sin(qJ(4));
t772 = cos(qJ(4));
t669 = -t680 * t768 + t772 * t685;
t797 = qJD(1) * t773;
t737 = -qJD(2) * t768 - t772 * t797;
t705 = qJD(4) * t737 + qJDD(2) * t772 - t743 * t768;
t736 = qJDD(4) + t742;
t738 = qJD(2) * t772 - t768 * t797;
t756 = t759 + qJD(4);
t666 = (t737 * t756 - t705) * pkin(8) + (t737 * t738 + t736) * pkin(4) + t669;
t670 = t772 * t680 + t768 * t685;
t704 = -qJD(4) * t738 - qJDD(2) * t768 - t743 * t772;
t713 = pkin(4) * t756 - pkin(8) * t738;
t735 = t737 ^ 2;
t667 = -pkin(4) * t735 + pkin(8) * t704 - t713 * t756 + t670;
t767 = sin(qJ(5));
t771 = cos(qJ(5));
t665 = t666 * t767 + t667 * t771;
t712 = -g(3) * t769 + t773 * t727;
t693 = pkin(2) * t775 - qJDD(2) * qJ(3) + qJD(2) * t814 - t739 * t797 - t712;
t684 = -pkin(7) * t766 * t776 + pkin(3) * t743 + qJD(2) * t751 - t693;
t672 = -pkin(4) * t704 - pkin(8) * t735 + t713 * t738 + t684;
t707 = t737 * t767 + t738 * t771;
t677 = -qJD(5) * t707 + t704 * t771 - t705 * t767;
t706 = t737 * t771 - t738 * t767;
t678 = qJD(5) * t706 + t704 * t767 + t705 * t771;
t754 = qJD(5) + t756;
t686 = Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t754;
t688 = Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t754;
t729 = qJDD(5) + t736;
t652 = -mrSges(6,1) * t672 + mrSges(6,3) * t665 + Ifges(6,4) * t678 + Ifges(6,2) * t677 + Ifges(6,6) * t729 - t686 * t707 + t688 * t754;
t664 = t666 * t771 - t667 * t767;
t687 = Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t754;
t653 = mrSges(6,2) * t672 - mrSges(6,3) * t664 + Ifges(6,1) * t678 + Ifges(6,4) * t677 + Ifges(6,5) * t729 + t686 * t706 - t687 * t754;
t697 = Ifges(5,5) * t738 + Ifges(5,6) * t737 + Ifges(5,3) * t756;
t699 = Ifges(5,1) * t738 + Ifges(5,4) * t737 + Ifges(5,5) * t756;
t695 = -mrSges(6,2) * t754 + mrSges(6,3) * t706;
t696 = mrSges(6,1) * t754 - mrSges(6,3) * t707;
t786 = m(6) * t672 - t677 * mrSges(6,1) + t678 * mrSges(6,2) - t695 * t706 + t707 * t696;
t690 = -mrSges(6,1) * t706 + mrSges(6,2) * t707;
t661 = m(6) * t664 + mrSges(6,1) * t729 - t678 * mrSges(6,3) - t690 * t707 + t695 * t754;
t662 = m(6) * t665 - mrSges(6,2) * t729 + t677 * mrSges(6,3) + t690 * t706 - t696 * t754;
t790 = -t661 * t767 + t771 * t662;
t636 = -mrSges(5,1) * t684 + mrSges(5,3) * t670 + Ifges(5,4) * t705 + Ifges(5,2) * t704 + Ifges(5,6) * t736 - pkin(4) * t786 + pkin(8) * t790 + t771 * t652 + t767 * t653 - t738 * t697 + t756 * t699;
t651 = t771 * t661 + t767 * t662;
t698 = Ifges(5,4) * t738 + Ifges(5,2) * t737 + Ifges(5,6) * t756;
t637 = mrSges(5,2) * t684 - mrSges(5,3) * t669 + Ifges(5,1) * t705 + Ifges(5,4) * t704 + Ifges(5,5) * t736 - pkin(8) * t651 - t652 * t767 + t653 * t771 + t697 * t737 - t698 * t756;
t740 = (mrSges(4,2) * t773 - mrSges(4,3) * t769) * qJD(1);
t749 = -mrSges(4,1) * t797 - qJD(2) * mrSges(4,3);
t708 = -mrSges(5,1) * t737 + mrSges(5,2) * t738;
t709 = -mrSges(5,2) * t756 + mrSges(5,3) * t737;
t648 = m(5) * t669 + mrSges(5,1) * t736 - mrSges(5,3) * t705 - t708 * t738 + t709 * t756 + t651;
t710 = mrSges(5,1) * t756 - mrSges(5,3) * t738;
t649 = m(5) * t670 - mrSges(5,2) * t736 + mrSges(5,3) * t704 + t708 * t737 - t710 * t756 + t790;
t645 = t648 * t772 + t649 * t768;
t785 = -m(4) * t694 - t742 * mrSges(4,1) - t645;
t644 = qJDD(2) * mrSges(4,2) + qJD(2) * t749 + t740 * t759 - t785;
t750 = mrSges(4,1) * t759 + qJD(2) * mrSges(4,2);
t781 = -m(5) * t684 + mrSges(5,1) * t704 - t705 * mrSges(5,2) + t709 * t737 - t738 * t710 - t786;
t779 = -m(4) * t693 + qJDD(2) * mrSges(4,3) + qJD(2) * t750 + t740 * t797 - t781;
t798 = t805 * qJD(2) + (t813 * t769 + t806 * t773) * qJD(1);
t799 = t804 * qJD(2) + (t806 * t769 + t812 * t773) * qJD(1);
t810 = (t799 * t769 - t798 * t773) * qJD(1) + t811 * qJDD(2) + t805 * t742 + t804 * t743 + mrSges(3,1) * t711 - mrSges(3,2) * t712 + mrSges(4,2) * t694 - mrSges(4,3) * t693 - pkin(2) * t644 - pkin(7) * t645 + qJ(3) * (mrSges(4,1) * t743 + t779) - t768 * t636 + t772 * t637;
t807 = pkin(6) * t776;
t741 = (-mrSges(3,1) * t773 + mrSges(3,2) * t769) * qJD(1);
t748 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t797;
t642 = m(3) * t711 - mrSges(3,3) * t742 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t748 - t749) * qJD(2) + (-t740 - t741) * t759 + t785;
t747 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t759;
t656 = t779 + t741 * t797 - qJDD(2) * mrSges(3,2) + m(3) * t712 - qJD(2) * t747 + (mrSges(3,3) + mrSges(4,1)) * t743;
t791 = -t642 * t769 + t773 * t656;
t633 = m(2) * t753 - mrSges(2,1) * t776 - qJDD(1) * mrSges(2,2) + t791;
t726 = t789 - t807;
t691 = -pkin(2) * t743 + t782 - t807;
t801 = -t768 * t648 + t772 * t649;
t788 = -m(4) * t691 - t743 * mrSges(4,2) + t750 * t759 - t801;
t778 = -m(3) * t726 + t748 * t797 + t743 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t742 + (-t747 * t769 - t749 * t773) * qJD(1) + t788;
t639 = m(2) * t752 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t776 + t778;
t802 = t770 * t633 + t774 * t639;
t635 = t773 * t642 + t769 * t656;
t800 = t811 * qJD(2) + (t805 * t769 + t804 * t773) * qJD(1);
t792 = t774 * t633 - t639 * t770;
t643 = -mrSges(4,3) * t742 + t749 * t797 - t788;
t628 = -mrSges(3,1) * t726 - mrSges(4,1) * t693 + mrSges(4,2) * t691 + mrSges(3,3) * t712 - pkin(2) * t643 - pkin(3) * t781 - pkin(7) * t801 + t798 * qJD(2) + t804 * qJDD(2) - t772 * t636 - t768 * t637 + t806 * t742 + t812 * t743 - t800 * t759;
t783 = mrSges(6,1) * t664 - mrSges(6,2) * t665 + Ifges(6,5) * t678 + Ifges(6,6) * t677 + Ifges(6,3) * t729 + t707 * t687 - t706 * t688;
t780 = mrSges(5,1) * t669 - mrSges(5,2) * t670 + Ifges(5,5) * t705 + Ifges(5,6) * t704 + Ifges(5,3) * t736 + pkin(4) * t651 + t738 * t698 - t737 * t699 + t783;
t630 = mrSges(4,1) * t694 + mrSges(3,2) * t726 - mrSges(3,3) * t711 - mrSges(4,3) * t691 + pkin(3) * t645 - qJ(3) * t643 - t799 * qJD(2) + t805 * qJDD(2) + t813 * t742 + t806 * t743 + t800 * t797 + t780;
t784 = mrSges(2,1) * t752 - mrSges(2,2) * t753 + Ifges(2,3) * qJDD(1) + pkin(1) * t778 + pkin(6) * t791 + t773 * t628 + t769 * t630;
t626 = mrSges(2,1) * g(3) + mrSges(2,3) * t753 + t776 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t635 - t810;
t625 = -mrSges(2,2) * g(3) - mrSges(2,3) * t752 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t776 - pkin(6) * t635 - t628 * t769 + t630 * t773;
t1 = [-m(1) * g(1) + t792; -m(1) * g(2) + t802; (-m(1) - m(2)) * g(3) + t635; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t802 + t774 * t625 - t770 * t626; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t792 + t770 * t625 + t774 * t626; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t784; t784; t810; t644; t780; t783;];
tauJB = t1;
