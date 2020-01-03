% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:37
% EndTime: 2019-12-31 20:12:41
% DurationCPUTime: 2.41s
% Computational Cost: add. (16941->271), mult. (34229->314), div. (0->0), fcn. (17731->6), ass. (0->111)
t811 = Ifges(5,1) + Ifges(6,1);
t794 = Ifges(5,4) - Ifges(6,5);
t807 = Ifges(6,4) + Ifges(5,5);
t810 = Ifges(5,2) + Ifges(6,3);
t805 = Ifges(5,6) - Ifges(6,6);
t809 = -2 * qJD(3);
t808 = Ifges(3,1) + Ifges(4,2);
t795 = Ifges(3,4) + Ifges(4,6);
t793 = Ifges(3,5) - Ifges(4,4);
t806 = Ifges(3,2) + Ifges(4,3);
t792 = Ifges(3,6) - Ifges(4,5);
t804 = Ifges(3,3) + Ifges(4,1);
t803 = Ifges(5,3) + Ifges(6,2);
t755 = sin(qJ(4));
t758 = cos(qJ(4));
t759 = cos(qJ(2));
t782 = qJD(1) * t759;
t726 = qJD(2) * t755 + t758 * t782;
t727 = qJD(2) * t758 - t755 * t782;
t756 = sin(qJ(2));
t783 = qJD(1) * t756;
t745 = qJD(4) + t783;
t802 = t810 * t726 - t794 * t727 - t805 * t745;
t801 = -t794 * t726 + t811 * t727 + t807 * t745;
t729 = (mrSges(4,2) * t759 - mrSges(4,3) * t756) * qJD(1);
t738 = -mrSges(4,1) * t782 - qJD(2) * mrSges(4,3);
t781 = qJD(1) * qJD(2);
t777 = t756 * t781;
t732 = qJDD(1) * t759 - t777;
t740 = pkin(3) * t783 - qJD(2) * pkin(7);
t754 = t759 ^ 2;
t762 = qJD(1) ^ 2;
t778 = t759 * t781;
t731 = qJDD(1) * t756 + t778;
t757 = sin(qJ(1));
t760 = cos(qJ(1));
t741 = t757 * g(1) - t760 * g(2);
t773 = -qJDD(1) * pkin(1) - t741;
t768 = pkin(2) * t777 + t783 * t809 + (-t731 - t778) * qJ(3) + t773;
t665 = -t740 * t783 + (-pkin(3) * t754 - pkin(6)) * t762 + (-pkin(2) - pkin(7)) * t732 + t768;
t742 = -g(1) * t760 - g(2) * t757;
t714 = -pkin(1) * t762 + qJDD(1) * pkin(6) + t742;
t700 = -t759 * g(3) - t756 * t714;
t728 = (-pkin(2) * t759 - qJ(3) * t756) * qJD(1);
t761 = qJD(2) ^ 2;
t674 = -qJDD(2) * pkin(2) - t761 * qJ(3) + t728 * t783 + qJDD(3) - t700;
t669 = (-t756 * t759 * t762 - qJDD(2)) * pkin(7) + (t731 - t778) * pkin(3) + t674;
t663 = t758 * t665 + t755 * t669;
t689 = qJD(4) * t727 + qJDD(2) * t755 + t758 * t732;
t698 = mrSges(5,1) * t745 - mrSges(5,3) * t727;
t725 = qJDD(4) + t731;
t693 = pkin(4) * t726 - qJ(5) * t727;
t743 = t745 ^ 2;
t657 = -pkin(4) * t743 + qJ(5) * t725 + 0.2e1 * qJD(5) * t745 - t693 * t726 + t663;
t699 = -mrSges(6,1) * t745 + mrSges(6,2) * t727;
t779 = m(6) * t657 + t725 * mrSges(6,3) + t745 * t699;
t694 = mrSges(6,1) * t726 - mrSges(6,3) * t727;
t787 = -mrSges(5,1) * t726 - mrSges(5,2) * t727 - t694;
t796 = -mrSges(5,3) - mrSges(6,2);
t647 = m(5) * t663 - t725 * mrSges(5,2) + t796 * t689 - t745 * t698 + t787 * t726 + t779;
t662 = -t665 * t755 + t669 * t758;
t690 = -qJD(4) * t726 + qJDD(2) * t758 - t732 * t755;
t696 = -mrSges(5,2) * t745 - mrSges(5,3) * t726;
t658 = -pkin(4) * t725 - qJ(5) * t743 + t693 * t727 + qJDD(5) - t662;
t697 = -mrSges(6,2) * t726 + mrSges(6,3) * t745;
t774 = -m(6) * t658 + t725 * mrSges(6,1) + t745 * t697;
t650 = m(5) * t662 + t725 * mrSges(5,1) + t796 * t690 + t745 * t696 + t787 * t727 + t774;
t642 = t755 * t647 + t758 * t650;
t770 = -m(4) * t674 - t731 * mrSges(4,1) - t642;
t639 = qJDD(2) * mrSges(4,2) + qJD(2) * t738 + t729 * t783 - t770;
t701 = -t756 * g(3) + t759 * t714;
t673 = t761 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t809 - t728 * t782 - t701;
t668 = -t754 * t762 * pkin(7) + t732 * pkin(3) + qJD(2) * t740 - t673;
t660 = -0.2e1 * qJD(5) * t727 + (t726 * t745 - t690) * qJ(5) + (t727 * t745 + t689) * pkin(4) + t668;
t654 = m(6) * t660 + t689 * mrSges(6,1) - t690 * mrSges(6,3) + t726 * t697 - t727 * t699;
t788 = t805 * t726 - t807 * t727 - t803 * t745;
t640 = -mrSges(5,1) * t668 - mrSges(6,1) * t660 + mrSges(6,2) * t657 + mrSges(5,3) * t663 - pkin(4) * t654 - t810 * t689 + t794 * t690 + t805 * t725 + t788 * t727 + t801 * t745;
t641 = mrSges(5,2) * t668 + mrSges(6,2) * t658 - mrSges(5,3) * t662 - mrSges(6,3) * t660 - qJ(5) * t654 - t794 * t689 + t811 * t690 + t807 * t725 + t788 * t726 + t802 * t745;
t739 = mrSges(4,1) * t783 + qJD(2) * mrSges(4,2);
t767 = m(5) * t668 + t689 * mrSges(5,1) + t690 * mrSges(5,2) + t726 * t696 + t727 * t698 + t654;
t766 = -m(4) * t673 + qJDD(2) * mrSges(4,3) + qJD(2) * t739 + t729 * t782 + t767;
t784 = t793 * qJD(2) + (t808 * t756 + t795 * t759) * qJD(1);
t785 = t792 * qJD(2) + (t795 * t756 + t806 * t759) * qJD(1);
t800 = (t785 * t756 - t784 * t759) * qJD(1) + t804 * qJDD(2) + t793 * t731 + t792 * t732 + mrSges(3,1) * t700 - mrSges(3,2) * t701 + mrSges(4,2) * t674 - mrSges(4,3) * t673 - pkin(2) * t639 - pkin(7) * t642 + qJ(3) * (mrSges(4,1) * t732 + t766) - t755 * t640 + t758 * t641;
t797 = t762 * pkin(6);
t730 = (-mrSges(3,1) * t759 + mrSges(3,2) * t756) * qJD(1);
t737 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t782;
t637 = m(3) * t700 - t731 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t737 - t738) * qJD(2) + (-t729 - t730) * t783 + t770;
t736 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t783;
t645 = (mrSges(3,3) + mrSges(4,1)) * t732 + t766 - qJDD(2) * mrSges(3,2) + m(3) * t701 - qJD(2) * t736 + t730 * t782;
t775 = -t637 * t756 + t759 * t645;
t630 = m(2) * t742 - mrSges(2,1) * t762 - qJDD(1) * mrSges(2,2) + t775;
t713 = t773 - t797;
t670 = -t732 * pkin(2) + t768 - t797;
t789 = t758 * t647 - t755 * t650;
t772 = -m(4) * t670 - t732 * mrSges(4,2) + t739 * t783 - t789;
t765 = -m(3) * t713 + t737 * t782 + t732 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t731 + (-t736 * t756 - t738 * t759) * qJD(1) + t772;
t634 = m(2) * t741 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t765;
t790 = t757 * t630 + t760 * t634;
t632 = t759 * t637 + t756 * t645;
t786 = t804 * qJD(2) + (t793 * t756 + t792 * t759) * qJD(1);
t776 = t760 * t630 - t634 * t757;
t638 = -t731 * mrSges(4,3) + t738 * t782 - t772;
t625 = -mrSges(3,1) * t713 - mrSges(4,1) * t673 + mrSges(4,2) * t670 + mrSges(3,3) * t701 - pkin(2) * t638 + pkin(3) * t767 - pkin(7) * t789 + t784 * qJD(2) + t792 * qJDD(2) - t758 * t640 - t755 * t641 + t795 * t731 + t806 * t732 - t786 * t783;
t653 = t690 * mrSges(6,2) + t727 * t694 - t774;
t763 = mrSges(5,1) * t662 - mrSges(6,1) * t658 - mrSges(5,2) * t663 + mrSges(6,3) * t657 - pkin(4) * t653 + qJ(5) * t779 - t802 * t727 + (-qJ(5) * t694 + t801) * t726 + t803 * t725 + t807 * t690 + (-qJ(5) * mrSges(6,2) - t805) * t689;
t627 = mrSges(4,1) * t674 + mrSges(3,2) * t713 - mrSges(3,3) * t700 - mrSges(4,3) * t670 + pkin(3) * t642 - qJ(3) * t638 - t785 * qJD(2) + t793 * qJDD(2) + t808 * t731 + t795 * t732 + t786 * t782 + t763;
t769 = mrSges(2,1) * t741 - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) + pkin(1) * t765 + pkin(6) * t775 + t759 * t625 + t756 * t627;
t623 = mrSges(2,1) * g(3) + mrSges(2,3) * t742 + t762 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t632 - t800;
t622 = -mrSges(2,2) * g(3) - mrSges(2,3) * t741 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t762 - pkin(6) * t632 - t625 * t756 + t627 * t759;
t1 = [-m(1) * g(1) + t776; -m(1) * g(2) + t790; (-m(1) - m(2)) * g(3) + t632; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t790 + t760 * t622 - t757 * t623; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t776 + t757 * t622 + t760 * t623; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t769; t769; t800; t639; t763; t653;];
tauJB = t1;
