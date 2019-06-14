% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:15
% EndTime: 2019-05-05 15:50:19
% DurationCPUTime: 3.37s
% Computational Cost: add. (38919->290), mult. (75362->344), div. (0->0), fcn. (44252->8), ass. (0->115)
t812 = 2 * qJD(1);
t775 = sin(qJ(1));
t779 = cos(qJ(1));
t747 = -t779 * g(1) - t775 * g(2);
t811 = qJDD(1) * qJ(2) + (qJD(2) * t812) + t747;
t746 = t775 * g(1) - t779 * g(2);
t780 = qJD(1) ^ 2;
t725 = -qJDD(1) * pkin(1) - t780 * qJ(2) + qJDD(2) - t746;
t789 = qJDD(1) * qJ(3) + (qJD(3) * t812) - t725;
t773 = sin(qJ(5));
t774 = sin(qJ(4));
t777 = cos(qJ(5));
t778 = cos(qJ(4));
t730 = (t773 * t778 + t774 * t777) * qJD(1);
t810 = -m(3) - m(4);
t809 = mrSges(2,1) - mrSges(3,2);
t808 = t780 * mrSges(4,3);
t723 = t780 * pkin(1) - t811;
t718 = qJDD(3) + (-pkin(1) - qJ(3)) * t780 + t811;
t712 = -qJDD(1) * pkin(7) + t718;
t704 = t774 * g(3) + t778 * t712;
t804 = qJD(1) * qJD(4);
t799 = t774 * t804;
t741 = t778 * qJDD(1) - t799;
t683 = (-t741 - t799) * pkin(8) + (-t774 * t778 * t780 + qJDD(4)) * pkin(4) + t704;
t705 = -t778 * g(3) + t774 * t712;
t740 = -t774 * qJDD(1) - t778 * t804;
t805 = qJD(1) * t778;
t745 = qJD(4) * pkin(4) - pkin(8) * t805;
t767 = t774 ^ 2;
t684 = -t767 * t780 * pkin(4) + t740 * pkin(8) - qJD(4) * t745 + t705;
t673 = t773 * t683 + t777 * t684;
t731 = (-t773 * t774 + t777 * t778) * qJD(1);
t694 = -t731 * qJD(5) + t777 * t740 - t773 * t741;
t706 = t730 * mrSges(6,1) + t731 * mrSges(6,2);
t756 = qJD(4) + qJD(5);
t720 = t756 * mrSges(6,1) - t731 * mrSges(6,3);
t755 = qJDD(4) + qJDD(5);
t707 = t730 * pkin(5) - t731 * pkin(9);
t754 = t756 ^ 2;
t669 = -t754 * pkin(5) + t755 * pkin(9) - t730 * t707 + t673;
t686 = -t740 * pkin(4) + t745 * t805 + (-pkin(8) * t767 - pkin(7)) * t780 + t789;
t695 = -t730 * qJD(5) + t773 * t740 + t777 * t741;
t670 = (t730 * t756 - t695) * pkin(9) + (t731 * t756 - t694) * pkin(5) + t686;
t772 = sin(qJ(6));
t776 = cos(qJ(6));
t666 = -t772 * t669 + t776 * t670;
t713 = -t772 * t731 + t776 * t756;
t676 = t713 * qJD(6) + t776 * t695 + t772 * t755;
t714 = t776 * t731 + t772 * t756;
t687 = -t713 * mrSges(7,1) + t714 * mrSges(7,2);
t693 = qJDD(6) - t694;
t726 = qJD(6) + t730;
t696 = -t726 * mrSges(7,2) + t713 * mrSges(7,3);
t663 = m(7) * t666 + t693 * mrSges(7,1) - t676 * mrSges(7,3) - t714 * t687 + t726 * t696;
t667 = t776 * t669 + t772 * t670;
t675 = -t714 * qJD(6) - t772 * t695 + t776 * t755;
t697 = t726 * mrSges(7,1) - t714 * mrSges(7,3);
t664 = m(7) * t667 - t693 * mrSges(7,2) + t675 * mrSges(7,3) + t713 * t687 - t726 * t697;
t795 = -t772 * t663 + t776 * t664;
t651 = m(6) * t673 - t755 * mrSges(6,2) + t694 * mrSges(6,3) - t730 * t706 - t756 * t720 + t795;
t672 = t777 * t683 - t773 * t684;
t719 = -t756 * mrSges(6,2) - t730 * mrSges(6,3);
t668 = -t755 * pkin(5) - t754 * pkin(9) + t731 * t707 - t672;
t787 = -m(7) * t668 + t675 * mrSges(7,1) - t676 * mrSges(7,2) + t713 * t696 - t714 * t697;
t659 = m(6) * t672 + t755 * mrSges(6,1) - t695 * mrSges(6,3) - t731 * t706 + t756 * t719 + t787;
t643 = t773 * t651 + t777 * t659;
t739 = (mrSges(5,1) * t774 + mrSges(5,2) * t778) * qJD(1);
t806 = qJD(1) * t774;
t743 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t806;
t640 = m(5) * t704 + qJDD(4) * mrSges(5,1) - t741 * mrSges(5,3) + qJD(4) * t743 - t739 * t805 + t643;
t744 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t805;
t796 = t777 * t651 - t773 * t659;
t641 = m(5) * t705 - qJDD(4) * mrSges(5,2) + t740 * mrSges(5,3) - qJD(4) * t744 - t739 * t806 + t796;
t634 = t778 * t640 + t774 * t641;
t794 = m(4) * t718 + qJDD(1) * mrSges(4,2) + t634;
t788 = -m(3) * t723 + (t780 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t794;
t630 = m(2) * t747 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t780) + t788;
t711 = -t780 * pkin(7) + t789;
t653 = t776 * t663 + t772 * t664;
t790 = m(6) * t686 - t694 * mrSges(6,1) + t695 * mrSges(6,2) + t730 * t719 + t731 * t720 + t653;
t786 = -m(5) * t711 + t740 * mrSges(5,1) - t741 * mrSges(5,2) - t743 * t806 - t744 * t805 - t790;
t648 = -m(4) * t789 - t780 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t786;
t783 = -m(3) * t725 + t780 * mrSges(3,3) - t648;
t645 = m(2) * t746 - t780 * mrSges(2,2) + t809 * qJDD(1) + t783;
t807 = t775 * t630 + t779 * t645;
t801 = (Ifges(2,5) - Ifges(3,4) + Ifges(4,5));
t800 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t798 = t779 * t630 - t775 * t645;
t797 = -t774 * t640 + t778 * t641;
t677 = Ifges(7,5) * t714 + Ifges(7,6) * t713 + Ifges(7,3) * t726;
t679 = Ifges(7,1) * t714 + Ifges(7,4) * t713 + Ifges(7,5) * t726;
t656 = -mrSges(7,1) * t668 + mrSges(7,3) * t667 + Ifges(7,4) * t676 + Ifges(7,2) * t675 + Ifges(7,6) * t693 - t714 * t677 + t726 * t679;
t678 = Ifges(7,4) * t714 + Ifges(7,2) * t713 + Ifges(7,6) * t726;
t657 = mrSges(7,2) * t668 - mrSges(7,3) * t666 + Ifges(7,1) * t676 + Ifges(7,4) * t675 + Ifges(7,5) * t693 + t713 * t677 - t726 * t678;
t699 = Ifges(6,4) * t731 - Ifges(6,2) * t730 + Ifges(6,6) * t756;
t700 = Ifges(6,1) * t731 - Ifges(6,4) * t730 + Ifges(6,5) * t756;
t785 = mrSges(6,1) * t672 - mrSges(6,2) * t673 + Ifges(6,5) * t695 + Ifges(6,6) * t694 + Ifges(6,3) * t755 + pkin(5) * t787 + pkin(9) * t795 + t776 * t656 + t772 * t657 + t731 * t699 + t730 * t700;
t784 = mrSges(7,1) * t666 - mrSges(7,2) * t667 + Ifges(7,5) * t676 + Ifges(7,6) * t675 + Ifges(7,3) * t693 + t714 * t678 - t713 * t679;
t728 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t778 - Ifges(5,2) * t774) * qJD(1);
t729 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t778 - Ifges(5,4) * t774) * qJD(1);
t782 = mrSges(5,1) * t704 - mrSges(5,2) * t705 + Ifges(5,5) * t741 + Ifges(5,6) * t740 + Ifges(5,3) * qJDD(4) + pkin(4) * t643 + t728 * t805 + t729 * t806 + t785;
t698 = Ifges(6,5) * t731 - Ifges(6,6) * t730 + Ifges(6,3) * t756;
t635 = mrSges(6,2) * t686 - mrSges(6,3) * t672 + Ifges(6,1) * t695 + Ifges(6,4) * t694 + Ifges(6,5) * t755 - pkin(9) * t653 - t772 * t656 + t776 * t657 - t730 * t698 - t756 * t699;
t636 = -mrSges(6,1) * t686 + mrSges(6,3) * t673 + Ifges(6,4) * t695 + Ifges(6,2) * t694 + Ifges(6,6) * t755 - pkin(5) * t653 - t731 * t698 + t756 * t700 - t784;
t727 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t778 - Ifges(5,6) * t774) * qJD(1);
t625 = -mrSges(5,1) * t711 + mrSges(5,3) * t705 + Ifges(5,4) * t741 + Ifges(5,2) * t740 + Ifges(5,6) * qJDD(4) - pkin(4) * t790 + pkin(8) * t796 + qJD(4) * t729 + t773 * t635 + t777 * t636 - t727 * t805;
t627 = mrSges(5,2) * t711 - mrSges(5,3) * t704 + Ifges(5,1) * t741 + Ifges(5,4) * t740 + Ifges(5,5) * qJDD(4) - pkin(8) * t643 - qJD(4) * t728 + t777 * t635 - t773 * t636 - t727 * t806;
t647 = qJDD(1) * mrSges(3,2) - t783;
t781 = -mrSges(2,2) * t747 - mrSges(3,3) * t723 + mrSges(4,3) * t789 - pkin(7) * t634 - qJ(3) * t648 - t774 * t625 + t778 * t627 + qJ(2) * (t788 - t808) - pkin(1) * t647 + mrSges(4,2) * t718 + mrSges(3,2) * t725 + mrSges(2,1) * t746 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t633 = t810 * g(3) + t797;
t632 = t794 - t808;
t624 = t782 + (t801 * t780) - qJ(3) * t797 + t800 * qJDD(1) + mrSges(2,3) * t747 + mrSges(4,1) * t718 - mrSges(3,1) * t723 + (qJ(3) * m(4) + mrSges(4,3) + t809) * g(3) + pkin(3) * t634 - pkin(1) * t633 + pkin(2) * t632;
t623 = -qJ(2) * t633 - mrSges(2,3) * t746 + pkin(2) * t648 + mrSges(3,1) * t725 + t774 * t627 + t778 * t625 + pkin(3) * t786 + pkin(7) * t797 - mrSges(4,1) * t789 - t800 * t780 + t801 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t798; -m(1) * g(2) + t807; (-m(1) - m(2) + t810) * g(3) + t797; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t807 + t779 * t623 - t775 * t624; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t798 + t775 * t623 + t779 * t624; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t781; t781; t647; t632; t782; t785; t784;];
tauJB  = t1;
