% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:15
% EndTime: 2019-12-31 18:32:19
% DurationCPUTime: 3.40s
% Computational Cost: add. (31451->272), mult. (75469->327), div. (0->0), fcn. (50910->8), ass. (0->122)
t794 = Ifges(4,1) + Ifges(5,2);
t787 = Ifges(4,4) + Ifges(5,6);
t786 = Ifges(4,5) - Ifges(5,4);
t793 = -Ifges(4,2) - Ifges(5,3);
t785 = Ifges(4,6) - Ifges(5,5);
t792 = Ifges(4,3) + Ifges(5,1);
t749 = qJD(1) ^ 2;
t744 = sin(qJ(3));
t742 = cos(pkin(8));
t789 = cos(qJ(3));
t770 = t742 * t789;
t741 = sin(pkin(8));
t774 = t741 * qJD(1);
t721 = -qJD(1) * t770 + t744 * t774;
t712 = t721 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t772 = qJDD(1) * t741;
t760 = t789 * t741 + t742 * t744;
t722 = t760 * qJD(1);
t775 = t722 * qJD(3);
t703 = -qJDD(1) * t770 + t744 * t772 + t775;
t714 = t722 * pkin(4) - qJD(3) * pkin(7);
t720 = t721 ^ 2;
t737 = t741 ^ 2;
t738 = t742 ^ 2;
t745 = sin(qJ(1));
t747 = cos(qJ(1));
t727 = t745 * g(1) - t747 * g(2);
t765 = qJDD(2) - t727;
t788 = pkin(2) * t742;
t702 = (-pkin(1) - t788) * qJDD(1) + (-qJ(2) + (-t737 - t738) * pkin(6)) * t749 + t765;
t776 = t721 * qJD(3);
t704 = t760 * qJDD(1) - t776;
t790 = -2 * qJD(4);
t750 = pkin(3) * t775 + t722 * t790 + (-t704 + t776) * qJ(4) + t702;
t654 = -t720 * pkin(4) - t722 * t714 + (pkin(3) + pkin(7)) * t703 + t750;
t728 = -t747 * g(1) - t745 * g(2);
t723 = -t749 * pkin(1) + qJDD(1) * qJ(2) + t728;
t773 = qJD(1) * qJD(2);
t769 = -t742 * g(3) - 0.2e1 * t741 * t773;
t682 = (-pkin(6) * qJDD(1) + t749 * t788 - t723) * t741 + t769;
t706 = -t741 * g(3) + (t723 + 0.2e1 * t773) * t742;
t771 = qJDD(1) * t742;
t782 = t738 * t749;
t683 = -pkin(2) * t782 + pkin(6) * t771 + t706;
t663 = t789 * t682 - t744 * t683;
t694 = t721 * pkin(3) - t722 * qJ(4);
t748 = qJD(3) ^ 2;
t661 = -qJDD(3) * pkin(3) - t748 * qJ(4) + t722 * t694 + qJDD(4) - t663;
t655 = (t721 * t722 - qJDD(3)) * pkin(7) + (t704 + t776) * pkin(4) + t661;
t743 = sin(qJ(5));
t746 = cos(qJ(5));
t652 = -t743 * t654 + t746 * t655;
t707 = -t743 * qJD(3) + t746 * t721;
t673 = t707 * qJD(5) + t746 * qJDD(3) + t743 * t703;
t708 = t746 * qJD(3) + t743 * t721;
t674 = -t707 * mrSges(6,1) + t708 * mrSges(6,2);
t718 = qJD(5) + t722;
t678 = -t718 * mrSges(6,2) + t707 * mrSges(6,3);
t701 = qJDD(5) + t704;
t649 = m(6) * t652 + t701 * mrSges(6,1) - t673 * mrSges(6,3) - t708 * t674 + t718 * t678;
t653 = t746 * t654 + t743 * t655;
t672 = -t708 * qJD(5) - t743 * qJDD(3) + t746 * t703;
t679 = t718 * mrSges(6,1) - t708 * mrSges(6,3);
t650 = m(6) * t653 - t701 * mrSges(6,2) + t672 * mrSges(6,3) + t707 * t674 - t718 * t679;
t640 = t746 * t649 + t743 * t650;
t696 = -t721 * mrSges(5,2) - t722 * mrSges(5,3);
t757 = -m(5) * t661 - t704 * mrSges(5,1) - t722 * t696 - t640;
t639 = qJDD(3) * mrSges(5,2) + qJD(3) * t712 - t757;
t664 = t744 * t682 + t789 * t683;
t756 = -t748 * pkin(3) + qJDD(3) * qJ(4) - t721 * t694 + t664;
t657 = -t703 * pkin(4) - t720 * pkin(7) + ((2 * qJD(4)) + t714) * qJD(3) + t756;
t665 = Ifges(6,5) * t708 + Ifges(6,6) * t707 + Ifges(6,3) * t718;
t667 = Ifges(6,1) * t708 + Ifges(6,4) * t707 + Ifges(6,5) * t718;
t641 = -mrSges(6,1) * t657 + mrSges(6,3) * t653 + Ifges(6,4) * t673 + Ifges(6,2) * t672 + Ifges(6,6) * t701 - t708 * t665 + t718 * t667;
t666 = Ifges(6,4) * t708 + Ifges(6,2) * t707 + Ifges(6,6) * t718;
t642 = mrSges(6,2) * t657 - mrSges(6,3) * t652 + Ifges(6,1) * t673 + Ifges(6,4) * t672 + Ifges(6,5) * t701 + t707 * t665 - t718 * t666;
t660 = qJD(3) * t790 - t756;
t713 = t722 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t759 = -m(6) * t657 + t672 * mrSges(6,1) - t673 * mrSges(6,2) + t707 * t678 - t708 * t679;
t754 = -m(5) * t660 + qJDD(3) * mrSges(5,3) + qJD(3) * t713 - t759;
t777 = t786 * qJD(3) - t787 * t721 + t794 * t722;
t778 = t785 * qJD(3) + t793 * t721 + t787 * t722;
t791 = t792 * qJDD(3) - t785 * t703 + t786 * t704 + t777 * t721 + t778 * t722 + mrSges(4,1) * t663 - mrSges(4,2) * t664 + mrSges(5,2) * t661 - mrSges(5,3) * t660 - pkin(3) * t639 - pkin(7) * t640 + qJ(4) * (-t703 * mrSges(5,1) - t721 * t696 + t754) - t743 * t641 + t746 * t642;
t783 = mrSges(3,2) * t741;
t695 = t721 * mrSges(4,1) + t722 * mrSges(4,2);
t710 = -qJD(3) * mrSges(4,2) - t721 * mrSges(4,3);
t637 = m(4) * t663 - t704 * mrSges(4,3) - t722 * t695 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t710 - t712) * qJD(3) + t757;
t711 = qJD(3) * mrSges(4,1) - t722 * mrSges(4,3);
t645 = m(4) * t664 - qJDD(3) * mrSges(4,2) - qJD(3) * t711 + (-t695 - t696) * t721 + (-mrSges(4,3) - mrSges(5,1)) * t703 + t754;
t631 = t789 * t637 + t744 * t645;
t705 = -t741 * t723 + t769;
t761 = mrSges(3,3) * qJDD(1) + t749 * (-mrSges(3,1) * t742 + t783);
t629 = m(3) * t705 - t761 * t741 + t631;
t766 = -t744 * t637 + t789 * t645;
t630 = m(3) * t706 + t761 * t742 + t766;
t767 = -t741 * t629 + t742 * t630;
t622 = m(2) * t728 - t749 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t767;
t719 = -qJDD(1) * pkin(1) - t749 * qJ(2) + t765;
t659 = t703 * pkin(3) + t750;
t780 = -t743 * t649 + t746 * t650;
t638 = m(5) * t659 - t703 * mrSges(5,2) - t704 * mrSges(5,3) - t721 * t712 - t722 * t713 + t780;
t753 = m(4) * t702 + t703 * mrSges(4,1) + t704 * mrSges(4,2) + t721 * t710 + t722 * t711 + t638;
t752 = -m(3) * t719 + mrSges(3,1) * t771 - t753 + (t737 * t749 + t782) * mrSges(3,3);
t633 = t752 + (mrSges(2,1) - t783) * qJDD(1) + m(2) * t727 - t749 * mrSges(2,2);
t781 = t745 * t622 + t747 * t633;
t624 = t742 * t629 + t741 * t630;
t779 = -t792 * qJD(3) + t785 * t721 - t786 * t722;
t768 = t747 * t622 - t745 * t633;
t764 = Ifges(3,1) * t741 + Ifges(3,4) * t742;
t763 = Ifges(3,4) * t741 + Ifges(3,2) * t742;
t762 = Ifges(3,5) * t741 + Ifges(3,6) * t742;
t619 = -mrSges(4,1) * t702 - mrSges(5,1) * t660 + mrSges(5,2) * t659 + mrSges(4,3) * t664 - pkin(3) * t638 - pkin(4) * t759 - pkin(7) * t780 + t777 * qJD(3) + t785 * qJDD(3) - t746 * t641 - t743 * t642 + t793 * t703 + t787 * t704 + t779 * t722;
t755 = mrSges(6,1) * t652 - mrSges(6,2) * t653 + Ifges(6,5) * t673 + Ifges(6,6) * t672 + Ifges(6,3) * t701 + t708 * t666 - t707 * t667;
t625 = mrSges(5,1) * t661 + mrSges(4,2) * t702 - mrSges(4,3) * t663 - mrSges(5,3) * t659 + pkin(4) * t640 - qJ(4) * t638 - t778 * qJD(3) + t786 * qJDD(3) - t787 * t703 + t794 * t704 + t779 * t721 + t755;
t725 = t762 * qJD(1);
t615 = -mrSges(3,1) * t719 + mrSges(3,3) * t706 - pkin(2) * t753 + pkin(6) * t766 + t763 * qJDD(1) + t789 * t619 + t744 * t625 - t725 * t774;
t618 = t742 * qJD(1) * t725 + mrSges(3,2) * t719 - mrSges(3,3) * t705 - pkin(6) * t631 + t764 * qJDD(1) - t744 * t619 + t789 * t625;
t635 = mrSges(3,2) * t772 - t752;
t758 = mrSges(2,1) * t727 - mrSges(2,2) * t728 + Ifges(2,3) * qJDD(1) - pkin(1) * t635 + qJ(2) * t767 + t742 * t615 + t741 * t618;
t616 = -pkin(2) * t631 - pkin(1) * t624 + mrSges(2,1) * g(3) + (Ifges(2,6) - t762) * qJDD(1) - mrSges(3,1) * t705 + mrSges(3,2) * t706 + mrSges(2,3) * t728 + (-t741 * t763 + t742 * t764 + Ifges(2,5)) * t749 - t791;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t727 + Ifges(2,5) * qJDD(1) - t749 * Ifges(2,6) - qJ(2) * t624 - t741 * t615 + t742 * t618;
t1 = [-m(1) * g(1) + t768; -m(1) * g(2) + t781; (-m(1) - m(2)) * g(3) + t624; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t781 + t747 * t613 - t745 * t616; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t768 + t745 * t613 + t747 * t616; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t758; t758; t635; t791; t639; t755;];
tauJB = t1;
