% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:32
% EndTime: 2019-12-31 18:12:34
% DurationCPUTime: 2.10s
% Computational Cost: add. (16137->245), mult. (38752->285), div. (0->0), fcn. (24311->6), ass. (0->108)
t774 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t754 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t753 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t773 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t752 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t772 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t726 = qJD(1) ^ 2;
t723 = sin(qJ(1));
t724 = cos(qJ(1));
t705 = -t724 * g(1) - t723 * g(2);
t700 = -t726 * pkin(1) + qJDD(1) * qJ(2) + t705;
t720 = sin(pkin(7));
t721 = cos(pkin(7));
t757 = qJD(1) * qJD(2);
t745 = -t721 * g(3) - 0.2e1 * t720 * t757;
t766 = pkin(2) * t721;
t652 = (-pkin(6) * qJDD(1) + t726 * t766 - t700) * t720 + t745;
t682 = -t720 * g(3) + (t700 + 0.2e1 * t757) * t721;
t755 = qJDD(1) * t721;
t715 = t721 ^ 2;
t763 = t715 * t726;
t653 = -pkin(2) * t763 + pkin(6) * t755 + t682;
t722 = sin(qJ(3));
t768 = cos(qJ(3));
t645 = t768 * t652 - t722 * t653;
t746 = t721 * t768;
t758 = t720 * qJD(1);
t698 = -qJD(1) * t746 + t722 * t758;
t734 = t768 * t720 + t721 * t722;
t699 = t734 * qJD(1);
t669 = t698 * pkin(3) - t699 * qJ(4);
t725 = qJD(3) ^ 2;
t644 = -qJDD(3) * pkin(3) - t725 * qJ(4) + t699 * t669 + qJDD(4) - t645;
t671 = -t698 * mrSges(5,2) - t699 * mrSges(5,3);
t760 = t698 * qJD(3);
t680 = t734 * qJDD(1) - t760;
t771 = -m(5) * t644 - t680 * mrSges(5,1) - t699 * t671;
t668 = -t699 * mrSges(6,2) + t698 * mrSges(6,3);
t638 = -0.2e1 * qJD(5) * qJD(3) + (t698 * t699 - qJDD(3)) * qJ(5) + (t680 + t760) * pkin(4) + t644;
t691 = -t698 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t741 = -m(6) * t638 + qJDD(3) * mrSges(6,3) + qJD(3) * t691;
t634 = t680 * mrSges(6,1) + t699 * t668 - t741;
t690 = t698 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t632 = qJDD(3) * mrSges(5,2) + qJD(3) * t690 + t634 - t771;
t756 = qJDD(1) * t720;
t759 = t699 * qJD(3);
t679 = -qJDD(1) * t746 + t722 * t756 + t759;
t688 = t699 * pkin(4) - qJD(3) * qJ(5);
t697 = t698 ^ 2;
t646 = t722 * t652 + t768 * t653;
t731 = -t725 * pkin(3) + qJDD(3) * qJ(4) - t698 * t669 + t646;
t640 = -t679 * pkin(4) - t697 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t688) * qJD(3) + t731;
t769 = -2 * qJD(4);
t643 = qJD(3) * t769 - t731;
t692 = t699 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t689 = t699 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t750 = m(6) * t640 + qJDD(3) * mrSges(6,2) + qJD(3) * t689;
t733 = -m(5) * t643 + qJDD(3) * mrSges(5,3) + qJD(3) * t692 + t750;
t747 = t753 * qJD(3) - t754 * t698 + t774 * t699;
t748 = t752 * qJD(3) + t773 * t698 + t754 * t699;
t761 = -t668 - t671;
t770 = t772 * qJDD(3) - t752 * t679 + t753 * t680 + t747 * t698 + t748 * t699 + mrSges(4,1) * t645 - mrSges(4,2) * t646 + mrSges(5,2) * t644 + mrSges(6,2) * t640 - mrSges(5,3) * t643 - mrSges(6,3) * t638 - pkin(3) * t632 + qJ(4) * (t761 * t698 + (-mrSges(5,1) - mrSges(6,1)) * t679 + t733) - qJ(5) * t634;
t765 = -mrSges(6,1) - mrSges(4,3);
t764 = mrSges(3,2) * t720;
t670 = t698 * mrSges(4,1) + t699 * mrSges(4,2);
t686 = -qJD(3) * mrSges(4,2) - t698 * mrSges(4,3);
t627 = m(4) * t645 + (-t668 - t670) * t699 + t765 * t680 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t686 - t690) * qJD(3) + t741 + t771;
t687 = qJD(3) * mrSges(4,1) - t699 * mrSges(4,3);
t630 = m(4) * t646 - qJDD(3) * mrSges(4,2) - qJD(3) * t687 + (-t670 + t761) * t698 + (-mrSges(5,1) + t765) * t679 + t733;
t621 = t768 * t627 + t722 * t630;
t681 = -t720 * t700 + t745;
t735 = mrSges(3,3) * qJDD(1) + t726 * (-mrSges(3,1) * t721 + t764);
t619 = m(3) * t681 - t735 * t720 + t621;
t742 = -t722 * t627 + t768 * t630;
t620 = m(3) * t682 + t735 * t721 + t742;
t743 = -t720 * t619 + t721 * t620;
t611 = m(2) * t705 - t726 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t743;
t704 = t723 * g(1) - t724 * g(2);
t740 = qJDD(2) - t704;
t696 = -qJDD(1) * pkin(1) - t726 * qJ(2) + t740;
t714 = t720 ^ 2;
t678 = (-pkin(1) - t766) * qJDD(1) + (-qJ(2) + (-t714 - t715) * pkin(6)) * t726 + t740;
t727 = pkin(3) * t759 + t699 * t769 + (-t680 + t760) * qJ(4) + t678;
t642 = t679 * pkin(3) + t727;
t637 = -t697 * pkin(4) + 0.2e1 * qJD(5) * t698 - t699 * t688 + (pkin(3) + qJ(5)) * t679 + t727;
t736 = m(6) * t637 - t680 * mrSges(6,2) + t679 * mrSges(6,3) - t699 * t689 + t698 * t691;
t631 = m(5) * t642 - t679 * mrSges(5,2) - t680 * mrSges(5,3) - t698 * t690 - t699 * t692 + t736;
t730 = m(4) * t678 + t679 * mrSges(4,1) + t680 * mrSges(4,2) + t698 * t686 + t699 * t687 + t631;
t728 = -m(3) * t696 + mrSges(3,1) * t755 - t730 + (t714 * t726 + t763) * mrSges(3,3);
t623 = t728 + (mrSges(2,1) - t764) * qJDD(1) - t726 * mrSges(2,2) + m(2) * t704;
t762 = t723 * t611 + t724 * t623;
t613 = t721 * t619 + t720 * t620;
t749 = -t772 * qJD(3) + t752 * t698 - t753 * t699;
t744 = t724 * t611 - t723 * t623;
t739 = Ifges(3,1) * t720 + Ifges(3,4) * t721;
t738 = Ifges(3,4) * t720 + Ifges(3,2) * t721;
t737 = Ifges(3,5) * t720 + Ifges(3,6) * t721;
t635 = -t679 * mrSges(6,1) - t698 * t668 + t750;
t614 = -mrSges(4,1) * t678 - mrSges(5,1) * t643 + mrSges(6,1) * t640 + mrSges(5,2) * t642 + mrSges(4,3) * t646 - mrSges(6,3) * t637 - pkin(3) * t631 + pkin(4) * t635 - qJ(5) * t736 + t747 * qJD(3) + t752 * qJDD(3) + t773 * t679 + t754 * t680 + t749 * t699;
t615 = mrSges(5,1) * t644 + mrSges(6,1) * t638 + mrSges(4,2) * t678 - mrSges(6,2) * t637 - mrSges(4,3) * t645 - mrSges(5,3) * t642 + pkin(4) * t634 - qJ(4) * t631 - t748 * qJD(3) + t753 * qJDD(3) - t754 * t679 + t774 * t680 + t749 * t698;
t702 = t737 * qJD(1);
t606 = -mrSges(3,1) * t696 + mrSges(3,3) * t682 - pkin(2) * t730 + pkin(6) * t742 + t738 * qJDD(1) + t768 * t614 + t722 * t615 - t702 * t758;
t608 = t721 * qJD(1) * t702 + mrSges(3,2) * t696 - mrSges(3,3) * t681 - pkin(6) * t621 + t739 * qJDD(1) - t722 * t614 + t768 * t615;
t625 = mrSges(3,2) * t756 - t728;
t732 = mrSges(2,1) * t704 - mrSges(2,2) * t705 + Ifges(2,3) * qJDD(1) - pkin(1) * t625 + qJ(2) * t743 + t721 * t606 + t720 * t608;
t604 = mrSges(2,1) * g(3) + (Ifges(2,6) - t737) * qJDD(1) + mrSges(2,3) * t705 - mrSges(3,1) * t681 + mrSges(3,2) * t682 - pkin(1) * t613 - pkin(2) * t621 + (-t720 * t738 + t721 * t739 + Ifges(2,5)) * t726 - t770;
t603 = -mrSges(2,2) * g(3) - mrSges(2,3) * t704 + Ifges(2,5) * qJDD(1) - t726 * Ifges(2,6) - qJ(2) * t613 - t720 * t606 + t721 * t608;
t1 = [-m(1) * g(1) + t744; -m(1) * g(2) + t762; (-m(1) - m(2)) * g(3) + t613; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t762 + t724 * t603 - t723 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t744 + t723 * t603 + t724 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t732; t732; t625; t770; t632; t635;];
tauJB = t1;
