% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:51:17
% EndTime: 2019-12-31 18:51:22
% DurationCPUTime: 4.10s
% Computational Cost: add. (41255->269), mult. (97020->323), div. (0->0), fcn. (67871->8), ass. (0->119)
t771 = Ifges(5,1) + Ifges(6,1);
t765 = Ifges(5,4) + Ifges(6,4);
t764 = Ifges(5,5) + Ifges(6,5);
t770 = Ifges(5,2) + Ifges(6,2);
t763 = Ifges(5,6) + Ifges(6,6);
t769 = Ifges(5,3) + Ifges(6,3);
t730 = qJD(1) ^ 2;
t721 = sin(pkin(8));
t722 = cos(pkin(8));
t724 = sin(qJ(3));
t727 = cos(qJ(3));
t737 = t721 * t724 - t722 * t727;
t704 = t737 * qJD(1);
t738 = t721 * t727 + t722 * t724;
t705 = t738 * qJD(1);
t752 = qJD(3) * t705;
t692 = -t737 * qJDD(1) - t752;
t753 = qJD(3) * t704;
t693 = t738 * qJDD(1) - t753;
t723 = sin(qJ(4));
t726 = cos(qJ(4));
t697 = qJD(3) * t726 - t705 * t723;
t667 = qJD(4) * t697 + qJDD(3) * t723 + t693 * t726;
t698 = qJD(3) * t723 + t705 * t726;
t669 = -mrSges(6,1) * t697 + mrSges(6,2) * t698;
t725 = sin(qJ(1));
t728 = cos(qJ(1));
t711 = -g(1) * t728 - g(2) * t725;
t706 = -pkin(1) * t730 + qJDD(1) * qJ(2) + t711;
t751 = qJD(1) * qJD(2);
t747 = -g(3) * t722 - 0.2e1 * t721 * t751;
t767 = pkin(2) * t722;
t680 = (-pkin(6) * qJDD(1) + t730 * t767 - t706) * t721 + t747;
t695 = -g(3) * t721 + (t706 + 0.2e1 * t751) * t722;
t750 = qJDD(1) * t722;
t719 = t722 ^ 2;
t760 = t719 * t730;
t681 = -pkin(2) * t760 + pkin(6) * t750 + t695;
t653 = t724 * t680 + t727 * t681;
t690 = pkin(3) * t704 - pkin(7) * t705;
t729 = qJD(3) ^ 2;
t648 = -pkin(3) * t729 + qJDD(3) * pkin(7) - t690 * t704 + t653;
t718 = t721 ^ 2;
t710 = g(1) * t725 - t728 * g(2);
t742 = qJDD(2) - t710;
t691 = (-pkin(1) - t767) * qJDD(1) + (-qJ(2) + (-t718 - t719) * pkin(6)) * t730 + t742;
t651 = (-t693 + t753) * pkin(7) + (-t692 + t752) * pkin(3) + t691;
t644 = -t648 * t723 + t726 * t651;
t689 = qJDD(4) - t692;
t702 = qJD(4) + t704;
t640 = -0.2e1 * qJD(5) * t698 + (t697 * t702 - t667) * qJ(5) + (t697 * t698 + t689) * pkin(4) + t644;
t673 = -mrSges(6,2) * t702 + mrSges(6,3) * t697;
t749 = m(6) * t640 + t689 * mrSges(6,1) + t702 * t673;
t637 = -t667 * mrSges(6,3) - t669 * t698 + t749;
t645 = t726 * t648 + t723 * t651;
t666 = -qJD(4) * t698 + qJDD(3) * t726 - t693 * t723;
t675 = pkin(4) * t702 - qJ(5) * t698;
t696 = t697 ^ 2;
t642 = -pkin(4) * t696 + t666 * qJ(5) + 0.2e1 * qJD(5) * t697 - t675 * t702 + t645;
t756 = t765 * t697 + t771 * t698 + t764 * t702;
t757 = -t770 * t697 - t765 * t698 - t763 * t702;
t768 = mrSges(5,1) * t644 + mrSges(6,1) * t640 - mrSges(5,2) * t645 - mrSges(6,2) * t642 + pkin(4) * t637 + t763 * t666 + t764 * t667 + t769 * t689 - t756 * t697 - t757 * t698;
t766 = -mrSges(5,2) - mrSges(6,2);
t761 = mrSges(3,2) * t721;
t670 = -mrSges(5,1) * t697 + mrSges(5,2) * t698;
t674 = -mrSges(5,2) * t702 + mrSges(5,3) * t697;
t630 = m(5) * t644 + mrSges(5,1) * t689 + t674 * t702 + (-t669 - t670) * t698 + (-mrSges(5,3) - mrSges(6,3)) * t667 + t749;
t748 = m(6) * t642 + t666 * mrSges(6,3) + t697 * t669;
t676 = mrSges(6,1) * t702 - mrSges(6,3) * t698;
t755 = -mrSges(5,1) * t702 + mrSges(5,3) * t698 - t676;
t633 = m(5) * t645 + t666 * mrSges(5,3) + t670 * t697 + t766 * t689 + t755 * t702 + t748;
t628 = -t630 * t723 + t726 * t633;
t687 = mrSges(4,1) * t704 + mrSges(4,2) * t705;
t700 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t705;
t625 = m(4) * t653 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t692 - qJD(3) * t700 - t687 * t704 + t628;
t652 = t680 * t727 - t724 * t681;
t647 = -qJDD(3) * pkin(3) - pkin(7) * t729 + t705 * t690 - t652;
t643 = -t666 * pkin(4) - qJ(5) * t696 + t675 * t698 + qJDD(5) + t647;
t743 = -m(6) * t643 + t666 * mrSges(6,1) + t697 * t673;
t636 = -m(5) * t647 + t666 * mrSges(5,1) + t766 * t667 + t697 * t674 + t755 * t698 + t743;
t699 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t704;
t635 = m(4) * t652 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t693 + qJD(3) * t699 - t687 * t705 + t636;
t617 = t724 * t625 + t727 * t635;
t694 = -t706 * t721 + t747;
t736 = mrSges(3,3) * qJDD(1) + t730 * (-mrSges(3,1) * t722 + t761);
t615 = m(3) * t694 - t736 * t721 + t617;
t744 = t727 * t625 - t635 * t724;
t616 = m(3) * t695 + t736 * t722 + t744;
t745 = -t615 * t721 + t722 * t616;
t608 = m(2) * t711 - mrSges(2,1) * t730 - qJDD(1) * mrSges(2,2) + t745;
t703 = -qJDD(1) * pkin(1) - qJ(2) * t730 + t742;
t627 = t726 * t630 + t723 * t633;
t734 = m(4) * t691 - t692 * mrSges(4,1) + mrSges(4,2) * t693 + t704 * t699 + t700 * t705 + t627;
t732 = -m(3) * t703 + mrSges(3,1) * t750 - t734 + (t718 * t730 + t760) * mrSges(3,3);
t620 = t732 + m(2) * t710 - mrSges(2,2) * t730 + (mrSges(2,1) - t761) * qJDD(1);
t759 = t725 * t608 + t728 * t620;
t610 = t722 * t615 + t721 * t616;
t758 = -t763 * t697 - t764 * t698 - t769 * t702;
t739 = Ifges(3,5) * t721 + Ifges(3,6) * t722;
t754 = t730 * t739;
t746 = t728 * t608 - t620 * t725;
t741 = Ifges(3,1) * t721 + Ifges(3,4) * t722;
t740 = Ifges(3,4) * t721 + Ifges(3,2) * t722;
t638 = t667 * mrSges(6,2) + t676 * t698 - t743;
t618 = -mrSges(5,1) * t647 + mrSges(5,3) * t645 - mrSges(6,1) * t643 + mrSges(6,3) * t642 - pkin(4) * t638 + qJ(5) * t748 + (-qJ(5) * t676 + t756) * t702 + t758 * t698 + (-mrSges(6,2) * qJ(5) + t763) * t689 + t765 * t667 + t770 * t666;
t626 = mrSges(5,2) * t647 + mrSges(6,2) * t643 - mrSges(5,3) * t644 - mrSges(6,3) * t640 - qJ(5) * t637 + t765 * t666 + t771 * t667 + t764 * t689 - t758 * t697 + t757 * t702;
t682 = Ifges(4,5) * t705 - Ifges(4,6) * t704 + Ifges(4,3) * qJD(3);
t683 = Ifges(4,4) * t705 - Ifges(4,2) * t704 + Ifges(4,6) * qJD(3);
t605 = mrSges(4,2) * t691 - mrSges(4,3) * t652 + Ifges(4,1) * t693 + Ifges(4,4) * t692 + Ifges(4,5) * qJDD(3) - pkin(7) * t627 - qJD(3) * t683 - t618 * t723 + t626 * t726 - t682 * t704;
t684 = Ifges(4,1) * t705 - Ifges(4,4) * t704 + Ifges(4,5) * qJD(3);
t611 = -mrSges(4,1) * t691 + mrSges(4,3) * t653 + Ifges(4,4) * t693 + Ifges(4,2) * t692 + Ifges(4,6) * qJDD(3) - pkin(3) * t627 + qJD(3) * t684 - t705 * t682 - t768;
t602 = -mrSges(3,1) * t703 + mrSges(3,3) * t695 - pkin(2) * t734 + pkin(6) * t744 + t740 * qJDD(1) + t724 * t605 + t727 * t611 - t721 * t754;
t604 = mrSges(3,2) * t703 - mrSges(3,3) * t694 - pkin(6) * t617 + t741 * qJDD(1) + t605 * t727 - t611 * t724 + t722 * t754;
t622 = qJDD(1) * t761 - t732;
t735 = mrSges(2,1) * t710 - mrSges(2,2) * t711 + Ifges(2,3) * qJDD(1) - pkin(1) * t622 + qJ(2) * t745 + t722 * t602 + t721 * t604;
t731 = mrSges(4,1) * t652 - mrSges(4,2) * t653 + Ifges(4,5) * t693 + Ifges(4,6) * t692 + Ifges(4,3) * qJDD(3) + pkin(3) * t636 + pkin(7) * t628 + t726 * t618 + t723 * t626 + t705 * t683 + t704 * t684;
t600 = -t731 - pkin(1) * t610 + mrSges(2,1) * g(3) + (Ifges(2,6) - t739) * qJDD(1) - pkin(2) * t617 - mrSges(3,1) * t694 + mrSges(3,2) * t695 + mrSges(2,3) * t711 + (-t721 * t740 + t722 * t741 + Ifges(2,5)) * t730;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t710 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t730 - qJ(2) * t610 - t602 * t721 + t604 * t722;
t1 = [-m(1) * g(1) + t746; -m(1) * g(2) + t759; (-m(1) - m(2)) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t759 + t728 * t599 - t725 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t746 + t725 * t599 + t728 * t600; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t735; t735; t622; t731; t768; t638;];
tauJB = t1;
