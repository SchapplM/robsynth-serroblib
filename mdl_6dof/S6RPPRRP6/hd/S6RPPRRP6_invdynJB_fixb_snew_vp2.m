% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:00:23
% EndTime: 2019-05-05 15:00:26
% DurationCPUTime: 2.00s
% Computational Cost: add. (17677->265), mult. (32377->298), div. (0->0), fcn. (16320->6), ass. (0->105)
t781 = Ifges(6,1) + Ifges(7,1);
t772 = Ifges(6,4) - Ifges(7,5);
t771 = -Ifges(6,5) - Ifges(7,4);
t780 = Ifges(6,2) + Ifges(7,3);
t770 = Ifges(6,6) - Ifges(7,6);
t779 = -Ifges(6,3) - Ifges(7,2);
t734 = sin(qJ(1));
t736 = cos(qJ(1));
t708 = g(1) * t734 - g(2) * t736;
t738 = qJD(1) ^ 2;
t685 = -qJDD(1) * pkin(1) - qJ(2) * t738 + qJDD(2) - t708;
t675 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t685;
t709 = -t736 * g(1) - t734 * g(2);
t778 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t709;
t732 = sin(qJ(5));
t735 = cos(qJ(4));
t761 = qJD(1) * t735;
t775 = cos(qJ(5));
t699 = -qJD(4) * t775 + t732 * t761;
t733 = sin(qJ(4));
t760 = qJD(1) * qJD(4);
t754 = t733 * t760;
t704 = qJDD(1) * t735 - t754;
t661 = -qJD(5) * t699 + qJDD(4) * t732 + t704 * t775;
t700 = qJD(4) * t732 + t761 * t775;
t668 = mrSges(7,1) * t699 - mrSges(7,3) * t700;
t671 = -t738 * pkin(7) - t675;
t753 = t735 * t760;
t703 = -qJDD(1) * t733 - t753;
t645 = (-t704 + t754) * pkin(8) + (-t703 + t753) * pkin(4) + t671;
t676 = qJDD(3) + (-pkin(1) - qJ(3)) * t738 + t778;
t672 = -qJDD(1) * pkin(7) + t676;
t663 = -g(3) * t735 + t672 * t733;
t702 = (pkin(4) * t733 - pkin(8) * t735) * qJD(1);
t737 = qJD(4) ^ 2;
t762 = qJD(1) * t733;
t648 = -pkin(4) * t737 + qJDD(4) * pkin(8) - t702 * t762 + t663;
t642 = t645 * t775 - t648 * t732;
t667 = pkin(5) * t699 - qJ(6) * t700;
t698 = qJDD(5) - t703;
t711 = qJD(5) + t762;
t710 = t711 ^ 2;
t640 = -pkin(5) * t698 - qJ(6) * t710 + t667 * t700 + qJDD(6) - t642;
t680 = -mrSges(7,2) * t699 + mrSges(7,3) * t711;
t748 = -m(7) * t640 + mrSges(7,1) * t698 + t680 * t711;
t637 = t661 * mrSges(7,2) + t700 * t668 - t748;
t643 = t645 * t732 + t648 * t775;
t639 = -pkin(5) * t710 + qJ(6) * t698 + 0.2e1 * qJD(6) * t711 - t667 * t699 + t643;
t660 = qJD(5) * t700 - qJDD(4) * t775 + t704 * t732;
t679 = -mrSges(7,1) * t711 + mrSges(7,2) * t700;
t755 = m(7) * t639 + mrSges(7,3) * t698 + t679 * t711;
t764 = t699 * t772 - t700 * t781 + t711 * t771;
t765 = t699 * t780 - t700 * t772 - t711 * t770;
t777 = -t770 * t660 - t771 * t661 - t779 * t698 - t764 * t699 - t765 * t700 + mrSges(6,1) * t642 - mrSges(7,1) * t640 - mrSges(6,2) * t643 + mrSges(7,3) * t639 - pkin(5) * t637 + qJ(6) * (-t660 * mrSges(7,2) - t699 * t668 + t755);
t776 = -m(3) - m(4);
t774 = mrSges(2,1) - mrSges(3,2);
t773 = -mrSges(6,3) - mrSges(7,2);
t768 = mrSges(4,3) * t738;
t683 = pkin(1) * t738 - t778;
t701 = (mrSges(5,1) * t733 + mrSges(5,2) * t735) * qJD(1);
t707 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t761;
t678 = mrSges(6,1) * t711 - mrSges(6,3) * t700;
t763 = -mrSges(6,1) * t699 - mrSges(6,2) * t700 - t668;
t631 = m(6) * t643 - t698 * mrSges(6,2) + t660 * t773 - t711 * t678 + t699 * t763 + t755;
t677 = -mrSges(6,2) * t711 - mrSges(6,3) * t699;
t633 = m(6) * t642 + t698 * mrSges(6,1) + t661 * t773 + t711 * t677 + t700 * t763 + t748;
t750 = t631 * t775 - t633 * t732;
t624 = m(5) * t663 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t703 - qJD(4) * t707 - t701 * t762 + t750;
t662 = t733 * g(3) + t735 * t672;
t706 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t762;
t647 = -qJDD(4) * pkin(4) - t737 * pkin(8) + t702 * t761 - t662;
t641 = -0.2e1 * qJD(6) * t700 + (t699 * t711 - t661) * qJ(6) + (t700 * t711 + t660) * pkin(5) + t647;
t635 = m(7) * t641 + mrSges(7,1) * t660 - mrSges(7,3) * t661 - t679 * t700 + t680 * t699;
t740 = -m(6) * t647 - mrSges(6,1) * t660 - mrSges(6,2) * t661 - t677 * t699 - t678 * t700 - t635;
t628 = m(5) * t662 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t704 + qJD(4) * t706 - t701 * t761 + t740;
t612 = t624 * t733 + t628 * t735;
t749 = m(4) * t676 + qJDD(1) * mrSges(4,2) + t612;
t744 = -m(3) * t683 + mrSges(3,2) * t738 + qJDD(1) * mrSges(3,3) + t749;
t608 = m(2) * t709 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t738 + t744;
t626 = t631 * t732 + t633 * t775;
t745 = -m(5) * t671 + mrSges(5,1) * t703 - mrSges(5,2) * t704 - t706 * t762 - t707 * t761 - t626;
t618 = m(4) * t675 - mrSges(4,2) * t738 - qJDD(1) * mrSges(4,3) + t745;
t741 = -m(3) * t685 + mrSges(3,3) * t738 - t618;
t614 = m(2) * t708 - t738 * mrSges(2,2) + qJDD(1) * t774 + t741;
t767 = t608 * t734 + t614 * t736;
t766 = t699 * t770 + t700 * t771 + t711 * t779;
t757 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t756 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t752 = t608 * t736 - t614 * t734;
t751 = t624 * t735 - t733 * t628;
t619 = -mrSges(6,1) * t647 - mrSges(7,1) * t641 + mrSges(7,2) * t639 + mrSges(6,3) * t643 - pkin(5) * t635 - t660 * t780 + t661 * t772 + t698 * t770 + t700 * t766 - t711 * t764;
t621 = mrSges(6,2) * t647 + mrSges(7,2) * t640 - mrSges(6,3) * t642 - mrSges(7,3) * t641 - qJ(6) * t635 - t660 * t772 + t661 * t781 - t698 * t771 + t699 * t766 + t711 * t765;
t688 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t735 - Ifges(5,2) * t733) * qJD(1);
t689 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t735 - Ifges(5,4) * t733) * qJD(1);
t743 = mrSges(5,1) * t662 - mrSges(5,2) * t663 + Ifges(5,5) * t704 + Ifges(5,6) * t703 + Ifges(5,3) * qJDD(4) + pkin(4) * t740 + pkin(8) * t750 + t619 * t775 + t621 * t732 + t688 * t761 + t689 * t762;
t687 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t735 - Ifges(5,6) * t733) * qJD(1);
t604 = mrSges(5,2) * t671 - mrSges(5,3) * t662 + Ifges(5,1) * t704 + Ifges(5,4) * t703 + Ifges(5,5) * qJDD(4) - pkin(8) * t626 - qJD(4) * t688 - t619 * t732 + t621 * t775 - t687 * t762;
t605 = -mrSges(5,1) * t671 + mrSges(5,3) * t663 + Ifges(5,4) * t704 + Ifges(5,2) * t703 + Ifges(5,6) * qJDD(4) - pkin(4) * t626 + qJD(4) * t689 - t687 * t761 - t777;
t616 = qJDD(1) * mrSges(3,2) - t741;
t739 = -mrSges(2,2) * t709 - mrSges(3,3) * t683 - mrSges(4,3) * t675 - pkin(7) * t612 - qJ(3) * t618 + t735 * t604 - t605 * t733 + qJ(2) * (t744 - t768) - pkin(1) * t616 + mrSges(4,2) * t676 + mrSges(3,2) * t685 + mrSges(2,1) * t708 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t611 = g(3) * t776 + t751;
t610 = t749 - t768;
t602 = t757 * t738 - qJ(3) * t751 + pkin(3) * t612 - pkin(1) * t611 + mrSges(2,3) * t709 + mrSges(4,1) * t676 - mrSges(3,1) * t683 + t743 + t756 * qJDD(1) + pkin(2) * t610 + (m(4) * qJ(3) + mrSges(4,3) + t774) * g(3);
t601 = -qJ(2) * t611 - mrSges(2,3) * t708 + pkin(2) * t618 + mrSges(3,1) * t685 + t735 * t605 + pkin(3) * t745 + pkin(7) * t751 + t733 * t604 + mrSges(4,1) * t675 - t756 * t738 + t757 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t752; -m(1) * g(2) + t767; (-m(1) - m(2) + t776) * g(3) + t751; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t767 + t736 * t601 - t734 * t602; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t752 + t734 * t601 + t736 * t602; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t739; t739; t616; t610; t743; t777; t637;];
tauJB  = t1;
