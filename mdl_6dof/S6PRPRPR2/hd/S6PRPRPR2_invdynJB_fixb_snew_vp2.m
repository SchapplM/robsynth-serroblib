% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:20:58
% EndTime: 2019-05-04 22:21:10
% DurationCPUTime: 12.00s
% Computational Cost: add. (213144->297), mult. (414753->383), div. (0->0), fcn. (292808->14), ass. (0->131)
t758 = sin(pkin(10));
t762 = cos(pkin(10));
t745 = t758 * g(1) - t762 * g(2);
t746 = -t762 * g(1) - t758 * g(2);
t755 = -g(3) + qJDD(1);
t766 = sin(qJ(2));
t763 = cos(pkin(6));
t769 = cos(qJ(2));
t789 = t763 * t769;
t759 = sin(pkin(6));
t791 = t759 * t769;
t708 = t745 * t789 - t766 * t746 + t755 * t791;
t706 = qJDD(2) * pkin(2) + t708;
t790 = t763 * t766;
t792 = t759 * t766;
t709 = t745 * t790 + t769 * t746 + t755 * t792;
t771 = qJD(2) ^ 2;
t707 = -t771 * pkin(2) + t709;
t757 = sin(pkin(11));
t761 = cos(pkin(11));
t692 = t757 * t706 + t761 * t707;
t690 = -t771 * pkin(3) + qJDD(2) * pkin(8) + t692;
t725 = -t759 * t745 + t763 * t755;
t720 = qJDD(3) + t725;
t765 = sin(qJ(4));
t768 = cos(qJ(4));
t686 = t768 * t690 + t765 * t720;
t741 = (-pkin(4) * t768 - qJ(5) * t765) * qJD(2);
t770 = qJD(4) ^ 2;
t786 = t768 * qJD(2);
t681 = -t770 * pkin(4) + qJDD(4) * qJ(5) + t741 * t786 + t686;
t691 = t761 * t706 - t757 * t707;
t689 = -qJDD(2) * pkin(3) - t771 * pkin(8) - t691;
t785 = qJD(2) * qJD(4);
t784 = t768 * t785;
t743 = t765 * qJDD(2) + t784;
t753 = t765 * t785;
t744 = t768 * qJDD(2) - t753;
t684 = (-t743 - t784) * qJ(5) + (-t744 + t753) * pkin(4) + t689;
t756 = sin(pkin(12));
t760 = cos(pkin(12));
t787 = qJD(2) * t765;
t737 = t756 * qJD(4) + t760 * t787;
t676 = -0.2e1 * qJD(5) * t737 - t756 * t681 + t760 * t684;
t723 = t756 * qJDD(4) + t760 * t743;
t736 = t760 * qJD(4) - t756 * t787;
t674 = (-t736 * t786 - t723) * pkin(9) + (t736 * t737 - t744) * pkin(5) + t676;
t677 = 0.2e1 * qJD(5) * t736 + t760 * t681 + t756 * t684;
t722 = t760 * qJDD(4) - t756 * t743;
t724 = -pkin(5) * t786 - t737 * pkin(9);
t735 = t736 ^ 2;
t675 = -t735 * pkin(5) + t722 * pkin(9) + t724 * t786 + t677;
t764 = sin(qJ(6));
t767 = cos(qJ(6));
t673 = t764 * t674 + t767 * t675;
t685 = -t765 * t690 + t768 * t720;
t680 = -qJDD(4) * pkin(4) - t770 * qJ(5) + t741 * t787 + qJDD(5) - t685;
t678 = -t722 * pkin(5) - t735 * pkin(9) + t737 * t724 + t680;
t715 = t764 * t736 + t767 * t737;
t694 = -t715 * qJD(6) + t767 * t722 - t764 * t723;
t714 = t767 * t736 - t764 * t737;
t695 = t714 * qJD(6) + t764 * t722 + t767 * t723;
t751 = qJD(6) - t786;
t696 = Ifges(7,5) * t715 + Ifges(7,6) * t714 + Ifges(7,3) * t751;
t698 = Ifges(7,1) * t715 + Ifges(7,4) * t714 + Ifges(7,5) * t751;
t739 = qJDD(6) - t744;
t662 = -mrSges(7,1) * t678 + mrSges(7,3) * t673 + Ifges(7,4) * t695 + Ifges(7,2) * t694 + Ifges(7,6) * t739 - t715 * t696 + t751 * t698;
t672 = t767 * t674 - t764 * t675;
t697 = Ifges(7,4) * t715 + Ifges(7,2) * t714 + Ifges(7,6) * t751;
t663 = mrSges(7,2) * t678 - mrSges(7,3) * t672 + Ifges(7,1) * t695 + Ifges(7,4) * t694 + Ifges(7,5) * t739 + t714 * t696 - t751 * t697;
t710 = Ifges(6,5) * t737 + Ifges(6,6) * t736 - Ifges(6,3) * t786;
t712 = Ifges(6,1) * t737 + Ifges(6,4) * t736 - Ifges(6,5) * t786;
t704 = -t751 * mrSges(7,2) + t714 * mrSges(7,3);
t705 = t751 * mrSges(7,1) - t715 * mrSges(7,3);
t775 = m(7) * t678 - t694 * mrSges(7,1) + t695 * mrSges(7,2) - t714 * t704 + t715 * t705;
t700 = -t714 * mrSges(7,1) + t715 * mrSges(7,2);
t667 = m(7) * t672 + t739 * mrSges(7,1) - t695 * mrSges(7,3) - t715 * t700 + t751 * t704;
t668 = m(7) * t673 - t739 * mrSges(7,2) + t694 * mrSges(7,3) + t714 * t700 - t751 * t705;
t780 = -t764 * t667 + t767 * t668;
t649 = -mrSges(6,1) * t680 + mrSges(6,3) * t677 + Ifges(6,4) * t723 + Ifges(6,2) * t722 - Ifges(6,6) * t744 - pkin(5) * t775 + pkin(9) * t780 + t767 * t662 + t764 * t663 - t737 * t710 - t712 * t786;
t661 = t767 * t667 + t764 * t668;
t711 = Ifges(6,4) * t737 + Ifges(6,2) * t736 - Ifges(6,6) * t786;
t650 = mrSges(6,2) * t680 - mrSges(6,3) * t676 + Ifges(6,1) * t723 + Ifges(6,4) * t722 - Ifges(6,5) * t744 - pkin(9) * t661 - t764 * t662 + t767 * t663 + t736 * t710 + t711 * t786;
t716 = -t736 * mrSges(6,1) + t737 * mrSges(6,2);
t777 = mrSges(6,2) * t786 + t736 * mrSges(6,3);
t659 = m(6) * t676 - t744 * mrSges(6,1) - t723 * mrSges(6,3) - t737 * t716 - t777 * t786 + t661;
t721 = -mrSges(6,1) * t786 - t737 * mrSges(6,3);
t660 = m(6) * t677 + t744 * mrSges(6,2) + t722 * mrSges(6,3) + t736 * t716 + t721 * t786 + t780;
t657 = -t756 * t659 + t760 * t660;
t671 = m(6) * t680 - t722 * mrSges(6,1) + t723 * mrSges(6,2) + t737 * t721 - t736 * t777 + t775;
t732 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t765 + Ifges(5,2) * t768) * qJD(2);
t733 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t765 + Ifges(5,4) * t768) * qJD(2);
t793 = mrSges(5,1) * t685 - mrSges(5,2) * t686 + Ifges(5,5) * t743 + Ifges(5,6) * t744 + Ifges(5,3) * qJDD(4) - pkin(4) * t671 + qJ(5) * t657 + t760 * t649 + t756 * t650 + (t765 * t732 - t768 * t733) * qJD(2);
t742 = (-mrSges(5,1) * t768 + mrSges(5,2) * t765) * qJD(2);
t747 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t787;
t655 = m(5) * t686 - qJDD(4) * mrSges(5,2) + t744 * mrSges(5,3) - qJD(4) * t747 + t742 * t786 + t657;
t748 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t786;
t670 = m(5) * t685 + qJDD(4) * mrSges(5,1) - t743 * mrSges(5,3) + qJD(4) * t748 - t742 * t787 - t671;
t781 = t768 * t655 - t765 * t670;
t644 = m(4) * t692 - t771 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t781;
t656 = t760 * t659 + t756 * t660;
t774 = -m(5) * t689 + t744 * mrSges(5,1) - t743 * mrSges(5,2) - t747 * t787 + t748 * t786 - t656;
t652 = m(4) * t691 + qJDD(2) * mrSges(4,1) - t771 * mrSges(4,2) + t774;
t640 = t757 * t644 + t761 * t652;
t638 = m(3) * t708 + qJDD(2) * mrSges(3,1) - t771 * mrSges(3,2) + t640;
t782 = t761 * t644 - t757 * t652;
t639 = m(3) * t709 - t771 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t782;
t648 = t765 * t655 + t768 * t670;
t647 = m(4) * t720 + t648;
t646 = m(3) * t725 + t647;
t626 = t638 * t789 + t639 * t790 - t759 * t646;
t624 = m(2) * t745 + t626;
t631 = -t766 * t638 + t769 * t639;
t630 = m(2) * t746 + t631;
t788 = t762 * t624 + t758 * t630;
t625 = t638 * t791 + t639 * t792 + t763 * t646;
t783 = -t758 * t624 + t762 * t630;
t779 = m(2) * t755 + t625;
t731 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t765 + Ifges(5,6) * t768) * qJD(2);
t632 = mrSges(5,2) * t689 - mrSges(5,3) * t685 + Ifges(5,1) * t743 + Ifges(5,4) * t744 + Ifges(5,5) * qJDD(4) - qJ(5) * t656 - qJD(4) * t732 - t756 * t649 + t760 * t650 + t731 * t786;
t773 = mrSges(7,1) * t672 - mrSges(7,2) * t673 + Ifges(7,5) * t695 + Ifges(7,6) * t694 + Ifges(7,3) * t739 + t715 * t697 - t714 * t698;
t641 = -t731 * t787 - t773 + Ifges(5,6) * qJDD(4) + (Ifges(5,2) + Ifges(6,3)) * t744 + t736 * t712 - t737 * t711 + Ifges(5,4) * t743 - Ifges(6,6) * t722 - Ifges(6,5) * t723 + qJD(4) * t733 + mrSges(5,3) * t686 - mrSges(5,1) * t689 - mrSges(6,1) * t676 + mrSges(6,2) * t677 - pkin(5) * t661 - pkin(4) * t656;
t622 = mrSges(4,2) * t720 - mrSges(4,3) * t691 + Ifges(4,5) * qJDD(2) - t771 * Ifges(4,6) - pkin(8) * t648 + t768 * t632 - t765 * t641;
t627 = -mrSges(4,1) * t720 + mrSges(4,3) * t692 + t771 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t648 - t793;
t619 = -mrSges(3,1) * t725 + mrSges(3,3) * t709 + t771 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t647 + qJ(3) * t782 + t757 * t622 + t761 * t627;
t620 = mrSges(3,2) * t725 - mrSges(3,3) * t708 + Ifges(3,5) * qJDD(2) - t771 * Ifges(3,6) - qJ(3) * t640 + t761 * t622 - t757 * t627;
t776 = pkin(7) * t631 + t619 * t769 + t620 * t766;
t621 = mrSges(3,1) * t708 - mrSges(3,2) * t709 + mrSges(4,1) * t691 - mrSges(4,2) * t692 + t765 * t632 + t768 * t641 + pkin(3) * t774 + pkin(8) * t781 + pkin(2) * t640 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t618 = mrSges(2,2) * t755 - mrSges(2,3) * t745 - t766 * t619 + t769 * t620 + (-t625 * t759 - t626 * t763) * pkin(7);
t617 = -mrSges(2,1) * t755 + mrSges(2,3) * t746 - pkin(1) * t625 - t759 * t621 + t763 * t776;
t1 = [-m(1) * g(1) + t783; -m(1) * g(2) + t788; -m(1) * g(3) + t779; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t788 - t758 * t617 + t762 * t618; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t783 + t762 * t617 + t758 * t618; -mrSges(1,1) * g(2) + mrSges(2,1) * t745 + mrSges(1,2) * g(1) - mrSges(2,2) * t746 + pkin(1) * t626 + t763 * t621 + t759 * t776; t779; t621; t647; t793; t671; t773;];
tauJB  = t1;
