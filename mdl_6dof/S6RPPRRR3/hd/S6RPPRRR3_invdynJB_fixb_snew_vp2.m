% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-05-05 15:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:31:16
% EndTime: 2019-05-05 15:31:22
% DurationCPUTime: 5.63s
% Computational Cost: add. (75035->295), mult. (140349->356), div. (0->0), fcn. (84378->10), ass. (0->124)
t759 = sin(qJ(1));
t763 = cos(qJ(1));
t732 = t759 * g(1) - g(2) * t763;
t723 = qJDD(1) * pkin(1) + t732;
t733 = -g(1) * t763 - g(2) * t759;
t765 = qJD(1) ^ 2;
t725 = -pkin(1) * t765 + t733;
t754 = sin(pkin(10));
t755 = cos(pkin(10));
t698 = t754 * t723 + t755 * t725;
t791 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t698;
t790 = -pkin(2) - pkin(7);
t789 = mrSges(3,1) - mrSges(4,2);
t788 = -Ifges(4,4) + Ifges(3,5);
t787 = Ifges(4,5) - Ifges(3,6);
t697 = t755 * t723 - t754 * t725;
t775 = -t765 * qJ(3) + qJDD(3) - t697;
t683 = qJDD(1) * t790 + t775;
t751 = -g(3) + qJDD(2);
t758 = sin(qJ(4));
t762 = cos(qJ(4));
t677 = t758 * t683 + t762 * t751;
t724 = (mrSges(5,1) * t758 + mrSges(5,2) * t762) * qJD(1);
t784 = qJD(1) * qJD(4);
t738 = t762 * t784;
t727 = -t758 * qJDD(1) - t738;
t785 = qJD(1) * t762;
t731 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t785;
t740 = t758 * qJD(1);
t680 = t765 * t790 - t791;
t781 = t758 * t784;
t728 = qJDD(1) * t762 - t781;
t666 = (-t728 + t781) * pkin(8) + (-t727 + t738) * pkin(4) + t680;
t726 = (pkin(4) * t758 - pkin(8) * t762) * qJD(1);
t764 = qJD(4) ^ 2;
t673 = -pkin(4) * t764 + qJDD(4) * pkin(8) - t726 * t740 + t677;
t757 = sin(qJ(5));
t761 = cos(qJ(5));
t655 = t761 * t666 - t757 * t673;
t721 = qJD(4) * t761 - t757 * t785;
t693 = qJD(5) * t721 + qJDD(4) * t757 + t728 * t761;
t720 = qJDD(5) - t727;
t722 = qJD(4) * t757 + t761 * t785;
t735 = t740 + qJD(5);
t653 = (t721 * t735 - t693) * pkin(9) + (t721 * t722 + t720) * pkin(5) + t655;
t656 = t757 * t666 + t761 * t673;
t692 = -qJD(5) * t722 + qJDD(4) * t761 - t728 * t757;
t702 = pkin(5) * t735 - pkin(9) * t722;
t719 = t721 ^ 2;
t654 = -pkin(5) * t719 + pkin(9) * t692 - t702 * t735 + t656;
t756 = sin(qJ(6));
t760 = cos(qJ(6));
t651 = t653 * t760 - t654 * t756;
t695 = t721 * t760 - t722 * t756;
t663 = qJD(6) * t695 + t692 * t756 + t693 * t760;
t696 = t721 * t756 + t722 * t760;
t674 = -mrSges(7,1) * t695 + mrSges(7,2) * t696;
t734 = qJD(6) + t735;
t681 = -mrSges(7,2) * t734 + mrSges(7,3) * t695;
t713 = qJDD(6) + t720;
t647 = m(7) * t651 + mrSges(7,1) * t713 - t663 * mrSges(7,3) - t674 * t696 + t681 * t734;
t652 = t653 * t756 + t654 * t760;
t662 = -qJD(6) * t696 + t692 * t760 - t693 * t756;
t682 = mrSges(7,1) * t734 - mrSges(7,3) * t696;
t648 = m(7) * t652 - mrSges(7,2) * t713 + t662 * mrSges(7,3) + t674 * t695 - t682 * t734;
t640 = t760 * t647 + t756 * t648;
t699 = -mrSges(6,1) * t721 + mrSges(6,2) * t722;
t700 = -mrSges(6,2) * t735 + mrSges(6,3) * t721;
t638 = m(6) * t655 + mrSges(6,1) * t720 - mrSges(6,3) * t693 - t699 * t722 + t700 * t735 + t640;
t701 = mrSges(6,1) * t735 - mrSges(6,3) * t722;
t776 = -t647 * t756 + t760 * t648;
t639 = m(6) * t656 - mrSges(6,2) * t720 + mrSges(6,3) * t692 + t699 * t721 - t701 * t735 + t776;
t777 = -t638 * t757 + t761 * t639;
t632 = m(5) * t677 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t727 - qJD(4) * t731 - t724 * t740 + t777;
t676 = t683 * t762 - t758 * t751;
t730 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t740;
t672 = -qJDD(4) * pkin(4) - pkin(8) * t764 + t726 * t785 - t676;
t657 = -pkin(5) * t692 - pkin(9) * t719 + t702 * t722 + t672;
t772 = m(7) * t657 - t662 * mrSges(7,1) + t663 * mrSges(7,2) - t695 * t681 + t682 * t696;
t768 = -m(6) * t672 + t692 * mrSges(6,1) - mrSges(6,2) * t693 + t721 * t700 - t701 * t722 - t772;
t643 = m(5) * t676 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t728 + qJD(4) * t730 - t724 * t785 + t768;
t624 = t758 * t632 + t762 * t643;
t686 = -qJDD(1) * pkin(2) + t775;
t774 = -m(4) * t686 + t765 * mrSges(4,3) - t624;
t617 = m(3) * t697 - t765 * mrSges(3,2) + qJDD(1) * t789 + t774;
t684 = t765 * pkin(2) + t791;
t634 = t761 * t638 + t757 * t639;
t773 = -m(5) * t680 + mrSges(5,1) * t727 - t728 * mrSges(5,2) - t730 * t740 - t731 * t785 - t634;
t769 = -m(4) * t684 + t765 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t773;
t629 = m(3) * t698 - mrSges(3,1) * t765 - qJDD(1) * mrSges(3,2) + t769;
t614 = t755 * t617 + t754 * t629;
t611 = m(2) * t732 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t765 + t614;
t779 = -t617 * t754 + t755 * t629;
t612 = m(2) * t733 - mrSges(2,1) * t765 - qJDD(1) * mrSges(2,2) + t779;
t786 = t763 * t611 + t759 * t612;
t780 = -t611 * t759 + t763 * t612;
t778 = t762 * t632 - t758 * t643;
t623 = m(4) * t751 + t778;
t622 = m(3) * t751 + t623;
t669 = Ifges(7,4) * t696 + Ifges(7,2) * t695 + Ifges(7,6) * t734;
t670 = Ifges(7,1) * t696 + Ifges(7,4) * t695 + Ifges(7,5) * t734;
t771 = -mrSges(7,1) * t651 + mrSges(7,2) * t652 - Ifges(7,5) * t663 - Ifges(7,6) * t662 - Ifges(7,3) * t713 - t696 * t669 + t695 * t670;
t668 = Ifges(7,5) * t696 + Ifges(7,6) * t695 + Ifges(7,3) * t734;
t641 = -mrSges(7,1) * t657 + mrSges(7,3) * t652 + Ifges(7,4) * t663 + Ifges(7,2) * t662 + Ifges(7,6) * t713 - t668 * t696 + t670 * t734;
t642 = mrSges(7,2) * t657 - mrSges(7,3) * t651 + Ifges(7,1) * t663 + Ifges(7,4) * t662 + Ifges(7,5) * t713 + t668 * t695 - t669 * t734;
t687 = Ifges(6,5) * t722 + Ifges(6,6) * t721 + Ifges(6,3) * t735;
t689 = Ifges(6,1) * t722 + Ifges(6,4) * t721 + Ifges(6,5) * t735;
t619 = -mrSges(6,1) * t672 + mrSges(6,3) * t656 + Ifges(6,4) * t693 + Ifges(6,2) * t692 + Ifges(6,6) * t720 - pkin(5) * t772 + pkin(9) * t776 + t760 * t641 + t756 * t642 - t722 * t687 + t735 * t689;
t688 = Ifges(6,4) * t722 + Ifges(6,2) * t721 + Ifges(6,6) * t735;
t626 = mrSges(6,2) * t672 - mrSges(6,3) * t655 + Ifges(6,1) * t693 + Ifges(6,4) * t692 + Ifges(6,5) * t720 - pkin(9) * t640 - t641 * t756 + t642 * t760 + t687 * t721 - t688 * t735;
t711 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t762 - Ifges(5,2) * t758) * qJD(1);
t712 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t762 - Ifges(5,4) * t758) * qJD(1);
t770 = mrSges(5,1) * t676 - mrSges(5,2) * t677 + Ifges(5,5) * t728 + Ifges(5,6) * t727 + Ifges(5,3) * qJDD(4) + pkin(4) * t768 + pkin(8) * t777 + t761 * t619 + t757 * t626 + t711 * t785 + t712 * t740;
t710 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t762 - Ifges(5,6) * t758) * qJD(1);
t607 = mrSges(5,2) * t680 - mrSges(5,3) * t676 + Ifges(5,1) * t728 + Ifges(5,4) * t727 + Ifges(5,5) * qJDD(4) - pkin(8) * t634 - qJD(4) * t711 - t619 * t757 + t626 * t761 - t710 * t740;
t766 = mrSges(6,1) * t655 - mrSges(6,2) * t656 + Ifges(6,5) * t693 + Ifges(6,6) * t692 + Ifges(6,3) * t720 + pkin(5) * t640 + t722 * t688 - t721 * t689 - t771;
t615 = -mrSges(5,1) * t680 + mrSges(5,3) * t677 + Ifges(5,4) * t728 + Ifges(5,2) * t727 + Ifges(5,6) * qJDD(4) - pkin(4) * t634 + qJD(4) * t712 - t710 * t785 - t766;
t621 = qJDD(1) * mrSges(4,2) - t774;
t767 = mrSges(2,1) * t732 + mrSges(3,1) * t697 - mrSges(2,2) * t733 - mrSges(3,2) * t698 + mrSges(4,2) * t686 - mrSges(4,3) * t684 + pkin(1) * t614 - pkin(2) * t621 - pkin(7) * t624 + qJ(3) * t769 + t762 * t607 - t615 * t758 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1);
t605 = (mrSges(3,2) - mrSges(4,3)) * t751 + t787 * t765 + t788 * qJDD(1) - mrSges(3,3) * t697 + mrSges(4,1) * t686 + t770 + pkin(3) * t624 - qJ(3) * t623;
t604 = -mrSges(4,1) * t684 + mrSges(3,3) * t698 - pkin(2) * t623 - pkin(3) * t773 - pkin(7) * t778 - qJDD(1) * t787 - t758 * t607 - t762 * t615 - t751 * t789 + t765 * t788;
t603 = -mrSges(2,2) * g(3) - mrSges(2,3) * t732 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t765 - qJ(2) * t614 - t604 * t754 + t605 * t755;
t602 = mrSges(2,1) * g(3) + mrSges(2,3) * t733 + t765 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t622 + qJ(2) * t779 + t755 * t604 + t754 * t605;
t1 = [-m(1) * g(1) + t780; -m(1) * g(2) + t786; (-m(1) - m(2)) * g(3) + t622; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t786 - t759 * t602 + t763 * t603; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t780 + t763 * t602 + t759 * t603; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t767; t767; t622; t621; t770; t766; -t771;];
tauJB  = t1;
