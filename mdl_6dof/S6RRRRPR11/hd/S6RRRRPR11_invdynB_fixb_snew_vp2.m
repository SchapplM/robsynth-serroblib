% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 23:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 23:26:23
% EndTime: 2019-05-07 23:27:51
% DurationCPUTime: 66.71s
% Computational Cost: add. (1110213->398), mult. (2388904->517), div. (0->0), fcn. (1934336->14), ass. (0->163)
t758 = sin(pkin(6));
t764 = sin(qJ(2));
t769 = cos(qJ(2));
t785 = qJD(1) * qJD(2);
t746 = (-qJDD(1) * t769 + t764 * t785) * t758;
t795 = pkin(8) * t758;
t760 = cos(pkin(6));
t794 = g(3) * t760;
t793 = t758 * t764;
t792 = t758 * t769;
t791 = t760 * t764;
t790 = t760 * t769;
t765 = sin(qJ(1));
t770 = cos(qJ(1));
t750 = t765 * g(1) - g(2) * t770;
t771 = qJD(1) ^ 2;
t741 = qJDD(1) * pkin(1) + t771 * t795 + t750;
t751 = -g(1) * t770 - g(2) * t765;
t742 = -pkin(1) * t771 + qJDD(1) * t795 + t751;
t788 = t741 * t791 + t769 * t742;
t716 = -g(3) * t793 + t788;
t754 = qJD(1) * t760 + qJD(2);
t787 = qJD(1) * t758;
t784 = t764 * t787;
t739 = mrSges(3,1) * t754 - mrSges(3,3) * t784;
t743 = (-mrSges(3,1) * t769 + mrSges(3,2) * t764) * t787;
t753 = qJDD(1) * t760 + qJDD(2);
t744 = (-pkin(2) * t769 - pkin(9) * t764) * t787;
t752 = t754 ^ 2;
t786 = qJD(1) * t769;
t694 = -pkin(2) * t752 + pkin(9) * t753 + (-g(3) * t764 + t744 * t786) * t758 + t788;
t745 = (qJDD(1) * t764 + t769 * t785) * t758;
t695 = pkin(2) * t746 - pkin(9) * t745 - t794 + (-t741 + (pkin(2) * t764 - pkin(9) * t769) * t754 * qJD(1)) * t758;
t763 = sin(qJ(3));
t768 = cos(qJ(3));
t671 = t768 * t694 + t763 * t695;
t734 = t754 * t763 + t768 * t784;
t713 = -t734 * qJD(3) - t763 * t745 + t753 * t768;
t733 = t754 * t768 - t763 * t784;
t717 = -mrSges(4,1) * t733 + mrSges(4,2) * t734;
t783 = t758 * t786;
t749 = qJD(3) - t783;
t723 = mrSges(4,1) * t749 - mrSges(4,3) * t734;
t738 = qJDD(3) + t746;
t718 = -pkin(3) * t733 - pkin(10) * t734;
t747 = t749 ^ 2;
t661 = -pkin(3) * t747 + pkin(10) * t738 + t718 * t733 + t671;
t715 = -g(3) * t792 + t741 * t790 - t764 * t742;
t693 = -pkin(2) * t753 - pkin(9) * t752 + t744 * t784 - t715;
t714 = qJD(3) * t733 + t745 * t768 + t753 * t763;
t664 = (-t733 * t749 - t714) * pkin(10) + (t734 * t749 - t713) * pkin(3) + t693;
t762 = sin(qJ(4));
t767 = cos(qJ(4));
t650 = -t661 * t762 + t767 * t664;
t720 = -t734 * t762 + t749 * t767;
t681 = qJD(4) * t720 + t714 * t767 + t738 * t762;
t711 = qJDD(4) - t713;
t721 = t734 * t767 + t749 * t762;
t732 = qJD(4) - t733;
t643 = (t720 * t732 - t681) * qJ(5) + (t720 * t721 + t711) * pkin(4) + t650;
t651 = t767 * t661 + t762 * t664;
t680 = -qJD(4) * t721 - t714 * t762 + t738 * t767;
t704 = pkin(4) * t732 - qJ(5) * t721;
t719 = t720 ^ 2;
t645 = -pkin(4) * t719 + qJ(5) * t680 - t704 * t732 + t651;
t757 = sin(pkin(12));
t759 = cos(pkin(12));
t700 = t720 * t757 + t721 * t759;
t637 = -0.2e1 * qJD(5) * t700 + t759 * t643 - t645 * t757;
t667 = t680 * t757 + t681 * t759;
t699 = t720 * t759 - t721 * t757;
t635 = (t699 * t732 - t667) * pkin(11) + (t699 * t700 + t711) * pkin(5) + t637;
t638 = 0.2e1 * qJD(5) * t699 + t757 * t643 + t759 * t645;
t666 = t680 * t759 - t681 * t757;
t684 = pkin(5) * t732 - pkin(11) * t700;
t698 = t699 ^ 2;
t636 = -pkin(5) * t698 + pkin(11) * t666 - t684 * t732 + t638;
t761 = sin(qJ(6));
t766 = cos(qJ(6));
t633 = t635 * t766 - t636 * t761;
t676 = t699 * t766 - t700 * t761;
t649 = qJD(6) * t676 + t666 * t761 + t667 * t766;
t677 = t699 * t761 + t700 * t766;
t658 = -mrSges(7,1) * t676 + mrSges(7,2) * t677;
t730 = qJD(6) + t732;
t668 = -mrSges(7,2) * t730 + mrSges(7,3) * t676;
t706 = qJDD(6) + t711;
t629 = m(7) * t633 + mrSges(7,1) * t706 - mrSges(7,3) * t649 - t658 * t677 + t668 * t730;
t634 = t635 * t761 + t636 * t766;
t648 = -qJD(6) * t677 + t666 * t766 - t667 * t761;
t669 = mrSges(7,1) * t730 - mrSges(7,3) * t677;
t630 = m(7) * t634 - mrSges(7,2) * t706 + mrSges(7,3) * t648 + t658 * t676 - t669 * t730;
t623 = t766 * t629 + t761 * t630;
t678 = -mrSges(6,1) * t699 + mrSges(6,2) * t700;
t682 = -mrSges(6,2) * t732 + mrSges(6,3) * t699;
t621 = m(6) * t637 + mrSges(6,1) * t711 - mrSges(6,3) * t667 - t678 * t700 + t682 * t732 + t623;
t683 = mrSges(6,1) * t732 - mrSges(6,3) * t700;
t778 = -t629 * t761 + t766 * t630;
t622 = m(6) * t638 - mrSges(6,2) * t711 + mrSges(6,3) * t666 + t678 * t699 - t683 * t732 + t778;
t617 = t759 * t621 + t757 * t622;
t701 = -mrSges(5,1) * t720 + mrSges(5,2) * t721;
t703 = -mrSges(5,2) * t732 + mrSges(5,3) * t720;
t615 = m(5) * t650 + mrSges(5,1) * t711 - mrSges(5,3) * t681 - t701 * t721 + t703 * t732 + t617;
t705 = mrSges(5,1) * t732 - mrSges(5,3) * t721;
t779 = -t621 * t757 + t759 * t622;
t616 = m(5) * t651 - mrSges(5,2) * t711 + mrSges(5,3) * t680 + t701 * t720 - t705 * t732 + t779;
t780 = -t615 * t762 + t767 * t616;
t610 = m(4) * t671 - mrSges(4,2) * t738 + mrSges(4,3) * t713 + t717 * t733 - t723 * t749 + t780;
t670 = -t763 * t694 + t695 * t768;
t722 = -mrSges(4,2) * t749 + mrSges(4,3) * t733;
t660 = -pkin(3) * t738 - pkin(10) * t747 + t734 * t718 - t670;
t652 = -pkin(4) * t680 - qJ(5) * t719 + t721 * t704 + qJDD(5) + t660;
t640 = -pkin(5) * t666 - pkin(11) * t698 + t684 * t700 + t652;
t777 = m(7) * t640 - t648 * mrSges(7,1) + t649 * mrSges(7,2) - t676 * t668 + t677 * t669;
t774 = m(6) * t652 - t666 * mrSges(6,1) + mrSges(6,2) * t667 - t699 * t682 + t683 * t700 + t777;
t772 = -m(5) * t660 + t680 * mrSges(5,1) - mrSges(5,2) * t681 + t720 * t703 - t705 * t721 - t774;
t632 = m(4) * t670 + mrSges(4,1) * t738 - mrSges(4,3) * t714 - t717 * t734 + t722 * t749 + t772;
t781 = t768 * t610 - t632 * t763;
t601 = m(3) * t716 - mrSges(3,2) * t753 - mrSges(3,3) * t746 - t739 * t754 + t743 * t783 + t781;
t604 = t763 * t610 + t768 * t632;
t727 = -t741 * t758 - t794;
t740 = -mrSges(3,2) * t754 + mrSges(3,3) * t783;
t603 = m(3) * t727 + mrSges(3,1) * t746 + mrSges(3,2) * t745 + (t739 * t764 - t740 * t769) * t787 + t604;
t611 = t615 * t767 + t616 * t762;
t773 = -m(4) * t693 + t713 * mrSges(4,1) - mrSges(4,2) * t714 + t733 * t722 - t723 * t734 - t611;
t607 = m(3) * t715 + mrSges(3,1) * t753 - mrSges(3,3) * t745 + t740 * t754 - t743 * t784 + t773;
t590 = t601 * t791 - t603 * t758 + t607 * t790;
t588 = m(2) * t750 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t771 + t590;
t594 = t769 * t601 - t607 * t764;
t593 = m(2) * t751 - mrSges(2,1) * t771 - qJDD(1) * mrSges(2,2) + t594;
t789 = t770 * t588 + t765 * t593;
t589 = t601 * t793 + t760 * t603 + t607 * t792;
t782 = -t588 * t765 + t770 * t593;
t653 = Ifges(7,5) * t677 + Ifges(7,6) * t676 + Ifges(7,3) * t730;
t655 = Ifges(7,1) * t677 + Ifges(7,4) * t676 + Ifges(7,5) * t730;
t624 = -mrSges(7,1) * t640 + mrSges(7,3) * t634 + Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t706 - t653 * t677 + t655 * t730;
t654 = Ifges(7,4) * t677 + Ifges(7,2) * t676 + Ifges(7,6) * t730;
t625 = mrSges(7,2) * t640 - mrSges(7,3) * t633 + Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t706 + t653 * t676 - t654 * t730;
t672 = Ifges(6,5) * t700 + Ifges(6,6) * t699 + Ifges(6,3) * t732;
t674 = Ifges(6,1) * t700 + Ifges(6,4) * t699 + Ifges(6,5) * t732;
t612 = -mrSges(6,1) * t652 + mrSges(6,3) * t638 + Ifges(6,4) * t667 + Ifges(6,2) * t666 + Ifges(6,6) * t711 - pkin(5) * t777 + pkin(11) * t778 + t766 * t624 + t761 * t625 - t700 * t672 + t732 * t674;
t673 = Ifges(6,4) * t700 + Ifges(6,2) * t699 + Ifges(6,6) * t732;
t613 = mrSges(6,2) * t652 - mrSges(6,3) * t637 + Ifges(6,1) * t667 + Ifges(6,4) * t666 + Ifges(6,5) * t711 - pkin(11) * t623 - t624 * t761 + t625 * t766 + t672 * t699 - t673 * t732;
t685 = Ifges(5,5) * t721 + Ifges(5,6) * t720 + Ifges(5,3) * t732;
t687 = Ifges(5,1) * t721 + Ifges(5,4) * t720 + Ifges(5,5) * t732;
t596 = -mrSges(5,1) * t660 + mrSges(5,3) * t651 + Ifges(5,4) * t681 + Ifges(5,2) * t680 + Ifges(5,6) * t711 - pkin(4) * t774 + qJ(5) * t779 + t759 * t612 + t757 * t613 - t721 * t685 + t732 * t687;
t686 = Ifges(5,4) * t721 + Ifges(5,2) * t720 + Ifges(5,6) * t732;
t597 = mrSges(5,2) * t660 - mrSges(5,3) * t650 + Ifges(5,1) * t681 + Ifges(5,4) * t680 + Ifges(5,5) * t711 - qJ(5) * t617 - t612 * t757 + t613 * t759 + t685 * t720 - t686 * t732;
t707 = Ifges(4,5) * t734 + Ifges(4,6) * t733 + Ifges(4,3) * t749;
t708 = Ifges(4,4) * t734 + Ifges(4,2) * t733 + Ifges(4,6) * t749;
t586 = mrSges(4,2) * t693 - mrSges(4,3) * t670 + Ifges(4,1) * t714 + Ifges(4,4) * t713 + Ifges(4,5) * t738 - pkin(10) * t611 - t596 * t762 + t597 * t767 + t707 * t733 - t708 * t749;
t709 = Ifges(4,1) * t734 + Ifges(4,4) * t733 + Ifges(4,5) * t749;
t595 = (-Ifges(5,3) - Ifges(6,3)) * t711 + Ifges(4,6) * t738 + t749 * t709 - t734 * t707 + t720 * t687 - t721 * t686 + Ifges(4,4) * t714 - Ifges(7,3) * t706 + Ifges(4,2) * t713 + t699 * t674 - t700 * t673 - mrSges(4,1) * t693 - Ifges(5,6) * t680 - Ifges(5,5) * t681 + t676 * t655 - t677 * t654 + mrSges(4,3) * t671 - Ifges(6,6) * t666 - Ifges(6,5) * t667 + mrSges(5,2) * t651 - Ifges(7,6) * t648 - Ifges(7,5) * t649 - mrSges(5,1) * t650 - mrSges(6,1) * t637 + mrSges(6,2) * t638 - pkin(3) * t611 + mrSges(7,2) * t634 - mrSges(7,1) * t633 - pkin(5) * t623 - pkin(4) * t617;
t724 = Ifges(3,3) * t754 + (Ifges(3,5) * t764 + Ifges(3,6) * t769) * t787;
t725 = Ifges(3,6) * t754 + (Ifges(3,4) * t764 + Ifges(3,2) * t769) * t787;
t584 = mrSges(3,2) * t727 - mrSges(3,3) * t715 + Ifges(3,1) * t745 - Ifges(3,4) * t746 + Ifges(3,5) * t753 - pkin(9) * t604 + t586 * t768 - t595 * t763 + t724 * t783 - t725 * t754;
t726 = Ifges(3,5) * t754 + (Ifges(3,1) * t764 + Ifges(3,4) * t769) * t787;
t585 = Ifges(3,4) * t745 - Ifges(3,2) * t746 + Ifges(3,6) * t753 - t724 * t784 + t754 * t726 - mrSges(3,1) * t727 + mrSges(3,3) * t716 - Ifges(4,5) * t714 - Ifges(4,6) * t713 - Ifges(4,3) * t738 - t734 * t708 + t733 * t709 - mrSges(4,1) * t670 + mrSges(4,2) * t671 - t762 * t597 - t767 * t596 - pkin(3) * t772 - pkin(10) * t780 - pkin(2) * t604;
t775 = pkin(8) * t594 + t584 * t764 + t585 * t769;
t583 = Ifges(3,5) * t745 - Ifges(3,6) * t746 + Ifges(3,3) * t753 + mrSges(3,1) * t715 - mrSges(3,2) * t716 + t763 * t586 + t768 * t595 + pkin(2) * t773 + pkin(9) * t781 + (t725 * t764 - t726 * t769) * t787;
t582 = -mrSges(2,2) * g(3) - mrSges(2,3) * t750 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t771 + t584 * t769 - t585 * t764 + (-t589 * t758 - t590 * t760) * pkin(8);
t581 = mrSges(2,1) * g(3) + mrSges(2,3) * t751 + Ifges(2,5) * t771 + Ifges(2,6) * qJDD(1) - pkin(1) * t589 - t583 * t758 + t760 * t775;
t1 = [-m(1) * g(1) + t782; -m(1) * g(2) + t789; (-m(1) - m(2)) * g(3) + t589; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t789 - t765 * t581 + t770 * t582; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t782 + t770 * t581 + t765 * t582; -mrSges(1,1) * g(2) + mrSges(2,1) * t750 + mrSges(1,2) * g(1) - mrSges(2,2) * t751 + Ifges(2,3) * qJDD(1) + pkin(1) * t590 + t583 * t760 + t758 * t775;];
tauB  = t1;
