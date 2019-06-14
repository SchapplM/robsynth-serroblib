% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:25:32
% EndTime: 2019-05-06 09:25:41
% DurationCPUTime: 5.62s
% Computational Cost: add. (57371->344), mult. (124905->405), div. (0->0), fcn. (74667->8), ass. (0->132)
t774 = -2 * qJD(3);
t773 = Ifges(3,1) + Ifges(4,2);
t772 = Ifges(6,1) + Ifges(7,1);
t764 = Ifges(3,4) + Ifges(4,6);
t763 = Ifges(6,4) - Ifges(7,5);
t762 = Ifges(7,4) + Ifges(6,5);
t761 = Ifges(3,5) - Ifges(4,4);
t771 = Ifges(3,2) + Ifges(4,3);
t770 = Ifges(6,2) + Ifges(7,3);
t760 = Ifges(3,6) - Ifges(4,5);
t759 = Ifges(6,6) - Ifges(7,6);
t769 = (Ifges(3,3) + Ifges(4,1));
t768 = Ifges(6,3) + Ifges(7,2);
t725 = sin(qJ(1));
t727 = cos(qJ(1));
t709 = -g(1) * t727 - g(2) * t725;
t729 = qJD(1) ^ 2;
t684 = -pkin(1) * t729 + qJDD(1) * pkin(7) + t709;
t724 = sin(qJ(2));
t726 = cos(qJ(2));
t662 = -t724 * g(3) + t726 * t684;
t696 = (-pkin(2) * t726 - qJ(3) * t724) * qJD(1);
t728 = qJD(2) ^ 2;
t748 = qJD(1) * t726;
t646 = t728 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t774) - t696 * t748 - t662;
t767 = cos(qJ(5));
t766 = t729 * pkin(7);
t765 = -mrSges(6,3) - mrSges(7,2);
t661 = -t726 * g(3) - t724 * t684;
t697 = (mrSges(4,2) * t726 - mrSges(4,3) * t724) * qJD(1);
t698 = (-mrSges(3,1) * t726 + mrSges(3,2) * t724) * qJD(1);
t747 = qJD(1) * qJD(2);
t744 = t726 * t747;
t699 = qJDD(1) * t724 + t744;
t704 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t748;
t706 = -mrSges(4,1) * t748 - qJD(2) * mrSges(4,3);
t743 = t724 * t747;
t700 = qJDD(1) * t726 - t743;
t749 = qJD(1) * t724;
t705 = pkin(3) * t749 - (qJD(2) * qJ(4));
t720 = t726 ^ 2;
t708 = t725 * g(1) - t727 * g(2);
t738 = -qJDD(1) * pkin(1) - t708;
t734 = pkin(2) * t743 + t749 * t774 + (-t699 - t744) * qJ(3) + t738;
t622 = -t705 * t749 + (-pkin(3) * t720 - pkin(7)) * t729 + (-pkin(2) - qJ(4)) * t700 + t734;
t647 = -qJDD(2) * pkin(2) - t728 * qJ(3) + t696 * t749 + qJDD(3) - t661;
t637 = (-t724 * t726 * t729 - qJDD(2)) * qJ(4) + (t699 - t744) * pkin(3) + t647;
t721 = sin(pkin(9));
t722 = cos(pkin(9));
t691 = qJD(2) * t722 - t721 * t748;
t617 = -0.2e1 * qJD(4) * t691 - t721 * t622 + t722 * t637;
t667 = qJDD(2) * t722 - t700 * t721;
t690 = -qJD(2) * t721 - t722 * t748;
t614 = (t690 * t749 - t667) * pkin(8) + (t690 * t691 + t699) * pkin(4) + t617;
t618 = 0.2e1 * qJD(4) * t690 + t722 * t622 + t721 * t637;
t666 = -qJDD(2) * t721 - t700 * t722;
t668 = pkin(4) * t749 - pkin(8) * t691;
t689 = t690 ^ 2;
t616 = -pkin(4) * t689 + pkin(8) * t666 - t668 * t749 + t618;
t723 = sin(qJ(5));
t610 = t723 * t614 + t616 * t767;
t658 = t723 * t690 + t691 * t767;
t625 = t658 * qJD(5) - t666 * t767 + t723 * t667;
t712 = qJD(5) + t749;
t650 = mrSges(6,1) * t712 - mrSges(6,3) * t658;
t657 = -t690 * t767 + t723 * t691;
t695 = qJDD(5) + t699;
t640 = pkin(5) * t657 - qJ(6) * t658;
t710 = t712 ^ 2;
t607 = -pkin(5) * t710 + qJ(6) * t695 + 0.2e1 * qJD(6) * t712 - t640 * t657 + t610;
t651 = -mrSges(7,1) * t712 + mrSges(7,2) * t658;
t745 = m(7) * t607 + t695 * mrSges(7,3) + t712 * t651;
t641 = mrSges(7,1) * t657 - mrSges(7,3) * t658;
t753 = -mrSges(6,1) * t657 - mrSges(6,2) * t658 - t641;
t599 = m(6) * t610 - t695 * mrSges(6,2) + t625 * t765 - t712 * t650 + t657 * t753 + t745;
t609 = t614 * t767 - t723 * t616;
t626 = -t657 * qJD(5) + t723 * t666 + t667 * t767;
t648 = -mrSges(6,2) * t712 - mrSges(6,3) * t657;
t608 = -t695 * pkin(5) - t710 * qJ(6) + t658 * t640 + qJDD(6) - t609;
t649 = -mrSges(7,2) * t657 + mrSges(7,3) * t712;
t739 = -m(7) * t608 + t695 * mrSges(7,1) + t712 * t649;
t601 = m(6) * t609 + t695 * mrSges(6,1) + t626 * t765 + t712 * t648 + t658 * t753 + t739;
t594 = t723 * t599 + t601 * t767;
t659 = -mrSges(5,1) * t690 + mrSges(5,2) * t691;
t664 = -mrSges(5,2) * t749 + mrSges(5,3) * t690;
t592 = m(5) * t617 + mrSges(5,1) * t699 - mrSges(5,3) * t667 - t659 * t691 + t664 * t749 + t594;
t665 = mrSges(5,1) * t749 - mrSges(5,3) * t691;
t740 = t599 * t767 - t601 * t723;
t593 = m(5) * t618 - mrSges(5,2) * t699 + mrSges(5,3) * t666 + t659 * t690 - t665 * t749 + t740;
t589 = t722 * t592 + t721 * t593;
t735 = -m(4) * t647 - t699 * mrSges(4,1) - t589;
t587 = m(3) * t661 - t699 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t704 - t706) * qJD(2) + (-t697 - t698) * t749 + t735;
t703 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t749;
t707 = mrSges(4,1) * t749 + qJD(2) * mrSges(4,2);
t630 = -t720 * t729 * qJ(4) + t700 * pkin(3) + qJD(2) * t705 + qJDD(4) - t646;
t620 = -t666 * pkin(4) - t689 * pkin(8) + t691 * t668 + t630;
t612 = (t658 * t712 + t625) * pkin(5) + (t657 * t712 - t626) * qJ(6) + t620 - 0.2e1 * qJD(6) * t658;
t605 = m(7) * t612 + t625 * mrSges(7,1) - t626 * mrSges(7,3) + t657 * t649 - t658 * t651;
t733 = m(6) * t620 + t625 * mrSges(6,1) + t626 * mrSges(6,2) + t657 * t648 + t658 * t650 + t605;
t731 = -m(5) * t630 + t666 * mrSges(5,1) - t667 * mrSges(5,2) + t690 * t664 - t691 * t665 - t733;
t730 = -m(4) * t646 + qJDD(2) * mrSges(4,3) + qJD(2) * t707 + t697 * t748 - t731;
t604 = -qJD(2) * t703 + t730 + t698 * t748 - qJDD(2) * mrSges(3,2) + m(3) * t662 + (mrSges(3,3) + mrSges(4,1)) * t700;
t741 = -t587 * t724 + t726 * t604;
t582 = m(2) * t709 - mrSges(2,1) * t729 - qJDD(1) * mrSges(2,2) + t741;
t683 = t738 - t766;
t643 = -t700 * pkin(2) + t734 - t766;
t757 = -t721 * t592 + t722 * t593;
t737 = -m(4) * t643 - t700 * mrSges(4,2) + t707 * t749 - t757;
t732 = -m(3) * t683 + t704 * t748 + t700 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t699 + (-t703 * t724 - t706 * t726) * qJD(1) + t737;
t585 = m(2) * t708 + qJDD(1) * mrSges(2,1) - t729 * mrSges(2,2) + t732;
t758 = t725 * t582 + t727 * t585;
t583 = t726 * t587 + t724 * t604;
t756 = t770 * t657 - t763 * t658 - t759 * t712;
t755 = t759 * t657 - t762 * t658 - t768 * t712;
t754 = -t763 * t657 + t772 * t658 + t762 * t712;
t752 = (t769 * qJD(2)) + (t761 * t724 + t760 * t726) * qJD(1);
t751 = -t760 * qJD(2) + (-t764 * t724 - t771 * t726) * qJD(1);
t750 = t761 * qJD(2) + (t773 * t724 + t764 * t726) * qJD(1);
t742 = t727 * t582 - t585 * t725;
t654 = Ifges(5,1) * t691 + Ifges(5,4) * t690 + Ifges(5,5) * t749;
t653 = Ifges(5,4) * t691 + Ifges(5,2) * t690 + Ifges(5,6) * t749;
t652 = Ifges(5,5) * t691 + Ifges(5,6) * t690 + Ifges(5,3) * t749;
t596 = mrSges(6,2) * t620 + mrSges(7,2) * t608 - mrSges(6,3) * t609 - mrSges(7,3) * t612 - qJ(6) * t605 - t763 * t625 + t772 * t626 + t755 * t657 + t762 * t695 + t756 * t712;
t595 = -mrSges(6,1) * t620 - mrSges(7,1) * t612 + mrSges(7,2) * t607 + mrSges(6,3) * t610 - pkin(5) * t605 - t770 * t625 + t763 * t626 + t755 * t658 + t759 * t695 + t754 * t712;
t588 = -t699 * mrSges(4,3) + t706 * t748 - t737;
t579 = mrSges(5,2) * t630 - mrSges(5,3) * t617 + Ifges(5,1) * t667 + Ifges(5,4) * t666 + Ifges(5,5) * t699 - pkin(8) * t594 - t723 * t595 + t596 * t767 + t690 * t652 - t653 * t749;
t578 = -mrSges(5,1) * t630 + mrSges(5,3) * t618 + Ifges(5,4) * t667 + Ifges(5,2) * t666 + Ifges(5,6) * t699 - pkin(4) * t733 + pkin(8) * t740 + t595 * t767 + t723 * t596 - t691 * t652 + t654 * t749;
t577 = mrSges(5,1) * t617 - mrSges(5,2) * t618 + pkin(5) * t739 - mrSges(7,1) * t608 + mrSges(6,1) * t609 - mrSges(6,2) * t610 - t690 * t654 + t691 * t653 - qJ(3) * t588 + t751 * qJD(2) + t752 * t748 + (-qJ(6) * t641 + t754) * t657 + (-pkin(5) * t641 - t756) * t658 + mrSges(7,3) * t607 - mrSges(4,3) * t643 + mrSges(4,1) * t647 + qJ(6) * t745 + Ifges(5,6) * t666 + Ifges(5,5) * t667 + mrSges(3,2) * t683 + pkin(3) * t589 + pkin(4) * t594 - mrSges(3,3) * t661 + (-mrSges(7,2) * qJ(6) - t759) * t625 + t761 * qJDD(2) + (-mrSges(7,2) * pkin(5) + t762) * t626 + t764 * t700 + (Ifges(5,3) + t773) * t699 + t768 * t695;
t576 = -mrSges(3,1) * t683 - mrSges(4,1) * t646 + mrSges(4,2) * t643 + mrSges(3,3) * t662 - pkin(2) * t588 - pkin(3) * t731 - qJ(4) * t757 + t750 * qJD(2) + t760 * qJDD(2) - t722 * t578 - t721 * t579 + t764 * t699 + t771 * t700 - t752 * t749;
t575 = -pkin(1) * t583 + mrSges(2,3) * t709 - pkin(2) * (-qJD(2) * t706 + t735) - qJ(3) * t730 - t722 * t579 + t721 * t578 + qJ(4) * t589 - mrSges(3,1) * t661 + mrSges(3,2) * t662 - mrSges(4,2) * t647 + mrSges(4,3) * t646 + t729 * Ifges(2,5) + mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (-mrSges(4,1) * qJ(3) - t760) * t700 - t761 * t699 + (mrSges(4,2) * pkin(2) - t769) * qJDD(2) + (t750 * t726 + (pkin(2) * t697 + t751) * t724) * qJD(1);
t574 = -mrSges(2,2) * g(3) - mrSges(2,3) * t708 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t729 - pkin(7) * t583 - t576 * t724 + t577 * t726;
t1 = [-m(1) * g(1) + t742; -m(1) * g(2) + t758; (-m(1) - m(2)) * g(3) + t583; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t758 + t727 * t574 - t725 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t742 + t725 * t574 + t727 * t575; -mrSges(1,1) * g(2) + mrSges(2,1) * t708 + mrSges(1,2) * g(1) - mrSges(2,2) * t709 + Ifges(2,3) * qJDD(1) + pkin(1) * t732 + pkin(7) * t741 + t726 * t576 + t724 * t577;];
tauB  = t1;
