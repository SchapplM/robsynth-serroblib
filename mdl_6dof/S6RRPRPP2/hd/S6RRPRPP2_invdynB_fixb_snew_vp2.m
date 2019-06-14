% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:23:22
% EndTime: 2019-05-06 12:23:34
% DurationCPUTime: 6.11s
% Computational Cost: add. (62127->345), mult. (140248->408), div. (0->0), fcn. (94776->8), ass. (0->132)
t781 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t761 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t760 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t780 = Ifges(5,2) + Ifges(6,3) + Ifges(7,2);
t759 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t779 = -Ifges(5,3) - Ifges(6,2) - Ifges(7,3);
t735 = sin(qJ(2));
t737 = cos(qJ(2));
t762 = qJD(1) * qJD(2);
t721 = qJDD(1) * t735 + t737 * t762;
t722 = qJDD(1) * t737 - t735 * t762;
t733 = sin(pkin(9));
t772 = cos(pkin(9));
t694 = t721 * t772 + t733 * t722;
t710 = (t733 * t737 + t735 * t772) * qJD(1);
t734 = sin(qJ(4));
t775 = cos(qJ(4));
t696 = -qJD(2) * t775 + t710 * t734;
t660 = -t696 * qJD(4) + t734 * qJDD(2) + t694 * t775;
t764 = qJD(1) * t737;
t765 = qJD(1) * t735;
t709 = -t733 * t765 + t772 * t764;
t708 = qJD(4) - t709;
t771 = t696 * t708;
t778 = (-t660 + t771) * qJ(5);
t697 = t734 * qJD(2) + t710 * t775;
t777 = -0.2e1 * t697;
t776 = 2 * qJD(5);
t740 = qJD(1) ^ 2;
t774 = pkin(2) * t740;
t773 = -mrSges(5,3) - mrSges(6,2);
t676 = mrSges(7,2) * t708 + mrSges(7,3) * t696;
t770 = t708 * t676;
t736 = sin(qJ(1));
t738 = cos(qJ(1));
t727 = -g(1) * t738 - g(2) * t736;
t716 = -pkin(1) * t740 + qJDD(1) * pkin(7) + t727;
t769 = t735 * t716;
t667 = qJDD(2) * pkin(2) - t721 * qJ(3) - t769 + (qJ(3) * t762 + t735 * t774 - g(3)) * t737;
t699 = -g(3) * t735 + t737 * t716;
t723 = qJD(2) * pkin(2) - qJ(3) * t765;
t732 = t737 ^ 2;
t668 = qJ(3) * t722 - qJD(2) * t723 - t732 * t774 + t699;
t635 = 0.2e1 * qJD(3) * t709 + t733 * t667 + t772 * t668;
t686 = -mrSges(4,1) * t709 + mrSges(4,2) * t710;
t693 = -t721 * t733 + t772 * t722;
t701 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t710;
t687 = -pkin(3) * t709 - pkin(8) * t710;
t739 = qJD(2) ^ 2;
t631 = -pkin(3) * t739 + qJDD(2) * pkin(8) + t687 * t709 + t635;
t726 = t736 * g(1) - t738 * g(2);
t747 = -qJDD(1) * pkin(1) - t726;
t674 = -t722 * pkin(2) + qJDD(3) + t723 * t765 + (-qJ(3) * t732 - pkin(7)) * t740 + t747;
t633 = (-qJD(2) * t709 - t694) * pkin(8) + (qJD(2) * t710 - t693) * pkin(3) + t674;
t626 = -t734 * t631 + t633 * t775;
t677 = -mrSges(5,2) * t708 - mrSges(5,3) * t696;
t692 = qJDD(4) - t693;
t669 = pkin(4) * t696 - qJ(5) * t697;
t707 = t708 ^ 2;
t624 = -t692 * pkin(4) - t707 * qJ(5) + t697 * t669 + qJDD(5) - t626;
t675 = -mrSges(6,2) * t696 + mrSges(6,3) * t708;
t617 = qJD(6) * t777 + (-t660 - t771) * qJ(6) + (t696 * t697 - t692) * pkin(5) + t624;
t671 = -mrSges(7,1) * t696 + mrSges(7,2) * t697;
t749 = -m(7) * t617 + t660 * mrSges(7,3) + t697 * t671;
t744 = -m(6) * t624 + t692 * mrSges(6,1) + t708 * t675 + t749;
t670 = mrSges(6,1) * t696 - mrSges(6,3) * t697;
t766 = -mrSges(5,1) * t696 - mrSges(5,2) * t697 - t670;
t613 = m(5) * t626 + (t676 + t677) * t708 + t766 * t697 + (mrSges(5,1) + mrSges(7,1)) * t692 + t773 * t660 + t744;
t627 = t775 * t631 + t734 * t633;
t659 = qJD(4) * t697 - qJDD(2) * t775 + t694 * t734;
t679 = -mrSges(7,1) * t708 - mrSges(7,3) * t697;
t680 = mrSges(5,1) * t708 - mrSges(5,3) * t697;
t623 = -pkin(4) * t707 + t692 * qJ(5) - t696 * t669 + t708 * t776 + t627;
t681 = -mrSges(6,1) * t708 + mrSges(6,2) * t697;
t678 = -pkin(5) * t708 - qJ(6) * t697;
t695 = t696 ^ 2;
t619 = -pkin(5) * t695 + qJ(6) * t659 + 0.2e1 * qJD(6) * t696 + t678 * t708 + t623;
t758 = m(7) * t619 + t659 * mrSges(7,3) + t696 * t671;
t746 = m(6) * t623 + t692 * mrSges(6,3) + t708 * t681 + t758;
t614 = m(5) * t627 + (t679 - t680) * t708 + t766 * t696 + (-mrSges(5,2) + mrSges(7,2)) * t692 + t773 * t659 + t746;
t751 = -t613 * t734 + t775 * t614;
t606 = m(4) * t635 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t693 - qJD(2) * t701 + t686 * t709 + t751;
t763 = qJD(3) * t710;
t704 = -0.2e1 * t763;
t767 = t772 * t667 - t733 * t668;
t634 = t704 + t767;
t700 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t709;
t745 = qJDD(2) * pkin(3) + t739 * pkin(8) - t710 * t687 + t767;
t630 = 0.2e1 * t763 - t745;
t625 = qJD(5) * t777 + t778 + (t697 * t708 + t659) * pkin(4) + t630;
t621 = -t695 * qJ(6) + qJDD(6) + t704 + (-pkin(4) - pkin(5)) * t659 - t778 + (-pkin(4) * t708 + t678 + t776) * t697 + t745;
t748 = -m(7) * t621 + t659 * mrSges(7,1) - t660 * mrSges(7,2) + t696 * t676 - t697 * t679;
t615 = m(6) * t625 + mrSges(6,1) * t659 - t660 * mrSges(6,3) + t675 * t696 - t697 * t681 + t748;
t741 = -m(5) * t630 - t659 * mrSges(5,1) - mrSges(5,2) * t660 - t696 * t677 - t680 * t697 - t615;
t609 = m(4) * t634 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t694 + qJD(2) * t700 - t686 * t710 + t741;
t600 = t733 * t606 + t772 * t609;
t698 = -t737 * g(3) - t769;
t720 = (-mrSges(3,1) * t737 + mrSges(3,2) * t735) * qJD(1);
t725 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t764;
t598 = m(3) * t698 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t721 + qJD(2) * t725 - t720 * t765 + t600;
t724 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t765;
t752 = t772 * t606 - t609 * t733;
t599 = m(3) * t699 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t722 - qJD(2) * t724 + t720 * t764 + t752;
t753 = -t598 * t735 + t737 * t599;
t591 = m(2) * t727 - mrSges(2,1) * t740 - qJDD(1) * mrSges(2,2) + t753;
t715 = -t740 * pkin(7) + t747;
t607 = t775 * t613 + t734 * t614;
t743 = m(4) * t674 - t693 * mrSges(4,1) + mrSges(4,2) * t694 - t709 * t700 + t701 * t710 + t607;
t742 = -m(3) * t715 + t722 * mrSges(3,1) - mrSges(3,2) * t721 - t724 * t765 + t725 * t764 - t743;
t603 = m(2) * t726 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t740 + t742;
t768 = t736 * t591 + t738 * t603;
t592 = t737 * t598 + t735 * t599;
t757 = t759 * t696 - t760 * t697 + t779 * t708;
t756 = t780 * t696 - t761 * t697 - t759 * t708;
t755 = t761 * t696 - t781 * t697 - t760 * t708;
t754 = t738 * t591 - t603 * t736;
t713 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t735 + Ifges(3,4) * t737) * qJD(1);
t712 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t735 + Ifges(3,2) * t737) * qJD(1);
t711 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t735 + Ifges(3,6) * t737) * qJD(1);
t684 = Ifges(4,1) * t710 + Ifges(4,4) * t709 + Ifges(4,5) * qJD(2);
t683 = Ifges(4,4) * t710 + Ifges(4,2) * t709 + Ifges(4,6) * qJD(2);
t682 = Ifges(4,5) * t710 + Ifges(4,6) * t709 + Ifges(4,3) * qJD(2);
t616 = -t692 * mrSges(7,1) - t749 - t770;
t601 = mrSges(5,2) * t630 + mrSges(6,2) * t624 + mrSges(7,2) * t621 - mrSges(5,3) * t626 - mrSges(6,3) * t625 - mrSges(7,3) * t617 - qJ(5) * t615 - qJ(6) * t616 - t761 * t659 + t781 * t660 + t760 * t692 + t757 * t696 + t756 * t708;
t594 = -mrSges(5,1) * t630 + mrSges(5,3) * t627 - mrSges(6,1) * t625 + mrSges(6,2) * t623 + mrSges(7,1) * t621 - mrSges(7,3) * t619 - pkin(5) * t748 - qJ(6) * t758 - pkin(4) * t615 + (-qJ(6) * t679 - t755) * t708 + t757 * t697 + (-mrSges(7,2) * qJ(6) + t759) * t692 + t761 * t660 - t780 * t659;
t593 = Ifges(4,6) * qJDD(2) - qJ(5) * (t708 * t679 + t746) - pkin(4) * (t744 + t770) - t710 * t682 + pkin(5) * t616 + mrSges(7,1) * t617 - mrSges(7,2) * t619 - mrSges(6,3) * t623 + mrSges(6,1) * t624 - mrSges(5,1) * t626 + mrSges(5,2) * t627 + mrSges(4,3) * t635 - mrSges(4,1) * t674 + qJD(2) * t684 + Ifges(4,2) * t693 + Ifges(4,4) * t694 - pkin(3) * t607 + (pkin(4) * t670 + t756) * t697 + (qJ(5) * t670 + t755) * t696 + (-mrSges(7,1) * pkin(4) - mrSges(7,2) * qJ(5) + t779) * t692 + (mrSges(6,2) * pkin(4) - t760) * t660 + (mrSges(6,2) * qJ(5) + t759) * t659;
t588 = mrSges(4,2) * t674 - mrSges(4,3) * t634 + Ifges(4,1) * t694 + Ifges(4,4) * t693 + Ifges(4,5) * qJDD(2) - pkin(8) * t607 - qJD(2) * t683 - t734 * t594 + t601 * t775 + t709 * t682;
t587 = mrSges(3,2) * t715 - mrSges(3,3) * t698 + Ifges(3,1) * t721 + Ifges(3,4) * t722 + Ifges(3,5) * qJDD(2) - qJ(3) * t600 - qJD(2) * t712 + t588 * t772 - t733 * t593 + t711 * t764;
t586 = -pkin(1) * t592 + mrSges(2,3) * t727 - pkin(2) * t600 - Ifges(3,5) * t721 - Ifges(3,6) * t722 - mrSges(3,1) * t698 + mrSges(3,2) * t699 - pkin(8) * t751 - mrSges(4,1) * t634 + mrSges(4,2) * t635 - t734 * t601 - t775 * t594 - pkin(3) * t741 - Ifges(4,5) * t694 - Ifges(4,6) * t693 + mrSges(2,1) * g(3) + t740 * Ifges(2,5) - t710 * t683 + t709 * t684 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t712 * t735 + t713 * t737) * qJD(1);
t585 = -mrSges(3,1) * t715 + mrSges(3,3) * t699 + Ifges(3,4) * t721 + Ifges(3,2) * t722 + Ifges(3,6) * qJDD(2) - pkin(2) * t743 + qJ(3) * t752 + qJD(2) * t713 + t733 * t588 + t593 * t772 - t711 * t765;
t584 = -mrSges(2,2) * g(3) - mrSges(2,3) * t726 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t740 - pkin(7) * t592 - t585 * t735 + t587 * t737;
t1 = [-m(1) * g(1) + t754; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t592; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t768 + t738 * t584 - t736 * t586; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t754 + t736 * t584 + t738 * t586; -mrSges(1,1) * g(2) + mrSges(2,1) * t726 + mrSges(1,2) * g(1) - mrSges(2,2) * t727 + Ifges(2,3) * qJDD(1) + pkin(1) * t742 + pkin(7) * t753 + t737 * t585 + t735 * t587;];
tauB  = t1;
