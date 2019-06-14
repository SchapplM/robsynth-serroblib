% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:48:12
% EndTime: 2019-05-08 04:48:35
% DurationCPUTime: 13.14s
% Computational Cost: add. (207729->364), mult. (413290->446), div. (0->0), fcn. (299208->10), ass. (0->141)
t768 = Ifges(6,1) + Ifges(7,1);
t762 = Ifges(6,4) - Ifges(7,5);
t761 = Ifges(7,4) + Ifges(6,5);
t767 = Ifges(6,2) + Ifges(7,3);
t760 = Ifges(6,6) - Ifges(7,6);
t766 = -Ifges(6,3) - Ifges(7,2);
t765 = cos(qJ(5));
t738 = qJD(1) ^ 2;
t764 = pkin(2) * t738;
t763 = -mrSges(6,3) - mrSges(7,2);
t733 = sin(qJ(1));
t737 = cos(qJ(1));
t722 = -g(1) * t737 - g(2) * t733;
t710 = -pkin(1) * t738 + qJDD(1) * pkin(7) + t722;
t732 = sin(qJ(2));
t759 = t710 * t732;
t736 = cos(qJ(2));
t751 = qJD(1) * qJD(2);
t716 = qJDD(1) * t732 + t736 * t751;
t671 = qJDD(2) * pkin(2) - pkin(8) * t716 - t759 + (pkin(8) * t751 + t732 * t764 - g(3)) * t736;
t696 = -g(3) * t732 + t736 * t710;
t717 = qJDD(1) * t736 - t732 * t751;
t753 = qJD(1) * t732;
t720 = qJD(2) * pkin(2) - pkin(8) * t753;
t728 = t736 ^ 2;
t672 = pkin(8) * t717 - qJD(2) * t720 - t728 * t764 + t696;
t731 = sin(qJ(3));
t735 = cos(qJ(3));
t649 = t731 * t671 + t735 * t672;
t708 = (t731 * t736 + t732 * t735) * qJD(1);
t679 = -t708 * qJD(3) - t731 * t716 + t717 * t735;
t752 = qJD(1) * t736;
t707 = -t731 * t753 + t735 * t752;
t689 = -mrSges(4,1) * t707 + mrSges(4,2) * t708;
t727 = qJD(2) + qJD(3);
t698 = mrSges(4,1) * t727 - mrSges(4,3) * t708;
t726 = qJDD(2) + qJDD(3);
t680 = qJD(3) * t707 + t716 * t735 + t717 * t731;
t721 = g(1) * t733 - t737 * g(2);
t743 = -qJDD(1) * pkin(1) - t721;
t681 = -pkin(2) * t717 + t720 * t753 + (-pkin(8) * t728 - pkin(7)) * t738 + t743;
t632 = (-t707 * t727 - t680) * pkin(9) + (t708 * t727 - t679) * pkin(3) + t681;
t690 = -pkin(3) * t707 - pkin(9) * t708;
t725 = t727 ^ 2;
t635 = -pkin(3) * t725 + pkin(9) * t726 + t690 * t707 + t649;
t730 = sin(qJ(4));
t734 = cos(qJ(4));
t622 = t734 * t632 - t635 * t730;
t693 = -t708 * t730 + t727 * t734;
t655 = qJD(4) * t693 + t680 * t734 + t726 * t730;
t678 = qJDD(4) - t679;
t694 = t708 * t734 + t727 * t730;
t703 = qJD(4) - t707;
t619 = (t693 * t703 - t655) * pkin(10) + (t693 * t694 + t678) * pkin(4) + t622;
t623 = t730 * t632 + t734 * t635;
t654 = -qJD(4) * t694 - t680 * t730 + t726 * t734;
t684 = pkin(4) * t703 - pkin(10) * t694;
t692 = t693 ^ 2;
t621 = -pkin(4) * t692 + pkin(10) * t654 - t684 * t703 + t623;
t729 = sin(qJ(5));
t615 = t729 * t619 + t765 * t621;
t666 = t729 * t693 + t765 * t694;
t628 = qJD(5) * t666 - t765 * t654 + t655 * t729;
t701 = qJD(5) + t703;
t658 = mrSges(6,1) * t701 - mrSges(6,3) * t666;
t665 = -t765 * t693 + t694 * t729;
t675 = qJDD(5) + t678;
t645 = pkin(5) * t665 - qJ(6) * t666;
t699 = t701 ^ 2;
t612 = -pkin(5) * t699 + qJ(6) * t675 + 0.2e1 * qJD(6) * t701 - t645 * t665 + t615;
t659 = -mrSges(7,1) * t701 + mrSges(7,2) * t666;
t750 = m(7) * t612 + t675 * mrSges(7,3) + t701 * t659;
t646 = mrSges(7,1) * t665 - mrSges(7,3) * t666;
t754 = -mrSges(6,1) * t665 - mrSges(6,2) * t666 - t646;
t605 = m(6) * t615 - mrSges(6,2) * t675 + t763 * t628 - t658 * t701 + t754 * t665 + t750;
t614 = t765 * t619 - t729 * t621;
t629 = -t665 * qJD(5) + t729 * t654 + t765 * t655;
t657 = -mrSges(6,2) * t701 - mrSges(6,3) * t665;
t613 = -t675 * pkin(5) - t699 * qJ(6) + t666 * t645 + qJDD(6) - t614;
t656 = -mrSges(7,2) * t665 + mrSges(7,3) * t701;
t744 = -m(7) * t613 + t675 * mrSges(7,1) + t701 * t656;
t607 = m(6) * t614 + mrSges(6,1) * t675 + t763 * t629 + t657 * t701 + t754 * t666 + t744;
t602 = t729 * t605 + t765 * t607;
t670 = -mrSges(5,1) * t693 + mrSges(5,2) * t694;
t682 = -mrSges(5,2) * t703 + mrSges(5,3) * t693;
t598 = m(5) * t622 + mrSges(5,1) * t678 - mrSges(5,3) * t655 - t670 * t694 + t682 * t703 + t602;
t683 = mrSges(5,1) * t703 - mrSges(5,3) * t694;
t745 = t765 * t605 - t607 * t729;
t599 = m(5) * t623 - mrSges(5,2) * t678 + mrSges(5,3) * t654 + t670 * t693 - t683 * t703 + t745;
t746 = -t598 * t730 + t734 * t599;
t593 = m(4) * t649 - mrSges(4,2) * t726 + mrSges(4,3) * t679 + t689 * t707 - t698 * t727 + t746;
t648 = t671 * t735 - t731 * t672;
t697 = -mrSges(4,2) * t727 + mrSges(4,3) * t707;
t634 = -pkin(3) * t726 - pkin(9) * t725 + t708 * t690 - t648;
t624 = -pkin(4) * t654 - pkin(10) * t692 + t694 * t684 + t634;
t617 = -0.2e1 * qJD(6) * t666 + (t665 * t701 - t629) * qJ(6) + (t666 * t701 + t628) * pkin(5) + t624;
t610 = m(7) * t617 + t628 * mrSges(7,1) - t629 * mrSges(7,3) + t665 * t656 - t666 * t659;
t741 = m(6) * t624 + t628 * mrSges(6,1) + mrSges(6,2) * t629 + t665 * t657 + t658 * t666 + t610;
t739 = -m(5) * t634 + t654 * mrSges(5,1) - mrSges(5,2) * t655 + t693 * t682 - t683 * t694 - t741;
t609 = m(4) * t648 + mrSges(4,1) * t726 - mrSges(4,3) * t680 - t689 * t708 + t697 * t727 + t739;
t588 = t731 * t593 + t735 * t609;
t695 = -g(3) * t736 - t759;
t715 = (-mrSges(3,1) * t736 + mrSges(3,2) * t732) * qJD(1);
t719 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t752;
t586 = m(3) * t695 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t716 + qJD(2) * t719 - t715 * t753 + t588;
t718 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t753;
t747 = t735 * t593 - t609 * t731;
t587 = m(3) * t696 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t717 - qJD(2) * t718 + t715 * t752 + t747;
t748 = -t586 * t732 + t736 * t587;
t578 = m(2) * t722 - mrSges(2,1) * t738 - qJDD(1) * mrSges(2,2) + t748;
t709 = -pkin(7) * t738 + t743;
t594 = t734 * t598 + t730 * t599;
t742 = m(4) * t681 - t679 * mrSges(4,1) + mrSges(4,2) * t680 - t707 * t697 + t698 * t708 + t594;
t740 = -m(3) * t709 + t717 * mrSges(3,1) - mrSges(3,2) * t716 - t718 * t753 + t719 * t752 - t742;
t590 = m(2) * t721 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t738 + t740;
t758 = t733 * t578 + t737 * t590;
t579 = t736 * t586 + t732 * t587;
t757 = t767 * t665 - t762 * t666 - t760 * t701;
t756 = t760 * t665 - t761 * t666 + t766 * t701;
t755 = -t762 * t665 + t768 * t666 + t761 * t701;
t749 = t737 * t578 - t590 * t733;
t706 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t732 + Ifges(3,4) * t736) * qJD(1);
t705 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t732 + Ifges(3,2) * t736) * qJD(1);
t704 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t732 + Ifges(3,6) * t736) * qJD(1);
t687 = Ifges(4,1) * t708 + Ifges(4,4) * t707 + Ifges(4,5) * t727;
t686 = Ifges(4,4) * t708 + Ifges(4,2) * t707 + Ifges(4,6) * t727;
t685 = Ifges(4,5) * t708 + Ifges(4,6) * t707 + Ifges(4,3) * t727;
t662 = Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t703;
t661 = Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t703;
t660 = Ifges(5,5) * t694 + Ifges(5,6) * t693 + Ifges(5,3) * t703;
t601 = mrSges(6,2) * t624 + mrSges(7,2) * t613 - mrSges(6,3) * t614 - mrSges(7,3) * t617 - qJ(6) * t610 - t762 * t628 + t768 * t629 + t756 * t665 + t761 * t675 + t757 * t701;
t600 = -mrSges(6,1) * t624 - mrSges(7,1) * t617 + mrSges(7,2) * t612 + mrSges(6,3) * t615 - pkin(5) * t610 - t767 * t628 + t762 * t629 + t756 * t666 + t760 * t675 + t755 * t701;
t582 = mrSges(5,2) * t634 - mrSges(5,3) * t622 + Ifges(5,1) * t655 + Ifges(5,4) * t654 + Ifges(5,5) * t678 - pkin(10) * t602 - t729 * t600 + t765 * t601 + t693 * t660 - t703 * t661;
t581 = -mrSges(5,1) * t634 + mrSges(5,3) * t623 + Ifges(5,4) * t655 + Ifges(5,2) * t654 + Ifges(5,6) * t678 - pkin(4) * t741 + pkin(10) * t745 + t765 * t600 + t729 * t601 - t694 * t660 + t703 * t662;
t580 = t766 * t675 + (mrSges(7,2) * qJ(6) + t760) * t628 + (mrSges(7,2) * pkin(5) - t761) * t629 + (qJ(6) * t646 - t755) * t665 + (pkin(5) * t646 + t757) * t666 - qJ(6) * t750 - pkin(5) * t744 - pkin(3) * t594 + t693 * t662 - t694 * t661 - pkin(4) * t602 - mrSges(7,3) * t612 + mrSges(7,1) * t613 + mrSges(4,3) * t649 - Ifges(5,3) * t678 + Ifges(4,2) * t679 + Ifges(4,4) * t680 + mrSges(6,2) * t615 - t708 * t685 + Ifges(4,6) * t726 + t727 * t687 - Ifges(5,6) * t654 - Ifges(5,5) * t655 - mrSges(4,1) * t681 - mrSges(5,1) * t622 + mrSges(5,2) * t623 - mrSges(6,1) * t614;
t575 = mrSges(4,2) * t681 - mrSges(4,3) * t648 + Ifges(4,1) * t680 + Ifges(4,4) * t679 + Ifges(4,5) * t726 - pkin(9) * t594 - t581 * t730 + t582 * t734 + t685 * t707 - t686 * t727;
t574 = mrSges(3,2) * t709 - mrSges(3,3) * t695 + Ifges(3,1) * t716 + Ifges(3,4) * t717 + Ifges(3,5) * qJDD(2) - pkin(8) * t588 - qJD(2) * t705 + t575 * t735 - t580 * t731 + t704 * t752;
t573 = mrSges(2,1) * g(3) + t738 * Ifges(2,5) - t708 * t686 + t707 * t687 + Ifges(2,6) * qJDD(1) - pkin(1) * t579 + mrSges(2,3) * t722 - Ifges(4,5) * t680 - Ifges(4,6) * t679 - Ifges(4,3) * t726 - mrSges(4,1) * t648 + mrSges(4,2) * t649 - t730 * t582 - t734 * t581 - pkin(3) * t739 - pkin(9) * t746 - Ifges(3,5) * t716 - Ifges(3,6) * t717 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t695 + mrSges(3,2) * t696 - pkin(2) * t588 + (-t705 * t732 + t706 * t736) * qJD(1);
t572 = -mrSges(3,1) * t709 + mrSges(3,3) * t696 + Ifges(3,4) * t716 + Ifges(3,2) * t717 + Ifges(3,6) * qJDD(2) - pkin(2) * t742 + pkin(8) * t747 + qJD(2) * t706 + t731 * t575 + t735 * t580 - t704 * t753;
t571 = -mrSges(2,2) * g(3) - mrSges(2,3) * t721 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t738 - pkin(7) * t579 - t572 * t732 + t574 * t736;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t758; (-m(1) - m(2)) * g(3) + t579; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t758 + t737 * t571 - t733 * t573; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t749 + t733 * t571 + t737 * t573; -mrSges(1,1) * g(2) + mrSges(2,1) * t721 + mrSges(1,2) * g(1) - mrSges(2,2) * t722 + Ifges(2,3) * qJDD(1) + pkin(1) * t740 + pkin(7) * t748 + t736 * t572 + t732 * t574;];
tauB  = t1;
