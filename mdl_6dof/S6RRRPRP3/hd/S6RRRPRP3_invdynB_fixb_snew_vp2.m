% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:37:59
% EndTime: 2019-05-07 07:38:15
% DurationCPUTime: 12.69s
% Computational Cost: add. (192451->363), mult. (397850->447), div. (0->0), fcn. (284856->10), ass. (0->139)
t754 = Ifges(6,1) + Ifges(7,1);
t746 = Ifges(6,4) - Ifges(7,5);
t753 = -Ifges(6,5) - Ifges(7,4);
t752 = Ifges(6,2) + Ifges(7,3);
t744 = Ifges(6,6) - Ifges(7,6);
t751 = -Ifges(6,3) - Ifges(7,2);
t750 = cos(qJ(3));
t749 = cos(qJ(5));
t722 = qJD(1) ^ 2;
t748 = pkin(2) * t722;
t747 = -mrSges(6,3) - mrSges(7,2);
t719 = sin(qJ(1));
t721 = cos(qJ(1));
t707 = -g(1) * t721 - g(2) * t719;
t695 = -pkin(1) * t722 + qJDD(1) * pkin(7) + t707;
t718 = sin(qJ(2));
t743 = t718 * t695;
t720 = cos(qJ(2));
t735 = qJD(1) * qJD(2);
t701 = qJDD(1) * t718 + t720 * t735;
t657 = qJDD(2) * pkin(2) - t701 * pkin(8) - t743 + (pkin(8) * t735 + t718 * t748 - g(3)) * t720;
t683 = -g(3) * t718 + t720 * t695;
t702 = qJDD(1) * t720 - t718 * t735;
t737 = qJD(1) * t718;
t705 = qJD(2) * pkin(2) - pkin(8) * t737;
t713 = t720 ^ 2;
t658 = pkin(8) * t702 - qJD(2) * t705 - t713 * t748 + t683;
t717 = sin(qJ(3));
t635 = t717 * t657 + t750 * t658;
t693 = (t717 * t720 + t718 * t750) * qJD(1);
t664 = qJD(3) * t693 + t701 * t717 - t750 * t702;
t736 = qJD(1) * t720;
t692 = t717 * t737 - t750 * t736;
t675 = mrSges(4,1) * t692 + mrSges(4,2) * t693;
t712 = qJD(2) + qJD(3);
t685 = mrSges(4,1) * t712 - mrSges(4,3) * t693;
t711 = qJDD(2) + qJDD(3);
t665 = -t692 * qJD(3) + t701 * t750 + t717 * t702;
t706 = t719 * g(1) - t721 * g(2);
t727 = -qJDD(1) * pkin(1) - t706;
t666 = -t702 * pkin(2) + t705 * t737 + (-pkin(8) * t713 - pkin(7)) * t722 + t727;
t618 = (t692 * t712 - t665) * qJ(4) + (t693 * t712 + t664) * pkin(3) + t666;
t674 = pkin(3) * t692 - qJ(4) * t693;
t710 = t712 ^ 2;
t621 = -pkin(3) * t710 + qJ(4) * t711 - t674 * t692 + t635;
t714 = sin(pkin(10));
t715 = cos(pkin(10));
t681 = t693 * t715 + t712 * t714;
t608 = -0.2e1 * qJD(4) * t681 + t715 * t618 - t714 * t621;
t652 = t665 * t715 + t711 * t714;
t680 = -t693 * t714 + t712 * t715;
t605 = (t680 * t692 - t652) * pkin(9) + (t680 * t681 + t664) * pkin(4) + t608;
t609 = 0.2e1 * qJD(4) * t680 + t714 * t618 + t715 * t621;
t651 = -t665 * t714 + t711 * t715;
t669 = pkin(4) * t692 - pkin(9) * t681;
t679 = t680 ^ 2;
t607 = -pkin(4) * t679 + pkin(9) * t651 - t669 * t692 + t609;
t716 = sin(qJ(5));
t601 = t716 * t605 + t749 * t607;
t650 = t716 * t680 + t681 * t749;
t616 = qJD(5) * t650 - t651 * t749 + t652 * t716;
t688 = qJD(5) + t692;
t640 = mrSges(6,1) * t688 - mrSges(6,3) * t650;
t649 = -t680 * t749 + t681 * t716;
t663 = qJDD(5) + t664;
t631 = pkin(5) * t649 - qJ(6) * t650;
t687 = t688 ^ 2;
t598 = -pkin(5) * t687 + qJ(6) * t663 + 0.2e1 * qJD(6) * t688 - t631 * t649 + t601;
t641 = -mrSges(7,1) * t688 + mrSges(7,2) * t650;
t734 = m(7) * t598 + t663 * mrSges(7,3) + t688 * t641;
t632 = mrSges(7,1) * t649 - mrSges(7,3) * t650;
t738 = -mrSges(6,1) * t649 - mrSges(6,2) * t650 - t632;
t591 = m(6) * t601 - t663 * mrSges(6,2) + t616 * t747 - t688 * t640 + t649 * t738 + t734;
t600 = t605 * t749 - t716 * t607;
t617 = -t649 * qJD(5) + t716 * t651 + t652 * t749;
t639 = -mrSges(6,2) * t688 - mrSges(6,3) * t649;
t599 = -t663 * pkin(5) - t687 * qJ(6) + t650 * t631 + qJDD(6) - t600;
t638 = -mrSges(7,2) * t649 + mrSges(7,3) * t688;
t728 = -m(7) * t599 + t663 * mrSges(7,1) + t688 * t638;
t593 = m(6) * t600 + t663 * mrSges(6,1) + t617 * t747 + t688 * t639 + t650 * t738 + t728;
t586 = t716 * t591 + t749 * t593;
t654 = -mrSges(5,1) * t680 + mrSges(5,2) * t681;
t667 = -mrSges(5,2) * t692 + mrSges(5,3) * t680;
t584 = m(5) * t608 + mrSges(5,1) * t664 - mrSges(5,3) * t652 - t654 * t681 + t667 * t692 + t586;
t668 = mrSges(5,1) * t692 - mrSges(5,3) * t681;
t729 = t749 * t591 - t593 * t716;
t585 = m(5) * t609 - mrSges(5,2) * t664 + mrSges(5,3) * t651 + t654 * t680 - t668 * t692 + t729;
t730 = -t584 * t714 + t715 * t585;
t579 = m(4) * t635 - mrSges(4,2) * t711 - mrSges(4,3) * t664 - t675 * t692 - t685 * t712 + t730;
t634 = t657 * t750 - t717 * t658;
t684 = -mrSges(4,2) * t712 - mrSges(4,3) * t692;
t620 = -t711 * pkin(3) - t710 * qJ(4) + t693 * t674 + qJDD(4) - t634;
t610 = -t651 * pkin(4) - t679 * pkin(9) + t681 * t669 + t620;
t603 = -0.2e1 * qJD(6) * t650 + (t649 * t688 - t617) * qJ(6) + (t650 * t688 + t616) * pkin(5) + t610;
t596 = m(7) * t603 + t616 * mrSges(7,1) - t617 * mrSges(7,3) + t649 * t638 - t650 * t641;
t725 = m(6) * t610 + t616 * mrSges(6,1) + mrSges(6,2) * t617 + t649 * t639 + t640 * t650 + t596;
t723 = -m(5) * t620 + t651 * mrSges(5,1) - mrSges(5,2) * t652 + t680 * t667 - t668 * t681 - t725;
t595 = m(4) * t634 + mrSges(4,1) * t711 - mrSges(4,3) * t665 - t675 * t693 + t684 * t712 + t723;
t574 = t717 * t579 + t750 * t595;
t682 = -t720 * g(3) - t743;
t700 = (-mrSges(3,1) * t720 + mrSges(3,2) * t718) * qJD(1);
t704 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t736;
t572 = m(3) * t682 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t701 + qJD(2) * t704 - t700 * t737 + t574;
t703 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t737;
t731 = t750 * t579 - t595 * t717;
t573 = m(3) * t683 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t702 - qJD(2) * t703 + t700 * t736 + t731;
t732 = -t572 * t718 + t720 * t573;
t564 = m(2) * t707 - mrSges(2,1) * t722 - qJDD(1) * mrSges(2,2) + t732;
t694 = -pkin(7) * t722 + t727;
t580 = t715 * t584 + t714 * t585;
t726 = m(4) * t666 + t664 * mrSges(4,1) + mrSges(4,2) * t665 + t692 * t684 + t685 * t693 + t580;
t724 = -m(3) * t694 + t702 * mrSges(3,1) - mrSges(3,2) * t701 - t703 * t737 + t704 * t736 - t726;
t576 = m(2) * t706 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t722 + t724;
t742 = t719 * t564 + t721 * t576;
t565 = t720 * t572 + t718 * t573;
t741 = t752 * t649 - t746 * t650 - t744 * t688;
t740 = t744 * t649 + t753 * t650 + t751 * t688;
t739 = -t746 * t649 + t754 * t650 - t753 * t688;
t733 = t721 * t564 - t576 * t719;
t691 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t718 + Ifges(3,4) * t720) * qJD(1);
t690 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t718 + Ifges(3,2) * t720) * qJD(1);
t689 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t718 + Ifges(3,6) * t720) * qJD(1);
t672 = Ifges(4,1) * t693 - Ifges(4,4) * t692 + Ifges(4,5) * t712;
t671 = Ifges(4,4) * t693 - Ifges(4,2) * t692 + Ifges(4,6) * t712;
t670 = Ifges(4,5) * t693 - Ifges(4,6) * t692 + Ifges(4,3) * t712;
t644 = Ifges(5,1) * t681 + Ifges(5,4) * t680 + Ifges(5,5) * t692;
t643 = Ifges(5,4) * t681 + Ifges(5,2) * t680 + Ifges(5,6) * t692;
t642 = Ifges(5,5) * t681 + Ifges(5,6) * t680 + Ifges(5,3) * t692;
t588 = mrSges(6,2) * t610 + mrSges(7,2) * t599 - mrSges(6,3) * t600 - mrSges(7,3) * t603 - qJ(6) * t596 - t746 * t616 + t754 * t617 + t740 * t649 - t663 * t753 + t741 * t688;
t587 = -mrSges(6,1) * t610 - mrSges(7,1) * t603 + mrSges(7,2) * t598 + mrSges(6,3) * t601 - pkin(5) * t596 - t752 * t616 + t746 * t617 + t740 * t650 + t744 * t663 + t739 * t688;
t568 = mrSges(5,2) * t620 - mrSges(5,3) * t608 + Ifges(5,1) * t652 + Ifges(5,4) * t651 + Ifges(5,5) * t664 - pkin(9) * t586 - t716 * t587 + t588 * t749 + t680 * t642 - t692 * t643;
t567 = -mrSges(5,1) * t620 + mrSges(5,3) * t609 + Ifges(5,4) * t652 + Ifges(5,2) * t651 + Ifges(5,6) * t664 - pkin(4) * t725 + pkin(9) * t729 + t587 * t749 + t716 * t588 - t681 * t642 + t692 * t644;
t566 = (mrSges(7,2) * qJ(6) + t744) * t616 + (qJ(6) * t632 - t739) * t649 + (pkin(5) * t632 + t741) * t650 - qJ(6) * t734 - pkin(5) * t728 + (mrSges(7,2) * pkin(5) + t753) * t617 + t751 * t663 + t680 * t644 - t681 * t643 + mrSges(7,1) * t599 - mrSges(6,1) * t600 + mrSges(4,3) * t635 - pkin(4) * t586 - mrSges(5,1) * t608 + mrSges(5,2) * t609 - t693 * t670 + Ifges(4,4) * t665 - mrSges(4,1) * t666 + Ifges(4,6) * t711 + t712 * t672 - mrSges(7,3) * t598 - Ifges(5,6) * t651 - Ifges(5,5) * t652 - pkin(3) * t580 + mrSges(6,2) * t601 + (-Ifges(4,2) - Ifges(5,3)) * t664;
t561 = mrSges(4,2) * t666 - mrSges(4,3) * t634 + Ifges(4,1) * t665 - Ifges(4,4) * t664 + Ifges(4,5) * t711 - qJ(4) * t580 - t714 * t567 + t568 * t715 - t670 * t692 - t712 * t671;
t560 = mrSges(3,2) * t694 - mrSges(3,3) * t682 + Ifges(3,1) * t701 + Ifges(3,4) * t702 + Ifges(3,5) * qJDD(2) - pkin(8) * t574 - qJD(2) * t690 + t561 * t750 - t717 * t566 + t689 * t736;
t559 = -pkin(1) * t565 + mrSges(2,3) * t707 - pkin(2) * t574 - mrSges(3,1) * t682 + mrSges(3,2) * t683 - Ifges(3,5) * t701 - Ifges(3,6) * t702 - Ifges(3,3) * qJDD(2) - t714 * t568 - t715 * t567 - pkin(3) * t723 - qJ(4) * t730 - Ifges(4,5) * t665 + Ifges(4,6) * t664 - Ifges(4,3) * t711 - mrSges(4,1) * t634 + mrSges(4,2) * t635 + mrSges(2,1) * g(3) - t693 * t671 - t692 * t672 + t722 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-t690 * t718 + t691 * t720) * qJD(1);
t558 = -mrSges(3,1) * t694 + mrSges(3,3) * t683 + Ifges(3,4) * t701 + Ifges(3,2) * t702 + Ifges(3,6) * qJDD(2) - pkin(2) * t726 + pkin(8) * t731 + qJD(2) * t691 + t717 * t561 + t566 * t750 - t689 * t737;
t557 = -mrSges(2,2) * g(3) - mrSges(2,3) * t706 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t722 - pkin(7) * t565 - t558 * t718 + t560 * t720;
t1 = [-m(1) * g(1) + t733; -m(1) * g(2) + t742; (-m(1) - m(2)) * g(3) + t565; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t742 + t721 * t557 - t719 * t559; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t733 + t719 * t557 + t721 * t559; -mrSges(1,1) * g(2) + mrSges(2,1) * t706 + mrSges(1,2) * g(1) - mrSges(2,2) * t707 + Ifges(2,3) * qJDD(1) + pkin(1) * t724 + pkin(7) * t732 + t720 * t558 + t718 * t560;];
tauB  = t1;
