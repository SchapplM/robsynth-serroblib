% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:47:39
% EndTime: 2019-05-05 17:47:50
% DurationCPUTime: 10.29s
% Computational Cost: add. (143474->341), mult. (352188->418), div. (0->0), fcn. (262982->10), ass. (0->141)
t759 = Ifges(6,1) + Ifges(7,1);
t752 = Ifges(6,4) - Ifges(7,5);
t751 = Ifges(7,4) + Ifges(6,5);
t758 = Ifges(6,2) + Ifges(7,3);
t757 = -Ifges(7,2) - Ifges(6,3);
t750 = Ifges(6,6) - Ifges(7,6);
t717 = qJD(1) ^ 2;
t756 = cos(qJ(3));
t755 = cos(qJ(5));
t711 = cos(pkin(9));
t754 = pkin(2) * t711;
t753 = -mrSges(6,3) - mrSges(7,2);
t709 = sin(pkin(9));
t749 = mrSges(3,2) * t709;
t707 = t711 ^ 2;
t748 = t707 * t717;
t714 = sin(qJ(1));
t715 = cos(qJ(1));
t697 = -g(1) * t715 - g(2) * t714;
t693 = -pkin(1) * t717 + qJDD(1) * qJ(2) + t697;
t739 = qJD(1) * qJD(2);
t734 = -t711 * g(3) - 0.2e1 * t709 * t739;
t656 = (-pkin(7) * qJDD(1) + t717 * t754 - t693) * t709 + t734;
t678 = -g(3) * t709 + (t693 + 0.2e1 * t739) * t711;
t737 = qJDD(1) * t711;
t663 = -pkin(2) * t748 + pkin(7) * t737 + t678;
t713 = sin(qJ(3));
t638 = t713 * t656 + t756 * t663;
t735 = t711 * t756;
t742 = qJD(1) * t709;
t691 = -qJD(1) * t735 + t713 * t742;
t722 = t756 * t709 + t711 * t713;
t692 = t722 * qJD(1);
t671 = mrSges(4,1) * t691 + mrSges(4,2) * t692;
t738 = qJDD(1) * t709;
t741 = qJD(3) * t692;
t675 = -qJDD(1) * t735 + t713 * t738 + t741;
t686 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t692;
t670 = pkin(3) * t691 - qJ(4) * t692;
t716 = qJD(3) ^ 2;
t619 = -pkin(3) * t716 + qJDD(3) * qJ(4) - t670 * t691 + t638;
t706 = t709 ^ 2;
t696 = t714 * g(1) - t715 * g(2);
t727 = qJDD(2) - t696;
t674 = (-pkin(1) - t754) * qJDD(1) + (-qJ(2) + (-t706 - t707) * pkin(7)) * t717 + t727;
t740 = t691 * qJD(3);
t676 = t722 * qJDD(1) - t740;
t622 = (-t676 + t740) * qJ(4) + (t675 + t741) * pkin(3) + t674;
t708 = sin(pkin(10));
t710 = cos(pkin(10));
t684 = qJD(3) * t708 + t692 * t710;
t609 = -0.2e1 * qJD(4) * t684 - t708 * t619 + t710 * t622;
t662 = qJDD(3) * t708 + t676 * t710;
t683 = qJD(3) * t710 - t692 * t708;
t606 = (t683 * t691 - t662) * pkin(8) + (t683 * t684 + t675) * pkin(4) + t609;
t610 = 0.2e1 * qJD(4) * t683 + t710 * t619 + t708 * t622;
t660 = pkin(4) * t691 - pkin(8) * t684;
t661 = qJDD(3) * t710 - t676 * t708;
t682 = t683 ^ 2;
t608 = -pkin(4) * t682 + pkin(8) * t661 - t660 * t691 + t610;
t712 = sin(qJ(5));
t602 = t712 * t606 + t755 * t608;
t649 = t712 * t683 + t755 * t684;
t617 = t649 * qJD(5) - t755 * t661 + t712 * t662;
t689 = qJD(5) + t691;
t641 = mrSges(6,1) * t689 - mrSges(6,3) * t649;
t648 = -t755 * t683 + t712 * t684;
t673 = qJDD(5) + t675;
t632 = pkin(5) * t648 - qJ(6) * t649;
t688 = t689 ^ 2;
t599 = -pkin(5) * t688 + qJ(6) * t673 + 0.2e1 * qJD(6) * t689 - t632 * t648 + t602;
t642 = -mrSges(7,1) * t689 + mrSges(7,2) * t649;
t736 = m(7) * t599 + t673 * mrSges(7,3) + t689 * t642;
t633 = mrSges(7,1) * t648 - mrSges(7,3) * t649;
t743 = -mrSges(6,1) * t648 - mrSges(6,2) * t649 - t633;
t592 = m(6) * t602 - t673 * mrSges(6,2) + t753 * t617 - t689 * t641 + t743 * t648 + t736;
t601 = t755 * t606 - t712 * t608;
t618 = -t648 * qJD(5) + t712 * t661 + t755 * t662;
t640 = -mrSges(6,2) * t689 - mrSges(6,3) * t648;
t600 = -t673 * pkin(5) - t688 * qJ(6) + t649 * t632 + qJDD(6) - t601;
t639 = -mrSges(7,2) * t648 + mrSges(7,3) * t689;
t728 = -m(7) * t600 + t673 * mrSges(7,1) + t689 * t639;
t594 = m(6) * t601 + t673 * mrSges(6,1) + t753 * t618 + t689 * t640 + t743 * t649 + t728;
t587 = t712 * t592 + t755 * t594;
t650 = -mrSges(5,1) * t683 + mrSges(5,2) * t684;
t658 = -mrSges(5,2) * t691 + mrSges(5,3) * t683;
t585 = m(5) * t609 + mrSges(5,1) * t675 - mrSges(5,3) * t662 - t650 * t684 + t658 * t691 + t587;
t659 = mrSges(5,1) * t691 - mrSges(5,3) * t684;
t729 = t755 * t592 - t594 * t712;
t586 = m(5) * t610 - mrSges(5,2) * t675 + mrSges(5,3) * t661 + t650 * t683 - t659 * t691 + t729;
t730 = -t585 * t708 + t710 * t586;
t580 = m(4) * t638 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t675 - qJD(3) * t686 - t671 * t691 + t730;
t637 = t756 * t656 - t713 * t663;
t685 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t691;
t616 = -qJDD(3) * pkin(3) - t716 * qJ(4) + t692 * t670 + qJDD(4) - t637;
t611 = -t661 * pkin(4) - t682 * pkin(8) + t684 * t660 + t616;
t604 = -0.2e1 * qJD(6) * t649 + (t648 * t689 - t618) * qJ(6) + (t649 * t689 + t617) * pkin(5) + t611;
t597 = m(7) * t604 + t617 * mrSges(7,1) - t618 * mrSges(7,3) + t648 * t639 - t649 * t642;
t720 = m(6) * t611 + t617 * mrSges(6,1) + mrSges(6,2) * t618 + t648 * t640 + t641 * t649 + t597;
t718 = -m(5) * t616 + t661 * mrSges(5,1) - mrSges(5,2) * t662 + t683 * t658 - t659 * t684 - t720;
t596 = m(4) * t637 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t676 + qJD(3) * t685 - t671 * t692 + t718;
t575 = t713 * t580 + t756 * t596;
t677 = -t709 * t693 + t734;
t723 = mrSges(3,3) * qJDD(1) + t717 * (-mrSges(3,1) * t711 + t749);
t573 = m(3) * t677 - t723 * t709 + t575;
t731 = t756 * t580 - t713 * t596;
t574 = m(3) * t678 + t723 * t711 + t731;
t732 = -t573 * t709 + t711 * t574;
t565 = m(2) * t697 - mrSges(2,1) * t717 - qJDD(1) * mrSges(2,2) + t732;
t690 = -qJDD(1) * pkin(1) - t717 * qJ(2) + t727;
t581 = t710 * t585 + t708 * t586;
t721 = m(4) * t674 + t675 * mrSges(4,1) + t676 * mrSges(4,2) + t691 * t685 + t692 * t686 + t581;
t719 = -m(3) * t690 + mrSges(3,1) * t737 - t721 + (t706 * t717 + t748) * mrSges(3,3);
t577 = t719 - t717 * mrSges(2,2) + m(2) * t696 + (mrSges(2,1) - t749) * qJDD(1);
t747 = t714 * t565 + t715 * t577;
t566 = t711 * t573 + t709 * t574;
t746 = t648 * t758 - t649 * t752 - t689 * t750;
t745 = t648 * t750 - t649 * t751 + t689 * t757;
t744 = -t752 * t648 + t649 * t759 + t751 * t689;
t733 = t715 * t565 - t577 * t714;
t726 = Ifges(3,1) * t709 + Ifges(3,4) * t711;
t725 = Ifges(3,4) * t709 + Ifges(3,2) * t711;
t724 = Ifges(3,5) * t709 + Ifges(3,6) * t711;
t695 = t724 * qJD(1);
t666 = Ifges(4,1) * t692 - Ifges(4,4) * t691 + Ifges(4,5) * qJD(3);
t665 = Ifges(4,4) * t692 - Ifges(4,2) * t691 + Ifges(4,6) * qJD(3);
t664 = Ifges(4,5) * t692 - Ifges(4,6) * t691 + Ifges(4,3) * qJD(3);
t645 = Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t691;
t644 = Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t691;
t643 = Ifges(5,5) * t684 + Ifges(5,6) * t683 + Ifges(5,3) * t691;
t589 = mrSges(6,2) * t611 + mrSges(7,2) * t600 - mrSges(6,3) * t601 - mrSges(7,3) * t604 - qJ(6) * t597 - t752 * t617 + t618 * t759 + t745 * t648 + t751 * t673 + t746 * t689;
t588 = -mrSges(6,1) * t611 - mrSges(7,1) * t604 + mrSges(7,2) * t599 + mrSges(6,3) * t602 - pkin(5) * t597 - t617 * t758 + t752 * t618 + t745 * t649 + t750 * t673 + t744 * t689;
t569 = mrSges(5,2) * t616 - mrSges(5,3) * t609 + Ifges(5,1) * t662 + Ifges(5,4) * t661 + Ifges(5,5) * t675 - pkin(8) * t587 - t712 * t588 + t755 * t589 + t683 * t643 - t691 * t644;
t568 = -mrSges(5,1) * t616 + mrSges(5,3) * t610 + Ifges(5,4) * t662 + Ifges(5,2) * t661 + Ifges(5,6) * t675 - pkin(4) * t720 + pkin(8) * t729 + t755 * t588 + t712 * t589 - t684 * t643 + t691 * t645;
t567 = (mrSges(7,2) * qJ(6) + t750) * t617 + (mrSges(7,2) * pkin(5) - t751) * t618 + (qJ(6) * t633 - t744) * t648 + (pkin(5) * t633 + t746) * t649 - qJ(6) * t736 - pkin(5) * t728 + Ifges(4,6) * qJDD(3) - t692 * t664 + t683 * t645 - t684 * t644 + Ifges(4,4) * t676 - mrSges(4,1) * t674 - Ifges(5,6) * t661 - Ifges(5,5) * t662 + qJD(3) * t666 + mrSges(4,3) * t638 - mrSges(5,1) * t609 + mrSges(5,2) * t610 - mrSges(6,1) * t601 + mrSges(6,2) * t602 - mrSges(7,3) * t599 + mrSges(7,1) * t600 - pkin(4) * t587 + (-Ifges(5,3) - Ifges(4,2)) * t675 - pkin(3) * t581 + t757 * t673;
t562 = mrSges(4,2) * t674 - mrSges(4,3) * t637 + Ifges(4,1) * t676 - Ifges(4,4) * t675 + Ifges(4,5) * qJDD(3) - qJ(4) * t581 - qJD(3) * t665 - t568 * t708 + t569 * t710 - t664 * t691;
t561 = t711 * qJD(1) * t695 + mrSges(3,2) * t690 - mrSges(3,3) * t677 - pkin(7) * t575 + t726 * qJDD(1) + t756 * t562 - t713 * t567;
t560 = mrSges(2,1) * g(3) - pkin(1) * t566 + mrSges(2,3) * t697 - pkin(2) * t575 + mrSges(3,2) * t678 - mrSges(3,1) * t677 - t710 * t568 - pkin(3) * t718 - qJ(4) * t730 - t708 * t569 - Ifges(4,5) * t676 + Ifges(4,6) * t675 - Ifges(4,3) * qJDD(3) - t692 * t665 - t691 * t666 - mrSges(4,1) * t637 + mrSges(4,2) * t638 + (Ifges(2,6) - t724) * qJDD(1) + (-t709 * t725 + t711 * t726 + Ifges(2,5)) * t717;
t559 = -mrSges(3,1) * t690 + mrSges(3,3) * t678 - pkin(2) * t721 + pkin(7) * t731 + t725 * qJDD(1) + t713 * t562 + t756 * t567 - t695 * t742;
t558 = -mrSges(2,2) * g(3) - mrSges(2,3) * t696 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t717 - qJ(2) * t566 - t559 * t709 + t561 * t711;
t1 = [-m(1) * g(1) + t733; -m(1) * g(2) + t747; (-m(1) - m(2)) * g(3) + t566; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t747 + t715 * t558 - t714 * t560; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t733 + t714 * t558 + t715 * t560; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t696 - mrSges(2,2) * t697 + t709 * t561 + t711 * t559 + pkin(1) * (-mrSges(3,2) * t738 + t719) + qJ(2) * t732;];
tauB  = t1;
