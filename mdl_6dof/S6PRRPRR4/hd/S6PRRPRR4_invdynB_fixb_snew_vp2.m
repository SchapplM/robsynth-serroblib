% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 05:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:06:55
% EndTime: 2019-05-05 05:07:04
% DurationCPUTime: 6.69s
% Computational Cost: add. (104486->322), mult. (206743->402), div. (0->0), fcn. (132990->12), ass. (0->139)
t753 = Ifges(4,1) + Ifges(5,1);
t748 = Ifges(4,4) - Ifges(5,5);
t747 = Ifges(4,5) + Ifges(5,4);
t752 = Ifges(4,2) + Ifges(5,3);
t746 = Ifges(4,6) - Ifges(5,6);
t751 = Ifges(4,3) + Ifges(5,2);
t709 = sin(qJ(5));
t710 = sin(qJ(3));
t713 = cos(qJ(5));
t714 = cos(qJ(3));
t662 = (t709 * t710 + t713 * t714) * qJD(2);
t750 = 2 * qJD(4);
t749 = mrSges(4,3) + mrSges(5,2);
t717 = qJD(2) ^ 2;
t745 = t714 ^ 2 * t717;
t704 = sin(pkin(6));
t711 = sin(qJ(2));
t744 = t704 * t711;
t715 = cos(qJ(2));
t743 = t704 * t715;
t706 = cos(pkin(6));
t742 = t706 * t711;
t741 = t706 * t715;
t703 = sin(pkin(11));
t705 = cos(pkin(11));
t681 = t703 * g(1) - t705 * g(2);
t701 = -g(3) + qJDD(1);
t649 = -t704 * t681 + t706 * t701;
t740 = t714 * t649;
t682 = -t705 * g(1) - t703 * g(2);
t639 = t681 * t742 + t715 * t682 + t701 * t744;
t634 = -t717 * pkin(2) + qJDD(2) * pkin(8) + t639;
t622 = t714 * t634 + t710 * t649;
t676 = (-mrSges(4,1) * t714 + mrSges(4,2) * t710) * qJD(2);
t733 = qJD(2) * qJD(3);
t732 = t710 * t733;
t678 = t714 * qJDD(2) - t732;
t735 = qJD(2) * t710;
t683 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t735;
t674 = (-pkin(3) * t714 - qJ(4) * t710) * qJD(2);
t716 = qJD(3) ^ 2;
t734 = qJD(2) * t714;
t611 = -t716 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t750 + t674 * t734 + t622;
t675 = (-mrSges(5,1) * t714 - mrSges(5,3) * t710) * qJD(2);
t684 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t735;
t689 = -qJD(3) * pkin(4) - pkin(9) * t735;
t606 = -pkin(4) * t745 - t678 * pkin(9) + qJD(3) * t689 + t611;
t731 = t714 * t733;
t677 = t710 * qJDD(2) + t731;
t627 = t710 * t634;
t725 = -t716 * qJ(4) + t674 * t735 + qJDD(4) + t627;
t607 = -t677 * pkin(9) + (-pkin(3) - pkin(4)) * qJDD(3) + (-pkin(4) * t710 * t717 + pkin(9) * t733 - t649) * t714 + t725;
t602 = t713 * t606 + t709 * t607;
t663 = (-t709 * t714 + t710 * t713) * qJD(2);
t629 = -t663 * qJD(5) - t709 * t677 - t713 * t678;
t642 = t662 * mrSges(6,1) + t663 * mrSges(6,2);
t696 = -qJD(3) + qJD(5);
t648 = t696 * mrSges(6,1) - t663 * mrSges(6,3);
t695 = -qJDD(3) + qJDD(5);
t643 = t662 * pkin(5) - t663 * pkin(10);
t694 = t696 ^ 2;
t600 = -t694 * pkin(5) + t695 * pkin(10) - t662 * t643 + t602;
t638 = t681 * t741 - t711 * t682 + t701 * t743;
t633 = -qJDD(2) * pkin(2) - t717 * pkin(8) - t638;
t721 = -t678 * pkin(3) + t633 + (-t677 - t731) * qJ(4);
t609 = -pkin(3) * t732 + t678 * pkin(4) - pkin(9) * t745 - t721 + (t689 + t750) * t735;
t630 = -t662 * qJD(5) + t713 * t677 - t709 * t678;
t603 = t609 + (t663 * t696 - t629) * pkin(5) + (t662 * t696 - t630) * pkin(10);
t708 = sin(qJ(6));
t712 = cos(qJ(6));
t597 = -t708 * t600 + t712 * t603;
t644 = -t708 * t663 + t712 * t696;
t616 = t644 * qJD(6) + t712 * t630 + t708 * t695;
t645 = t712 * t663 + t708 * t696;
t623 = -t644 * mrSges(7,1) + t645 * mrSges(7,2);
t626 = qJDD(6) - t629;
t653 = qJD(6) + t662;
t631 = -t653 * mrSges(7,2) + t644 * mrSges(7,3);
t595 = m(7) * t597 + t626 * mrSges(7,1) - t616 * mrSges(7,3) - t645 * t623 + t653 * t631;
t598 = t712 * t600 + t708 * t603;
t615 = -t645 * qJD(6) - t708 * t630 + t712 * t695;
t632 = t653 * mrSges(7,1) - t645 * mrSges(7,3);
t596 = m(7) * t598 - t626 * mrSges(7,2) + t615 * mrSges(7,3) + t644 * t623 - t653 * t632;
t727 = -t708 * t595 + t712 * t596;
t587 = m(6) * t602 - t695 * mrSges(6,2) + t629 * mrSges(6,3) - t662 * t642 - t696 * t648 + t727;
t601 = -t709 * t606 + t713 * t607;
t647 = -t696 * mrSges(6,2) - t662 * mrSges(6,3);
t599 = -t695 * pkin(5) - t694 * pkin(10) + t663 * t643 - t601;
t720 = -m(7) * t599 + t615 * mrSges(7,1) - t616 * mrSges(7,2) + t644 * t631 - t645 * t632;
t591 = m(6) * t601 + t695 * mrSges(6,1) - t630 * mrSges(6,3) - t663 * t642 + t696 * t647 + t720;
t728 = t713 * t587 - t709 * t591;
t722 = m(5) * t611 + qJDD(3) * mrSges(5,3) + qJD(3) * t684 + t675 * t734 + t728;
t579 = m(4) * t622 - qJDD(3) * mrSges(4,2) - qJD(3) * t683 + t676 * t734 + t749 * t678 + t722;
t621 = -t627 + t740;
t685 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t734;
t581 = t709 * t587 + t713 * t591;
t612 = -qJDD(3) * pkin(3) + t725 - t740;
t686 = mrSges(5,2) * t734 + qJD(3) * mrSges(5,3);
t719 = -m(5) * t612 + qJDD(3) * mrSges(5,1) + qJD(3) * t686 - t581;
t580 = m(4) * t621 + qJDD(3) * mrSges(4,1) + qJD(3) * t685 - t749 * t677 + (-t675 - t676) * t735 + t719;
t729 = t714 * t579 - t710 * t580;
t570 = m(3) * t639 - t717 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t729;
t573 = t710 * t579 + t714 * t580;
t572 = m(3) * t649 + t573;
t613 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t735 + t721;
t588 = t712 * t595 + t708 * t596;
t724 = -m(6) * t609 + t629 * mrSges(6,1) - t630 * mrSges(6,2) - t662 * t647 - t663 * t648 - t588;
t585 = m(5) * t613 - t678 * mrSges(5,1) - t677 * mrSges(5,3) - t684 * t735 - t686 * t734 + t724;
t718 = -m(4) * t633 + t678 * mrSges(4,1) - t677 * mrSges(4,2) - t683 * t735 + t685 * t734 - t585;
t584 = m(3) * t638 + qJDD(2) * mrSges(3,1) - t717 * mrSges(3,2) + t718;
t561 = t570 * t742 - t704 * t572 + t584 * t741;
t559 = m(2) * t681 + t561;
t566 = t715 * t570 - t711 * t584;
t565 = m(2) * t682 + t566;
t739 = t705 * t559 + t703 * t565;
t738 = t751 * qJD(3) + (t747 * t710 + t746 * t714) * qJD(2);
t737 = -t746 * qJD(3) + (-t748 * t710 - t752 * t714) * qJD(2);
t736 = t747 * qJD(3) + (t753 * t710 + t748 * t714) * qJD(2);
t560 = t570 * t744 + t706 * t572 + t584 * t743;
t730 = -t703 * t559 + t705 * t565;
t617 = Ifges(7,5) * t645 + Ifges(7,6) * t644 + Ifges(7,3) * t653;
t619 = Ifges(7,1) * t645 + Ifges(7,4) * t644 + Ifges(7,5) * t653;
t589 = -mrSges(7,1) * t599 + mrSges(7,3) * t598 + Ifges(7,4) * t616 + Ifges(7,2) * t615 + Ifges(7,6) * t626 - t645 * t617 + t653 * t619;
t618 = Ifges(7,4) * t645 + Ifges(7,2) * t644 + Ifges(7,6) * t653;
t590 = mrSges(7,2) * t599 - mrSges(7,3) * t597 + Ifges(7,1) * t616 + Ifges(7,4) * t615 + Ifges(7,5) * t626 + t644 * t617 - t653 * t618;
t635 = Ifges(6,5) * t663 - Ifges(6,6) * t662 + Ifges(6,3) * t696;
t636 = Ifges(6,4) * t663 - Ifges(6,2) * t662 + Ifges(6,6) * t696;
t574 = mrSges(6,2) * t609 - mrSges(6,3) * t601 + Ifges(6,1) * t630 + Ifges(6,4) * t629 + Ifges(6,5) * t695 - pkin(10) * t588 - t708 * t589 + t712 * t590 - t662 * t635 - t696 * t636;
t637 = Ifges(6,1) * t663 - Ifges(6,4) * t662 + Ifges(6,5) * t696;
t575 = -mrSges(6,1) * t609 - mrSges(7,1) * t597 + mrSges(7,2) * t598 + mrSges(6,3) * t602 + Ifges(6,4) * t630 - Ifges(7,5) * t616 + Ifges(6,2) * t629 + Ifges(6,6) * t695 - Ifges(7,6) * t615 - Ifges(7,3) * t626 - pkin(5) * t588 - t645 * t618 + t644 * t619 - t663 * t635 + t696 * t637;
t557 = -mrSges(4,1) * t633 - mrSges(5,1) * t613 + mrSges(5,2) * t611 + mrSges(4,3) * t622 - pkin(3) * t585 - pkin(4) * t724 - pkin(9) * t728 + t736 * qJD(3) + t746 * qJDD(3) - t709 * t574 - t713 * t575 + t748 * t677 + t752 * t678 - t738 * t735;
t562 = mrSges(4,2) * t633 + mrSges(5,2) * t612 - mrSges(4,3) * t621 - mrSges(5,3) * t613 - pkin(9) * t581 - qJ(4) * t585 + t737 * qJD(3) + t747 * qJDD(3) + t713 * t574 - t709 * t575 + t753 * t677 + t748 * t678 + t738 * t734;
t555 = mrSges(3,2) * t649 - mrSges(3,3) * t638 + Ifges(3,5) * qJDD(2) - t717 * Ifges(3,6) - pkin(8) * t573 - t710 * t557 + t714 * t562;
t556 = t717 * Ifges(3,5) - mrSges(3,1) * t649 + t712 * t589 + t708 * t590 + Ifges(6,6) * t629 + Ifges(6,5) * t630 + mrSges(3,3) * t639 - mrSges(4,1) * t621 + mrSges(4,2) * t622 - mrSges(5,3) * t611 + mrSges(5,1) * t612 + mrSges(6,1) * t601 - mrSges(6,2) * t602 + Ifges(6,3) * t695 + t662 * t637 + t663 * t636 - pkin(3) * t719 + pkin(5) * t720 - qJ(4) * t722 + (-qJ(4) * mrSges(5,2) - t746) * t678 + (pkin(3) * mrSges(5,2) - t747) * t677 + pkin(10) * t727 + pkin(4) * t581 - pkin(2) * t573 + (t736 * t714 + (pkin(3) * t675 + t737) * t710) * qJD(2) + Ifges(3,6) * qJDD(2) - t751 * qJDD(3);
t723 = pkin(7) * t566 + t555 * t711 + t556 * t715;
t554 = mrSges(3,1) * t638 - mrSges(3,2) * t639 + Ifges(3,3) * qJDD(2) + pkin(2) * t718 + pkin(8) * t729 + t714 * t557 + t710 * t562;
t553 = mrSges(2,2) * t701 - mrSges(2,3) * t681 + t715 * t555 - t711 * t556 + (-t560 * t704 - t561 * t706) * pkin(7);
t552 = -mrSges(2,1) * t701 + mrSges(2,3) * t682 - pkin(1) * t560 - t704 * t554 + t723 * t706;
t1 = [-m(1) * g(1) + t730; -m(1) * g(2) + t739; -m(1) * g(3) + m(2) * t701 + t560; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t739 - t703 * t552 + t705 * t553; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t730 + t705 * t552 + t703 * t553; -mrSges(1,1) * g(2) + mrSges(2,1) * t681 + mrSges(1,2) * g(1) - mrSges(2,2) * t682 + pkin(1) * t561 + t706 * t554 + t723 * t704;];
tauB  = t1;
