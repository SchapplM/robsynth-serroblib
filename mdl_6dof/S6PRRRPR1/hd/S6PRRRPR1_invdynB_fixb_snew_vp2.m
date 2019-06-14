% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:58:58
% EndTime: 2019-05-05 06:59:21
% DurationCPUTime: 21.95s
% Computational Cost: add. (369765->342), mult. (789476->443), div. (0->0), fcn. (586950->14), ass. (0->143)
t679 = sin(pkin(11));
t682 = cos(pkin(11));
t667 = g(1) * t679 - g(2) * t682;
t668 = -g(1) * t682 - g(2) * t679;
t677 = -g(3) + qJDD(1);
t680 = sin(pkin(6));
t683 = cos(pkin(6));
t687 = sin(qJ(2));
t691 = cos(qJ(2));
t639 = -t687 * t668 + (t667 * t683 + t677 * t680) * t691;
t714 = 2 * qJD(5);
t692 = qJD(2) ^ 2;
t695 = -qJDD(2) * pkin(2) - t639;
t632 = -pkin(8) * t692 + t695;
t686 = sin(qJ(3));
t690 = cos(qJ(3));
t707 = qJD(2) * qJD(3);
t706 = t690 * t707;
t665 = qJDD(2) * t686 + t706;
t666 = qJDD(2) * t690 - t686 * t707;
t709 = qJD(2) * t686;
t669 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t709;
t708 = qJD(2) * t690;
t670 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t708;
t672 = qJD(3) * pkin(3) - pkin(9) * t709;
t676 = t690 ^ 2;
t621 = -pkin(3) * t666 + t672 * t709 + (-pkin(9) * t676 - pkin(8)) * t692 + t695;
t685 = sin(qJ(4));
t689 = cos(qJ(4));
t658 = (t685 * t690 + t686 * t689) * qJD(2);
t628 = -qJD(4) * t658 - t665 * t685 + t666 * t689;
t657 = (-t685 * t686 + t689 * t690) * qJD(2);
t629 = qJD(4) * t657 + t665 * t689 + t666 * t685;
t675 = qJD(3) + qJD(4);
t648 = -mrSges(5,2) * t675 + mrSges(5,3) * t657;
t650 = mrSges(5,1) * t675 - mrSges(5,3) * t658;
t711 = t683 * t687;
t712 = t680 * t687;
t640 = t667 * t711 + t691 * t668 + t677 * t712;
t633 = -pkin(2) * t692 + qJDD(2) * pkin(8) + t640;
t651 = -t667 * t680 + t677 * t683;
t622 = -t633 * t686 + t690 * t651;
t607 = (-t665 + t706) * pkin(9) + (t686 * t690 * t692 + qJDD(3)) * pkin(3) + t622;
t623 = t690 * t633 + t686 * t651;
t610 = -pkin(3) * t676 * t692 + pkin(9) * t666 - qJD(3) * t672 + t623;
t590 = t689 * t607 - t610 * t685;
t674 = qJDD(3) + qJDD(4);
t587 = (t657 * t675 - t629) * qJ(5) + (t657 * t658 + t674) * pkin(4) + t590;
t591 = t685 * t607 + t689 * t610;
t649 = pkin(4) * t675 - qJ(5) * t658;
t653 = t657 ^ 2;
t589 = -pkin(4) * t653 + qJ(5) * t628 - t649 * t675 + t591;
t678 = sin(pkin(12));
t681 = cos(pkin(12));
t643 = t657 * t681 - t658 * t678;
t584 = t678 * t587 + t681 * t589 + t643 * t714;
t644 = t657 * t678 + t658 * t681;
t620 = -pkin(5) * t643 - pkin(10) * t644;
t673 = t675 ^ 2;
t582 = -pkin(5) * t673 + pkin(10) * t674 + t620 * t643 + t584;
t593 = -pkin(4) * t628 - qJ(5) * t653 + t658 * t649 + qJDD(5) + t621;
t608 = t628 * t681 - t629 * t678;
t609 = t628 * t678 + t629 * t681;
t585 = (-t643 * t675 - t609) * pkin(10) + (t644 * t675 - t608) * pkin(5) + t593;
t684 = sin(qJ(6));
t688 = cos(qJ(6));
t579 = -t582 * t684 + t585 * t688;
t625 = -t644 * t684 + t675 * t688;
t596 = qJD(6) * t625 + t609 * t688 + t674 * t684;
t606 = qJDD(6) - t608;
t626 = t644 * t688 + t675 * t684;
t611 = -mrSges(7,1) * t625 + mrSges(7,2) * t626;
t635 = qJD(6) - t643;
t612 = -mrSges(7,2) * t635 + mrSges(7,3) * t625;
t577 = m(7) * t579 + mrSges(7,1) * t606 - mrSges(7,3) * t596 - t611 * t626 + t612 * t635;
t580 = t582 * t688 + t585 * t684;
t595 = -qJD(6) * t626 - t609 * t684 + t674 * t688;
t613 = mrSges(7,1) * t635 - mrSges(7,3) * t626;
t578 = m(7) * t580 - mrSges(7,2) * t606 + mrSges(7,3) * t595 + t611 * t625 - t613 * t635;
t569 = t688 * t577 + t684 * t578;
t630 = -mrSges(6,2) * t675 + mrSges(6,3) * t643;
t631 = mrSges(6,1) * t675 - mrSges(6,3) * t644;
t698 = m(6) * t593 - t608 * mrSges(6,1) + t609 * mrSges(6,2) - t643 * t630 + t644 * t631 + t569;
t694 = m(5) * t621 - t628 * mrSges(5,1) + mrSges(5,2) * t629 - t657 * t648 + t650 * t658 + t698;
t693 = -m(4) * t632 + t666 * mrSges(4,1) - mrSges(4,2) * t665 - t669 * t709 + t670 * t708 - t694;
t565 = m(3) * t639 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t692 + t693;
t713 = t565 * t691;
t619 = -mrSges(6,1) * t643 + mrSges(6,2) * t644;
t701 = -t577 * t684 + t688 * t578;
t568 = m(6) * t584 - mrSges(6,2) * t674 + mrSges(6,3) * t608 + t619 * t643 - t631 * t675 + t701;
t700 = -t587 * t681 + t589 * t678;
t583 = -0.2e1 * qJD(5) * t644 - t700;
t581 = -pkin(5) * t674 - pkin(10) * t673 + (t714 + t620) * t644 + t700;
t696 = -m(7) * t581 + t595 * mrSges(7,1) - mrSges(7,2) * t596 + t625 * t612 - t613 * t626;
t573 = m(6) * t583 + mrSges(6,1) * t674 - mrSges(6,3) * t609 - t619 * t644 + t630 * t675 + t696;
t562 = t678 * t568 + t681 * t573;
t645 = -mrSges(5,1) * t657 + mrSges(5,2) * t658;
t560 = m(5) * t590 + mrSges(5,1) * t674 - mrSges(5,3) * t629 - t645 * t658 + t648 * t675 + t562;
t702 = t681 * t568 - t573 * t678;
t561 = m(5) * t591 - mrSges(5,2) * t674 + mrSges(5,3) * t628 + t645 * t657 - t650 * t675 + t702;
t554 = t689 * t560 + t685 * t561;
t664 = (-mrSges(4,1) * t690 + mrSges(4,2) * t686) * qJD(2);
t552 = m(4) * t622 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t665 + qJD(3) * t670 - t664 * t709 + t554;
t703 = -t560 * t685 + t689 * t561;
t553 = m(4) * t623 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t666 - qJD(3) * t669 + t664 * t708 + t703;
t704 = -t686 * t552 + t690 * t553;
t544 = m(3) * t640 - mrSges(3,1) * t692 - qJDD(2) * mrSges(3,2) + t704;
t547 = t690 * t552 + t686 * t553;
t546 = m(3) * t651 + t547;
t535 = t544 * t711 - t546 * t680 + t683 * t713;
t533 = m(2) * t667 + t535;
t539 = t691 * t544 - t565 * t687;
t538 = m(2) * t668 + t539;
t710 = t682 * t533 + t679 * t538;
t534 = t544 * t712 + t683 * t546 + t680 * t713;
t705 = -t533 * t679 + t682 * t538;
t597 = Ifges(7,5) * t626 + Ifges(7,6) * t625 + Ifges(7,3) * t635;
t599 = Ifges(7,1) * t626 + Ifges(7,4) * t625 + Ifges(7,5) * t635;
t570 = -mrSges(7,1) * t581 + mrSges(7,3) * t580 + Ifges(7,4) * t596 + Ifges(7,2) * t595 + Ifges(7,6) * t606 - t597 * t626 + t599 * t635;
t598 = Ifges(7,4) * t626 + Ifges(7,2) * t625 + Ifges(7,6) * t635;
t571 = mrSges(7,2) * t581 - mrSges(7,3) * t579 + Ifges(7,1) * t596 + Ifges(7,4) * t595 + Ifges(7,5) * t606 + t597 * t625 - t598 * t635;
t614 = Ifges(6,5) * t644 + Ifges(6,6) * t643 + Ifges(6,3) * t675;
t615 = Ifges(6,4) * t644 + Ifges(6,2) * t643 + Ifges(6,6) * t675;
t555 = mrSges(6,2) * t593 - mrSges(6,3) * t583 + Ifges(6,1) * t609 + Ifges(6,4) * t608 + Ifges(6,5) * t674 - pkin(10) * t569 - t570 * t684 + t571 * t688 + t614 * t643 - t615 * t675;
t616 = Ifges(6,1) * t644 + Ifges(6,4) * t643 + Ifges(6,5) * t675;
t556 = -mrSges(6,1) * t593 - mrSges(7,1) * t579 + mrSges(7,2) * t580 + mrSges(6,3) * t584 + Ifges(6,4) * t609 - Ifges(7,5) * t596 + Ifges(6,2) * t608 + Ifges(6,6) * t674 - Ifges(7,6) * t595 - Ifges(7,3) * t606 - pkin(5) * t569 - t598 * t626 + t599 * t625 - t614 * t644 + t616 * t675;
t636 = Ifges(5,5) * t658 + Ifges(5,6) * t657 + Ifges(5,3) * t675;
t638 = Ifges(5,1) * t658 + Ifges(5,4) * t657 + Ifges(5,5) * t675;
t540 = -mrSges(5,1) * t621 + mrSges(5,3) * t591 + Ifges(5,4) * t629 + Ifges(5,2) * t628 + Ifges(5,6) * t674 - pkin(4) * t698 + qJ(5) * t702 + t678 * t555 + t681 * t556 - t658 * t636 + t675 * t638;
t637 = Ifges(5,4) * t658 + Ifges(5,2) * t657 + Ifges(5,6) * t675;
t548 = mrSges(5,2) * t621 - mrSges(5,3) * t590 + Ifges(5,1) * t629 + Ifges(5,4) * t628 + Ifges(5,5) * t674 - qJ(5) * t562 + t555 * t681 - t556 * t678 + t636 * t657 - t637 * t675;
t654 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t686 + Ifges(4,6) * t690) * qJD(2);
t656 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t686 + Ifges(4,4) * t690) * qJD(2);
t529 = -mrSges(4,1) * t632 + mrSges(4,3) * t623 + Ifges(4,4) * t665 + Ifges(4,2) * t666 + Ifges(4,6) * qJDD(3) - pkin(3) * t694 + pkin(9) * t703 + qJD(3) * t656 + t689 * t540 + t685 * t548 - t654 * t709;
t655 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t686 + Ifges(4,2) * t690) * qJD(2);
t530 = mrSges(4,2) * t632 - mrSges(4,3) * t622 + Ifges(4,1) * t665 + Ifges(4,4) * t666 + Ifges(4,5) * qJDD(3) - pkin(9) * t554 - qJD(3) * t655 - t540 * t685 + t548 * t689 + t654 * t708;
t528 = mrSges(3,2) * t651 - mrSges(3,3) * t639 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t692 - pkin(8) * t547 - t529 * t686 + t530 * t690;
t531 = -pkin(5) * t696 + Ifges(3,6) * qJDD(2) - Ifges(4,3) * qJDD(3) - pkin(10) * t701 + (-Ifges(5,3) - Ifges(6,3)) * t674 + (-t686 * t655 + t690 * t656) * qJD(2) + t692 * Ifges(3,5) - t684 * t571 - t688 * t570 - Ifges(4,5) * t665 - Ifges(4,6) * t666 + t657 * t638 - t658 * t637 + mrSges(3,3) * t640 + t643 * t616 - t644 * t615 - mrSges(3,1) * t651 - Ifges(5,6) * t628 - Ifges(5,5) * t629 - mrSges(4,1) * t622 + mrSges(4,2) * t623 - Ifges(6,6) * t608 - Ifges(6,5) * t609 - mrSges(5,1) * t590 + mrSges(5,2) * t591 + mrSges(6,2) * t584 - mrSges(6,1) * t583 - pkin(4) * t562 - pkin(3) * t554 - pkin(2) * t547;
t697 = pkin(7) * t539 + t528 * t687 + t531 * t691;
t527 = mrSges(3,1) * t639 - mrSges(3,2) * t640 + Ifges(3,3) * qJDD(2) + pkin(2) * t693 + pkin(8) * t704 + t690 * t529 + t686 * t530;
t526 = mrSges(2,2) * t677 - mrSges(2,3) * t667 + t528 * t691 - t531 * t687 + (-t534 * t680 - t535 * t683) * pkin(7);
t525 = -mrSges(2,1) * t677 + mrSges(2,3) * t668 - pkin(1) * t534 - t527 * t680 + t683 * t697;
t1 = [-m(1) * g(1) + t705; -m(1) * g(2) + t710; -m(1) * g(3) + m(2) * t677 + t534; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t710 - t679 * t525 + t682 * t526; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t705 + t682 * t525 + t679 * t526; -mrSges(1,1) * g(2) + mrSges(2,1) * t667 + mrSges(1,2) * g(1) - mrSges(2,2) * t668 + pkin(1) * t535 + t527 * t683 + t680 * t697;];
tauB  = t1;
