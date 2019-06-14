% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:33:03
% EndTime: 2019-05-05 16:33:10
% DurationCPUTime: 5.99s
% Computational Cost: add. (65136->326), mult. (141711->391), div. (0->0), fcn. (88614->10), ass. (0->130)
t700 = -2 * qJD(4);
t699 = Ifges(5,1) + Ifges(6,2);
t698 = -Ifges(6,1) - Ifges(5,3);
t695 = Ifges(5,4) + Ifges(6,6);
t694 = Ifges(5,5) - Ifges(6,4);
t697 = Ifges(5,2) + Ifges(6,3);
t693 = Ifges(5,6) - Ifges(6,5);
t660 = sin(qJ(1));
t663 = cos(qJ(1));
t645 = t660 * g(1) - g(2) * t663;
t637 = qJDD(1) * pkin(1) + t645;
t646 = -g(1) * t663 - g(2) * t660;
t665 = qJD(1) ^ 2;
t639 = -pkin(1) * t665 + t646;
t656 = sin(pkin(9));
t657 = cos(pkin(9));
t608 = t656 * t637 + t657 * t639;
t595 = -pkin(2) * t665 + qJDD(1) * pkin(7) + t608;
t654 = -g(3) + qJDD(2);
t659 = sin(qJ(3));
t662 = cos(qJ(3));
t583 = -t659 * t595 + t662 * t654;
t681 = qJD(1) * qJD(3);
t679 = t662 * t681;
t640 = qJDD(1) * t659 + t679;
t570 = (-t640 + t679) * qJ(4) + (t659 * t662 * t665 + qJDD(3)) * pkin(3) + t583;
t584 = t662 * t595 + t659 * t654;
t641 = qJDD(1) * t662 - t659 * t681;
t685 = qJD(1) * t659;
t642 = qJD(3) * pkin(3) - qJ(4) * t685;
t653 = t662 ^ 2;
t571 = -pkin(3) * t653 * t665 + qJ(4) * t641 - qJD(3) * t642 + t584;
t655 = sin(pkin(10));
t692 = cos(pkin(10));
t626 = (t655 * t662 + t692 * t659) * qJD(1);
t565 = t692 * t570 - t655 * t571 + t626 * t700;
t696 = -2 * qJD(5);
t684 = qJD(1) * t662;
t625 = t655 * t685 - t692 * t684;
t600 = mrSges(5,1) * t625 + mrSges(5,2) * t626;
t610 = t692 * t640 + t655 * t641;
t614 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t625;
t616 = mrSges(6,1) * t625 - qJD(3) * mrSges(6,3);
t599 = pkin(4) * t625 - qJ(5) * t626;
t664 = qJD(3) ^ 2;
t562 = -qJDD(3) * pkin(4) - t664 * qJ(5) + t626 * t599 + qJDD(5) - t565;
t683 = qJD(3) * t625;
t557 = (t625 * t626 - qJDD(3)) * pkin(8) + (t610 + t683) * pkin(5) + t562;
t609 = t655 * t640 - t692 * t641;
t618 = pkin(5) * t626 - qJD(3) * pkin(8);
t624 = t625 ^ 2;
t607 = t657 * t637 - t656 * t639;
t673 = -qJDD(1) * pkin(2) - t607;
t572 = -t641 * pkin(3) + qJDD(4) + t642 * t685 + (-qJ(4) * t653 - pkin(7)) * t665 + t673;
t667 = (-t610 + t683) * qJ(5) + t572 + (qJD(3) * pkin(4) + t696) * t626;
t560 = -t624 * pkin(5) - t626 * t618 + (pkin(4) + pkin(8)) * t609 + t667;
t658 = sin(qJ(6));
t661 = cos(qJ(6));
t555 = t557 * t661 - t560 * t658;
t611 = -qJD(3) * t658 + t625 * t661;
t579 = qJD(6) * t611 + qJDD(3) * t661 + t609 * t658;
t612 = qJD(3) * t661 + t625 * t658;
t580 = -mrSges(7,1) * t611 + mrSges(7,2) * t612;
t623 = qJD(6) + t626;
t585 = -mrSges(7,2) * t623 + mrSges(7,3) * t611;
t606 = qJDD(6) + t610;
t553 = m(7) * t555 + mrSges(7,1) * t606 - mrSges(7,3) * t579 - t580 * t612 + t585 * t623;
t556 = t557 * t658 + t560 * t661;
t578 = -qJD(6) * t612 - qJDD(3) * t658 + t609 * t661;
t586 = mrSges(7,1) * t623 - mrSges(7,3) * t612;
t554 = m(7) * t556 - mrSges(7,2) * t606 + mrSges(7,3) * t578 + t580 * t611 - t586 * t623;
t545 = t661 * t553 + t658 * t554;
t601 = -mrSges(6,2) * t625 - mrSges(6,3) * t626;
t670 = -m(6) * t562 - t610 * mrSges(6,1) - t626 * t601 - t545;
t543 = m(5) * t565 - t610 * mrSges(5,3) - t626 * t600 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t614 - t616) * qJD(3) + t670;
t621 = t625 * t700;
t689 = t655 * t570 + t692 * t571;
t566 = t621 + t689;
t615 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t626;
t672 = t664 * pkin(4) - qJDD(3) * qJ(5) - t689;
t561 = qJD(3) * t696 + ((2 * qJD(4)) + t599) * t625 + t672;
t617 = mrSges(6,1) * t626 + qJD(3) * mrSges(6,2);
t559 = -t609 * pkin(5) - t624 * pkin(8) - t625 * t599 + t621 + ((2 * qJD(5)) + t618) * qJD(3) - t672;
t671 = -m(7) * t559 + t578 * mrSges(7,1) - t579 * mrSges(7,2) + t611 * t585 - t612 * t586;
t669 = -m(6) * t561 + qJDD(3) * mrSges(6,3) + qJD(3) * t617 - t671;
t550 = m(5) * t566 - qJDD(3) * mrSges(5,2) - qJD(3) * t615 + (-t600 - t601) * t625 + (-mrSges(5,3) - mrSges(6,1)) * t609 + t669;
t539 = t692 * t543 + t655 * t550;
t638 = (-mrSges(4,1) * t662 + mrSges(4,2) * t659) * qJD(1);
t644 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t684;
t537 = m(4) * t583 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t640 + qJD(3) * t644 - t638 * t685 + t539;
t643 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t685;
t675 = -t543 * t655 + t692 * t550;
t538 = m(4) * t584 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t641 - qJD(3) * t643 + t638 * t684 + t675;
t676 = -t537 * t659 + t662 * t538;
t531 = m(3) * t608 - mrSges(3,1) * t665 - qJDD(1) * mrSges(3,2) + t676;
t594 = -t665 * pkin(7) + t673;
t564 = t609 * pkin(4) + t667;
t690 = -t658 * t553 + t661 * t554;
t544 = m(6) * t564 - t609 * mrSges(6,2) - t610 * mrSges(6,3) - t625 * t616 - t626 * t617 + t690;
t668 = m(5) * t572 + t609 * mrSges(5,1) + mrSges(5,2) * t610 + t625 * t614 + t615 * t626 + t544;
t666 = -m(4) * t594 + t641 * mrSges(4,1) - mrSges(4,2) * t640 - t643 * t685 + t644 * t684 - t668;
t541 = m(3) * t607 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t665 + t666;
t527 = t656 * t531 + t657 * t541;
t525 = m(2) * t645 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t665 + t527;
t677 = t657 * t531 - t541 * t656;
t526 = m(2) * t646 - mrSges(2,1) * t665 - qJDD(1) * mrSges(2,2) + t677;
t691 = t663 * t525 + t660 * t526;
t532 = t662 * t537 + t659 * t538;
t688 = qJD(3) * t698 + t625 * t693 - t626 * t694;
t687 = -qJD(3) * t693 + t625 * t697 - t626 * t695;
t686 = qJD(3) * t694 - t625 * t695 + t626 * t699;
t680 = m(3) * t654 + t532;
t678 = -t525 * t660 + t663 * t526;
t632 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t659 + Ifges(4,4) * t662) * qJD(1);
t631 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t659 + Ifges(4,2) * t662) * qJD(1);
t630 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t659 + Ifges(4,6) * t662) * qJD(1);
t575 = Ifges(7,1) * t612 + Ifges(7,4) * t611 + Ifges(7,5) * t623;
t574 = Ifges(7,4) * t612 + Ifges(7,2) * t611 + Ifges(7,6) * t623;
t573 = Ifges(7,5) * t612 + Ifges(7,6) * t611 + Ifges(7,3) * t623;
t547 = mrSges(7,2) * t559 - mrSges(7,3) * t555 + Ifges(7,1) * t579 + Ifges(7,4) * t578 + Ifges(7,5) * t606 + t573 * t611 - t574 * t623;
t546 = -mrSges(7,1) * t559 + mrSges(7,3) * t556 + Ifges(7,4) * t579 + Ifges(7,2) * t578 + Ifges(7,6) * t606 - t573 * t612 + t575 * t623;
t533 = mrSges(6,1) * t562 + mrSges(7,1) * t555 + mrSges(5,2) * t572 - mrSges(7,2) * t556 - mrSges(5,3) * t565 - mrSges(6,3) * t564 + Ifges(7,5) * t579 + Ifges(7,6) * t578 + Ifges(7,3) * t606 + pkin(5) * t545 - qJ(5) * t544 + t612 * t574 - t611 * t575 + t688 * t625 + t699 * t610 - t695 * t609 + t694 * qJDD(3) + t687 * qJD(3);
t528 = -mrSges(5,1) * t572 - mrSges(6,1) * t561 + mrSges(6,2) * t564 + mrSges(5,3) * t566 - pkin(4) * t544 - pkin(5) * t671 - pkin(8) * t690 + t686 * qJD(3) + t693 * qJDD(3) - t661 * t546 - t658 * t547 - t609 * t697 + t695 * t610 + t688 * t626;
t521 = mrSges(4,2) * t594 - mrSges(4,3) * t583 + Ifges(4,1) * t640 + Ifges(4,4) * t641 + Ifges(4,5) * qJDD(3) - qJ(4) * t539 - qJD(3) * t631 - t655 * t528 + t692 * t533 + t630 * t684;
t520 = (mrSges(6,2) * pkin(4) - Ifges(4,3) + t698) * qJDD(3) + Ifges(3,6) * qJDD(1) + (qJ(5) * t601 - t686) * t625 + t687 * t626 - pkin(2) * t532 + (-t631 * t659 + t632 * t662) * qJD(1) - qJ(5) * t669 - pkin(4) * (-qJD(3) * t616 + t670) + t665 * Ifges(3,5) + t658 * t546 - t661 * t547 - mrSges(3,1) * t654 - Ifges(4,5) * t640 - Ifges(4,6) * t641 + mrSges(3,3) * t608 - pkin(3) * t539 + pkin(8) * t545 + mrSges(6,3) * t561 - mrSges(6,2) * t562 - mrSges(5,1) * t565 + mrSges(5,2) * t566 - mrSges(4,1) * t583 + mrSges(4,2) * t584 + (mrSges(6,1) * qJ(5) + t693) * t609 - t694 * t610;
t519 = -mrSges(4,1) * t594 + mrSges(4,3) * t584 + Ifges(4,4) * t640 + Ifges(4,2) * t641 + Ifges(4,6) * qJDD(3) - pkin(3) * t668 + qJ(4) * t675 + qJD(3) * t632 + t692 * t528 + t655 * t533 - t630 * t685;
t518 = mrSges(3,2) * t654 - mrSges(3,3) * t607 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t665 - pkin(7) * t532 - t519 * t659 + t521 * t662;
t517 = -mrSges(2,2) * g(3) - mrSges(2,3) * t645 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t665 - qJ(2) * t527 + t518 * t657 - t520 * t656;
t516 = mrSges(2,1) * g(3) + mrSges(2,3) * t646 + t665 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t680 + qJ(2) * t677 + t656 * t518 + t657 * t520;
t1 = [-m(1) * g(1) + t678; -m(1) * g(2) + t691; (-m(1) - m(2)) * g(3) + t680; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t691 - t660 * t516 + t663 * t517; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t678 + t663 * t516 + t660 * t517; pkin(1) * t527 + mrSges(2,1) * t645 - mrSges(2,2) * t646 + t659 * t521 + t662 * t519 + pkin(2) * t666 + pkin(7) * t676 + mrSges(3,1) * t607 - mrSges(3,2) * t608 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
