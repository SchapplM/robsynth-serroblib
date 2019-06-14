% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-05-05 17:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:15:25
% EndTime: 2019-05-05 17:15:31
% DurationCPUTime: 3.73s
% Computational Cost: add. (34691->321), mult. (76267->375), div. (0->0), fcn. (46890->8), ass. (0->128)
t706 = -2 * qJD(4);
t705 = Ifges(5,1) + Ifges(6,2);
t704 = (-Ifges(6,1) - Ifges(5,3));
t696 = Ifges(5,4) + Ifges(6,6);
t694 = Ifges(5,5) - Ifges(6,4);
t703 = -Ifges(5,2) - Ifges(6,3);
t692 = Ifges(5,6) - Ifges(6,5);
t658 = qJD(1) ^ 2;
t653 = sin(qJ(1));
t656 = cos(qJ(1));
t639 = -t656 * g(1) - t653 * g(2);
t670 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t639;
t699 = -pkin(1) - pkin(7);
t607 = t699 * t658 + t670;
t652 = sin(qJ(3));
t655 = cos(qJ(3));
t679 = qJD(1) * qJD(3);
t633 = -qJDD(1) * t652 - t655 * t679;
t677 = t652 * t679;
t634 = qJDD(1) * t655 - t677;
t683 = qJD(1) * t652;
t635 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t683;
t682 = qJD(1) * t655;
t637 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t682;
t636 = (qJD(3) * pkin(3)) - qJ(4) * t682;
t647 = t652 ^ 2;
t576 = -t633 * pkin(3) + qJDD(4) + t636 * t682 + (-qJ(4) * t647 + t699) * t658 + t670;
t650 = cos(pkin(9));
t690 = sin(pkin(9));
t597 = -t650 * t633 + t690 * t634;
t598 = t690 * t633 + t650 * t634;
t623 = t650 * t682 - t690 * t683;
t609 = (qJD(3) * mrSges(5,1)) - mrSges(5,3) * t623;
t622 = (t650 * t652 + t690 * t655) * qJD(1);
t681 = qJD(3) * t622;
t700 = -2 * qJD(5);
t660 = (-t598 + t681) * qJ(5) + t576 + (qJD(3) * pkin(4) + t700) * t623;
t560 = t597 * pkin(4) + t660;
t611 = mrSges(6,1) * t623 + (qJD(3) * mrSges(6,2));
t638 = t653 * g(1) - t656 * g(2);
t667 = -t658 * qJ(2) + qJDD(2) - t638;
t613 = t699 * qJDD(1) + t667;
t599 = t652 * g(3) + t655 * t613;
t573 = (-t634 - t677) * qJ(4) + (-t652 * t655 * t658 + qJDD(3)) * pkin(3) + t599;
t600 = -g(3) * t655 + t652 * t613;
t574 = -pkin(3) * t647 * t658 + qJ(4) * t633 - qJD(3) * t636 + t600;
t561 = t650 * t573 - t690 * t574 + t623 * t706;
t590 = pkin(4) * t622 - qJ(5) * t623;
t657 = qJD(3) ^ 2;
t558 = -qJDD(3) * pkin(4) - t657 * qJ(5) + t623 * t590 + qJDD(5) - t561;
t553 = (t622 * t623 - qJDD(3)) * pkin(8) + (t598 + t681) * pkin(5) + t558;
t612 = pkin(5) * t623 - qJD(3) * pkin(8);
t621 = t622 ^ 2;
t556 = t660 + (pkin(4) + pkin(8)) * t597 - t623 * t612 - t621 * pkin(5);
t651 = sin(qJ(6));
t654 = cos(qJ(6));
t551 = t553 * t654 - t556 * t651;
t601 = -qJD(3) * t651 + t622 * t654;
t569 = qJD(6) * t601 + qJDD(3) * t654 + t597 * t651;
t602 = qJD(3) * t654 + t622 * t651;
t577 = -mrSges(7,1) * t601 + mrSges(7,2) * t602;
t619 = qJD(6) + t623;
t580 = -mrSges(7,2) * t619 + mrSges(7,3) * t601;
t596 = qJDD(6) + t598;
t549 = m(7) * t551 + mrSges(7,1) * t596 - mrSges(7,3) * t569 - t577 * t602 + t580 * t619;
t552 = t553 * t651 + t556 * t654;
t568 = -qJD(6) * t602 - qJDD(3) * t651 + t597 * t654;
t581 = mrSges(7,1) * t619 - mrSges(7,3) * t602;
t550 = m(7) * t552 - mrSges(7,2) * t596 + mrSges(7,3) * t568 + t577 * t601 - t581 * t619;
t674 = -t651 * t549 + t654 * t550;
t669 = -m(6) * t560 + t598 * mrSges(6,3) + t623 * t611 - t674;
t610 = mrSges(6,1) * t622 - qJD(3) * mrSges(6,3);
t684 = -(qJD(3) * mrSges(5,2)) - mrSges(5,3) * t622 - t610;
t697 = mrSges(5,1) - mrSges(6,2);
t661 = m(5) * t576 + t598 * mrSges(5,2) + t697 * t597 + t623 * t609 + t684 * t622 - t669;
t702 = -m(4) * t607 + t633 * mrSges(4,1) - t634 * mrSges(4,2) - t635 * t683 - t637 * t682 - t661;
t698 = mrSges(2,1) - mrSges(3,2);
t695 = Ifges(2,5) - Ifges(3,4);
t693 = -Ifges(2,6) + Ifges(3,5);
t591 = mrSges(5,1) * t622 + mrSges(5,2) * t623;
t542 = t654 * t549 + t651 * t550;
t592 = -mrSges(6,2) * t622 - mrSges(6,3) * t623;
t664 = -m(6) * t558 - t598 * mrSges(6,1) - t623 * t592 - t542;
t540 = m(5) * t561 - t598 * mrSges(5,3) + t684 * qJD(3) + t697 * qJDD(3) - t623 * t591 + t664;
t617 = t622 * t706;
t688 = t690 * t573 + t650 * t574;
t562 = t617 + t688;
t668 = (t657 * pkin(4)) - qJDD(3) * qJ(5) - t688;
t557 = (qJD(3) * t700) + ((2 * qJD(4)) + t590) * t622 + t668;
t555 = -t597 * pkin(5) - t621 * pkin(8) - t622 * t590 + t617 + ((2 * qJD(5)) + t612) * qJD(3) - t668;
t666 = -m(7) * t555 + t568 * mrSges(7,1) - t569 * mrSges(7,2) + t601 * t580 - t602 * t581;
t662 = -m(6) * t557 + qJDD(3) * mrSges(6,3) + qJD(3) * t611 - t666;
t547 = m(5) * t562 - qJDD(3) * mrSges(5,2) - qJD(3) * t609 + (-t591 - t592) * t622 + (-mrSges(5,3) - mrSges(6,1)) * t597 + t662;
t535 = t650 * t540 + t690 * t547;
t632 = (mrSges(4,1) * t652 + mrSges(4,2) * t655) * qJD(1);
t533 = m(4) * t599 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t634 + qJD(3) * t635 - t632 * t682 + t535;
t672 = -t690 * t540 + t650 * t547;
t534 = m(4) * t600 - qJDD(3) * mrSges(4,2) + t633 * mrSges(4,3) - qJD(3) * t637 - t632 * t683 + t672;
t530 = t655 * t533 + t652 * t534;
t620 = -qJDD(1) * pkin(1) + t667;
t665 = -m(3) * t620 + t658 * mrSges(3,3) - t530;
t528 = m(2) * t638 - t658 * mrSges(2,2) + t698 * qJDD(1) + t665;
t614 = t658 * pkin(1) - t670;
t659 = -m(3) * t614 + t658 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t702;
t538 = m(2) * t639 - t658 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t659;
t689 = t656 * t528 + t653 * t538;
t687 = (t704 * qJD(3)) + t692 * t622 - t694 * t623;
t686 = t692 * qJD(3) + t703 * t622 + t696 * t623;
t685 = t694 * qJD(3) - t696 * t622 + t705 * t623;
t676 = -t528 * t653 + t656 * t538;
t675 = -t652 * t533 + t655 * t534;
t626 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t655 - Ifges(4,4) * t652) * qJD(1);
t625 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t655 - Ifges(4,2) * t652) * qJD(1);
t624 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t655 - Ifges(4,6) * t652) * qJD(1);
t565 = Ifges(7,1) * t602 + Ifges(7,4) * t601 + Ifges(7,5) * t619;
t564 = Ifges(7,4) * t602 + Ifges(7,2) * t601 + Ifges(7,6) * t619;
t563 = Ifges(7,5) * t602 + Ifges(7,6) * t601 + Ifges(7,3) * t619;
t544 = mrSges(7,2) * t555 - mrSges(7,3) * t551 + Ifges(7,1) * t569 + Ifges(7,4) * t568 + Ifges(7,5) * t596 + t563 * t601 - t564 * t619;
t543 = -mrSges(7,1) * t555 + mrSges(7,3) * t552 + Ifges(7,4) * t569 + Ifges(7,2) * t568 + Ifges(7,6) * t596 - t563 * t602 + t565 * t619;
t541 = -t597 * mrSges(6,2) - t622 * t610 - t669;
t531 = mrSges(6,1) * t558 + mrSges(7,1) * t551 + mrSges(5,2) * t576 - mrSges(7,2) * t552 - mrSges(5,3) * t561 - mrSges(6,3) * t560 + Ifges(7,5) * t569 + Ifges(7,6) * t568 + Ifges(7,3) * t596 + pkin(5) * t542 - qJ(5) * t541 + t602 * t564 - t601 * t565 + t687 * t622 + t705 * t598 - t696 * t597 + t694 * qJDD(3) - t686 * qJD(3);
t529 = -m(3) * g(3) + t675;
t526 = -mrSges(5,1) * t576 - mrSges(6,1) * t557 + mrSges(6,2) * t560 + mrSges(5,3) * t562 - pkin(4) * t541 - pkin(5) * t666 - pkin(8) * t674 + t685 * qJD(3) + t692 * qJDD(3) - t654 * t543 - t651 * t544 + t703 * t597 + t696 * t598 + t687 * t623;
t525 = mrSges(4,2) * t607 - mrSges(4,3) * t599 + Ifges(4,1) * t634 + Ifges(4,4) * t633 + Ifges(4,5) * qJDD(3) - qJ(4) * t535 - qJD(3) * t625 - t690 * t526 + t650 * t531 - t624 * t683;
t524 = -mrSges(4,1) * t607 + mrSges(4,3) * t600 + Ifges(4,4) * t634 + Ifges(4,2) * t633 + Ifges(4,6) * qJDD(3) - pkin(3) * t661 + qJ(4) * t672 + qJD(3) * t626 + t650 * t526 + t690 * t531 - t624 * t682;
t523 = (-mrSges(6,1) * qJ(5) - t692) * t597 + t693 * t658 + t694 * t598 + t695 * qJDD(1) + qJ(5) * t662 + pkin(4) * (-qJD(3) * t610 + t664) + (t625 * t655 + t626 * t652) * qJD(1) + (-mrSges(6,2) * pkin(4) + Ifges(4,3) - t704) * qJDD(3) + (-qJ(5) * t592 + t685) * t622 + t686 * t623 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t654 * t544 - t651 * t543 + Ifges(4,6) * t633 + Ifges(4,5) * t634 - mrSges(2,3) * t638 + mrSges(3,1) * t620 + mrSges(4,1) * t599 - mrSges(4,2) * t600 + mrSges(6,2) * t558 + mrSges(5,1) * t561 - mrSges(5,2) * t562 - mrSges(6,3) * t557 - pkin(8) * t542 + pkin(3) * t535 + pkin(2) * t530 - qJ(2) * t529;
t522 = -mrSges(3,1) * t614 + mrSges(2,3) * t639 - pkin(1) * t529 - pkin(2) * t702 - pkin(7) * t675 + t698 * g(3) - t693 * qJDD(1) - t655 * t524 - t652 * t525 + t695 * t658;
t1 = [-m(1) * g(1) + t676; -m(1) * g(2) + t689; (-m(1) - m(2) - m(3)) * g(3) + t675; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t689 - t653 * t522 + t656 * t523; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t676 + t656 * t522 + t653 * t523; pkin(1) * t665 + qJ(2) * t659 + t655 * t525 - t652 * t524 - pkin(7) * t530 + mrSges(2,1) * t638 - mrSges(2,2) * t639 + mrSges(3,2) * t620 - mrSges(3,3) * t614 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
