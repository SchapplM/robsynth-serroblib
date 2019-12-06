% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:04
% EndTime: 2019-12-05 15:22:08
% DurationCPUTime: 3.68s
% Computational Cost: add. (31256->223), mult. (75804->315), div. (0->0), fcn. (49449->10), ass. (0->117)
t648 = sin(pkin(7));
t651 = cos(pkin(7));
t630 = t648 * g(1) - t651 * g(2);
t631 = -t651 * g(1) - t648 * g(2);
t653 = sin(qJ(2));
t655 = cos(qJ(2));
t611 = t653 * t630 + t655 * t631;
t656 = qJD(2) ^ 2;
t699 = -t656 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t611;
t610 = t655 * t630 - t653 * t631;
t664 = -t656 * qJ(3) + qJDD(3) - t610;
t647 = sin(pkin(8));
t650 = cos(pkin(8));
t671 = -pkin(3) * t650 - qJ(4) * t647;
t691 = t647 * qJD(2);
t698 = (-pkin(2) + t671) * qJDD(2) + t664 - 0.2e1 * qJD(4) * t691;
t645 = -g(3) + qJDD(1);
t596 = t650 * t645 - t699 * t647;
t697 = mrSges(4,2) * t647;
t643 = t647 ^ 2;
t696 = t643 * t656;
t646 = sin(pkin(9));
t695 = t646 * t647;
t649 = cos(pkin(9));
t694 = t647 * t649;
t597 = t647 * t645 + t699 * t650;
t626 = (-mrSges(4,1) * t650 + t697) * qJD(2);
t625 = t671 * qJD(2);
t690 = t650 * qJD(2);
t589 = t625 * t690 + t597;
t669 = -pkin(4) * t650 - pkin(6) * t694;
t692 = t698 * t649;
t582 = t669 * qJDD(2) + (-t589 + (-pkin(4) * t643 * t649 + pkin(6) * t647 * t650) * t656) * t646 + t692;
t585 = t649 * t589 + t698 * t646;
t621 = t669 * qJD(2);
t686 = t646 ^ 2 * t696;
t688 = qJDD(2) * t647;
t583 = -t646 * pkin(6) * t688 - pkin(4) * t686 + t621 * t690 + t585;
t652 = sin(qJ(5));
t654 = cos(qJ(5));
t580 = t654 * t582 - t652 * t583;
t666 = (-t646 * t654 - t649 * t652) * t647;
t615 = qJD(2) * t666;
t665 = (-t646 * t652 + t649 * t654) * t647;
t616 = qJD(2) * t665;
t598 = -t615 * mrSges(6,1) + t616 * mrSges(6,2);
t602 = t615 * qJD(5) + qJDD(2) * t665;
t633 = qJD(5) - t690;
t607 = -t633 * mrSges(6,2) + t615 * mrSges(6,3);
t687 = t650 * qJDD(2);
t632 = qJDD(5) - t687;
t577 = m(6) * t580 + t632 * mrSges(6,1) - t602 * mrSges(6,3) - t616 * t598 + t633 * t607;
t581 = t652 * t582 + t654 * t583;
t601 = -t616 * qJD(5) + qJDD(2) * t666;
t608 = t633 * mrSges(6,1) - t616 * mrSges(6,3);
t578 = m(6) * t581 - t632 * mrSges(6,2) + t601 * mrSges(6,3) + t615 * t598 - t633 * t608;
t569 = t654 * t577 + t652 * t578;
t584 = -t646 * t589 + t692;
t675 = mrSges(5,1) * t646 + mrSges(5,2) * t649;
t617 = t675 * t691;
t667 = mrSges(5,2) * t650 - mrSges(5,3) * t695;
t619 = t667 * qJD(2);
t668 = -mrSges(5,1) * t650 - mrSges(5,3) * t694;
t567 = m(5) * t584 + t668 * qJDD(2) + (-t617 * t694 - t619 * t650) * qJD(2) + t569;
t620 = t668 * qJD(2);
t679 = -t652 * t577 + t654 * t578;
t568 = m(5) * t585 + t667 * qJDD(2) + (-t617 * t695 + t620 * t650) * qJD(2) + t679;
t680 = -t646 * t567 + t649 * t568;
t562 = m(4) * t597 + (qJDD(2) * mrSges(4,3) + qJD(2) * t626) * t650 + t680;
t588 = t625 * t691 + qJDD(4) - t596;
t586 = -pkin(6) * t686 + (pkin(4) * qJDD(2) * t646 + qJD(2) * t621 * t649) * t647 + t588;
t662 = m(6) * t586 - t601 * mrSges(6,1) + t602 * mrSges(6,2) - t615 * t607 + t616 * t608;
t658 = m(5) * t588 + t662;
t670 = t619 * t646 + t620 * t649;
t573 = m(4) * t596 + ((-mrSges(4,3) - t675) * qJDD(2) + (-t626 - t670) * qJD(2)) * t647 - t658;
t681 = t650 * t562 - t647 * t573;
t553 = m(3) * t611 - t656 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t681;
t565 = t649 * t567 + t646 * t568;
t605 = -qJDD(2) * pkin(2) + t664;
t659 = -m(4) * t605 + mrSges(4,1) * t687 - t565 + (t650 ^ 2 * t656 + t696) * mrSges(4,3);
t559 = m(3) * t610 - t656 * mrSges(3,2) + (mrSges(3,1) - t697) * qJDD(2) + t659;
t548 = t653 * t553 + t655 * t559;
t546 = m(2) * t630 + t548;
t682 = t655 * t553 - t653 * t559;
t547 = m(2) * t631 + t682;
t693 = t651 * t546 + t648 * t547;
t555 = t647 * t562 + t650 * t573;
t685 = m(3) * t645 + t555;
t683 = -t648 * t546 + t651 * t547;
t676 = m(2) * t645 + t685;
t674 = Ifges(4,1) * t647 + Ifges(4,4) * t650;
t673 = Ifges(4,5) * t647 + Ifges(4,6) * t650;
t672 = Ifges(5,5) * t649 - Ifges(5,6) * t646;
t590 = Ifges(6,5) * t616 + Ifges(6,6) * t615 + Ifges(6,3) * t633;
t592 = Ifges(6,1) * t616 + Ifges(6,4) * t615 + Ifges(6,5) * t633;
t570 = -mrSges(6,1) * t586 + mrSges(6,3) * t581 + Ifges(6,4) * t602 + Ifges(6,2) * t601 + Ifges(6,6) * t632 - t616 * t590 + t633 * t592;
t591 = Ifges(6,4) * t616 + Ifges(6,2) * t615 + Ifges(6,6) * t633;
t571 = mrSges(6,2) * t586 - mrSges(6,3) * t580 + Ifges(6,1) * t602 + Ifges(6,4) * t601 + Ifges(6,5) * t632 + t615 * t590 - t633 * t591;
t612 = (-Ifges(5,3) * t650 + t672 * t647) * qJD(2);
t661 = -Ifges(5,5) * t650 + (Ifges(5,1) * t649 - Ifges(5,4) * t646) * t647;
t614 = t661 * qJD(2);
t660 = -Ifges(5,6) * t650 + (Ifges(5,4) * t649 - Ifges(5,2) * t646) * t647;
t556 = -mrSges(5,1) * t588 + mrSges(5,3) * t585 + t652 * t571 + t654 * t570 - pkin(4) * t662 + pkin(6) * t679 + (-t612 * t694 - t650 * t614) * qJD(2) + t660 * qJDD(2);
t613 = t660 * qJD(2);
t557 = mrSges(5,2) * t588 - mrSges(5,3) * t584 - pkin(6) * t569 - t652 * t570 + t654 * t571 + (-t612 * t695 + t613 * t650) * qJD(2) + t661 * qJDD(2);
t627 = t673 * qJD(2);
t542 = mrSges(4,2) * t605 - mrSges(4,3) * t596 - qJ(4) * t565 + t674 * qJDD(2) - t646 * t556 + t649 * t557 + t627 * t690;
t657 = mrSges(6,1) * t580 - mrSges(6,2) * t581 + Ifges(6,5) * t602 + Ifges(6,6) * t601 + Ifges(6,3) * t632 + t616 * t591 - t615 * t592;
t550 = -mrSges(4,1) * t605 - mrSges(5,1) * t584 + mrSges(5,2) * t585 + mrSges(4,3) * t597 - pkin(3) * t565 - pkin(4) * t569 + (Ifges(4,2) + Ifges(5,3)) * t687 + ((Ifges(4,4) - t672) * qJDD(2) + (-t613 * t649 - t614 * t646 - t627) * qJD(2)) * t647 - t657;
t564 = mrSges(4,2) * t688 - t659;
t663 = mrSges(3,1) * t610 - mrSges(3,2) * t611 + Ifges(3,3) * qJDD(2) - pkin(2) * t564 + qJ(3) * t681 + t647 * t542 + t650 * t550;
t579 = (t670 * qJD(2) + t675 * qJDD(2)) * t647 + t658;
t540 = -mrSges(3,1) * t645 + mrSges(3,3) * t611 - mrSges(4,1) * t596 + mrSges(4,2) * t597 - t646 * t557 - t649 * t556 + pkin(3) * t579 - qJ(4) * t680 - pkin(2) * t555 + (Ifges(3,6) - t673) * qJDD(2) + (Ifges(3,5) - t647 * (Ifges(4,4) * t647 + Ifges(4,2) * t650) + t650 * t674) * t656;
t539 = mrSges(3,2) * t645 - mrSges(3,3) * t610 + Ifges(3,5) * qJDD(2) - t656 * Ifges(3,6) - qJ(3) * t555 + t650 * t542 - t647 * t550;
t538 = mrSges(2,2) * t645 - mrSges(2,3) * t630 - pkin(5) * t548 + t655 * t539 - t653 * t540;
t537 = -mrSges(2,1) * t645 + mrSges(2,3) * t631 - pkin(1) * t685 + pkin(5) * t682 + t653 * t539 + t655 * t540;
t1 = [-m(1) * g(1) + t683; -m(1) * g(2) + t693; -m(1) * g(3) + t676; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t693 - t648 * t537 + t651 * t538; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t683 + t651 * t537 + t648 * t538; -mrSges(1,1) * g(2) + mrSges(2,1) * t630 + mrSges(1,2) * g(1) - mrSges(2,2) * t631 + pkin(1) * t548 + t663; t676; t663; t564; t579; t657;];
tauJB = t1;
