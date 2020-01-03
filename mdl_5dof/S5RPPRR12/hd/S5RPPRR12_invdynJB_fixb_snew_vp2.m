% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:09
% EndTime: 2019-12-31 18:07:12
% DurationCPUTime: 2.45s
% Computational Cost: add. (24205->240), mult. (53586->291), div. (0->0), fcn. (34580->8), ass. (0->109)
t655 = sin(qJ(1));
t658 = cos(qJ(1));
t632 = t655 * g(1) - t658 * g(2);
t660 = qJD(1) ^ 2;
t669 = -t660 * qJ(2) + qJDD(2) - t632;
t691 = -pkin(1) - qJ(3);
t699 = -(2 * qJD(1) * qJD(3)) + t691 * qJDD(1) + t669;
t633 = -t658 * g(1) - t655 * g(2);
t698 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t633;
t651 = sin(pkin(8));
t652 = cos(pkin(8));
t654 = sin(qJ(4));
t657 = cos(qJ(4));
t673 = t651 * t657 + t652 * t654;
t627 = t673 * qJD(1);
t623 = t660 * pkin(1) - t698;
t697 = -m(3) * t623 + t660 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t672 = -t651 * t654 + t652 * t657;
t628 = t672 * qJD(1);
t685 = t628 * qJD(4);
t612 = -t673 * qJDD(1) - t685;
t695 = pkin(3) * t660;
t694 = mrSges(2,1) - mrSges(3,2);
t693 = -Ifges(3,4) + Ifges(2,5);
t692 = -Ifges(2,6) + Ifges(3,5);
t690 = mrSges(4,2) * t652;
t609 = t651 * g(3) + t699 * t652;
t594 = (-pkin(6) * qJDD(1) - t651 * t695) * t652 + t609;
t610 = -t652 * g(3) + t699 * t651;
t642 = t651 ^ 2;
t683 = qJDD(1) * t651;
t595 = -pkin(6) * t683 - t642 * t695 + t610;
t582 = t654 * t594 + t657 * t595;
t604 = t627 * mrSges(5,1) + t628 * mrSges(5,2);
t621 = qJD(4) * mrSges(5,1) - t628 * mrSges(5,3);
t611 = t627 * pkin(4) - t628 * pkin(7);
t659 = qJD(4) ^ 2;
t578 = -t659 * pkin(4) + qJDD(4) * pkin(7) - t627 * t611 + t582;
t668 = qJDD(3) + t698;
t688 = -t652 ^ 2 - t642;
t599 = pkin(3) * t683 + (t688 * pkin(6) + t691) * t660 + t668;
t686 = t627 * qJD(4);
t613 = t672 * qJDD(1) - t686;
t579 = (-t613 + t686) * pkin(7) + (-t612 + t685) * pkin(4) + t599;
t653 = sin(qJ(5));
t656 = cos(qJ(5));
t575 = -t653 * t578 + t656 * t579;
t615 = t656 * qJD(4) - t653 * t628;
t589 = t615 * qJD(5) + t653 * qJDD(4) + t656 * t613;
t616 = t653 * qJD(4) + t656 * t628;
t592 = -t615 * mrSges(6,1) + t616 * mrSges(6,2);
t625 = qJD(5) + t627;
t596 = -t625 * mrSges(6,2) + t615 * mrSges(6,3);
t608 = qJDD(5) - t612;
t572 = m(6) * t575 + t608 * mrSges(6,1) - t589 * mrSges(6,3) - t616 * t592 + t625 * t596;
t576 = t656 * t578 + t653 * t579;
t588 = -t616 * qJD(5) + t656 * qJDD(4) - t653 * t613;
t597 = t625 * mrSges(6,1) - t616 * mrSges(6,3);
t573 = m(6) * t576 - t608 * mrSges(6,2) + t588 * mrSges(6,3) + t615 * t592 - t625 * t597;
t677 = -t653 * t572 + t656 * t573;
t560 = m(5) * t582 - qJDD(4) * mrSges(5,2) + t612 * mrSges(5,3) - qJD(4) * t621 - t627 * t604 + t677;
t581 = t657 * t594 - t654 * t595;
t620 = -qJD(4) * mrSges(5,2) - t627 * mrSges(5,3);
t577 = -qJDD(4) * pkin(4) - t659 * pkin(7) + t628 * t611 - t581;
t666 = -m(6) * t577 + t588 * mrSges(6,1) - t589 * mrSges(6,2) + t615 * t596 - t616 * t597;
t568 = m(5) * t581 + qJDD(4) * mrSges(5,1) - t613 * mrSges(5,3) + qJD(4) * t620 - t628 * t604 + t666;
t552 = t654 * t560 + t657 * t568;
t671 = -qJDD(1) * mrSges(4,3) - t660 * (mrSges(4,1) * t651 + t690);
t550 = m(4) * t609 + t671 * t652 + t552;
t678 = t657 * t560 - t654 * t568;
t551 = m(4) * t610 + t671 * t651 + t678;
t546 = t652 * t550 + t651 * t551;
t626 = -qJDD(1) * pkin(1) + t669;
t667 = -m(3) * t626 + t660 * mrSges(3,3) - t546;
t542 = m(2) * t632 - t660 * mrSges(2,2) + t694 * qJDD(1) + t667;
t619 = t691 * t660 + t668;
t562 = t656 * t572 + t653 * t573;
t665 = m(5) * t599 - t612 * mrSges(5,1) + t613 * mrSges(5,2) + t627 * t620 + t628 * t621 + t562;
t664 = m(4) * t619 + mrSges(4,1) * t683 + qJDD(1) * t690 + t665;
t681 = t688 * mrSges(4,3);
t555 = (-mrSges(2,1) + t681) * t660 + t664 + m(2) * t633 - qJDD(1) * mrSges(2,2) + t697;
t689 = t658 * t542 + t655 * t555;
t674 = Ifges(4,5) * t652 - Ifges(4,6) * t651;
t687 = t660 * t674;
t680 = -t655 * t542 + t658 * t555;
t679 = -t651 * t550 + t652 * t551;
t676 = Ifges(4,1) * t652 - Ifges(4,4) * t651;
t675 = Ifges(4,4) * t652 - Ifges(4,2) * t651;
t583 = Ifges(6,5) * t616 + Ifges(6,6) * t615 + Ifges(6,3) * t625;
t585 = Ifges(6,1) * t616 + Ifges(6,4) * t615 + Ifges(6,5) * t625;
t565 = -mrSges(6,1) * t577 + mrSges(6,3) * t576 + Ifges(6,4) * t589 + Ifges(6,2) * t588 + Ifges(6,6) * t608 - t616 * t583 + t625 * t585;
t584 = Ifges(6,4) * t616 + Ifges(6,2) * t615 + Ifges(6,6) * t625;
t566 = mrSges(6,2) * t577 - mrSges(6,3) * t575 + Ifges(6,1) * t589 + Ifges(6,4) * t588 + Ifges(6,5) * t608 + t615 * t583 - t625 * t584;
t601 = Ifges(5,4) * t628 - Ifges(5,2) * t627 + Ifges(5,6) * qJD(4);
t602 = Ifges(5,1) * t628 - Ifges(5,4) * t627 + Ifges(5,5) * qJD(4);
t663 = mrSges(5,1) * t581 - mrSges(5,2) * t582 + Ifges(5,5) * t613 + Ifges(5,6) * t612 + Ifges(5,3) * qJDD(4) + pkin(4) * t666 + pkin(7) * t677 + t656 * t565 + t653 * t566 + t628 * t601 + t627 * t602;
t600 = Ifges(5,5) * t628 - Ifges(5,6) * t627 + Ifges(5,3) * qJD(4);
t547 = mrSges(5,2) * t599 - mrSges(5,3) * t581 + Ifges(5,1) * t613 + Ifges(5,4) * t612 + Ifges(5,5) * qJDD(4) - pkin(7) * t562 - qJD(4) * t601 - t653 * t565 + t656 * t566 - t627 * t600;
t661 = mrSges(6,1) * t575 - mrSges(6,2) * t576 + Ifges(6,5) * t589 + Ifges(6,6) * t588 + Ifges(6,3) * t608 + t616 * t584 - t615 * t585;
t548 = -mrSges(5,1) * t599 + mrSges(5,3) * t582 + Ifges(5,4) * t613 + Ifges(5,2) * t612 + Ifges(5,6) * qJDD(4) - pkin(4) * t562 + qJD(4) * t602 - t628 * t600 - t661;
t538 = -mrSges(4,1) * t619 + mrSges(4,3) * t610 - pkin(3) * t665 + pkin(6) * t678 + t675 * qJDD(1) + t654 * t547 + t657 * t548 - t652 * t687;
t540 = mrSges(4,2) * t619 - mrSges(4,3) * t609 - pkin(6) * t552 + t676 * qJDD(1) + t657 * t547 - t654 * t548 - t651 * t687;
t544 = qJDD(1) * mrSges(3,2) - t667;
t557 = t660 * t681 + t664;
t662 = -mrSges(2,2) * t633 - mrSges(3,3) * t623 - pkin(1) * t544 - qJ(3) * t546 - t651 * t538 + t652 * t540 + qJ(2) * (t557 + t697) + mrSges(3,2) * t626 + mrSges(2,1) * t632 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t545 = -m(3) * g(3) + t679;
t537 = (t674 + t693) * qJDD(1) + t663 + mrSges(4,1) * t609 - mrSges(4,2) * t610 + mrSges(3,1) * t626 - mrSges(2,3) * t632 + pkin(3) * t552 + pkin(2) * t546 - qJ(2) * t545 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t651 * t676 + t652 * t675 + t692) * t660;
t536 = -mrSges(3,1) * t623 + mrSges(2,3) * t633 - pkin(1) * t545 + pkin(2) * t557 + t694 * g(3) - qJ(3) * t679 - t692 * qJDD(1) - t652 * t538 - t651 * t540 + t693 * t660;
t1 = [-m(1) * g(1) + t680; -m(1) * g(2) + t689; (-m(1) - m(2) - m(3)) * g(3) + t679; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t689 - t655 * t536 + t658 * t537; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t680 + t658 * t536 + t655 * t537; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t662; t662; t544; t557; t663; t661;];
tauJB = t1;
