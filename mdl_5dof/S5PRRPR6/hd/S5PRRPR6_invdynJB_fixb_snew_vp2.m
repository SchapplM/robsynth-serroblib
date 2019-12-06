% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:55
% EndTime: 2019-12-05 16:31:05
% DurationCPUTime: 7.33s
% Computational Cost: add. (84140->266), mult. (171872->348), div. (0->0), fcn. (117937->12), ass. (0->116)
t649 = sin(pkin(9));
t652 = cos(pkin(9));
t639 = t649 * g(1) - t652 * g(2);
t640 = -t652 * g(1) - t649 * g(2);
t647 = -g(3) + qJDD(1);
t650 = sin(pkin(5));
t653 = cos(pkin(5));
t656 = sin(qJ(2));
t659 = cos(qJ(2));
t604 = -t656 * t640 + (t639 * t653 + t647 * t650) * t659;
t679 = t653 * t656;
t680 = t650 * t656;
t605 = t639 * t679 + t659 * t640 + t647 * t680;
t661 = qJD(2) ^ 2;
t601 = -t661 * pkin(2) + qJDD(2) * pkin(7) + t605;
t619 = -t650 * t639 + t653 * t647;
t655 = sin(qJ(3));
t658 = cos(qJ(3));
t596 = t658 * t601 + t655 * t619;
t635 = (-pkin(3) * t658 - qJ(4) * t655) * qJD(2);
t660 = qJD(3) ^ 2;
t676 = t658 * qJD(2);
t584 = -t660 * pkin(3) + qJDD(3) * qJ(4) + t635 * t676 + t596;
t600 = -qJDD(2) * pkin(2) - t661 * pkin(7) - t604;
t675 = qJD(2) * qJD(3);
t674 = t658 * t675;
t637 = t655 * qJDD(2) + t674;
t645 = t655 * t675;
t638 = t658 * qJDD(2) - t645;
t588 = (-t637 - t674) * qJ(4) + (-t638 + t645) * pkin(3) + t600;
t648 = sin(pkin(10));
t651 = cos(pkin(10));
t677 = qJD(2) * t655;
t630 = t648 * qJD(3) + t651 * t677;
t579 = -0.2e1 * qJD(4) * t630 - t648 * t584 + t651 * t588;
t617 = t648 * qJDD(3) + t651 * t637;
t629 = t651 * qJD(3) - t648 * t677;
t577 = (-t629 * t676 - t617) * pkin(8) + (t629 * t630 - t638) * pkin(4) + t579;
t580 = 0.2e1 * qJD(4) * t629 + t651 * t584 + t648 * t588;
t616 = t651 * qJDD(3) - t648 * t637;
t618 = -pkin(4) * t676 - t630 * pkin(8);
t628 = t629 ^ 2;
t578 = -t628 * pkin(4) + t616 * pkin(8) + t618 * t676 + t580;
t654 = sin(qJ(5));
t657 = cos(qJ(5));
t576 = t654 * t577 + t657 * t578;
t595 = -t655 * t601 + t658 * t619;
t583 = -qJDD(3) * pkin(3) - t660 * qJ(4) + t635 * t677 + qJDD(4) - t595;
t581 = -t616 * pkin(4) - t628 * pkin(8) + t630 * t618 + t583;
t611 = t654 * t629 + t657 * t630;
t589 = -t611 * qJD(5) + t657 * t616 - t654 * t617;
t610 = t657 * t629 - t654 * t630;
t590 = t610 * qJD(5) + t654 * t616 + t657 * t617;
t644 = qJD(5) - t676;
t591 = Ifges(6,5) * t611 + Ifges(6,6) * t610 + Ifges(6,3) * t644;
t593 = Ifges(6,1) * t611 + Ifges(6,4) * t610 + Ifges(6,5) * t644;
t632 = qJDD(5) - t638;
t565 = -mrSges(6,1) * t581 + mrSges(6,3) * t576 + Ifges(6,4) * t590 + Ifges(6,2) * t589 + Ifges(6,6) * t632 - t611 * t591 + t644 * t593;
t575 = t657 * t577 - t654 * t578;
t592 = Ifges(6,4) * t611 + Ifges(6,2) * t610 + Ifges(6,6) * t644;
t566 = mrSges(6,2) * t581 - mrSges(6,3) * t575 + Ifges(6,1) * t590 + Ifges(6,4) * t589 + Ifges(6,5) * t632 + t610 * t591 - t644 * t592;
t606 = Ifges(5,5) * t630 + Ifges(5,6) * t629 - Ifges(5,3) * t676;
t608 = Ifges(5,1) * t630 + Ifges(5,4) * t629 - Ifges(5,5) * t676;
t602 = -t644 * mrSges(6,2) + t610 * mrSges(6,3);
t603 = t644 * mrSges(6,1) - t611 * mrSges(6,3);
t665 = m(6) * t581 - t589 * mrSges(6,1) + t590 * mrSges(6,2) - t610 * t602 + t611 * t603;
t597 = -t610 * mrSges(6,1) + t611 * mrSges(6,2);
t572 = m(6) * t575 + t632 * mrSges(6,1) - t590 * mrSges(6,3) - t611 * t597 + t644 * t602;
t573 = m(6) * t576 - t632 * mrSges(6,2) + t589 * mrSges(6,3) + t610 * t597 - t644 * t603;
t671 = -t654 * t572 + t657 * t573;
t551 = -mrSges(5,1) * t583 + mrSges(5,3) * t580 + Ifges(5,4) * t617 + Ifges(5,2) * t616 - Ifges(5,6) * t638 - pkin(4) * t665 + pkin(8) * t671 + t657 * t565 + t654 * t566 - t630 * t606 - t608 * t676;
t564 = t657 * t572 + t654 * t573;
t607 = Ifges(5,4) * t630 + Ifges(5,2) * t629 - Ifges(5,6) * t676;
t552 = mrSges(5,2) * t583 - mrSges(5,3) * t579 + Ifges(5,1) * t617 + Ifges(5,4) * t616 - Ifges(5,5) * t638 - pkin(8) * t564 - t654 * t565 + t657 * t566 + t629 * t606 + t607 * t676;
t612 = -t629 * mrSges(5,1) + t630 * mrSges(5,2);
t667 = mrSges(5,2) * t676 + t629 * mrSges(5,3);
t562 = m(5) * t579 - t638 * mrSges(5,1) - t617 * mrSges(5,3) - t630 * t612 - t667 * t676 + t564;
t615 = -mrSges(5,1) * t676 - t630 * mrSges(5,3);
t563 = m(5) * t580 + t638 * mrSges(5,2) + t616 * mrSges(5,3) + t629 * t612 + t615 * t676 + t671;
t560 = -t648 * t562 + t651 * t563;
t574 = m(5) * t583 - t616 * mrSges(5,1) + t617 * mrSges(5,2) + t630 * t615 - t629 * t667 + t665;
t625 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t655 + Ifges(4,2) * t658) * qJD(2);
t626 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t655 + Ifges(4,4) * t658) * qJD(2);
t682 = mrSges(4,1) * t595 - mrSges(4,2) * t596 + Ifges(4,5) * t637 + Ifges(4,6) * t638 + Ifges(4,3) * qJDD(3) - pkin(3) * t574 + qJ(4) * t560 + t651 * t551 + t648 * t552 + (t655 * t625 - t658 * t626) * qJD(2);
t559 = t651 * t562 + t648 * t563;
t641 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t677;
t642 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t676;
t664 = -m(4) * t600 + t638 * mrSges(4,1) - t637 * mrSges(4,2) - t641 * t677 + t642 * t676 - t559;
t555 = m(3) * t604 + qJDD(2) * mrSges(3,1) - t661 * mrSges(3,2) + t664;
t681 = t555 * t659;
t636 = (-mrSges(4,1) * t658 + mrSges(4,2) * t655) * qJD(2);
t558 = m(4) * t596 - qJDD(3) * mrSges(4,2) + t638 * mrSges(4,3) - qJD(3) * t641 + t636 * t676 + t560;
t568 = m(4) * t595 + qJDD(3) * mrSges(4,1) - t637 * mrSges(4,3) + qJD(3) * t642 - t636 * t677 - t574;
t672 = t658 * t558 - t655 * t568;
t547 = m(3) * t605 - t661 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t672;
t550 = t655 * t558 + t658 * t568;
t549 = m(3) * t619 + t550;
t537 = t547 * t679 - t650 * t549 + t653 * t681;
t535 = m(2) * t639 + t537;
t542 = t659 * t547 - t656 * t555;
t541 = m(2) * t640 + t542;
t678 = t652 * t535 + t649 * t541;
t536 = t547 * t680 + t653 * t549 + t650 * t681;
t673 = -t649 * t535 + t652 * t541;
t670 = m(2) * t647 + t536;
t624 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t655 + Ifges(4,6) * t658) * qJD(2);
t538 = mrSges(4,2) * t600 - mrSges(4,3) * t595 + Ifges(4,1) * t637 + Ifges(4,4) * t638 + Ifges(4,5) * qJDD(3) - qJ(4) * t559 - qJD(3) * t625 - t648 * t551 + t651 * t552 + t624 * t676;
t663 = mrSges(6,1) * t575 - mrSges(6,2) * t576 + Ifges(6,5) * t590 + Ifges(6,6) * t589 + Ifges(6,3) * t632 + t611 * t592 - t610 * t593;
t543 = -t624 * t677 + Ifges(4,6) * qJDD(3) + (Ifges(4,2) + Ifges(5,3)) * t638 + Ifges(4,4) * t637 + qJD(3) * t626 + t629 * t608 - t630 * t607 - t663 - Ifges(5,6) * t616 - Ifges(5,5) * t617 + mrSges(4,3) * t596 - mrSges(4,1) * t600 - mrSges(5,1) * t579 + mrSges(5,2) * t580 - pkin(4) * t564 - pkin(3) * t559;
t532 = mrSges(3,2) * t619 - mrSges(3,3) * t604 + Ifges(3,5) * qJDD(2) - t661 * Ifges(3,6) - pkin(7) * t550 + t658 * t538 - t655 * t543;
t533 = -mrSges(3,1) * t619 + mrSges(3,3) * t605 + t661 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t550 - t682;
t666 = pkin(6) * t542 + t532 * t656 + t533 * t659;
t531 = mrSges(3,1) * t604 - mrSges(3,2) * t605 + Ifges(3,3) * qJDD(2) + pkin(2) * t664 + pkin(7) * t672 + t655 * t538 + t658 * t543;
t530 = mrSges(2,2) * t647 - mrSges(2,3) * t639 + t659 * t532 - t656 * t533 + (-t536 * t650 - t537 * t653) * pkin(6);
t529 = -mrSges(2,1) * t647 + mrSges(2,3) * t640 - pkin(1) * t536 - t650 * t531 + t666 * t653;
t1 = [-m(1) * g(1) + t673; -m(1) * g(2) + t678; -m(1) * g(3) + t670; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t678 - t649 * t529 + t652 * t530; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t673 + t652 * t529 + t649 * t530; -mrSges(1,1) * g(2) + mrSges(2,1) * t639 + mrSges(1,2) * g(1) - mrSges(2,2) * t640 + pkin(1) * t537 + t653 * t531 + t666 * t650; t670; t531; t682; t574; t663;];
tauJB = t1;
