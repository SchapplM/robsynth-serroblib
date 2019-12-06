% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:50
% EndTime: 2019-12-05 16:43:54
% DurationCPUTime: 2.81s
% Computational Cost: add. (25489->237), mult. (51653->292), div. (0->0), fcn. (32333->8), ass. (0->100)
t682 = Ifges(5,4) + Ifges(6,4);
t691 = Ifges(5,2) + Ifges(6,2);
t687 = Ifges(5,6) + Ifges(6,6);
t652 = sin(qJ(4));
t653 = sin(qJ(3));
t655 = cos(qJ(4));
t656 = cos(qJ(3));
t622 = (-t652 * t653 + t655 * t656) * qJD(2);
t676 = qJD(2) * qJD(3);
t672 = t656 * t676;
t631 = t653 * qJDD(2) + t672;
t632 = t656 * qJDD(2) - t653 * t676;
t592 = t622 * qJD(4) + t655 * t631 + t652 * t632;
t623 = (t652 * t656 + t653 * t655) * qJD(2);
t605 = -t622 * mrSges(6,1) + t623 * mrSges(6,2);
t650 = sin(pkin(8));
t651 = cos(pkin(8));
t633 = t650 * g(1) - t651 * g(2);
t634 = -t651 * g(1) - t650 * g(2);
t654 = sin(qJ(2));
t657 = cos(qJ(2));
t612 = t654 * t633 + t657 * t634;
t658 = qJD(2) ^ 2;
t608 = -t658 * pkin(2) + qJDD(2) * pkin(6) + t612;
t649 = -g(3) + qJDD(1);
t593 = -t653 * t608 + t656 * t649;
t577 = (-t631 + t672) * pkin(7) + (t653 * t656 * t658 + qJDD(3)) * pkin(3) + t593;
t594 = t656 * t608 + t653 * t649;
t678 = qJD(2) * t653;
t637 = qJD(3) * pkin(3) - pkin(7) * t678;
t648 = t656 ^ 2;
t578 = -t648 * t658 * pkin(3) + t632 * pkin(7) - qJD(3) * t637 + t594;
t572 = t655 * t577 - t652 * t578;
t645 = qJDD(3) + qJDD(4);
t646 = qJD(3) + qJD(4);
t566 = -0.2e1 * qJD(5) * t623 + (t622 * t646 - t592) * qJ(5) + (t622 * t623 + t645) * pkin(4) + t572;
t613 = -t646 * mrSges(6,2) + t622 * mrSges(6,3);
t675 = m(6) * t566 + t645 * mrSges(6,1) + t646 * t613;
t562 = -t592 * mrSges(6,3) - t623 * t605 + t675;
t573 = t652 * t577 + t655 * t578;
t591 = -t623 * qJD(4) - t652 * t631 + t655 * t632;
t615 = t646 * pkin(4) - t623 * qJ(5);
t618 = t622 ^ 2;
t568 = -t618 * pkin(4) + t591 * qJ(5) + 0.2e1 * qJD(5) * t622 - t646 * t615 + t573;
t688 = Ifges(5,5) + Ifges(6,5);
t689 = Ifges(5,1) + Ifges(6,1);
t679 = t682 * t622 + t689 * t623 + t688 * t646;
t685 = t691 * t622 + t682 * t623 + t687 * t646;
t686 = Ifges(5,3) + Ifges(6,3);
t690 = mrSges(5,1) * t572 + mrSges(6,1) * t566 - mrSges(5,2) * t573 - mrSges(6,2) * t568 + pkin(4) * t562 + t687 * t591 + t688 * t592 - t679 * t622 + t685 * t623 + t686 * t645;
t606 = -t622 * mrSges(5,1) + t623 * mrSges(5,2);
t614 = -t646 * mrSges(5,2) + t622 * mrSges(5,3);
t556 = m(5) * t572 + t645 * mrSges(5,1) + t646 * t614 + (-t605 - t606) * t623 + (-mrSges(5,3) - mrSges(6,3)) * t592 + t675;
t616 = t646 * mrSges(6,1) - t623 * mrSges(6,3);
t617 = t646 * mrSges(5,1) - t623 * mrSges(5,3);
t674 = m(6) * t568 + t591 * mrSges(6,3) + t622 * t605;
t559 = m(5) * t573 + t591 * mrSges(5,3) + t622 * t606 + (-t616 - t617) * t646 + (-mrSges(5,2) - mrSges(6,2)) * t645 + t674;
t552 = t655 * t556 + t652 * t559;
t620 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t653 + Ifges(4,2) * t656) * qJD(2);
t621 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t653 + Ifges(4,4) * t656) * qJD(2);
t684 = mrSges(4,1) * t593 - mrSges(4,2) * t594 + Ifges(4,5) * t631 + Ifges(4,6) * t632 + Ifges(4,3) * qJDD(3) + pkin(3) * t552 + (t653 * t620 - t656 * t621) * qJD(2) + t690;
t630 = (-mrSges(4,1) * t656 + mrSges(4,2) * t653) * qJD(2);
t677 = qJD(2) * t656;
t636 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t677;
t550 = m(4) * t593 + qJDD(3) * mrSges(4,1) - t631 * mrSges(4,3) + qJD(3) * t636 - t630 * t678 + t552;
t635 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t678;
t667 = -t652 * t556 + t655 * t559;
t551 = m(4) * t594 - qJDD(3) * mrSges(4,2) + t632 * mrSges(4,3) - qJD(3) * t635 + t630 * t677 + t667;
t668 = -t653 * t550 + t656 * t551;
t542 = m(3) * t612 - t658 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t668;
t611 = t657 * t633 - t654 * t634;
t664 = -qJDD(2) * pkin(2) - t611;
t607 = -t658 * pkin(6) + t664;
t579 = -t632 * pkin(3) + t637 * t678 + (-pkin(7) * t648 - pkin(6)) * t658 + t664;
t570 = -t591 * pkin(4) - t618 * qJ(5) + t623 * t615 + qJDD(5) + t579;
t563 = m(6) * t570 - t591 * mrSges(6,1) + t592 * mrSges(6,2) - t622 * t613 + t623 * t616;
t662 = m(5) * t579 - t591 * mrSges(5,1) + t592 * mrSges(5,2) - t622 * t614 + t623 * t617 + t563;
t660 = -m(4) * t607 + t632 * mrSges(4,1) - t631 * mrSges(4,2) - t635 * t678 + t636 * t677 - t662;
t554 = m(3) * t611 + qJDD(2) * mrSges(3,1) - t658 * mrSges(3,2) + t660;
t539 = t654 * t542 + t657 * t554;
t537 = m(2) * t633 + t539;
t669 = t657 * t542 - t654 * t554;
t538 = m(2) * t634 + t669;
t681 = t651 * t537 + t650 * t538;
t544 = t656 * t550 + t653 * t551;
t680 = -t687 * t622 - t688 * t623 - t686 * t646;
t673 = m(3) * t649 + t544;
t670 = -t650 * t537 + t651 * t538;
t666 = m(2) * t649 + t673;
t545 = -mrSges(5,1) * t579 + mrSges(5,3) * t573 - mrSges(6,1) * t570 + mrSges(6,3) * t568 - pkin(4) * t563 + qJ(5) * t674 + (-qJ(5) * t616 + t679) * t646 + (-qJ(5) * mrSges(6,2) + t687) * t645 + t680 * t623 + t682 * t592 + t691 * t591;
t546 = mrSges(5,2) * t579 + mrSges(6,2) * t570 - mrSges(5,3) * t572 - mrSges(6,3) * t566 - qJ(5) * t562 + t682 * t591 + t689 * t592 - t680 * t622 + t688 * t645 - t685 * t646;
t619 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t653 + Ifges(4,6) * t656) * qJD(2);
t531 = -mrSges(4,1) * t607 + mrSges(4,3) * t594 + Ifges(4,4) * t631 + Ifges(4,2) * t632 + Ifges(4,6) * qJDD(3) - pkin(3) * t662 + pkin(7) * t667 + qJD(3) * t621 + t655 * t545 + t652 * t546 - t619 * t678;
t533 = mrSges(4,2) * t607 - mrSges(4,3) * t593 + Ifges(4,1) * t631 + Ifges(4,4) * t632 + Ifges(4,5) * qJDD(3) - pkin(7) * t552 - qJD(3) * t620 - t652 * t545 + t655 * t546 + t619 * t677;
t663 = mrSges(3,1) * t611 - mrSges(3,2) * t612 + Ifges(3,3) * qJDD(2) + pkin(2) * t660 + pkin(6) * t668 + t656 * t531 + t653 * t533;
t529 = -mrSges(3,1) * t649 + mrSges(3,3) * t612 + t658 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t544 - t684;
t528 = mrSges(3,2) * t649 - mrSges(3,3) * t611 + Ifges(3,5) * qJDD(2) - t658 * Ifges(3,6) - pkin(6) * t544 - t653 * t531 + t656 * t533;
t527 = mrSges(2,2) * t649 - mrSges(2,3) * t633 - pkin(5) * t539 + t657 * t528 - t654 * t529;
t526 = -mrSges(2,1) * t649 + mrSges(2,3) * t634 - pkin(1) * t673 + pkin(5) * t669 + t654 * t528 + t657 * t529;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t681; -m(1) * g(3) + t666; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t681 - t650 * t526 + t651 * t527; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t670 + t651 * t526 + t650 * t527; -mrSges(1,1) * g(2) + mrSges(2,1) * t633 + mrSges(1,2) * g(1) - mrSges(2,2) * t634 + pkin(1) * t539 + t663; t666; t663; t684; t690; t563;];
tauJB = t1;
