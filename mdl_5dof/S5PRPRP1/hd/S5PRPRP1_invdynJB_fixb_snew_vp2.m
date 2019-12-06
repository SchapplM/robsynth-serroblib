% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:31
% EndTime: 2019-12-05 15:28:34
% DurationCPUTime: 2.16s
% Computational Cost: add. (19003->214), mult. (42262->263), div. (0->0), fcn. (27225->8), ass. (0->104)
t686 = Ifges(5,1) + Ifges(6,1);
t679 = Ifges(5,4) - Ifges(6,5);
t678 = Ifges(5,5) + Ifges(6,4);
t685 = Ifges(5,2) + Ifges(6,3);
t677 = Ifges(5,6) - Ifges(6,6);
t684 = Ifges(5,3) + Ifges(6,2);
t642 = qJD(2) ^ 2;
t638 = sin(qJ(4));
t636 = cos(pkin(8));
t682 = cos(qJ(4));
t659 = t636 * t682;
t634 = sin(pkin(8));
t665 = t634 * qJD(2);
t608 = -qJD(2) * t659 + t638 * t665;
t647 = t682 * t634 + t636 * t638;
t609 = t647 * qJD(2);
t590 = t608 * mrSges(6,1) - t609 * mrSges(6,3);
t667 = t608 * qJD(4);
t597 = t647 * qJDD(2) - t667;
t635 = sin(pkin(7));
t637 = cos(pkin(7));
t615 = t635 * g(1) - t637 * g(2);
t616 = -t637 * g(1) - t635 * g(2);
t639 = sin(qJ(2));
t640 = cos(qJ(2));
t601 = t639 * t615 + t640 * t616;
t598 = -t642 * pkin(2) + qJDD(2) * qJ(3) + t601;
t633 = -g(3) + qJDD(1);
t664 = qJD(2) * qJD(3);
t668 = t636 * t633 - 0.2e1 * t634 * t664;
t681 = pkin(3) * t636;
t574 = (-pkin(6) * qJDD(2) + t642 * t681 - t598) * t634 + t668;
t578 = t634 * t633 + (t598 + 0.2e1 * t664) * t636;
t662 = qJDD(2) * t636;
t629 = t636 ^ 2;
t674 = t629 * t642;
t575 = -pkin(3) * t674 + pkin(6) * t662 + t578;
t570 = t682 * t574 - t638 * t575;
t589 = t608 * pkin(4) - t609 * qJ(5);
t641 = qJD(4) ^ 2;
t567 = -qJDD(4) * pkin(4) - t641 * qJ(5) + t609 * t589 + qJDD(5) - t570;
t607 = -t608 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t653 = -m(6) * t567 + qJDD(4) * mrSges(6,1) + qJD(4) * t607;
t564 = t597 * mrSges(6,2) + t609 * t590 - t653;
t571 = t638 * t574 + t682 * t575;
t566 = -t641 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t608 * t589 + t571;
t663 = qJDD(2) * t634;
t666 = t609 * qJD(4);
t596 = -qJDD(2) * t659 + t638 * t663 + t666;
t606 = -qJD(4) * mrSges(6,1) + t609 * mrSges(6,2);
t661 = m(6) * t566 + qJDD(4) * mrSges(6,3) + qJD(4) * t606;
t670 = t678 * qJD(4) - t679 * t608 + t686 * t609;
t672 = -t677 * qJD(4) + t685 * t608 - t679 * t609;
t683 = t684 * qJDD(4) - t677 * t596 + t678 * t597 + t670 * t608 - t672 * t609 + mrSges(5,1) * t570 - mrSges(6,1) * t567 - mrSges(5,2) * t571 + mrSges(6,3) * t566 - pkin(4) * t564 + qJ(5) * (-t596 * mrSges(6,2) - t608 * t590 + t661);
t680 = -mrSges(5,3) - mrSges(6,2);
t675 = mrSges(4,2) * t634;
t605 = qJD(4) * mrSges(5,1) - t609 * mrSges(5,3);
t669 = -t608 * mrSges(5,1) - t609 * mrSges(5,2) - t590;
t560 = m(5) * t571 - qJDD(4) * mrSges(5,2) - qJD(4) * t605 + t680 * t596 + t669 * t608 + t661;
t604 = -qJD(4) * mrSges(5,2) - t608 * mrSges(5,3);
t561 = m(5) * t570 + qJDD(4) * mrSges(5,1) + qJD(4) * t604 + t680 * t597 + t669 * t609 + t653;
t552 = t638 * t560 + t682 * t561;
t577 = -t634 * t598 + t668;
t648 = mrSges(4,3) * qJDD(2) + t642 * (-mrSges(4,1) * t636 + t675);
t550 = m(4) * t577 - t648 * t634 + t552;
t655 = t682 * t560 - t638 * t561;
t551 = m(4) * t578 + t648 * t636 + t655;
t656 = -t634 * t550 + t636 * t551;
t542 = m(3) * t601 - t642 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t656;
t600 = t640 * t615 - t639 * t616;
t649 = qJDD(3) - t600;
t595 = -qJDD(2) * pkin(2) - t642 * qJ(3) + t649;
t628 = t634 ^ 2;
t576 = (-pkin(2) - t681) * qJDD(2) + (-qJ(3) + (-t628 - t629) * pkin(6)) * t642 + t649;
t569 = -0.2e1 * qJD(5) * t609 + (-t597 + t667) * qJ(5) + (t596 + t666) * pkin(4) + t576;
t562 = m(6) * t569 + t596 * mrSges(6,1) - t597 * mrSges(6,3) - t609 * t606 + t608 * t607;
t645 = m(5) * t576 + t596 * mrSges(5,1) + t597 * mrSges(5,2) + t608 * t604 + t609 * t605 + t562;
t643 = -m(4) * t595 + mrSges(4,1) * t662 - t645 + (t628 * t642 + t674) * mrSges(4,3);
t554 = t643 + (mrSges(3,1) - t675) * qJDD(2) - t642 * mrSges(3,2) + m(3) * t600;
t539 = t639 * t542 + t640 * t554;
t537 = m(2) * t615 + t539;
t657 = t640 * t542 - t639 * t554;
t538 = m(2) * t616 + t657;
t673 = t637 * t537 + t635 * t538;
t544 = t636 * t550 + t634 * t551;
t671 = -t684 * qJD(4) + t677 * t608 - t678 * t609;
t660 = m(3) * t633 + t544;
t658 = -t635 * t537 + t637 * t538;
t654 = m(2) * t633 + t660;
t652 = Ifges(4,1) * t634 + Ifges(4,4) * t636;
t651 = Ifges(4,4) * t634 + Ifges(4,2) * t636;
t650 = Ifges(4,5) * t634 + Ifges(4,6) * t636;
t545 = -mrSges(5,1) * t576 - mrSges(6,1) * t569 + mrSges(6,2) * t566 + mrSges(5,3) * t571 - pkin(4) * t562 + t670 * qJD(4) + t677 * qJDD(4) - t685 * t596 + t679 * t597 + t671 * t609;
t546 = mrSges(5,2) * t576 + mrSges(6,2) * t567 - mrSges(5,3) * t570 - mrSges(6,3) * t569 - qJ(5) * t562 + t672 * qJD(4) + t678 * qJDD(4) - t679 * t596 + t686 * t597 + t671 * t608;
t614 = t650 * qJD(2);
t531 = -mrSges(4,1) * t595 + mrSges(4,3) * t578 - pkin(3) * t645 + pkin(6) * t655 + t651 * qJDD(2) + t682 * t545 + t638 * t546 - t614 * t665;
t533 = t636 * qJD(2) * t614 + mrSges(4,2) * t595 - mrSges(4,3) * t577 - pkin(6) * t552 + t652 * qJDD(2) - t638 * t545 + t682 * t546;
t556 = mrSges(4,2) * t663 - t643;
t646 = mrSges(3,1) * t600 - mrSges(3,2) * t601 + Ifges(3,3) * qJDD(2) - pkin(2) * t556 + qJ(3) * t656 + t636 * t531 + t634 * t533;
t529 = (Ifges(3,6) - t650) * qJDD(2) - mrSges(3,1) * t633 + mrSges(3,3) * t601 - mrSges(4,1) * t577 + mrSges(4,2) * t578 - pkin(3) * t552 - pkin(2) * t544 + (-t634 * t651 + t636 * t652 + Ifges(3,5)) * t642 - t683;
t528 = mrSges(3,2) * t633 - mrSges(3,3) * t600 + Ifges(3,5) * qJDD(2) - t642 * Ifges(3,6) - qJ(3) * t544 - t634 * t531 + t636 * t533;
t527 = mrSges(2,2) * t633 - mrSges(2,3) * t615 - pkin(5) * t539 + t640 * t528 - t639 * t529;
t526 = -mrSges(2,1) * t633 + mrSges(2,3) * t616 - pkin(1) * t660 + pkin(5) * t657 + t639 * t528 + t640 * t529;
t1 = [-m(1) * g(1) + t658; -m(1) * g(2) + t673; -m(1) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t673 - t635 * t526 + t637 * t527; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t658 + t637 * t526 + t635 * t527; -mrSges(1,1) * g(2) + mrSges(2,1) * t615 + mrSges(1,2) * g(1) - mrSges(2,2) * t616 + pkin(1) * t539 + t646; t654; t646; t556; t683; t564;];
tauJB = t1;
