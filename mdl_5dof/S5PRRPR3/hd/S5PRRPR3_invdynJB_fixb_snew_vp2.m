% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:27
% EndTime: 2019-12-05 16:19:33
% DurationCPUTime: 5.17s
% Computational Cost: add. (57040->259), mult. (125738->335), div. (0->0), fcn. (84148->10), ass. (0->108)
t638 = sin(pkin(8));
t640 = cos(pkin(8));
t622 = t638 * g(1) - t640 * g(2);
t623 = -t640 * g(1) - t638 * g(2);
t643 = sin(qJ(2));
t646 = cos(qJ(2));
t601 = t643 * t622 + t646 * t623;
t647 = qJD(2) ^ 2;
t596 = -t647 * pkin(2) + qJDD(2) * pkin(6) + t601;
t636 = -g(3) + qJDD(1);
t642 = sin(qJ(3));
t645 = cos(qJ(3));
t584 = -t642 * t596 + t645 * t636;
t663 = qJD(2) * qJD(3);
t661 = t645 * t663;
t620 = t642 * qJDD(2) + t661;
t579 = (-t620 + t661) * qJ(4) + (t642 * t645 * t647 + qJDD(3)) * pkin(3) + t584;
t585 = t645 * t596 + t642 * t636;
t621 = t645 * qJDD(2) - t642 * t663;
t665 = qJD(2) * t642;
t624 = qJD(3) * pkin(3) - qJ(4) * t665;
t635 = t645 ^ 2;
t580 = -t635 * t647 * pkin(3) + t621 * qJ(4) - qJD(3) * t624 + t585;
t637 = sin(pkin(9));
t639 = cos(pkin(9));
t609 = (t637 * t645 + t639 * t642) * qJD(2);
t559 = -0.2e1 * qJD(4) * t609 + t639 * t579 - t637 * t580;
t598 = t639 * t620 + t637 * t621;
t608 = (-t637 * t642 + t639 * t645) * qJD(2);
t557 = (qJD(3) * t608 - t598) * pkin(7) + (t608 * t609 + qJDD(3)) * pkin(4) + t559;
t560 = 0.2e1 * qJD(4) * t608 + t637 * t579 + t639 * t580;
t597 = -t637 * t620 + t639 * t621;
t604 = qJD(3) * pkin(4) - t609 * pkin(7);
t607 = t608 ^ 2;
t558 = -t607 * pkin(4) + t597 * pkin(7) - qJD(3) * t604 + t560;
t641 = sin(qJ(5));
t644 = cos(qJ(5));
t555 = t644 * t557 - t641 * t558;
t589 = t644 * t608 - t641 * t609;
t569 = t589 * qJD(5) + t641 * t597 + t644 * t598;
t590 = t641 * t608 + t644 * t609;
t575 = -t589 * mrSges(6,1) + t590 * mrSges(6,2);
t633 = qJD(3) + qJD(5);
t582 = -t633 * mrSges(6,2) + t589 * mrSges(6,3);
t632 = qJDD(3) + qJDD(5);
t551 = m(6) * t555 + t632 * mrSges(6,1) - t569 * mrSges(6,3) - t590 * t575 + t633 * t582;
t556 = t641 * t557 + t644 * t558;
t568 = -t590 * qJD(5) + t644 * t597 - t641 * t598;
t583 = t633 * mrSges(6,1) - t590 * mrSges(6,3);
t552 = m(6) * t556 - t632 * mrSges(6,2) + t568 * mrSges(6,3) + t589 * t575 - t633 * t583;
t542 = t644 * t551 + t641 * t552;
t592 = -t608 * mrSges(5,1) + t609 * mrSges(5,2);
t602 = -qJD(3) * mrSges(5,2) + t608 * mrSges(5,3);
t540 = m(5) * t559 + qJDD(3) * mrSges(5,1) - t598 * mrSges(5,3) + qJD(3) * t602 - t609 * t592 + t542;
t603 = qJD(3) * mrSges(5,1) - t609 * mrSges(5,3);
t656 = -t641 * t551 + t644 * t552;
t541 = m(5) * t560 - qJDD(3) * mrSges(5,2) + t597 * mrSges(5,3) - qJD(3) * t603 + t608 * t592 + t656;
t536 = t639 * t540 + t637 * t541;
t587 = Ifges(5,4) * t609 + Ifges(5,2) * t608 + Ifges(5,6) * qJD(3);
t588 = Ifges(5,1) * t609 + Ifges(5,4) * t608 + Ifges(5,5) * qJD(3);
t611 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t642 + Ifges(4,2) * t645) * qJD(2);
t612 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t642 + Ifges(4,4) * t645) * qJD(2);
t571 = Ifges(6,4) * t590 + Ifges(6,2) * t589 + Ifges(6,6) * t633;
t572 = Ifges(6,1) * t590 + Ifges(6,4) * t589 + Ifges(6,5) * t633;
t650 = -mrSges(6,1) * t555 + mrSges(6,2) * t556 - Ifges(6,5) * t569 - Ifges(6,6) * t568 - Ifges(6,3) * t632 - t590 * t571 + t589 * t572;
t668 = mrSges(4,1) * t584 + mrSges(5,1) * t559 - mrSges(4,2) * t585 - mrSges(5,2) * t560 + Ifges(4,5) * t620 + Ifges(5,5) * t598 + Ifges(4,6) * t621 + Ifges(5,6) * t597 + pkin(3) * t536 + pkin(4) * t542 + (t642 * t611 - t645 * t612) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t609 * t587 - t608 * t588 - t650;
t619 = (-mrSges(4,1) * t645 + mrSges(4,2) * t642) * qJD(2);
t664 = qJD(2) * t645;
t626 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t664;
t534 = m(4) * t584 + qJDD(3) * mrSges(4,1) - t620 * mrSges(4,3) + qJD(3) * t626 - t619 * t665 + t536;
t625 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t665;
t657 = -t637 * t540 + t639 * t541;
t535 = m(4) * t585 - qJDD(3) * mrSges(4,2) + t621 * mrSges(4,3) - qJD(3) * t625 + t619 * t664 + t657;
t658 = -t642 * t534 + t645 * t535;
t526 = m(3) * t601 - t647 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t658;
t600 = t646 * t622 - t643 * t623;
t652 = -qJDD(2) * pkin(2) - t600;
t581 = -t621 * pkin(3) + qJDD(4) + t624 * t665 + (-qJ(4) * t635 - pkin(6)) * t647 + t652;
t562 = -t597 * pkin(4) - t607 * pkin(7) + t609 * t604 + t581;
t654 = m(6) * t562 - t568 * mrSges(6,1) + t569 * mrSges(6,2) - t589 * t582 + t590 * t583;
t553 = m(5) * t581 - t597 * mrSges(5,1) + t598 * mrSges(5,2) - t608 * t602 + t609 * t603 + t654;
t595 = -t647 * pkin(6) + t652;
t649 = -m(4) * t595 + t621 * mrSges(4,1) - t620 * mrSges(4,2) - t625 * t665 + t626 * t664 - t553;
t546 = m(3) * t600 + qJDD(2) * mrSges(3,1) - t647 * mrSges(3,2) + t649;
t523 = t643 * t526 + t646 * t546;
t521 = m(2) * t622 + t523;
t659 = t646 * t526 - t643 * t546;
t522 = m(2) * t623 + t659;
t666 = t640 * t521 + t638 * t522;
t528 = t645 * t534 + t642 * t535;
t662 = m(3) * t636 + t528;
t660 = -t638 * t521 + t640 * t522;
t655 = m(2) * t636 + t662;
t570 = Ifges(6,5) * t590 + Ifges(6,6) * t589 + Ifges(6,3) * t633;
t543 = -mrSges(6,1) * t562 + mrSges(6,3) * t556 + Ifges(6,4) * t569 + Ifges(6,2) * t568 + Ifges(6,6) * t632 - t590 * t570 + t633 * t572;
t544 = mrSges(6,2) * t562 - mrSges(6,3) * t555 + Ifges(6,1) * t569 + Ifges(6,4) * t568 + Ifges(6,5) * t632 + t589 * t570 - t633 * t571;
t586 = Ifges(5,5) * t609 + Ifges(5,6) * t608 + Ifges(5,3) * qJD(3);
t529 = -mrSges(5,1) * t581 + mrSges(5,3) * t560 + Ifges(5,4) * t598 + Ifges(5,2) * t597 + Ifges(5,6) * qJDD(3) - pkin(4) * t654 + pkin(7) * t656 + qJD(3) * t588 + t644 * t543 + t641 * t544 - t609 * t586;
t530 = mrSges(5,2) * t581 - mrSges(5,3) * t559 + Ifges(5,1) * t598 + Ifges(5,4) * t597 + Ifges(5,5) * qJDD(3) - pkin(7) * t542 - qJD(3) * t587 - t641 * t543 + t644 * t544 + t608 * t586;
t610 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t642 + Ifges(4,6) * t645) * qJD(2);
t515 = -mrSges(4,1) * t595 + mrSges(4,3) * t585 + Ifges(4,4) * t620 + Ifges(4,2) * t621 + Ifges(4,6) * qJDD(3) - pkin(3) * t553 + qJ(4) * t657 + qJD(3) * t612 + t639 * t529 + t637 * t530 - t610 * t665;
t517 = mrSges(4,2) * t595 - mrSges(4,3) * t584 + Ifges(4,1) * t620 + Ifges(4,4) * t621 + Ifges(4,5) * qJDD(3) - qJ(4) * t536 - qJD(3) * t611 - t637 * t529 + t639 * t530 + t610 * t664;
t651 = mrSges(3,1) * t600 - mrSges(3,2) * t601 + Ifges(3,3) * qJDD(2) + pkin(2) * t649 + pkin(6) * t658 + t645 * t515 + t642 * t517;
t513 = -mrSges(3,1) * t636 + mrSges(3,3) * t601 + t647 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t528 - t668;
t512 = mrSges(3,2) * t636 - mrSges(3,3) * t600 + Ifges(3,5) * qJDD(2) - t647 * Ifges(3,6) - pkin(6) * t528 - t642 * t515 + t645 * t517;
t511 = mrSges(2,2) * t636 - mrSges(2,3) * t622 - pkin(5) * t523 + t646 * t512 - t643 * t513;
t510 = -mrSges(2,1) * t636 + mrSges(2,3) * t623 - pkin(1) * t662 + pkin(5) * t659 + t643 * t512 + t646 * t513;
t1 = [-m(1) * g(1) + t660; -m(1) * g(2) + t666; -m(1) * g(3) + t655; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t666 - t638 * t510 + t640 * t511; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t660 + t640 * t510 + t638 * t511; -mrSges(1,1) * g(2) + mrSges(2,1) * t622 + mrSges(1,2) * g(1) - mrSges(2,2) * t623 + pkin(1) * t523 + t651; t655; t651; t668; t553; -t650;];
tauJB = t1;
