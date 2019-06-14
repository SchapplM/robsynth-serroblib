% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:23:01
% EndTime: 2019-05-05 18:23:14
% DurationCPUTime: 11.86s
% Computational Cost: add. (189660->343), mult. (411627->433), div. (0->0), fcn. (280021->12), ass. (0->134)
t656 = sin(qJ(1));
t660 = cos(qJ(1));
t640 = g(1) * t656 - g(2) * t660;
t632 = qJDD(1) * pkin(1) + t640;
t641 = -g(1) * t660 - g(2) * t656;
t662 = qJD(1) ^ 2;
t634 = -pkin(1) * t662 + t641;
t650 = sin(pkin(10));
t652 = cos(pkin(10));
t608 = t632 * t650 + t634 * t652;
t600 = -pkin(2) * t662 + qJDD(1) * pkin(7) + t608;
t648 = -g(3) + qJDD(2);
t655 = sin(qJ(3));
t659 = cos(qJ(3));
t590 = -t655 * t600 + t648 * t659;
t677 = qJD(1) * qJD(3);
t675 = t659 * t677;
t635 = qJDD(1) * t655 + t675;
t574 = (-t635 + t675) * qJ(4) + (t655 * t659 * t662 + qJDD(3)) * pkin(3) + t590;
t591 = t600 * t659 + t648 * t655;
t636 = qJDD(1) * t659 - t655 * t677;
t680 = qJD(1) * t655;
t637 = qJD(3) * pkin(3) - qJ(4) * t680;
t647 = t659 ^ 2;
t577 = -pkin(3) * t647 * t662 + qJ(4) * t636 - qJD(3) * t637 + t591;
t649 = sin(pkin(11));
t651 = cos(pkin(11));
t621 = (t649 * t659 + t651 * t655) * qJD(1);
t561 = -0.2e1 * qJD(4) * t621 + t574 * t651 - t577 * t649;
t679 = qJD(1) * t659;
t620 = -t649 * t680 + t651 * t679;
t562 = 0.2e1 * qJD(4) * t620 + t574 * t649 + t577 * t651;
t602 = -mrSges(5,1) * t620 + mrSges(5,2) * t621;
t609 = -t635 * t649 + t636 * t651;
t615 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t621;
t603 = -pkin(4) * t620 - pkin(8) * t621;
t661 = qJD(3) ^ 2;
t557 = -pkin(4) * t661 + qJDD(3) * pkin(8) + t603 * t620 + t562;
t607 = t652 * t632 - t634 * t650;
t667 = -qJDD(1) * pkin(2) - t607;
t578 = -t636 * pkin(3) + qJDD(4) + t637 * t680 + (-qJ(4) * t647 - pkin(7)) * t662 + t667;
t610 = t635 * t651 + t636 * t649;
t565 = (-qJD(3) * t620 - t610) * pkin(8) + (qJD(3) * t621 - t609) * pkin(4) + t578;
t654 = sin(qJ(5));
t658 = cos(qJ(5));
t552 = -t654 * t557 + t565 * t658;
t612 = qJD(3) * t658 - t621 * t654;
t585 = qJD(5) * t612 + qJDD(3) * t654 + t610 * t658;
t606 = qJDD(5) - t609;
t613 = qJD(3) * t654 + t621 * t658;
t619 = qJD(5) - t620;
t550 = (t612 * t619 - t585) * pkin(9) + (t612 * t613 + t606) * pkin(5) + t552;
t553 = t557 * t658 + t565 * t654;
t584 = -qJD(5) * t613 + qJDD(3) * t658 - t610 * t654;
t594 = pkin(5) * t619 - pkin(9) * t613;
t611 = t612 ^ 2;
t551 = -pkin(5) * t611 + pkin(9) * t584 - t594 * t619 + t553;
t653 = sin(qJ(6));
t657 = cos(qJ(6));
t548 = t550 * t657 - t551 * t653;
t586 = t612 * t657 - t613 * t653;
t560 = qJD(6) * t586 + t584 * t653 + t585 * t657;
t587 = t612 * t653 + t613 * t657;
t570 = -mrSges(7,1) * t586 + mrSges(7,2) * t587;
t616 = qJD(6) + t619;
t575 = -mrSges(7,2) * t616 + mrSges(7,3) * t586;
t604 = qJDD(6) + t606;
t546 = m(7) * t548 + mrSges(7,1) * t604 - mrSges(7,3) * t560 - t570 * t587 + t575 * t616;
t549 = t550 * t653 + t551 * t657;
t559 = -qJD(6) * t587 + t584 * t657 - t585 * t653;
t576 = mrSges(7,1) * t616 - mrSges(7,3) * t587;
t547 = m(7) * t549 - mrSges(7,2) * t604 + mrSges(7,3) * t559 + t570 * t586 - t576 * t616;
t538 = t546 * t657 + t547 * t653;
t588 = -mrSges(6,1) * t612 + mrSges(6,2) * t613;
t592 = -mrSges(6,2) * t619 + mrSges(6,3) * t612;
t536 = m(6) * t552 + mrSges(6,1) * t606 - mrSges(6,3) * t585 - t588 * t613 + t592 * t619 + t538;
t593 = mrSges(6,1) * t619 - mrSges(6,3) * t613;
t669 = -t546 * t653 + t547 * t657;
t537 = m(6) * t553 - mrSges(6,2) * t606 + mrSges(6,3) * t584 + t588 * t612 - t593 * t619 + t669;
t670 = -t536 * t654 + t537 * t658;
t531 = m(5) * t562 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t609 - qJD(3) * t615 + t602 * t620 + t670;
t614 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t620;
t556 = -qJDD(3) * pkin(4) - pkin(8) * t661 + t603 * t621 - t561;
t554 = -pkin(5) * t584 - pkin(9) * t611 + t594 * t613 + t556;
t666 = m(7) * t554 - mrSges(7,1) * t559 + mrSges(7,2) * t560 - t575 * t586 + t576 * t587;
t664 = -m(6) * t556 + mrSges(6,1) * t584 - mrSges(6,2) * t585 + t592 * t612 - t593 * t613 - t666;
t542 = m(5) * t561 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t610 + qJD(3) * t614 - t602 * t621 + t664;
t524 = t531 * t649 + t542 * t651;
t633 = (-mrSges(4,1) * t659 + mrSges(4,2) * t655) * qJD(1);
t639 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t679;
t522 = m(4) * t590 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t635 + qJD(3) * t639 - t633 * t680 + t524;
t638 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t680;
t671 = t531 * t651 - t542 * t649;
t523 = m(4) * t591 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t636 - qJD(3) * t638 + t633 * t679 + t671;
t672 = -t522 * t655 + t523 * t659;
t516 = m(3) * t608 - mrSges(3,1) * t662 - qJDD(1) * mrSges(3,2) + t672;
t599 = -t662 * pkin(7) + t667;
t532 = t536 * t658 + t537 * t654;
t665 = m(5) * t578 - mrSges(5,1) * t609 + mrSges(5,2) * t610 - t614 * t620 + t615 * t621 + t532;
t663 = -m(4) * t599 + mrSges(4,1) * t636 - mrSges(4,2) * t635 - t638 * t680 + t639 * t679 - t665;
t528 = m(3) * t607 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t662 + t663;
t512 = t516 * t650 + t528 * t652;
t510 = m(2) * t640 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t662 + t512;
t673 = t516 * t652 - t528 * t650;
t511 = m(2) * t641 - mrSges(2,1) * t662 - qJDD(1) * mrSges(2,2) + t673;
t681 = t510 * t660 + t511 * t656;
t517 = t522 * t659 + t523 * t655;
t676 = m(3) * t648 + t517;
t674 = -t510 * t656 + t511 * t660;
t627 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t655 + Ifges(4,4) * t659) * qJD(1);
t626 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t655 + Ifges(4,2) * t659) * qJD(1);
t625 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t655 + Ifges(4,6) * t659) * qJD(1);
t598 = Ifges(5,1) * t621 + Ifges(5,4) * t620 + Ifges(5,5) * qJD(3);
t597 = Ifges(5,4) * t621 + Ifges(5,2) * t620 + Ifges(5,6) * qJD(3);
t596 = Ifges(5,5) * t621 + Ifges(5,6) * t620 + Ifges(5,3) * qJD(3);
t581 = Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t619;
t580 = Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t619;
t579 = Ifges(6,5) * t613 + Ifges(6,6) * t612 + Ifges(6,3) * t619;
t568 = Ifges(7,1) * t587 + Ifges(7,4) * t586 + Ifges(7,5) * t616;
t567 = Ifges(7,4) * t587 + Ifges(7,2) * t586 + Ifges(7,6) * t616;
t566 = Ifges(7,5) * t587 + Ifges(7,6) * t586 + Ifges(7,3) * t616;
t540 = mrSges(7,2) * t554 - mrSges(7,3) * t548 + Ifges(7,1) * t560 + Ifges(7,4) * t559 + Ifges(7,5) * t604 + t566 * t586 - t567 * t616;
t539 = -mrSges(7,1) * t554 + mrSges(7,3) * t549 + Ifges(7,4) * t560 + Ifges(7,2) * t559 + Ifges(7,6) * t604 - t566 * t587 + t568 * t616;
t526 = mrSges(6,2) * t556 - mrSges(6,3) * t552 + Ifges(6,1) * t585 + Ifges(6,4) * t584 + Ifges(6,5) * t606 - pkin(9) * t538 - t539 * t653 + t540 * t657 + t579 * t612 - t580 * t619;
t525 = -mrSges(6,1) * t556 + mrSges(6,3) * t553 + Ifges(6,4) * t585 + Ifges(6,2) * t584 + Ifges(6,6) * t606 - pkin(5) * t666 + pkin(9) * t669 + t657 * t539 + t653 * t540 - t613 * t579 + t619 * t581;
t518 = Ifges(5,4) * t610 + Ifges(5,2) * t609 + Ifges(5,6) * qJDD(3) - t621 * t596 + qJD(3) * t598 - mrSges(5,1) * t578 + mrSges(5,3) * t562 - Ifges(6,5) * t585 - Ifges(6,6) * t584 - Ifges(6,3) * t606 - t613 * t580 + t612 * t581 - mrSges(6,1) * t552 + mrSges(6,2) * t553 - Ifges(7,5) * t560 - Ifges(7,6) * t559 - Ifges(7,3) * t604 - t587 * t567 + t586 * t568 - mrSges(7,1) * t548 + mrSges(7,2) * t549 - pkin(5) * t538 - pkin(4) * t532;
t513 = mrSges(5,2) * t578 - mrSges(5,3) * t561 + Ifges(5,1) * t610 + Ifges(5,4) * t609 + Ifges(5,5) * qJDD(3) - pkin(8) * t532 - qJD(3) * t597 - t525 * t654 + t526 * t658 + t596 * t620;
t506 = mrSges(4,2) * t599 - mrSges(4,3) * t590 + Ifges(4,1) * t635 + Ifges(4,4) * t636 + Ifges(4,5) * qJDD(3) - qJ(4) * t524 - qJD(3) * t626 + t513 * t651 - t518 * t649 + t625 * t679;
t505 = -pkin(2) * t517 - mrSges(3,1) * t648 + mrSges(3,3) * t608 - pkin(3) * t524 - Ifges(4,5) * t635 - Ifges(4,6) * t636 - mrSges(4,1) * t590 + mrSges(4,2) * t591 - t654 * t526 - t658 * t525 - pkin(4) * t664 - pkin(8) * t670 - Ifges(5,5) * t610 - Ifges(5,6) * t609 - mrSges(5,1) * t561 + mrSges(5,2) * t562 + t662 * Ifges(3,5) - t621 * t597 + t620 * t598 + Ifges(3,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t626 * t655 + t627 * t659) * qJD(1);
t504 = -mrSges(4,1) * t599 + mrSges(4,3) * t591 + Ifges(4,4) * t635 + Ifges(4,2) * t636 + Ifges(4,6) * qJDD(3) - pkin(3) * t665 + qJ(4) * t671 + qJD(3) * t627 + t649 * t513 + t651 * t518 - t625 * t680;
t503 = mrSges(3,2) * t648 - mrSges(3,3) * t607 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t662 - pkin(7) * t517 - t504 * t655 + t506 * t659;
t502 = -mrSges(2,2) * g(3) - mrSges(2,3) * t640 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t662 - qJ(2) * t512 + t503 * t652 - t505 * t650;
t501 = mrSges(2,1) * g(3) + mrSges(2,3) * t641 + t662 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t676 + qJ(2) * t673 + t650 * t503 + t652 * t505;
t1 = [-m(1) * g(1) + t674; -m(1) * g(2) + t681; (-m(1) - m(2)) * g(3) + t676; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t681 - t501 * t656 + t502 * t660; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t674 + t660 * t501 + t656 * t502; pkin(1) * t512 + mrSges(2,1) * t640 - mrSges(2,2) * t641 + pkin(2) * t663 + pkin(7) * t672 + t655 * t506 + t659 * t504 + mrSges(3,1) * t607 - mrSges(3,2) * t608 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
