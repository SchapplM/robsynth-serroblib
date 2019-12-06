% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:49
% EndTime: 2019-12-05 16:35:58
% DurationCPUTime: 6.25s
% Computational Cost: add. (69716->269), mult. (142327->350), div. (0->0), fcn. (96873->12), ass. (0->121)
t645 = sin(pkin(9));
t648 = cos(pkin(9));
t635 = t645 * g(1) - t648 * g(2);
t636 = -t648 * g(1) - t645 * g(2);
t643 = -g(3) + qJDD(1);
t646 = sin(pkin(5));
t649 = cos(pkin(5));
t652 = sin(qJ(2));
t655 = cos(qJ(2));
t601 = -t652 * t636 + (t635 * t649 + t643 * t646) * t655;
t679 = t649 * t652;
t680 = t646 * t652;
t602 = t635 * t679 + t655 * t636 + t643 * t680;
t657 = qJD(2) ^ 2;
t598 = -t657 * pkin(2) + qJDD(2) * pkin(7) + t602;
t617 = -t646 * t635 + t649 * t643;
t651 = sin(qJ(3));
t654 = cos(qJ(3));
t586 = t654 * t598 + t651 * t617;
t631 = (-pkin(3) * t654 - qJ(4) * t651) * qJD(2);
t656 = qJD(3) ^ 2;
t675 = t654 * qJD(2);
t582 = -t656 * pkin(3) + qJDD(3) * qJ(4) + t631 * t675 + t586;
t597 = -qJDD(2) * pkin(2) - t657 * pkin(7) - t601;
t674 = qJD(2) * qJD(3);
t672 = t654 * t674;
t633 = t651 * qJDD(2) + t672;
t673 = t651 * t674;
t634 = t654 * qJDD(2) - t673;
t584 = (-t633 - t672) * qJ(4) + (-t634 + t673) * pkin(3) + t597;
t644 = sin(pkin(10));
t647 = cos(pkin(10));
t676 = qJD(2) * t651;
t626 = t647 * qJD(3) - t644 * t676;
t683 = 2 * qJD(4);
t578 = t647 * t582 + t644 * t584 + t626 * t683;
t627 = t644 * qJD(3) + t647 * t676;
t608 = -t626 * pkin(4) - t627 * pkin(8);
t681 = t654 ^ 2 * t657;
t576 = -pkin(4) * t681 - t634 * pkin(8) + t626 * t608 + t578;
t615 = t647 * qJDD(3) - t644 * t633;
t616 = t644 * qJDD(3) + t647 * t633;
t595 = t651 * t598;
t662 = -qJDD(3) * pkin(3) - t656 * qJ(4) + t631 * t676 + qJDD(4) + t595;
t579 = -t615 * pkin(4) - t616 * pkin(8) + (-t617 + (-pkin(4) * t627 + pkin(8) * t626) * qJD(2)) * t654 + t662;
t650 = sin(qJ(5));
t653 = cos(qJ(5));
t573 = -t650 * t576 + t653 * t579;
t609 = -t650 * t627 - t653 * t675;
t593 = t609 * qJD(5) + t653 * t616 - t650 * t634;
t610 = t653 * t627 - t650 * t675;
t594 = -t609 * mrSges(6,1) + t610 * mrSges(6,2);
t624 = qJD(5) - t626;
t599 = -t624 * mrSges(6,2) + t609 * mrSges(6,3);
t613 = qJDD(5) - t615;
t571 = m(6) * t573 + t613 * mrSges(6,1) - t593 * mrSges(6,3) - t610 * t594 + t624 * t599;
t574 = t653 * t576 + t650 * t579;
t592 = -t610 * qJD(5) - t650 * t616 - t653 * t634;
t600 = t624 * mrSges(6,1) - t610 * mrSges(6,3);
t572 = m(6) * t574 - t613 * mrSges(6,2) + t592 * mrSges(6,3) + t609 * t594 - t624 * t600;
t566 = t653 * t571 + t650 * t572;
t667 = t644 * t582 - t647 * t584;
t575 = -pkin(8) * t681 + t634 * pkin(4) + (t683 + t608) * t627 + t667;
t587 = Ifges(6,5) * t610 + Ifges(6,6) * t609 + Ifges(6,3) * t624;
t589 = Ifges(6,1) * t610 + Ifges(6,4) * t609 + Ifges(6,5) * t624;
t567 = -mrSges(6,1) * t575 + mrSges(6,3) * t574 + Ifges(6,4) * t593 + Ifges(6,2) * t592 + Ifges(6,6) * t613 - t610 * t587 + t624 * t589;
t588 = Ifges(6,4) * t610 + Ifges(6,2) * t609 + Ifges(6,6) * t624;
t568 = mrSges(6,2) * t575 - mrSges(6,3) * t573 + Ifges(6,1) * t593 + Ifges(6,4) * t592 + Ifges(6,5) * t613 + t609 * t587 - t624 * t588;
t577 = -0.2e1 * qJD(4) * t627 - t667;
t678 = t654 * t617;
t581 = t662 - t678;
t603 = Ifges(5,5) * t627 + Ifges(5,6) * t626 - Ifges(5,3) * t675;
t604 = Ifges(5,4) * t627 + Ifges(5,2) * t626 - Ifges(5,6) * t675;
t551 = mrSges(5,2) * t581 - mrSges(5,3) * t577 + Ifges(5,1) * t616 + Ifges(5,4) * t615 - Ifges(5,5) * t634 - pkin(8) * t566 - t650 * t567 + t653 * t568 + t626 * t603 + t604 * t675;
t605 = Ifges(5,1) * t627 + Ifges(5,4) * t626 - Ifges(5,5) * t675;
t659 = mrSges(6,1) * t573 - mrSges(6,2) * t574 + Ifges(6,5) * t593 + Ifges(6,6) * t592 + Ifges(6,3) * t613 + t610 * t588 - t609 * t589;
t552 = -mrSges(5,1) * t581 + mrSges(5,3) * t578 + Ifges(5,4) * t616 + Ifges(5,2) * t615 - Ifges(5,6) * t634 - pkin(4) * t566 - t627 * t603 - t605 * t675 - t659;
t607 = -t626 * mrSges(5,1) + t627 * mrSges(5,2);
t614 = -mrSges(5,1) * t675 - t627 * mrSges(5,3);
t669 = -t650 * t571 + t653 * t572;
t564 = m(5) * t578 + t634 * mrSges(5,2) + t615 * mrSges(5,3) + t626 * t607 + t614 * t675 + t669;
t661 = -m(6) * t575 + t592 * mrSges(6,1) - t593 * mrSges(6,2) + t609 * t599 - t610 * t600;
t664 = mrSges(5,2) * t675 + t626 * mrSges(5,3);
t569 = m(5) * t577 - t634 * mrSges(5,1) - t616 * mrSges(5,3) - t627 * t607 - t664 * t675 + t661;
t560 = t647 * t564 - t644 * t569;
t565 = m(5) * t581 - t615 * mrSges(5,1) + t616 * mrSges(5,2) + t627 * t614 - t626 * t664 + t566;
t585 = -t595 + t678;
t622 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t651 + Ifges(4,2) * t654) * qJD(2);
t623 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t651 + Ifges(4,4) * t654) * qJD(2);
t684 = mrSges(4,1) * t585 - mrSges(4,2) * t586 + Ifges(4,5) * t633 + Ifges(4,6) * t634 + Ifges(4,3) * qJDD(3) - pkin(3) * t565 + qJ(4) * t560 + t644 * t551 + t647 * t552 + (t651 * t622 - t654 * t623) * qJD(2);
t559 = t644 * t564 + t647 * t569;
t637 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t676;
t638 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t675;
t660 = -m(4) * t597 + t634 * mrSges(4,1) - t633 * mrSges(4,2) - t637 * t676 + t638 * t675 - t559;
t555 = m(3) * t601 + qJDD(2) * mrSges(3,1) - t657 * mrSges(3,2) + t660;
t682 = t555 * t655;
t632 = (-mrSges(4,1) * t654 + mrSges(4,2) * t651) * qJD(2);
t558 = m(4) * t586 - qJDD(3) * mrSges(4,2) + t634 * mrSges(4,3) - qJD(3) * t637 + t632 * t675 + t560;
t562 = m(4) * t585 + qJDD(3) * mrSges(4,1) - t633 * mrSges(4,3) + qJD(3) * t638 - t632 * t676 - t565;
t670 = t654 * t558 - t651 * t562;
t547 = m(3) * t602 - t657 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t670;
t550 = t651 * t558 + t654 * t562;
t549 = m(3) * t617 + t550;
t537 = t547 * t679 - t646 * t549 + t649 * t682;
t535 = m(2) * t635 + t537;
t543 = t655 * t547 - t652 * t555;
t542 = m(2) * t636 + t543;
t677 = t648 * t535 + t645 * t542;
t536 = t547 * t680 + t649 * t549 + t646 * t682;
t671 = -t645 * t535 + t648 * t542;
t668 = m(2) * t643 + t536;
t621 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t651 + Ifges(4,6) * t654) * qJD(2);
t538 = mrSges(4,2) * t597 - mrSges(4,3) * t585 + Ifges(4,1) * t633 + Ifges(4,4) * t634 + Ifges(4,5) * qJDD(3) - qJ(4) * t559 - qJD(3) * t622 + t647 * t551 - t644 * t552 + t621 * t675;
t539 = Ifges(4,4) * t633 + Ifges(4,6) * qJDD(3) - t621 * t676 + qJD(3) * t623 - mrSges(4,1) * t597 + mrSges(4,3) * t586 - Ifges(5,5) * t616 - Ifges(5,6) * t615 - t627 * t604 + t626 * t605 - mrSges(5,1) * t577 + mrSges(5,2) * t578 - t650 * t568 - t653 * t567 - pkin(4) * t661 - pkin(8) * t669 - pkin(3) * t559 + (Ifges(4,2) + Ifges(5,3)) * t634;
t532 = mrSges(3,2) * t617 - mrSges(3,3) * t601 + Ifges(3,5) * qJDD(2) - t657 * Ifges(3,6) - pkin(7) * t550 + t654 * t538 - t651 * t539;
t533 = -mrSges(3,1) * t617 + mrSges(3,3) * t602 + t657 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t550 - t684;
t663 = pkin(6) * t543 + t532 * t652 + t533 * t655;
t531 = mrSges(3,1) * t601 - mrSges(3,2) * t602 + Ifges(3,3) * qJDD(2) + pkin(2) * t660 + pkin(7) * t670 + t651 * t538 + t654 * t539;
t530 = mrSges(2,2) * t643 - mrSges(2,3) * t635 + t655 * t532 - t652 * t533 + (-t536 * t646 - t537 * t649) * pkin(6);
t529 = -mrSges(2,1) * t643 + mrSges(2,3) * t636 - pkin(1) * t536 - t646 * t531 + t663 * t649;
t1 = [-m(1) * g(1) + t671; -m(1) * g(2) + t677; -m(1) * g(3) + t668; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t677 - t645 * t529 + t648 * t530; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t671 + t648 * t529 + t645 * t530; -mrSges(1,1) * g(2) + mrSges(2,1) * t635 + mrSges(1,2) * g(1) - mrSges(2,2) * t636 + pkin(1) * t537 + t649 * t531 + t663 * t646; t668; t531; t684; t565; t659;];
tauJB = t1;
