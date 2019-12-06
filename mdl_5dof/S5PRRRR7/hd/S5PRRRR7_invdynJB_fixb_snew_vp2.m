% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:16
% EndTime: 2019-12-05 17:12:22
% DurationCPUTime: 5.20s
% Computational Cost: add. (59711->262), mult. (123825->334), div. (0->0), fcn. (84128->10), ass. (0->110)
t642 = sin(pkin(9));
t670 = cos(pkin(9));
t627 = -t670 * g(1) - t642 * g(2);
t641 = -g(3) + qJDD(1);
t646 = sin(qJ(2));
t650 = cos(qJ(2));
t609 = t650 * t627 + t646 * t641;
t651 = qJD(2) ^ 2;
t604 = -t651 * pkin(2) + qJDD(2) * pkin(6) + t609;
t626 = t642 * g(1) - t670 * g(2);
t645 = sin(qJ(3));
t649 = cos(qJ(3));
t593 = -t645 * t604 - t649 * t626;
t666 = qJD(2) * qJD(3);
t664 = t649 * t666;
t624 = t645 * qJDD(2) + t664;
t582 = (-t624 + t664) * pkin(7) + (t645 * t649 * t651 + qJDD(3)) * pkin(3) + t593;
t594 = t649 * t604 - t645 * t626;
t625 = t649 * qJDD(2) - t645 * t666;
t668 = qJD(2) * t645;
t630 = qJD(3) * pkin(3) - pkin(7) * t668;
t640 = t649 ^ 2;
t583 = -t640 * t651 * pkin(3) + t625 * pkin(7) - qJD(3) * t630 + t594;
t644 = sin(qJ(4));
t648 = cos(qJ(4));
t564 = t648 * t582 - t644 * t583;
t614 = (-t644 * t645 + t648 * t649) * qJD(2);
t590 = t614 * qJD(4) + t648 * t624 + t644 * t625;
t615 = (t644 * t649 + t645 * t648) * qJD(2);
t638 = qJDD(3) + qJDD(4);
t639 = qJD(3) + qJD(4);
t559 = (t614 * t639 - t590) * pkin(8) + (t614 * t615 + t638) * pkin(4) + t564;
t565 = t644 * t582 + t648 * t583;
t589 = -t615 * qJD(4) - t644 * t624 + t648 * t625;
t607 = t639 * pkin(4) - t615 * pkin(8);
t610 = t614 ^ 2;
t560 = -t610 * pkin(4) + t589 * pkin(8) - t639 * t607 + t565;
t643 = sin(qJ(5));
t647 = cos(qJ(5));
t557 = t647 * t559 - t643 * t560;
t599 = t647 * t614 - t643 * t615;
t571 = t599 * qJD(5) + t643 * t589 + t647 * t590;
t600 = t643 * t614 + t647 * t615;
t578 = -t599 * mrSges(6,1) + t600 * mrSges(6,2);
t635 = qJD(5) + t639;
t591 = -t635 * mrSges(6,2) + t599 * mrSges(6,3);
t634 = qJDD(5) + t638;
t554 = m(6) * t557 + t634 * mrSges(6,1) - t571 * mrSges(6,3) - t600 * t578 + t635 * t591;
t558 = t643 * t559 + t647 * t560;
t570 = -t600 * qJD(5) + t647 * t589 - t643 * t590;
t592 = t635 * mrSges(6,1) - t600 * mrSges(6,3);
t555 = m(6) * t558 - t634 * mrSges(6,2) + t570 * mrSges(6,3) + t599 * t578 - t635 * t592;
t545 = t647 * t554 + t643 * t555;
t601 = -t614 * mrSges(5,1) + t615 * mrSges(5,2);
t605 = -t639 * mrSges(5,2) + t614 * mrSges(5,3);
t542 = m(5) * t564 + t638 * mrSges(5,1) - t590 * mrSges(5,3) - t615 * t601 + t639 * t605 + t545;
t606 = t639 * mrSges(5,1) - t615 * mrSges(5,3);
t660 = -t643 * t554 + t647 * t555;
t543 = m(5) * t565 - t638 * mrSges(5,2) + t589 * mrSges(5,3) + t614 * t601 - t639 * t606 + t660;
t538 = t648 * t542 + t644 * t543;
t612 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t645 + Ifges(4,2) * t649) * qJD(2);
t613 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t645 + Ifges(4,4) * t649) * qJD(2);
t596 = Ifges(5,4) * t615 + Ifges(5,2) * t614 + Ifges(5,6) * t639;
t597 = Ifges(5,1) * t615 + Ifges(5,4) * t614 + Ifges(5,5) * t639;
t574 = Ifges(6,4) * t600 + Ifges(6,2) * t599 + Ifges(6,6) * t635;
t575 = Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t635;
t656 = -mrSges(6,1) * t557 + mrSges(6,2) * t558 - Ifges(6,5) * t571 - Ifges(6,6) * t570 - Ifges(6,3) * t634 - t600 * t574 + t599 * t575;
t653 = -mrSges(5,1) * t564 + mrSges(5,2) * t565 - Ifges(5,5) * t590 - Ifges(5,6) * t589 - Ifges(5,3) * t638 - pkin(4) * t545 - t615 * t596 + t614 * t597 + t656;
t671 = mrSges(4,1) * t593 - mrSges(4,2) * t594 + Ifges(4,5) * t624 + Ifges(4,6) * t625 + Ifges(4,3) * qJDD(3) + pkin(3) * t538 + (t645 * t612 - t649 * t613) * qJD(2) - t653;
t623 = (-mrSges(4,1) * t649 + mrSges(4,2) * t645) * qJD(2);
t667 = qJD(2) * t649;
t629 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t667;
t536 = m(4) * t593 + qJDD(3) * mrSges(4,1) - t624 * mrSges(4,3) + qJD(3) * t629 - t623 * t668 + t538;
t628 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t668;
t661 = -t644 * t542 + t648 * t543;
t537 = m(4) * t594 - qJDD(3) * mrSges(4,2) + t625 * mrSges(4,3) - qJD(3) * t628 + t623 * t667 + t661;
t532 = -t645 * t536 + t649 * t537;
t528 = m(3) * t609 - t651 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t532;
t608 = -t646 * t627 + t650 * t641;
t657 = -qJDD(2) * pkin(2) - t608;
t603 = -t651 * pkin(6) + t657;
t584 = -t625 * pkin(3) + t630 * t668 + (-pkin(7) * t640 - pkin(6)) * t651 + t657;
t562 = -t589 * pkin(4) - t610 * pkin(8) + t615 * t607 + t584;
t659 = m(6) * t562 - t570 * mrSges(6,1) + t571 * mrSges(6,2) - t599 * t591 + t600 * t592;
t655 = m(5) * t584 - t589 * mrSges(5,1) + t590 * mrSges(5,2) - t614 * t605 + t615 * t606 + t659;
t550 = -m(4) * t603 + t625 * mrSges(4,1) - t624 * mrSges(4,2) - t628 * t668 + t629 * t667 - t655;
t549 = m(3) * t608 + qJDD(2) * mrSges(3,1) - t651 * mrSges(3,2) + t550;
t662 = t650 * t528 - t646 * t549;
t524 = m(2) * t627 + t662;
t531 = t649 * t536 + t645 * t537;
t530 = (m(2) + m(3)) * t626 - t531;
t669 = t642 * t524 + t670 * t530;
t525 = t646 * t528 + t650 * t549;
t665 = m(2) * t641 + t525;
t663 = t670 * t524 - t642 * t530;
t573 = Ifges(6,5) * t600 + Ifges(6,6) * t599 + Ifges(6,3) * t635;
t546 = -mrSges(6,1) * t562 + mrSges(6,3) * t558 + Ifges(6,4) * t571 + Ifges(6,2) * t570 + Ifges(6,6) * t634 - t600 * t573 + t635 * t575;
t547 = mrSges(6,2) * t562 - mrSges(6,3) * t557 + Ifges(6,1) * t571 + Ifges(6,4) * t570 + Ifges(6,5) * t634 + t599 * t573 - t635 * t574;
t595 = Ifges(5,5) * t615 + Ifges(5,6) * t614 + Ifges(5,3) * t639;
t533 = -mrSges(5,1) * t584 + mrSges(5,3) * t565 + Ifges(5,4) * t590 + Ifges(5,2) * t589 + Ifges(5,6) * t638 - pkin(4) * t659 + pkin(8) * t660 + t647 * t546 + t643 * t547 - t615 * t595 + t639 * t597;
t534 = mrSges(5,2) * t584 - mrSges(5,3) * t564 + Ifges(5,1) * t590 + Ifges(5,4) * t589 + Ifges(5,5) * t638 - pkin(8) * t545 - t643 * t546 + t647 * t547 + t614 * t595 - t639 * t596;
t611 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t645 + Ifges(4,6) * t649) * qJD(2);
t520 = -mrSges(4,1) * t603 + mrSges(4,3) * t594 + Ifges(4,4) * t624 + Ifges(4,2) * t625 + Ifges(4,6) * qJDD(3) - pkin(3) * t655 + pkin(7) * t661 + qJD(3) * t613 + t648 * t533 + t644 * t534 - t611 * t668;
t521 = mrSges(4,2) * t603 - mrSges(4,3) * t593 + Ifges(4,1) * t624 + Ifges(4,4) * t625 + Ifges(4,5) * qJDD(3) - pkin(7) * t538 - qJD(3) * t612 - t644 * t533 + t648 * t534 + t611 * t667;
t654 = mrSges(3,1) * t608 - mrSges(3,2) * t609 + Ifges(3,3) * qJDD(2) + pkin(2) * t550 + pkin(6) * t532 + t649 * t520 + t645 * t521;
t519 = mrSges(3,1) * t626 + mrSges(3,3) * t609 + t651 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t531 - t671;
t518 = -mrSges(3,2) * t626 - mrSges(3,3) * t608 + Ifges(3,5) * qJDD(2) - t651 * Ifges(3,6) - pkin(6) * t531 - t645 * t520 + t649 * t521;
t517 = -mrSges(2,1) * t641 + mrSges(2,3) * t627 - pkin(1) * t525 - t654;
t516 = mrSges(2,2) * t641 - mrSges(2,3) * t626 - pkin(5) * t525 + t650 * t518 - t646 * t519;
t1 = [-m(1) * g(1) + t663; -m(1) * g(2) + t669; -m(1) * g(3) + t665; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t669 + t670 * t516 - t642 * t517; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t663 + t642 * t516 + t670 * t517; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t626 - mrSges(2,2) * t627 + t646 * t518 + t650 * t519 + pkin(1) * (m(3) * t626 - t531) + pkin(5) * t662; t665; t654; t671; -t653; -t656;];
tauJB = t1;
