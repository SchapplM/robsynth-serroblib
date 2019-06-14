% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:51:01
% EndTime: 2019-05-05 21:51:05
% DurationCPUTime: 2.35s
% Computational Cost: add. (19565->294), mult. (36827->329), div. (0->0), fcn. (20600->6), ass. (0->113)
t661 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t643 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t642 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t660 = Ifges(5,2) + Ifges(7,2) + Ifges(6,3);
t641 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t659 = -Ifges(5,3) - Ifges(6,1) - Ifges(7,1);
t614 = sin(qJ(1));
t616 = cos(qJ(1));
t601 = -t616 * g(1) - t614 * g(2);
t658 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t601;
t657 = -2 * qJD(5);
t656 = -pkin(1) - pkin(7);
t655 = cos(qJ(4));
t654 = mrSges(2,1) - mrSges(3,2);
t653 = -mrSges(7,1) - mrSges(5,3);
t652 = -Ifges(3,4) + Ifges(2,5);
t651 = (Ifges(3,5) - Ifges(2,6));
t612 = sin(qJ(4));
t615 = cos(qJ(3));
t647 = qJD(1) * t615;
t592 = -t655 * qJD(3) + t612 * t647;
t613 = sin(qJ(3));
t646 = t613 * qJD(1);
t603 = qJD(4) + t646;
t650 = t592 * t603;
t600 = t614 * g(1) - t616 * g(2);
t618 = qJD(1) ^ 2;
t628 = -t618 * qJ(2) + qJDD(2) - t600;
t576 = t656 * qJDD(1) + t628;
t564 = -t615 * g(3) + t613 * t576;
t594 = (mrSges(4,1) * t613 + mrSges(4,2) * t615) * qJD(1);
t645 = qJD(1) * qJD(3);
t635 = t615 * t645;
t596 = -t613 * qJDD(1) - t635;
t599 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t647;
t574 = t656 * t618 - t658;
t636 = t613 * t645;
t597 = t615 * qJDD(1) - t636;
t529 = (-t597 + t636) * pkin(8) + (-t596 + t635) * pkin(3) + t574;
t595 = (pkin(3) * t613 - pkin(8) * t615) * qJD(1);
t617 = qJD(3) ^ 2;
t533 = -t617 * pkin(3) + qJDD(3) * pkin(8) - t595 * t646 + t564;
t526 = t655 * t529 - t612 * t533;
t555 = -t592 * qJD(4) + t612 * qJDD(3) + t655 * t597;
t593 = t612 * qJD(3) + t655 * t647;
t559 = -t593 * mrSges(7,2) + t592 * mrSges(7,3);
t561 = t592 * mrSges(5,1) + t593 * mrSges(5,2);
t565 = -t603 * mrSges(5,2) - t592 * mrSges(5,3);
t569 = t592 * mrSges(6,1) - t603 * mrSges(6,3);
t591 = qJDD(4) - t596;
t560 = t592 * pkin(4) - t593 * qJ(5);
t602 = t603 ^ 2;
t524 = -t591 * pkin(4) - t602 * qJ(5) + t593 * t560 + qJDD(5) - t526;
t562 = -t592 * mrSges(6,2) - t593 * mrSges(6,3);
t518 = -0.2e1 * qJD(6) * t603 + (t592 * t593 - t591) * qJ(6) + (t555 + t650) * pkin(5) + t524;
t570 = -t592 * mrSges(7,1) + t603 * mrSges(7,2);
t631 = m(7) * t518 - t591 * mrSges(7,3) - t603 * t570;
t624 = -m(6) * t524 - t555 * mrSges(6,1) - t593 * t562 - t631;
t512 = m(5) * t526 + (t565 - t569) * t603 + (-t559 - t561) * t593 + (mrSges(5,1) - mrSges(6,2)) * t591 + t653 * t555 + t624;
t527 = t612 * t529 + t655 * t533;
t554 = t593 * qJD(4) - t655 * qJDD(3) + t612 * t597;
t566 = t603 * mrSges(5,1) - t593 * mrSges(5,3);
t622 = -t602 * pkin(4) + t591 * qJ(5) - t592 * t560 + t527;
t523 = t603 * t657 - t622;
t571 = t593 * mrSges(6,1) + t603 * mrSges(6,2);
t567 = t593 * pkin(5) - t603 * qJ(6);
t590 = t592 ^ 2;
t522 = -t554 * pkin(5) - t590 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t567) * t603 + t622;
t568 = t593 * mrSges(7,1) - t603 * mrSges(7,3);
t640 = m(7) * t522 + t591 * mrSges(7,2) + t603 * t568;
t627 = -m(6) * t523 + t591 * mrSges(6,3) + t603 * t571 + t640;
t648 = -t559 - t562;
t515 = m(5) * t527 - t591 * mrSges(5,2) - t603 * t566 + (-t561 + t648) * t592 + (-mrSges(6,1) + t653) * t554 + t627;
t632 = -t612 * t512 + t655 * t515;
t508 = m(4) * t564 - qJDD(3) * mrSges(4,2) + t596 * mrSges(4,3) - qJD(3) * t599 - t594 * t646 + t632;
t563 = t613 * g(3) + t615 * t576;
t598 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t646;
t532 = -qJDD(3) * pkin(3) - t617 * pkin(8) + t595 * t647 - t563;
t620 = (-t555 + t650) * qJ(5) + t532 + (t603 * pkin(4) + t657) * t593;
t525 = t554 * pkin(4) + t620;
t520 = -t590 * pkin(5) + 0.2e1 * qJD(6) * t592 - t593 * t567 + (pkin(4) + qJ(6)) * t554 + t620;
t630 = m(7) * t520 - t555 * mrSges(7,2) + t554 * mrSges(7,3) - t593 * t568 + t592 * t570;
t623 = -m(6) * t525 + t554 * mrSges(6,2) + t592 * t569 - t630;
t619 = -m(5) * t532 - t554 * mrSges(5,1) - t592 * t565 + (-t566 + t571) * t593 + (-mrSges(5,2) + mrSges(6,3)) * t555 + t623;
t510 = m(4) * t563 + qJDD(3) * mrSges(4,1) - t597 * mrSges(4,3) + qJD(3) * t598 - t594 * t647 + t619;
t502 = t613 * t508 + t615 * t510;
t578 = -qJDD(1) * pkin(1) + t628;
t626 = -m(3) * t578 + (t618 * mrSges(3,3)) - t502;
t500 = m(2) * t600 - (t618 * mrSges(2,2)) + t654 * qJDD(1) + t626;
t577 = t618 * pkin(1) + t658;
t509 = t655 * t512 + t612 * t515;
t625 = -m(4) * t574 + t596 * mrSges(4,1) - t597 * mrSges(4,2) - t598 * t646 - t599 * t647 - t509;
t621 = -m(3) * t577 + (t618 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t625;
t506 = m(2) * t601 - (t618 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t621;
t649 = t616 * t500 + t614 * t506;
t639 = t641 * t592 - t642 * t593 + t659 * t603;
t638 = t660 * t592 - t643 * t593 - t641 * t603;
t637 = -t643 * t592 + t661 * t593 + t642 * t603;
t634 = -t614 * t500 + t616 * t506;
t633 = t615 * t508 - t613 * t510;
t583 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t615 - Ifges(4,4) * t613) * qJD(1);
t582 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t615 - Ifges(4,2) * t613) * qJD(1);
t581 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t615 - Ifges(4,6) * t613) * qJD(1);
t517 = t555 * mrSges(7,1) + t593 * t559 + t631;
t516 = -t555 * mrSges(6,3) - t593 * t571 - t623;
t503 = mrSges(6,1) * t524 + mrSges(7,1) * t518 + mrSges(5,2) * t532 - mrSges(7,2) * t520 - mrSges(5,3) * t526 - mrSges(6,3) * t525 + pkin(5) * t517 - qJ(5) * t516 - t643 * t554 + t661 * t555 + t642 * t591 + t639 * t592 + t638 * t603;
t501 = -m(3) * g(3) + t633;
t498 = -mrSges(5,1) * t532 + mrSges(5,3) * t527 - mrSges(6,1) * t523 + mrSges(6,2) * t525 + mrSges(7,1) * t522 - mrSges(7,3) * t520 - pkin(5) * (t592 * t559 - t640) - qJ(6) * t630 - pkin(4) * t516 + t637 * t603 + t639 * t593 + t641 * t591 + t643 * t555 + (-pkin(5) * mrSges(7,1) - t660) * t554;
t497 = -qJ(5) * t627 - t581 * t647 + Ifges(4,6) * qJDD(3) - pkin(4) * (-t603 * t569 + t624) + Ifges(4,4) * t597 + Ifges(4,2) * t596 - mrSges(4,1) * t574 + qJD(3) * t583 + mrSges(4,3) * t564 + mrSges(6,3) * t523 - mrSges(6,2) * t524 - mrSges(5,1) * t526 + mrSges(5,2) * t527 + qJ(6) * t517 + mrSges(7,3) * t518 - mrSges(7,2) * t522 - pkin(3) * t509 + (pkin(4) * t559 + t638) * t593 + (pkin(4) * mrSges(6,2) + t659) * t591 + (pkin(4) * mrSges(7,1) - t642) * t555 + (-qJ(5) * t648 - t637) * t592 + (-qJ(5) * (-mrSges(6,1) - mrSges(7,1)) + t641) * t554;
t496 = mrSges(4,2) * t574 - mrSges(4,3) * t563 + Ifges(4,1) * t597 + Ifges(4,4) * t596 + Ifges(4,5) * qJDD(3) - pkin(8) * t509 - qJD(3) * t582 - t612 * t498 + t655 * t503 - t581 * t646;
t495 = -qJ(2) * t501 - mrSges(2,3) * t600 + pkin(2) * t502 + mrSges(3,1) * t578 + pkin(8) * t632 + t612 * t503 + t655 * t498 + pkin(3) * t619 + mrSges(4,1) * t563 - mrSges(4,2) * t564 + Ifges(4,5) * t597 + Ifges(4,6) * t596 + Ifges(4,3) * qJDD(3) + (t651 * t618) + t652 * qJDD(1) + (t615 * t582 + t613 * t583) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t494 = -mrSges(3,1) * t577 + mrSges(2,3) * t601 - pkin(1) * t501 - pkin(2) * t625 - pkin(7) * t633 + t654 * g(3) - t651 * qJDD(1) - t613 * t496 - t615 * t497 + t652 * t618;
t1 = [-m(1) * g(1) + t634; -m(1) * g(2) + t649; (-m(1) - m(2) - m(3)) * g(3) + t633; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t649 - t614 * t494 + t616 * t495; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t634 + t616 * t494 + t614 * t495; pkin(1) * t626 + qJ(2) * t621 - mrSges(2,2) * t601 + t615 * t496 - t613 * t497 - pkin(7) * t502 + mrSges(2,1) * t600 + mrSges(3,2) * t578 - mrSges(3,3) * t577 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
