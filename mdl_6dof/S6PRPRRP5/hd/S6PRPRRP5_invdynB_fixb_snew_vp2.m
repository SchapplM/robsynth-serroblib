% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-04 23:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:54:33
% EndTime: 2019-05-04 23:54:40
% DurationCPUTime: 4.04s
% Computational Cost: add. (45019->270), mult. (82724->324), div. (0->0), fcn. (51636->10), ass. (0->117)
t637 = Ifges(6,1) + Ifges(7,1);
t628 = Ifges(6,4) + Ifges(7,4);
t626 = Ifges(6,5) + Ifges(7,5);
t636 = Ifges(6,2) + Ifges(7,2);
t635 = Ifges(6,6) + Ifges(7,6);
t634 = Ifges(6,3) + Ifges(7,3);
t585 = sin(pkin(10));
t587 = cos(pkin(10));
t573 = t585 * g(1) - t587 * g(2);
t574 = -t587 * g(1) - t585 * g(2);
t582 = -g(3) + qJDD(1);
t594 = cos(qJ(2));
t588 = cos(pkin(6));
t591 = sin(qJ(2));
t622 = t588 * t591;
t586 = sin(pkin(6));
t623 = t586 * t591;
t528 = t573 * t622 + t594 * t574 + t582 * t623;
t633 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t528;
t527 = -t591 * t574 + (t573 * t588 + t582 * t586) * t594;
t632 = -pkin(2) - pkin(8);
t631 = mrSges(3,1) - mrSges(4,2);
t630 = -mrSges(6,2) - mrSges(7,2);
t629 = (-Ifges(4,4) + Ifges(3,5));
t627 = Ifges(4,5) - Ifges(3,6);
t596 = qJD(2) ^ 2;
t599 = -t596 * qJ(3) + qJDD(3) - t527;
t523 = t632 * qJDD(2) + t599;
t553 = -t586 * t573 + t588 * t582;
t590 = sin(qJ(4));
t593 = cos(qJ(4));
t519 = t590 * t523 + t593 * t553;
t569 = (mrSges(5,1) * t590 + mrSges(5,2) * t593) * qJD(2);
t614 = qJD(2) * qJD(4);
t609 = t593 * t614;
t571 = -t590 * qJDD(2) - t609;
t616 = qJD(2) * t593;
t576 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t616;
t570 = (pkin(4) * t590 - pkin(9) * t593) * qJD(2);
t595 = qJD(4) ^ 2;
t615 = t590 * qJD(2);
t514 = -t595 * pkin(4) + qJDD(4) * pkin(9) - t570 * t615 + t519;
t522 = t632 * t596 - t633;
t610 = t590 * t614;
t572 = t593 * qJDD(2) - t610;
t517 = (-t572 + t610) * pkin(9) + (-t571 + t609) * pkin(4) + t522;
t589 = sin(qJ(5));
t592 = cos(qJ(5));
t509 = -t589 * t514 + t592 * t517;
t567 = t592 * qJD(4) - t589 * t616;
t541 = t567 * qJD(5) + t589 * qJDD(4) + t592 * t572;
t568 = t589 * qJD(4) + t592 * t616;
t543 = -t567 * mrSges(7,1) + t568 * mrSges(7,2);
t544 = -t567 * mrSges(6,1) + t568 * mrSges(6,2);
t578 = qJD(5) + t615;
t548 = -t578 * mrSges(6,2) + t567 * mrSges(6,3);
t564 = qJDD(5) - t571;
t506 = -0.2e1 * qJD(6) * t568 + (t567 * t578 - t541) * qJ(6) + (t567 * t568 + t564) * pkin(5) + t509;
t547 = -t578 * mrSges(7,2) + t567 * mrSges(7,3);
t612 = m(7) * t506 + t564 * mrSges(7,1) + t578 * t547;
t499 = m(6) * t509 + t564 * mrSges(6,1) + t578 * t548 + (-t543 - t544) * t568 + (-mrSges(6,3) - mrSges(7,3)) * t541 + t612;
t510 = t592 * t514 + t589 * t517;
t540 = -t568 * qJD(5) + t592 * qJDD(4) - t589 * t572;
t549 = t578 * pkin(5) - t568 * qJ(6);
t563 = t567 ^ 2;
t508 = -t563 * pkin(5) + t540 * qJ(6) + 0.2e1 * qJD(6) * t567 - t578 * t549 + t510;
t611 = m(7) * t508 + t540 * mrSges(7,3) + t567 * t543;
t550 = t578 * mrSges(7,1) - t568 * mrSges(7,3);
t617 = -t578 * mrSges(6,1) + t568 * mrSges(6,3) - t550;
t502 = m(6) * t510 + t540 * mrSges(6,3) + t567 * t544 + t630 * t564 + t617 * t578 + t611;
t606 = -t589 * t499 + t592 * t502;
t495 = m(5) * t519 - qJDD(4) * mrSges(5,2) + t571 * mrSges(5,3) - qJD(4) * t576 - t569 * t615 + t606;
t518 = t593 * t523 - t590 * t553;
t575 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t615;
t513 = -qJDD(4) * pkin(4) - t595 * pkin(9) + t570 * t616 - t518;
t511 = -t540 * pkin(5) - t563 * qJ(6) + t568 * t549 + qJDD(6) + t513;
t604 = m(7) * t511 - t540 * mrSges(7,1) - t567 * t547;
t597 = -m(6) * t513 + t540 * mrSges(6,1) + t630 * t541 + t567 * t548 + t617 * t568 - t604;
t503 = m(5) * t518 + qJDD(4) * mrSges(5,1) - t572 * mrSges(5,3) + qJD(4) * t575 - t569 * t616 + t597;
t488 = t590 * t495 + t593 * t503;
t525 = -qJDD(2) * pkin(2) + t599;
t601 = -m(4) * t525 + (t596 * mrSges(4,3)) - t488;
t484 = m(3) * t527 - (t596 * mrSges(3,2)) + t631 * qJDD(2) + t601;
t624 = t484 * t594;
t607 = t593 * t495 - t590 * t503;
t487 = m(4) * t553 + t607;
t486 = m(3) * t553 + t487;
t524 = t596 * pkin(2) + t633;
t497 = t592 * t499 + t589 * t502;
t600 = -m(5) * t522 + t571 * mrSges(5,1) - t572 * mrSges(5,2) - t575 * t615 - t576 * t616 - t497;
t598 = -m(4) * t524 + (t596 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t600;
t493 = m(3) * t528 - (t596 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t598;
t475 = -t586 * t486 + t493 * t622 + t588 * t624;
t473 = m(2) * t573 + t475;
t480 = -t591 * t484 + t594 * t493;
t479 = m(2) * t574 + t480;
t621 = t587 * t473 + t585 * t479;
t620 = t635 * t567 + t626 * t568 + t634 * t578;
t619 = -t636 * t567 - t628 * t568 - t635 * t578;
t618 = t628 * t567 + t637 * t568 + t626 * t578;
t474 = t588 * t486 + t493 * t623 + t586 * t624;
t608 = -t585 * t473 + t587 * t479;
t489 = -mrSges(6,1) * t513 + mrSges(6,3) * t510 - mrSges(7,1) * t511 + mrSges(7,3) * t508 - pkin(5) * t604 + qJ(6) * t611 + (-qJ(6) * t550 + t618) * t578 + (-pkin(5) * t550 - t620) * t568 + (-qJ(6) * mrSges(7,2) + t635) * t564 + (-pkin(5) * mrSges(7,2) + t628) * t541 + t636 * t540;
t504 = -t541 * mrSges(7,3) - t568 * t543 + t612;
t496 = mrSges(6,2) * t513 + mrSges(7,2) * t511 - mrSges(6,3) * t509 - mrSges(7,3) * t506 - qJ(6) * t504 + t628 * t540 + t637 * t541 + t626 * t564 + t620 * t567 + t619 * t578;
t556 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t593 - Ifges(5,6) * t590) * qJD(2);
t557 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t593 - Ifges(5,2) * t590) * qJD(2);
t476 = mrSges(5,2) * t522 - mrSges(5,3) * t518 + Ifges(5,1) * t572 + Ifges(5,4) * t571 + Ifges(5,5) * qJDD(4) - pkin(9) * t497 - qJD(4) * t557 - t589 * t489 + t592 * t496 - t556 * t615;
t558 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t593 - Ifges(5,4) * t590) * qJD(2);
t481 = -t556 * t616 - mrSges(5,1) * t522 - mrSges(6,1) * t509 - mrSges(7,1) * t506 + mrSges(6,2) * t510 + mrSges(7,2) * t508 + mrSges(5,3) * t519 + Ifges(5,4) * t572 + Ifges(5,2) * t571 + Ifges(5,6) * qJDD(4) - pkin(4) * t497 - pkin(5) * t504 + qJD(4) * t558 + t619 * t568 + t618 * t567 - t634 * t564 - t626 * t541 - t635 * t540;
t470 = -mrSges(4,1) * t524 + mrSges(3,3) * t528 - pkin(2) * t487 - pkin(3) * t600 - pkin(8) * t607 - t627 * qJDD(2) - t590 * t476 - t593 * t481 - t631 * t553 + (t629 * t596);
t471 = -qJ(3) * t487 - mrSges(3,3) * t527 + pkin(3) * t488 + mrSges(4,1) * t525 + pkin(9) * t606 + t589 * t496 + t592 * t489 + pkin(4) * t597 + Ifges(5,5) * t572 + Ifges(5,6) * t571 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t518 - mrSges(5,2) * t519 + t627 * t596 + (mrSges(3,2) - mrSges(4,3)) * t553 + t629 * qJDD(2) + (t593 * t557 + t590 * t558) * qJD(2);
t602 = pkin(7) * t480 + t470 * t594 + t471 * t591;
t469 = mrSges(3,1) * t527 - mrSges(3,2) * t528 + mrSges(4,2) * t525 - mrSges(4,3) * t524 + t593 * t476 - t590 * t481 - pkin(8) * t488 + pkin(2) * t601 + qJ(3) * t598 + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t468 = mrSges(2,2) * t582 - mrSges(2,3) * t573 - t591 * t470 + t594 * t471 + (-t474 * t586 - t475 * t588) * pkin(7);
t467 = -mrSges(2,1) * t582 + mrSges(2,3) * t574 - pkin(1) * t474 - t586 * t469 + t602 * t588;
t1 = [-m(1) * g(1) + t608; -m(1) * g(2) + t621; -m(1) * g(3) + m(2) * t582 + t474; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t621 - t585 * t467 + t587 * t468; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t608 + t587 * t467 + t585 * t468; -mrSges(1,1) * g(2) + mrSges(2,1) * t573 + mrSges(1,2) * g(1) - mrSges(2,2) * t574 + pkin(1) * t475 + t588 * t469 + t602 * t586;];
tauB  = t1;
