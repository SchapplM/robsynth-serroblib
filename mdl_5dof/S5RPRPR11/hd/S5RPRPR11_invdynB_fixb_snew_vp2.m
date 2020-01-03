% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:25
% EndTime: 2019-12-31 18:27:29
% DurationCPUTime: 3.19s
% Computational Cost: add. (27296->271), mult. (66002->327), div. (0->0), fcn. (44794->8), ass. (0->117)
t589 = sin(pkin(8));
t583 = t589 ^ 2;
t590 = cos(pkin(8));
t584 = t590 ^ 2;
t597 = qJD(1) ^ 2;
t593 = sin(qJ(1));
t595 = cos(qJ(1));
t570 = t593 * g(1) - t595 * g(2);
t614 = -qJDD(2) + t570;
t635 = pkin(2) * t590;
t550 = -(qJ(2) + (t583 + t584) * pkin(6)) * t597 - (pkin(1) + t635) * qJDD(1) - t614;
t592 = sin(qJ(3));
t636 = cos(qJ(3));
t602 = t636 * t589 + t590 * t592;
t615 = t590 * t636;
t620 = t589 * qJD(1);
t565 = -qJD(1) * t615 + t592 * t620;
t622 = t565 * qJD(3);
t552 = t602 * qJDD(1) - t622;
t644 = t550 + (-t552 + t622) * qJ(4);
t643 = Ifges(4,1) + Ifges(5,1);
t633 = Ifges(4,4) - Ifges(5,5);
t632 = Ifges(4,5) + Ifges(5,4);
t642 = Ifges(4,2) + Ifges(5,3);
t631 = Ifges(4,6) - Ifges(5,6);
t641 = -Ifges(4,3) - Ifges(5,2);
t637 = 2 * qJD(4);
t634 = -mrSges(4,3) - mrSges(5,2);
t630 = mrSges(3,2) * t589;
t629 = t584 * t597;
t571 = -t595 * g(1) - t593 * g(2);
t567 = -t597 * pkin(1) + qJDD(1) * qJ(2) + t571;
t619 = qJD(1) * qJD(2);
t613 = -t590 * g(3) - 0.2e1 * t589 * t619;
t531 = (-pkin(6) * qJDD(1) + t597 * t635 - t567) * t589 + t613;
t554 = -t589 * g(3) + (t567 + 0.2e1 * t619) * t590;
t617 = qJDD(1) * t590;
t532 = -pkin(2) * t629 + pkin(6) * t617 + t554;
t517 = t592 * t531 + t636 * t532;
t618 = qJDD(1) * t589;
t566 = t602 * qJD(1);
t621 = t566 * qJD(3);
t551 = -qJDD(1) * t615 + t592 * t618 + t621;
t558 = qJD(3) * mrSges(4,1) - t566 * mrSges(4,3);
t544 = t565 * pkin(3) - t566 * qJ(4);
t596 = qJD(3) ^ 2;
t510 = -t596 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t637 - t565 * t544 + t517;
t559 = -qJD(3) * mrSges(5,1) + t566 * mrSges(5,2);
t516 = t636 * t531 - t592 * t532;
t511 = -qJDD(3) * pkin(3) - t596 * qJ(4) + t566 * t544 + qJDD(4) - t516;
t505 = (-t552 - t622) * pkin(7) + (t565 * t566 - qJDD(3)) * pkin(4) + t511;
t561 = -qJD(3) * pkin(4) - t566 * pkin(7);
t564 = t565 ^ 2;
t506 = -t564 * pkin(4) + t551 * pkin(7) + qJD(3) * t561 + t510;
t591 = sin(qJ(5));
t594 = cos(qJ(5));
t501 = t594 * t505 - t591 * t506;
t539 = t594 * t565 - t591 * t566;
t515 = t539 * qJD(5) + t591 * t551 + t594 * t552;
t540 = t591 * t565 + t594 * t566;
t523 = -t539 * mrSges(6,1) + t540 * mrSges(6,2);
t585 = -qJD(3) + qJD(5);
t527 = -t585 * mrSges(6,2) + t539 * mrSges(6,3);
t582 = -qJDD(3) + qJDD(5);
t499 = m(6) * t501 + t582 * mrSges(6,1) - t515 * mrSges(6,3) - t540 * t523 + t585 * t527;
t502 = t591 * t505 + t594 * t506;
t514 = -t540 * qJD(5) + t594 * t551 - t591 * t552;
t528 = t585 * mrSges(6,1) - t540 * mrSges(6,3);
t500 = m(6) * t502 - t582 * mrSges(6,2) + t514 * mrSges(6,3) + t539 * t523 - t585 * t528;
t609 = -t591 * t499 + t594 * t500;
t601 = m(5) * t510 + qJDD(3) * mrSges(5,3) + qJD(3) * t559 + t609;
t545 = t565 * mrSges(5,1) - t566 * mrSges(5,3);
t624 = -t565 * mrSges(4,1) - t566 * mrSges(4,2) - t545;
t490 = m(4) * t517 - qJDD(3) * mrSges(4,2) - qJD(3) * t558 + t634 * t551 + t624 * t565 + t601;
t557 = -qJD(3) * mrSges(4,2) - t565 * mrSges(4,3);
t492 = t594 * t499 + t591 * t500;
t560 = -t565 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t600 = -m(5) * t511 + qJDD(3) * mrSges(5,1) + qJD(3) * t560 - t492;
t491 = m(4) * t516 + qJDD(3) * mrSges(4,1) + qJD(3) * t557 + t634 * t552 + t624 * t566 + t600;
t486 = t592 * t490 + t636 * t491;
t553 = -t589 * t567 + t613;
t603 = mrSges(3,3) * qJDD(1) + t597 * (-mrSges(3,1) * t590 + t630);
t484 = m(3) * t553 - t589 * t603 + t486;
t610 = t636 * t490 - t592 * t491;
t485 = m(3) * t554 + t590 * t603 + t610;
t611 = -t589 * t484 + t590 * t485;
t477 = m(2) * t571 - t597 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t611;
t563 = -qJDD(1) * pkin(1) - t597 * qJ(2) - t614;
t508 = -0.2e1 * qJD(4) * t566 + (t551 + t621) * pkin(3) + t644;
t504 = -t564 * pkin(7) + (-pkin(3) - pkin(4)) * t551 + (-pkin(3) * qJD(3) + t561 + t637) * t566 - t644;
t604 = -m(6) * t504 + t514 * mrSges(6,1) - t515 * mrSges(6,2) + t539 * t527 - t540 * t528;
t497 = m(5) * t508 + t551 * mrSges(5,1) - t552 * mrSges(5,3) - t566 * t559 + t565 * t560 + t604;
t599 = m(4) * t550 + t551 * mrSges(4,1) + t552 * mrSges(4,2) + t565 * t557 + t566 * t558 + t497;
t598 = -m(3) * t563 + mrSges(3,1) * t617 - t599 + (t583 * t597 + t629) * mrSges(3,3);
t494 = (mrSges(2,1) - t630) * qJDD(1) + m(2) * t570 - t597 * mrSges(2,2) + t598;
t628 = t593 * t477 + t595 * t494;
t478 = t590 * t484 + t589 * t485;
t627 = -t631 * qJD(3) + t642 * t565 - t633 * t566;
t626 = t641 * qJD(3) + t631 * t565 - t632 * t566;
t625 = t632 * qJD(3) - t633 * t565 + t643 * t566;
t612 = t595 * t477 - t593 * t494;
t607 = Ifges(3,1) * t589 + Ifges(3,4) * t590;
t606 = Ifges(3,4) * t589 + Ifges(3,2) * t590;
t605 = Ifges(3,5) * t589 + Ifges(3,6) * t590;
t569 = t605 * qJD(1);
t520 = Ifges(6,1) * t540 + Ifges(6,4) * t539 + Ifges(6,5) * t585;
t519 = Ifges(6,4) * t540 + Ifges(6,2) * t539 + Ifges(6,6) * t585;
t518 = Ifges(6,5) * t540 + Ifges(6,6) * t539 + Ifges(6,3) * t585;
t496 = mrSges(6,2) * t504 - mrSges(6,3) * t501 + Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t582 + t539 * t518 - t585 * t519;
t495 = -mrSges(6,1) * t504 + mrSges(6,3) * t502 + Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t582 - t540 * t518 + t585 * t520;
t480 = mrSges(4,2) * t550 + mrSges(5,2) * t511 - mrSges(4,3) * t516 - mrSges(5,3) * t508 - pkin(7) * t492 - qJ(4) * t497 + t627 * qJD(3) + t632 * qJDD(3) - t591 * t495 + t594 * t496 - t633 * t551 + t643 * t552 + t626 * t565;
t479 = -mrSges(4,1) * t550 - mrSges(5,1) * t508 + mrSges(5,2) * t510 + mrSges(4,3) * t517 - pkin(3) * t497 - pkin(4) * t604 - pkin(7) * t609 + t625 * qJD(3) + t631 * qJDD(3) - t594 * t495 - t591 * t496 - t642 * t551 + t633 * t552 + t626 * t566;
t474 = t590 * qJD(1) * t569 + mrSges(3,2) * t563 - mrSges(3,3) * t553 - pkin(6) * t486 + t607 * qJDD(1) - t592 * t479 + t636 * t480;
t473 = -mrSges(3,1) * t563 + mrSges(3,3) * t554 - pkin(2) * t599 + pkin(6) * t610 + t606 * qJDD(1) + t636 * t479 + t592 * t480 - t569 * t620;
t472 = t641 * qJDD(3) - mrSges(3,1) * t553 + mrSges(3,2) * t554 + mrSges(2,3) * t571 - t539 * t520 + t540 * t519 - mrSges(4,1) * t516 + mrSges(4,2) * t517 + Ifges(6,6) * t514 + Ifges(6,5) * t515 - mrSges(5,3) * t510 + mrSges(5,1) * t511 + mrSges(6,1) * t501 - mrSges(6,2) * t502 + Ifges(6,3) * t582 + (Ifges(2,6) - t605) * qJDD(1) - qJ(4) * t601 - pkin(3) * t600 + pkin(4) * t492 - pkin(2) * t486 - pkin(1) * t478 + (qJ(4) * mrSges(5,2) + t631) * t551 + (pkin(3) * mrSges(5,2) - t632) * t552 + (qJ(4) * t545 - t625) * t565 + (pkin(3) * t545 + t627) * t566 + mrSges(2,1) * g(3) + (-t589 * t606 + t590 * t607 + Ifges(2,5)) * t597;
t471 = -mrSges(2,2) * g(3) - mrSges(2,3) * t570 + Ifges(2,5) * qJDD(1) - t597 * Ifges(2,6) - qJ(2) * t478 - t589 * t473 + t590 * t474;
t1 = [-m(1) * g(1) + t612; -m(1) * g(2) + t628; (-m(1) - m(2)) * g(3) + t478; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t628 + t595 * t471 - t593 * t472; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t612 + t593 * t471 + t595 * t472; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t570 - mrSges(2,2) * t571 + t589 * t474 + t590 * t473 + pkin(1) * (-mrSges(3,2) * t618 + t598) + qJ(2) * t611;];
tauB = t1;
