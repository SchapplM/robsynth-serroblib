% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:42
% EndTime: 2019-12-31 19:46:47
% DurationCPUTime: 3.15s
% Computational Cost: add. (27915->294), mult. (60932->357), div. (0->0), fcn. (34590->8), ass. (0->114)
t617 = -2 * qJD(3);
t616 = Ifges(3,1) + Ifges(4,2);
t612 = Ifges(3,4) + Ifges(4,6);
t611 = Ifges(3,5) - Ifges(4,4);
t615 = Ifges(3,2) + Ifges(4,3);
t610 = Ifges(3,6) - Ifges(4,5);
t614 = (Ifges(3,3) + Ifges(4,1));
t581 = sin(qJ(1));
t584 = cos(qJ(1));
t566 = -g(1) * t584 - g(2) * t581;
t586 = qJD(1) ^ 2;
t543 = -pkin(1) * t586 + qJDD(1) * pkin(6) + t566;
t580 = sin(qJ(2));
t583 = cos(qJ(2));
t522 = -g(3) * t580 + t583 * t543;
t553 = (-pkin(2) * t583 - qJ(3) * t580) * qJD(1);
t585 = qJD(2) ^ 2;
t604 = qJD(1) * t583;
t510 = pkin(2) * t585 - qJDD(2) * qJ(3) + (qJD(2) * t617) - t553 * t604 - t522;
t613 = t586 * pkin(6);
t521 = -t583 * g(3) - t580 * t543;
t554 = (mrSges(4,2) * t583 - mrSges(4,3) * t580) * qJD(1);
t555 = (-mrSges(3,1) * t583 + mrSges(3,2) * t580) * qJD(1);
t602 = qJD(1) * qJD(2);
t599 = t583 * t602;
t556 = qJDD(1) * t580 + t599;
t561 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t604;
t563 = -mrSges(4,1) * t604 - qJD(2) * mrSges(4,3);
t600 = t580 * t602;
t557 = qJDD(1) * t583 - t600;
t603 = t580 * qJD(1);
t562 = pkin(3) * t603 - (qJD(2) * qJ(4));
t576 = t583 ^ 2;
t565 = t581 * g(1) - t584 * g(2);
t595 = -qJDD(1) * pkin(1) - t565;
t590 = pkin(2) * t600 + t603 * t617 + (-t556 - t599) * qJ(3) + t595;
t495 = -t562 * t603 + (-pkin(3) * t576 - pkin(6)) * t586 + (-pkin(2) - qJ(4)) * t557 + t590;
t511 = -qJDD(2) * pkin(2) - t585 * qJ(3) + t553 * t603 + qJDD(3) - t521;
t506 = (-t580 * t583 * t586 - qJDD(2)) * qJ(4) + (t556 - t599) * pkin(3) + t511;
t577 = sin(pkin(8));
t578 = cos(pkin(8));
t548 = qJD(2) * t578 - t577 * t604;
t490 = -0.2e1 * qJD(4) * t548 - t577 * t495 + t578 * t506;
t527 = qJDD(2) * t578 - t557 * t577;
t547 = -qJD(2) * t577 - t578 * t604;
t488 = (t547 * t603 - t527) * pkin(7) + (t547 * t548 + t556) * pkin(4) + t490;
t491 = 0.2e1 * qJD(4) * t547 + t495 * t578 + t577 * t506;
t526 = -qJDD(2) * t577 - t557 * t578;
t528 = pkin(4) * t603 - pkin(7) * t548;
t546 = t547 ^ 2;
t489 = -pkin(4) * t546 + pkin(7) * t526 - t528 * t603 + t491;
t579 = sin(qJ(5));
t582 = cos(qJ(5));
t486 = t488 * t582 - t489 * t579;
t518 = t547 * t582 - t548 * t579;
t498 = qJD(5) * t518 + t526 * t579 + t527 * t582;
t519 = t547 * t579 + t548 * t582;
t508 = -mrSges(6,1) * t518 + mrSges(6,2) * t519;
t568 = qJD(5) + t603;
t512 = -mrSges(6,2) * t568 + mrSges(6,3) * t518;
t552 = qJDD(5) + t556;
t484 = m(6) * t486 + mrSges(6,1) * t552 - mrSges(6,3) * t498 - t508 * t519 + t512 * t568;
t487 = t488 * t579 + t489 * t582;
t497 = -qJD(5) * t519 + t526 * t582 - t527 * t579;
t513 = mrSges(6,1) * t568 - mrSges(6,3) * t519;
t485 = m(6) * t487 - mrSges(6,2) * t552 + mrSges(6,3) * t497 + t508 * t518 - t513 * t568;
t475 = t484 * t582 + t485 * t579;
t520 = -mrSges(5,1) * t547 + mrSges(5,2) * t548;
t524 = -mrSges(5,2) * t603 + mrSges(5,3) * t547;
t473 = m(5) * t490 + mrSges(5,1) * t556 - mrSges(5,3) * t527 - t520 * t548 + t524 * t603 + t475;
t525 = mrSges(5,1) * t603 - mrSges(5,3) * t548;
t596 = -t484 * t579 + t485 * t582;
t474 = m(5) * t491 - mrSges(5,2) * t556 + mrSges(5,3) * t526 + t520 * t547 - t525 * t603 + t596;
t470 = t578 * t473 + t577 * t474;
t591 = -m(4) * t511 - t556 * mrSges(4,1) - t470;
t468 = m(3) * t521 - t556 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t561 - t563) * qJD(2) + (-t554 - t555) * t603 + t591;
t560 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t603;
t564 = mrSges(4,1) * t603 + qJD(2) * mrSges(4,2);
t502 = -qJ(4) * t576 * t586 + pkin(3) * t557 + qJD(2) * t562 + qJDD(4) - t510;
t493 = -pkin(4) * t526 - pkin(7) * t546 + t528 * t548 + t502;
t592 = m(6) * t493 - t497 * mrSges(6,1) + mrSges(6,2) * t498 - t518 * t512 + t519 * t513;
t589 = -m(5) * t502 + t526 * mrSges(5,1) - t527 * mrSges(5,2) + t547 * t524 - t548 * t525 - t592;
t588 = -m(4) * t510 + qJDD(2) * mrSges(4,3) + qJD(2) * t564 + t554 * t604 - t589;
t480 = t555 * t604 + t588 - qJD(2) * t560 + m(3) * t522 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t557;
t597 = -t468 * t580 + t480 * t583;
t461 = m(2) * t566 - mrSges(2,1) * t586 - qJDD(1) * mrSges(2,2) + t597;
t542 = t595 - t613;
t509 = -t557 * pkin(2) + t590 - t613;
t608 = -t473 * t577 + t474 * t578;
t594 = -m(4) * t509 - t557 * mrSges(4,2) + t564 * t603 - t608;
t587 = -m(3) * t542 + t561 * t604 + t557 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t556 + (-t560 * t580 - t563 * t583) * qJD(1) + t594;
t466 = m(2) * t565 + qJDD(1) * mrSges(2,1) - t586 * mrSges(2,2) + t587;
t609 = t461 * t581 + t466 * t584;
t462 = t468 * t583 + t480 * t580;
t607 = (t614 * qJD(2)) + (t580 * t611 + t583 * t610) * qJD(1);
t606 = -t610 * qJD(2) + (-t580 * t612 - t583 * t615) * qJD(1);
t605 = t611 * qJD(2) + (t580 * t616 + t583 * t612) * qJD(1);
t598 = t461 * t584 - t466 * t581;
t516 = Ifges(5,1) * t548 + Ifges(5,4) * t547 + Ifges(5,5) * t603;
t515 = Ifges(5,4) * t548 + Ifges(5,2) * t547 + Ifges(5,6) * t603;
t514 = Ifges(5,5) * t548 + Ifges(5,6) * t547 + Ifges(5,3) * t603;
t505 = Ifges(6,1) * t519 + Ifges(6,4) * t518 + Ifges(6,5) * t568;
t504 = Ifges(6,4) * t519 + Ifges(6,2) * t518 + Ifges(6,6) * t568;
t503 = Ifges(6,5) * t519 + Ifges(6,6) * t518 + Ifges(6,3) * t568;
t477 = mrSges(6,2) * t493 - mrSges(6,3) * t486 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t552 + t503 * t518 - t504 * t568;
t476 = -mrSges(6,1) * t493 + mrSges(6,3) * t487 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t552 - t503 * t519 + t505 * t568;
t469 = -t556 * mrSges(4,3) + t563 * t604 - t594;
t464 = mrSges(5,2) * t502 - mrSges(5,3) * t490 + Ifges(5,1) * t527 + Ifges(5,4) * t526 + Ifges(5,5) * t556 - pkin(7) * t475 - t476 * t579 + t477 * t582 + t514 * t547 - t515 * t603;
t463 = -mrSges(5,1) * t502 + mrSges(5,3) * t491 + Ifges(5,4) * t527 + Ifges(5,2) * t526 + Ifges(5,6) * t556 - pkin(4) * t592 + pkin(7) * t596 + t582 * t476 + t579 * t477 - t548 * t514 + t516 * t603;
t458 = t611 * qJDD(2) + t612 * t557 + t606 * qJD(2) + t607 * t604 + (Ifges(5,3) + t616) * t556 + Ifges(6,3) * t552 + mrSges(3,2) * t542 - t547 * t516 + t548 * t515 + Ifges(5,6) * t526 + Ifges(5,5) * t527 - t518 * t505 + t519 * t504 - mrSges(3,3) * t521 - mrSges(4,3) * t509 + mrSges(4,1) * t511 + Ifges(6,6) * t497 + Ifges(6,5) * t498 + mrSges(5,1) * t490 - mrSges(5,2) * t491 + mrSges(6,1) * t486 - mrSges(6,2) * t487 + pkin(4) * t475 + pkin(3) * t470 - qJ(3) * t469;
t457 = -mrSges(3,1) * t542 - mrSges(4,1) * t510 + mrSges(4,2) * t509 + mrSges(3,3) * t522 - pkin(2) * t469 - pkin(3) * t589 - qJ(4) * t608 + t605 * qJD(2) + t610 * qJDD(2) - t578 * t463 - t577 * t464 + t612 * t556 + t557 * t615 - t607 * t603;
t456 = -pkin(1) * t462 + mrSges(2,3) * t566 - pkin(2) * (-qJD(2) * t563 + t591) - qJ(3) * t588 - t578 * t464 + t577 * t463 + qJ(4) * t470 - mrSges(3,1) * t521 + mrSges(3,2) * t522 - mrSges(4,2) * t511 + mrSges(4,3) * t510 + mrSges(2,1) * g(3) + t586 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-mrSges(4,1) * qJ(3) - t610) * t557 - t611 * t556 + (mrSges(4,2) * pkin(2) - t614) * qJDD(2) + (t605 * t583 + (pkin(2) * t554 + t606) * t580) * qJD(1);
t455 = -mrSges(2,2) * g(3) - mrSges(2,3) * t565 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t586 - pkin(6) * t462 - t457 * t580 + t458 * t583;
t1 = [-m(1) * g(1) + t598; -m(1) * g(2) + t609; (-m(1) - m(2)) * g(3) + t462; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t609 + t455 * t584 - t456 * t581; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t598 + t581 * t455 + t584 * t456; -mrSges(1,1) * g(2) + mrSges(2,1) * t565 + mrSges(1,2) * g(1) - mrSges(2,2) * t566 + Ifges(2,3) * qJDD(1) + pkin(1) * t587 + pkin(6) * t597 + t583 * t457 + t580 * t458;];
tauB = t1;
