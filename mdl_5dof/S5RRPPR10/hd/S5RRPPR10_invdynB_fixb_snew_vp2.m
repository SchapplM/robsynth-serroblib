% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:19
% EndTime: 2019-12-31 19:43:24
% DurationCPUTime: 3.20s
% Computational Cost: add. (26807->292), mult. (58576->351), div. (0->0), fcn. (36026->8), ass. (0->114)
t619 = -2 * qJD(3);
t618 = Ifges(4,1) + Ifges(5,1);
t610 = Ifges(4,4) - Ifges(5,5);
t617 = Ifges(4,5) + Ifges(5,4);
t616 = -Ifges(4,2) - Ifges(5,3);
t615 = Ifges(5,2) + Ifges(4,3);
t614 = Ifges(4,6) - Ifges(5,6);
t577 = sin(qJ(1));
t580 = cos(qJ(1));
t564 = t577 * g(1) - t580 * g(2);
t582 = qJD(1) ^ 2;
t547 = -qJDD(1) * pkin(1) - t582 * pkin(6) - t564;
t576 = sin(qJ(2));
t579 = cos(qJ(2));
t597 = qJD(1) * qJD(2);
t595 = t579 * t597;
t560 = t576 * qJDD(1) + t595;
t568 = t576 * t597;
t561 = t579 * qJDD(1) - t568;
t506 = (-t560 - t595) * qJ(3) + (-t561 + t568) * pkin(2) + t547;
t565 = -t580 * g(1) - t577 * g(2);
t548 = -t582 * pkin(1) + qJDD(1) * pkin(6) + t565;
t529 = -t576 * g(3) + t579 * t548;
t558 = (-pkin(2) * t579 - qJ(3) * t576) * qJD(1);
t581 = qJD(2) ^ 2;
t598 = t579 * qJD(1);
t510 = -t581 * pkin(2) + qJDD(2) * qJ(3) + t558 * t598 + t529;
t574 = sin(pkin(8));
t600 = qJD(1) * t576;
t607 = cos(pkin(8));
t553 = t574 * qJD(2) + t607 * t600;
t493 = t607 * t506 - t574 * t510 + t553 * t619;
t537 = t574 * qJDD(2) + t607 * t560;
t528 = -t579 * g(3) - t576 * t548;
t585 = qJDD(2) * pkin(2) + t581 * qJ(3) - t558 * t600 - qJDD(3) + t528;
t552 = -t607 * qJD(2) + t574 * t600;
t596 = t552 * t598;
t613 = -(t537 + t596) * qJ(4) - t585;
t612 = -2 * qJD(4);
t611 = -mrSges(4,3) - mrSges(5,2);
t606 = t579 ^ 2 * t582;
t559 = (-mrSges(3,1) * t579 + mrSges(3,2) * t576) * qJD(1);
t562 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t600;
t494 = t574 * t506 + t607 * t510 + t552 * t619;
t534 = -mrSges(4,1) * t598 - t553 * mrSges(4,3);
t535 = mrSges(5,1) * t598 + t553 * mrSges(5,2);
t536 = -t607 * qJDD(2) + t574 * t560;
t525 = t552 * pkin(3) - t553 * qJ(4);
t490 = -pkin(3) * t606 - t561 * qJ(4) - t552 * t525 + t598 * t612 + t494;
t491 = t561 * pkin(3) - qJ(4) * t606 + t553 * t525 + qJDD(4) - t493;
t485 = (-t537 + t596) * pkin(7) + (t552 * t553 + t561) * pkin(4) + t491;
t538 = pkin(4) * t598 - t553 * pkin(7);
t550 = t552 ^ 2;
t486 = -t550 * pkin(4) + t536 * pkin(7) - t538 * t598 + t490;
t575 = sin(qJ(5));
t578 = cos(qJ(5));
t483 = t578 * t485 - t575 * t486;
t523 = t578 * t552 - t575 * t553;
t498 = t523 * qJD(5) + t575 * t536 + t578 * t537;
t524 = t575 * t552 + t578 * t553;
t504 = -t523 * mrSges(6,1) + t524 * mrSges(6,2);
t566 = qJD(5) + t598;
t511 = -t566 * mrSges(6,2) + t523 * mrSges(6,3);
t557 = qJDD(5) + t561;
t481 = m(6) * t483 + t557 * mrSges(6,1) - t498 * mrSges(6,3) - t524 * t504 + t566 * t511;
t484 = t575 * t485 + t578 * t486;
t497 = -t524 * qJD(5) + t578 * t536 - t575 * t537;
t512 = t566 * mrSges(6,1) - t524 * mrSges(6,3);
t482 = m(6) * t484 - t557 * mrSges(6,2) + t497 * mrSges(6,3) + t523 * t504 - t566 * t512;
t591 = -t575 * t481 + t578 * t482;
t587 = m(5) * t490 - t561 * mrSges(5,3) + t591;
t526 = t552 * mrSges(5,1) - t553 * mrSges(5,3);
t601 = -t552 * mrSges(4,1) - t553 * mrSges(4,2) - t526;
t472 = m(4) * t494 + t561 * mrSges(4,2) + t601 * t552 + t611 * t536 + (t534 - t535) * t598 + t587;
t532 = -t552 * mrSges(5,2) - mrSges(5,3) * t598;
t533 = mrSges(4,2) * t598 - t552 * mrSges(4,3);
t474 = t578 * t481 + t575 * t482;
t586 = -m(5) * t491 - t561 * mrSges(5,1) - t474;
t473 = m(4) * t493 - t561 * mrSges(4,1) + t601 * t553 + t611 * t537 + (-t532 - t533) * t598 + t586;
t592 = t607 * t472 - t574 * t473;
t469 = m(3) * t529 - qJDD(2) * mrSges(3,2) + t561 * mrSges(3,3) - qJD(2) * t562 + t559 * t598 + t592;
t563 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t598;
t492 = t553 * t612 + (-t553 * t598 + t536) * pkin(3) + t613;
t488 = -t550 * pkin(7) + (-pkin(3) - pkin(4)) * t536 + (pkin(3) * t598 + (2 * qJD(4)) + t538) * t553 - t613;
t588 = -m(6) * t488 + t497 * mrSges(6,1) - t498 * mrSges(6,2) + t523 * t511 - t524 * t512;
t479 = m(5) * t492 + t536 * mrSges(5,1) - t537 * mrSges(5,3) + t552 * t532 - t553 * t535 + t588;
t583 = m(4) * t585 - t536 * mrSges(4,1) - t537 * mrSges(4,2) - t552 * t533 - t553 * t534 - t479;
t478 = m(3) * t528 + qJDD(2) * mrSges(3,1) - t560 * mrSges(3,3) + qJD(2) * t563 - t559 * t600 + t583;
t593 = t579 * t469 - t576 * t478;
t463 = m(2) * t565 - t582 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t593;
t470 = t574 * t472 + t607 * t473;
t584 = -m(3) * t547 + t561 * mrSges(3,1) - t560 * mrSges(3,2) - t562 * t600 + t563 * t598 - t470;
t466 = m(2) * t564 + qJDD(1) * mrSges(2,1) - t582 * mrSges(2,2) + t584;
t605 = t577 * t463 + t580 * t466;
t464 = t576 * t469 + t579 * t478;
t604 = t616 * t552 + t610 * t553 - t614 * t598;
t603 = t614 * t552 - t617 * t553 + t615 * t598;
t602 = t610 * t552 - t618 * t553 + t617 * t598;
t594 = t580 * t463 - t577 * t466;
t546 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t576 + Ifges(3,4) * t579) * qJD(1);
t545 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t576 + Ifges(3,2) * t579) * qJD(1);
t544 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t576 + Ifges(3,6) * t579) * qJD(1);
t501 = Ifges(6,1) * t524 + Ifges(6,4) * t523 + Ifges(6,5) * t566;
t500 = Ifges(6,4) * t524 + Ifges(6,2) * t523 + Ifges(6,6) * t566;
t499 = Ifges(6,5) * t524 + Ifges(6,6) * t523 + Ifges(6,3) * t566;
t476 = mrSges(6,2) * t488 - mrSges(6,3) * t483 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t557 + t523 * t499 - t566 * t500;
t475 = -mrSges(6,1) * t488 + mrSges(6,3) * t484 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t557 - t524 * t499 + t566 * t501;
t460 = -mrSges(4,2) * t585 + mrSges(5,2) * t491 - mrSges(4,3) * t493 - mrSges(5,3) * t492 - pkin(7) * t474 - qJ(4) * t479 - t575 * t475 + t578 * t476 - t610 * t536 + t618 * t537 + t603 * t552 - t561 * t617 + t604 * t598;
t459 = mrSges(4,1) * t585 - mrSges(5,1) * t492 + mrSges(5,2) * t490 + mrSges(4,3) * t494 - pkin(3) * t479 - pkin(4) * t588 - pkin(7) * t591 - t578 * t475 - t575 * t476 + t616 * t536 + t610 * t537 + t603 * t553 - t561 * t614 + t602 * t598;
t458 = (Ifges(3,2) + t615) * t561 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t560 + Ifges(6,3) * t557 + qJD(2) * t546 - mrSges(3,1) * t547 - t523 * t501 + t524 * t500 + mrSges(3,3) * t529 - mrSges(4,1) * t493 + mrSges(4,2) * t494 + Ifges(6,6) * t497 + Ifges(6,5) * t498 - mrSges(5,3) * t490 + mrSges(5,1) * t491 + mrSges(6,1) * t483 - mrSges(6,2) * t484 + pkin(4) * t474 - pkin(2) * t470 + (qJ(4) * mrSges(5,2) + t614) * t536 + (pkin(3) * mrSges(5,2) - t617) * t537 - qJ(4) * (-t535 * t598 + t587) - pkin(3) * (-t532 * t598 + t586) - t544 * t600 + (qJ(4) * t526 + t602) * t552 + (pkin(3) * t526 - t604) * t553;
t457 = mrSges(3,2) * t547 - mrSges(3,3) * t528 + Ifges(3,1) * t560 + Ifges(3,4) * t561 + Ifges(3,5) * qJDD(2) - qJ(3) * t470 - qJD(2) * t545 - t574 * t459 + t607 * t460 + t544 * t598;
t456 = Ifges(2,6) * qJDD(1) + t582 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t565 - Ifges(3,5) * t560 - Ifges(3,6) * t561 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t528 + mrSges(3,2) * t529 - t574 * t460 - t607 * t459 - pkin(2) * t583 - qJ(3) * t592 - pkin(1) * t464 + (-t576 * t545 + t579 * t546) * qJD(1);
t455 = -mrSges(2,2) * g(3) - mrSges(2,3) * t564 + Ifges(2,5) * qJDD(1) - t582 * Ifges(2,6) - pkin(6) * t464 + t579 * t457 - t576 * t458;
t1 = [-m(1) * g(1) + t594; -m(1) * g(2) + t605; (-m(1) - m(2)) * g(3) + t464; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t605 + t580 * t455 - t577 * t456; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t594 + t577 * t455 + t580 * t456; -mrSges(1,1) * g(2) + mrSges(2,1) * t564 + mrSges(1,2) * g(1) - mrSges(2,2) * t565 + Ifges(2,3) * qJDD(1) + pkin(1) * t584 + pkin(6) * t593 + t576 * t457 + t579 * t458;];
tauB = t1;
