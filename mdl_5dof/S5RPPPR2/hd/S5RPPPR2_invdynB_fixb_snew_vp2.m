% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:20
% EndTime: 2019-12-05 17:31:24
% DurationCPUTime: 4.84s
% Computational Cost: add. (36210->282), mult. (101411->386), div. (0->0), fcn. (68373->10), ass. (0->128)
t561 = sin(pkin(7));
t564 = cos(pkin(7));
t566 = sin(qJ(1));
t568 = cos(qJ(1));
t545 = t566 * g(2) - t568 * g(3);
t569 = qJD(1) ^ 2;
t614 = -t569 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t545;
t517 = -t561 * g(1) + t614 * t564;
t581 = -pkin(2) * t564 - qJ(3) * t561;
t540 = t581 * qJD(1);
t600 = t564 * qJD(1);
t508 = t540 * t600 + t517;
t560 = sin(pkin(8));
t563 = cos(pkin(8));
t546 = t568 * g(2) + t566 * g(3);
t574 = -t569 * qJ(2) + qJDD(2) - t546;
t601 = t561 * qJD(1);
t613 = (-pkin(1) + t581) * qJDD(1) + t574 - 0.2e1 * qJD(3) * t601;
t488 = -t560 * t508 + t613 * t563;
t516 = -t564 * g(1) - t614 * t561;
t559 = sin(pkin(9));
t562 = cos(pkin(9));
t604 = t561 * t563;
t576 = t559 * t604 + t562 * t564;
t529 = t576 * qJD(1);
t527 = t576 * qJDD(1);
t612 = 2 * qJD(4);
t611 = mrSges(3,2) * t561;
t610 = Ifges(4,4) * t563;
t609 = Ifges(3,6) * t564;
t608 = Ifges(4,6) * t564;
t607 = t561 ^ 2 * t569;
t606 = t564 ^ 2 * t569;
t605 = t560 * t561;
t603 = t564 * t569;
t541 = (-mrSges(3,1) * t564 + t611) * qJD(1);
t489 = t563 * t508 + t613 * t560;
t584 = mrSges(4,1) * t560 + mrSges(4,2) * t563;
t534 = t584 * t601;
t578 = -mrSges(4,1) * t564 - mrSges(4,3) * t604;
t538 = t578 * qJD(1);
t577 = mrSges(4,2) * t564 - mrSges(4,3) * t605;
t533 = (pkin(3) * t560 - qJ(4) * t563) * t601;
t595 = t560 * t601;
t597 = qJDD(1) * t564;
t486 = -pkin(3) * t606 - qJ(4) * t597 - t533 * t595 + t489;
t507 = t540 * t601 + qJDD(3) - t516;
t494 = ((-qJDD(1) * t563 - t560 * t603) * qJ(4) + (qJDD(1) * t560 - t563 * t603) * pkin(3)) * t561 + t507;
t482 = t562 * t486 + t559 * t494 - t529 * t612;
t594 = t563 * t601;
t530 = -t559 * t600 + t562 * t594;
t509 = t529 * mrSges(5,1) + t530 * mrSges(5,2);
t515 = mrSges(5,1) * t595 - t530 * mrSges(5,3);
t510 = t529 * pkin(4) - t530 * pkin(6);
t598 = qJDD(1) * t561;
t591 = t560 * t598;
t596 = t560 ^ 2 * t607;
t480 = -pkin(4) * t596 + pkin(6) * t591 - t529 * t510 + t482;
t485 = pkin(3) * t597 - qJ(4) * t606 + t533 * t594 + qJDD(4) - t488;
t528 = (-t559 * t564 + t562 * t604) * qJDD(1);
t483 = (t529 * t595 - t528) * pkin(6) + (t530 * t595 + t527) * pkin(4) + t485;
t565 = sin(qJ(5));
t567 = cos(qJ(5));
t477 = -t565 * t480 + t567 * t483;
t511 = -t565 * t530 + t567 * t595;
t512 = t567 * t530 + t565 * t595;
t496 = -t511 * mrSges(6,1) + t512 * mrSges(6,2);
t498 = t511 * qJD(5) + t567 * t528 + t565 * t591;
t526 = qJD(5) + t529;
t499 = -t526 * mrSges(6,2) + t511 * mrSges(6,3);
t525 = qJDD(5) + t527;
t475 = m(6) * t477 + t525 * mrSges(6,1) - t498 * mrSges(6,3) - t512 * t496 + t526 * t499;
t478 = t567 * t480 + t565 * t483;
t497 = -t512 * qJD(5) - t565 * t528 + t567 * t591;
t500 = t526 * mrSges(6,1) - t512 * mrSges(6,3);
t476 = m(6) * t478 - t525 * mrSges(6,2) + t497 * mrSges(6,3) + t511 * t496 - t526 * t500;
t586 = -t565 * t475 + t567 * t476;
t469 = m(5) * t482 - t527 * mrSges(5,3) - t529 * t509 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t515) * t605 + t586;
t579 = t559 * t486 - t562 * t494;
t481 = -0.2e1 * qJD(4) * t530 - t579;
t514 = -mrSges(5,2) * t595 - t529 * mrSges(5,3);
t479 = -pkin(4) * t591 - pkin(6) * t596 + (t612 + t510) * t530 + t579;
t573 = -m(6) * t479 + t497 * mrSges(6,1) - t498 * mrSges(6,2) + t511 * t499 - t512 * t500;
t473 = m(5) * t481 - t528 * mrSges(5,3) - t530 * t509 + (mrSges(5,1) * qJDD(1) + qJD(1) * t514) * t605 + t573;
t587 = t562 * t469 - t559 * t473;
t465 = m(4) * t489 + t577 * qJDD(1) + (-t534 * t605 + t538 * t564) * qJD(1) + t587;
t537 = t577 * qJD(1);
t470 = t567 * t475 + t565 * t476;
t570 = -m(5) * t485 - t527 * mrSges(5,1) - t528 * mrSges(5,2) - t529 * t514 - t530 * t515 - t470;
t467 = m(4) * t488 + t578 * qJDD(1) + (-t534 * t604 - t537 * t564) * qJD(1) + t570;
t588 = t563 * t465 - t560 * t467;
t458 = m(3) * t517 + (qJDD(1) * mrSges(3,3) + qJD(1) * t541) * t564 + t588;
t466 = t559 * t469 + t562 * t473;
t575 = -m(4) * t507 - t466;
t463 = m(3) * t516 + ((-mrSges(3,3) - t584) * qJDD(1) + (-t537 * t560 - t538 * t563 - t541) * qJD(1)) * t561 + t575;
t454 = t561 * t458 + t564 * t463;
t589 = t564 * t458 - t561 * t463;
t453 = m(2) * t545 - t569 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t589;
t459 = t560 * t465 + t563 * t467;
t536 = -qJDD(1) * pkin(1) + t574;
t571 = -m(3) * t536 + mrSges(3,1) * t597 - t459 + (t606 + t607) * mrSges(3,3);
t455 = m(2) * t546 - t569 * mrSges(2,2) + (mrSges(2,1) - t611) * qJDD(1) + t571;
t590 = t568 * t453 - t566 * t455;
t583 = Ifges(3,1) * t561 + Ifges(3,4) * t564;
t582 = Ifges(4,5) * t563 - Ifges(4,6) * t560;
t580 = -t566 * t453 - t568 * t455;
t572 = -Ifges(4,5) * t564 + (Ifges(4,1) * t563 - Ifges(4,4) * t560) * t561;
t542 = (Ifges(3,5) * t561 + t609) * qJD(1);
t522 = t572 * qJD(1);
t521 = (-t608 + (-Ifges(4,2) * t560 + t610) * t561) * qJD(1);
t520 = (-Ifges(4,3) * t564 + t582 * t561) * qJD(1);
t503 = Ifges(5,1) * t530 - Ifges(5,4) * t529 + Ifges(5,5) * t595;
t502 = Ifges(5,4) * t530 - Ifges(5,2) * t529 + Ifges(5,6) * t595;
t501 = Ifges(5,5) * t530 - Ifges(5,6) * t529 + Ifges(5,3) * t595;
t492 = Ifges(6,1) * t512 + Ifges(6,4) * t511 + Ifges(6,5) * t526;
t491 = Ifges(6,4) * t512 + Ifges(6,2) * t511 + Ifges(6,6) * t526;
t490 = Ifges(6,5) * t512 + Ifges(6,6) * t511 + Ifges(6,3) * t526;
t472 = mrSges(6,2) * t479 - mrSges(6,3) * t477 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t525 + t511 * t490 - t526 * t491;
t471 = -mrSges(6,1) * t479 + mrSges(6,3) * t478 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t525 - t512 * t490 + t526 * t492;
t461 = -mrSges(5,1) * t485 - mrSges(6,1) * t477 + mrSges(6,2) * t478 + mrSges(5,3) * t482 + Ifges(5,4) * t528 - Ifges(6,5) * t498 - Ifges(5,2) * t527 - Ifges(6,6) * t497 - Ifges(6,3) * t525 - pkin(4) * t470 - t512 * t491 + t511 * t492 - t530 * t501 + (Ifges(5,6) * qJDD(1) + qJD(1) * t503) * t605;
t460 = mrSges(5,2) * t485 - mrSges(5,3) * t481 + Ifges(5,1) * t528 - Ifges(5,4) * t527 - pkin(6) * t470 - t565 * t471 + t567 * t472 - t529 * t501 + (Ifges(5,5) * qJDD(1) - qJD(1) * t502) * t605;
t451 = -mrSges(4,1) * t507 + mrSges(4,3) * t489 - Ifges(5,5) * t528 + Ifges(5,6) * t527 - t530 * t502 - t529 * t503 - mrSges(5,1) * t481 + mrSges(5,2) * t482 - t565 * t472 - t567 * t471 - pkin(4) * t573 - pkin(6) * t586 - pkin(3) * t466 + (-t520 * t604 - t564 * t522) * qJD(1) + (-t608 + (t610 + (-Ifges(4,2) - Ifges(5,3)) * t560) * t561) * qJDD(1);
t450 = mrSges(4,2) * t507 - mrSges(4,3) * t488 - qJ(4) * t466 + t562 * t460 - t559 * t461 + (-t520 * t605 + t521 * t564) * qJD(1) + t572 * qJDD(1);
t449 = -mrSges(3,1) * t536 + mrSges(3,3) * t517 - mrSges(4,1) * t488 + mrSges(4,2) * t489 - t559 * t460 - t562 * t461 - pkin(3) * t570 - qJ(4) * t587 - pkin(2) * t459 + (Ifges(3,2) + Ifges(4,3)) * t597 + ((Ifges(3,4) - t582) * qJDD(1) + (-t521 * t563 - t522 * t560 - t542) * qJD(1)) * t561;
t448 = mrSges(3,2) * t536 - mrSges(3,3) * t516 - qJ(3) * t459 + t583 * qJDD(1) + t563 * t450 - t560 * t451 + t542 * t600;
t447 = t569 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t545 - mrSges(3,1) * t516 + mrSges(3,2) * t517 - t560 * t450 - t563 * t451 - pkin(2) * t575 - qJ(3) * t588 - pkin(1) * t454 + (-t609 + Ifges(2,6) + (pkin(2) * t584 - Ifges(3,5)) * t561) * qJDD(1) + (-pkin(2) * (-t537 * t605 - t538 * t604) + (-t561 * (Ifges(3,4) * t561 + Ifges(3,2) * t564) + t564 * t583) * qJD(1)) * qJD(1);
t446 = -mrSges(2,2) * g(1) - mrSges(2,3) * t546 + Ifges(2,5) * qJDD(1) - t569 * Ifges(2,6) - qJ(2) * t454 + t564 * t448 - t561 * t449;
t1 = [(-m(1) - m(2)) * g(1) + t454; -m(1) * g(2) + t580; -m(1) * g(3) + t590; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t546 - mrSges(2,2) * t545 + t561 * t448 + t564 * t449 + pkin(1) * (-mrSges(3,2) * t598 + t571) + qJ(2) * t589; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t590 - t566 * t446 - t568 * t447; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t580 + t568 * t446 - t566 * t447;];
tauB = t1;
