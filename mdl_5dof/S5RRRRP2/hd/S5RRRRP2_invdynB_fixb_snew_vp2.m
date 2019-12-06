% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:47
% EndTime: 2019-12-05 18:47:50
% DurationCPUTime: 2.62s
% Computational Cost: add. (37574->251), mult. (47662->304), div. (0->0), fcn. (28134->8), ass. (0->100)
t553 = Ifges(5,1) + Ifges(6,1);
t549 = Ifges(5,4) + Ifges(6,4);
t548 = Ifges(5,5) + Ifges(6,5);
t552 = Ifges(5,2) + Ifges(6,2);
t547 = -Ifges(5,6) - Ifges(6,6);
t551 = -Ifges(5,3) - Ifges(6,3);
t518 = qJD(1) + qJD(2);
t514 = t518 ^ 2;
t550 = pkin(3) * t514;
t521 = sin(qJ(3));
t546 = t518 * t521;
t525 = cos(qJ(3));
t545 = t518 * t525;
t523 = sin(qJ(1));
t527 = cos(qJ(1));
t509 = t527 * g(2) + t523 * g(3);
t503 = qJDD(1) * pkin(1) + t509;
t508 = t523 * g(2) - t527 * g(3);
t528 = qJD(1) ^ 2;
t504 = -t528 * pkin(1) + t508;
t522 = sin(qJ(2));
t526 = cos(qJ(2));
t481 = t522 * t503 + t526 * t504;
t516 = qJDD(1) + qJDD(2);
t479 = -t514 * pkin(2) + t516 * pkin(7) + t481;
t544 = t521 * t479;
t540 = qJD(3) * t518;
t498 = t521 * t516 + t525 * t540;
t454 = qJDD(3) * pkin(3) - t498 * pkin(8) - t544 + (pkin(8) * t540 + t521 * t550 - g(1)) * t525;
t465 = -t521 * g(1) + t525 * t479;
t499 = t525 * t516 - t521 * t540;
t507 = qJD(3) * pkin(3) - pkin(8) * t546;
t519 = t525 ^ 2;
t455 = t499 * pkin(8) - qJD(3) * t507 - t519 * t550 + t465;
t520 = sin(qJ(4));
t524 = cos(qJ(4));
t449 = t524 * t454 - t520 * t455;
t492 = (-t520 * t521 + t524 * t525) * t518;
t463 = t492 * qJD(4) + t524 * t498 + t520 * t499;
t493 = (t520 * t525 + t521 * t524) * t518;
t476 = -t492 * mrSges(6,1) + t493 * mrSges(6,2);
t477 = -t492 * mrSges(5,1) + t493 * mrSges(5,2);
t517 = qJD(3) + qJD(4);
t484 = -t517 * mrSges(5,2) + t492 * mrSges(5,3);
t515 = qJDD(3) + qJDD(4);
t444 = -0.2e1 * qJD(5) * t493 + (t492 * t517 - t463) * qJ(5) + (t492 * t493 + t515) * pkin(4) + t449;
t483 = -t517 * mrSges(6,2) + t492 * mrSges(6,3);
t539 = m(6) * t444 + t515 * mrSges(6,1) + t517 * t483;
t438 = m(5) * t449 + t515 * mrSges(5,1) + t517 * t484 + (-t476 - t477) * t493 + (-mrSges(5,3) - mrSges(6,3)) * t463 + t539;
t450 = t520 * t454 + t524 * t455;
t462 = -t493 * qJD(4) - t520 * t498 + t524 * t499;
t486 = t517 * mrSges(6,1) - t493 * mrSges(6,3);
t487 = t517 * mrSges(5,1) - t493 * mrSges(5,3);
t485 = t517 * pkin(4) - t493 * qJ(5);
t488 = t492 ^ 2;
t446 = -t488 * pkin(4) + t462 * qJ(5) + 0.2e1 * qJD(5) * t492 - t517 * t485 + t450;
t538 = m(6) * t446 + t462 * mrSges(6,3) + t492 * t476;
t441 = m(5) * t450 + t462 * mrSges(5,3) + t492 * t477 + (-t486 - t487) * t517 + (-mrSges(5,2) - mrSges(6,2)) * t515 + t538;
t434 = t524 * t438 + t520 * t441;
t464 = -t525 * g(1) - t544;
t497 = (-mrSges(4,1) * t525 + mrSges(4,2) * t521) * t518;
t506 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t545;
t432 = m(4) * t464 + qJDD(3) * mrSges(4,1) - t498 * mrSges(4,3) + qJD(3) * t506 - t497 * t546 + t434;
t505 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t546;
t534 = -t520 * t438 + t524 * t441;
t433 = m(4) * t465 - qJDD(3) * mrSges(4,2) + t499 * mrSges(4,3) - qJD(3) * t505 + t497 * t545 + t534;
t535 = -t521 * t432 + t525 * t433;
t425 = m(3) * t481 - t514 * mrSges(3,1) - t516 * mrSges(3,2) + t535;
t480 = t526 * t503 - t522 * t504;
t531 = -t516 * pkin(2) - t480;
t478 = -t514 * pkin(7) + t531;
t456 = -t499 * pkin(3) + t507 * t546 + (-pkin(8) * t519 - pkin(7)) * t514 + t531;
t448 = -t462 * pkin(4) - t488 * qJ(5) + t493 * t485 + qJDD(5) + t456;
t533 = m(6) * t448 - t462 * mrSges(6,1) + t463 * mrSges(6,2) - t492 * t483 + t493 * t486;
t530 = m(5) * t456 - t462 * mrSges(5,1) + t463 * mrSges(5,2) - t492 * t484 + t493 * t487 + t533;
t529 = -m(4) * t478 + t499 * mrSges(4,1) - t498 * mrSges(4,2) - t505 * t546 + t506 * t545 - t530;
t436 = m(3) * t480 + t516 * mrSges(3,1) - t514 * mrSges(3,2) + t529;
t422 = t522 * t425 + t526 * t436;
t426 = t525 * t432 + t521 * t433;
t543 = t547 * t492 - t548 * t493 + t551 * t517;
t542 = -t552 * t492 - t549 * t493 + t547 * t517;
t541 = t549 * t492 + t553 * t493 + t548 * t517;
t536 = t526 * t425 - t522 * t436;
t420 = m(2) * t508 - t528 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t536;
t421 = m(2) * t509 + qJDD(1) * mrSges(2,1) - t528 * mrSges(2,2) + t422;
t537 = t527 * t420 - t523 * t421;
t532 = -t523 * t420 - t527 * t421;
t491 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t521 + Ifges(4,4) * t525) * t518;
t490 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t521 + Ifges(4,2) * t525) * t518;
t489 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t521 + Ifges(4,6) * t525) * t518;
t442 = -t463 * mrSges(6,3) - t493 * t476 + t539;
t428 = mrSges(5,2) * t456 + mrSges(6,2) * t448 - mrSges(5,3) * t449 - mrSges(6,3) * t444 - qJ(5) * t442 + t549 * t462 + t553 * t463 - t543 * t492 + t548 * t515 + t542 * t517;
t427 = -mrSges(5,1) * t456 + mrSges(5,3) * t450 - mrSges(6,1) * t448 + mrSges(6,3) * t446 - pkin(4) * t533 + qJ(5) * t538 + (-qJ(5) * t486 + t541) * t517 + (-qJ(5) * mrSges(6,2) - t547) * t515 + t543 * t493 + t549 * t463 + t552 * t462;
t418 = mrSges(4,2) * t478 - mrSges(4,3) * t464 + Ifges(4,1) * t498 + Ifges(4,4) * t499 + Ifges(4,5) * qJDD(3) - pkin(8) * t434 - qJD(3) * t490 - t520 * t427 + t524 * t428 + t489 * t545;
t417 = -mrSges(4,1) * t478 + mrSges(4,3) * t465 + Ifges(4,4) * t498 + Ifges(4,2) * t499 + Ifges(4,6) * qJDD(3) - pkin(3) * t530 + pkin(8) * t534 + qJD(3) * t491 + t524 * t427 + t520 * t428 - t489 * t546;
t416 = mrSges(3,1) * g(1) - Ifges(4,3) * qJDD(3) + Ifges(3,6) * t516 + t514 * Ifges(3,5) - Ifges(4,5) * t498 - Ifges(4,6) * t499 + mrSges(3,3) * t481 - mrSges(4,1) * t464 + mrSges(4,2) * t465 + mrSges(6,2) * t446 - mrSges(5,1) * t449 + mrSges(5,2) * t450 - pkin(4) * t442 - mrSges(6,1) * t444 - pkin(3) * t434 - pkin(2) * t426 + (-t521 * t490 + t525 * t491) * t518 + t551 * t515 + t542 * t493 + t541 * t492 - t548 * t463 + t547 * t462;
t415 = -mrSges(3,2) * g(1) - mrSges(3,3) * t480 + Ifges(3,5) * t516 - t514 * Ifges(3,6) - pkin(7) * t426 - t521 * t417 + t525 * t418;
t414 = -mrSges(2,2) * g(1) - mrSges(2,3) * t509 + Ifges(2,5) * qJDD(1) - t528 * Ifges(2,6) - pkin(6) * t422 + t526 * t415 - t522 * t416;
t413 = Ifges(2,6) * qJDD(1) + t528 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t508 + t522 * t415 + t526 * t416 - pkin(1) * (-m(3) * g(1) + t426) + pkin(6) * t536;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t426; -m(1) * g(2) + t532; -m(1) * g(3) + t537; mrSges(2,1) * t509 + mrSges(3,1) * t480 - mrSges(1,2) * g(3) - mrSges(2,2) * t508 - mrSges(3,2) * t481 + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t516 + pkin(1) * t422 + pkin(2) * t529 + pkin(7) * t535 + t525 * t417 + t521 * t418; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t537 - t527 * t413 - t523 * t414; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t532 - t523 * t413 + t527 * t414;];
tauB = t1;
