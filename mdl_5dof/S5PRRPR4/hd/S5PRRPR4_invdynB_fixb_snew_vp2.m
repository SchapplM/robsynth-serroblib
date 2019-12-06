% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:56
% EndTime: 2019-12-05 16:22:04
% DurationCPUTime: 4.28s
% Computational Cost: add. (46559->259), mult. (102327->335), div. (0->0), fcn. (68382->10), ass. (0->103)
t508 = cos(pkin(8));
t485 = sin(pkin(8));
t474 = -t508 * g(1) - t485 * g(2);
t483 = -g(3) + qJDD(1);
t489 = sin(qJ(2));
t492 = cos(qJ(2));
t454 = t492 * t474 + t489 * t483;
t493 = qJD(2) ^ 2;
t449 = -t493 * pkin(2) + qJDD(2) * pkin(6) + t454;
t473 = t485 * g(1) - t508 * g(2);
t488 = sin(qJ(3));
t491 = cos(qJ(3));
t435 = -t488 * t449 - t491 * t473;
t504 = qJD(2) * qJD(3);
t503 = t491 * t504;
t471 = t488 * qJDD(2) + t503;
t430 = (-t471 + t503) * qJ(4) + (t488 * t491 * t493 + qJDD(3)) * pkin(3) + t435;
t436 = t491 * t449 - t488 * t473;
t472 = t491 * qJDD(2) - t488 * t504;
t506 = qJD(2) * t488;
t475 = qJD(3) * pkin(3) - qJ(4) * t506;
t482 = t491 ^ 2;
t431 = -t482 * t493 * pkin(3) + t472 * qJ(4) - qJD(3) * t475 + t436;
t484 = sin(pkin(9));
t486 = cos(pkin(9));
t459 = (t484 * t491 + t486 * t488) * qJD(2);
t413 = -0.2e1 * qJD(4) * t459 + t486 * t430 - t484 * t431;
t446 = t486 * t471 + t484 * t472;
t458 = (-t484 * t488 + t486 * t491) * qJD(2);
t411 = (qJD(3) * t458 - t446) * pkin(7) + (t458 * t459 + qJDD(3)) * pkin(4) + t413;
t414 = 0.2e1 * qJD(4) * t458 + t484 * t430 + t486 * t431;
t445 = -t484 * t471 + t486 * t472;
t452 = qJD(3) * pkin(4) - t459 * pkin(7);
t457 = t458 ^ 2;
t412 = -t457 * pkin(4) + t445 * pkin(7) - qJD(3) * t452 + t414;
t487 = sin(qJ(5));
t490 = cos(qJ(5));
t409 = t490 * t411 - t487 * t412;
t440 = t490 * t458 - t487 * t459;
t420 = t440 * qJD(5) + t487 * t445 + t490 * t446;
t441 = t487 * t458 + t490 * t459;
t426 = -t440 * mrSges(6,1) + t441 * mrSges(6,2);
t481 = qJD(3) + qJD(5);
t433 = -t481 * mrSges(6,2) + t440 * mrSges(6,3);
t480 = qJDD(3) + qJDD(5);
t407 = m(6) * t409 + t480 * mrSges(6,1) - t420 * mrSges(6,3) - t441 * t426 + t481 * t433;
t410 = t487 * t411 + t490 * t412;
t419 = -t441 * qJD(5) + t490 * t445 - t487 * t446;
t434 = t481 * mrSges(6,1) - t441 * mrSges(6,3);
t408 = m(6) * t410 - t480 * mrSges(6,2) + t419 * mrSges(6,3) + t440 * t426 - t481 * t434;
t399 = t490 * t407 + t487 * t408;
t443 = -t458 * mrSges(5,1) + t459 * mrSges(5,2);
t450 = -qJD(3) * mrSges(5,2) + t458 * mrSges(5,3);
t397 = m(5) * t413 + qJDD(3) * mrSges(5,1) - t446 * mrSges(5,3) + qJD(3) * t450 - t459 * t443 + t399;
t451 = qJD(3) * mrSges(5,1) - t459 * mrSges(5,3);
t498 = -t487 * t407 + t490 * t408;
t398 = m(5) * t414 - qJDD(3) * mrSges(5,2) + t445 * mrSges(5,3) - qJD(3) * t451 + t458 * t443 + t498;
t393 = t486 * t397 + t484 * t398;
t470 = (-mrSges(4,1) * t491 + mrSges(4,2) * t488) * qJD(2);
t505 = qJD(2) * t491;
t477 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t505;
t391 = m(4) * t435 + qJDD(3) * mrSges(4,1) - t471 * mrSges(4,3) + qJD(3) * t477 - t470 * t506 + t393;
t476 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t506;
t499 = -t484 * t397 + t486 * t398;
t392 = m(4) * t436 - qJDD(3) * mrSges(4,2) + t472 * mrSges(4,3) - qJD(3) * t476 + t470 * t505 + t499;
t500 = -t488 * t391 + t491 * t392;
t384 = m(3) * t454 - t493 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t500;
t453 = -t489 * t474 + t492 * t483;
t496 = -qJDD(2) * pkin(2) - t453;
t448 = -t493 * pkin(6) + t496;
t432 = -t472 * pkin(3) + qJDD(4) + t475 * t506 + (-qJ(4) * t482 - pkin(6)) * t493 + t496;
t416 = -t445 * pkin(4) - t457 * pkin(7) + t459 * t452 + t432;
t497 = m(6) * t416 - t419 * mrSges(6,1) + t420 * mrSges(6,2) - t440 * t433 + t441 * t434;
t495 = m(5) * t432 - t445 * mrSges(5,1) + t446 * mrSges(5,2) - t458 * t450 + t459 * t451 + t497;
t494 = -m(4) * t448 + t472 * mrSges(4,1) - t471 * mrSges(4,2) - t476 * t506 + t477 * t505 - t495;
t403 = m(3) * t453 + qJDD(2) * mrSges(3,1) - t493 * mrSges(3,2) + t494;
t501 = t492 * t384 - t489 * t403;
t380 = m(2) * t474 + t501;
t387 = t491 * t391 + t488 * t392;
t386 = (m(2) + m(3)) * t473 - t387;
t507 = t485 * t380 + t508 * t386;
t381 = t489 * t384 + t492 * t403;
t502 = t508 * t380 - t485 * t386;
t462 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t488 + Ifges(4,4) * t491) * qJD(2);
t461 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t488 + Ifges(4,2) * t491) * qJD(2);
t460 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t488 + Ifges(4,6) * t491) * qJD(2);
t439 = Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * qJD(3);
t438 = Ifges(5,4) * t459 + Ifges(5,2) * t458 + Ifges(5,6) * qJD(3);
t437 = Ifges(5,5) * t459 + Ifges(5,6) * t458 + Ifges(5,3) * qJD(3);
t423 = Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t481;
t422 = Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t481;
t421 = Ifges(6,5) * t441 + Ifges(6,6) * t440 + Ifges(6,3) * t481;
t401 = mrSges(6,2) * t416 - mrSges(6,3) * t409 + Ifges(6,1) * t420 + Ifges(6,4) * t419 + Ifges(6,5) * t480 + t440 * t421 - t481 * t422;
t400 = -mrSges(6,1) * t416 + mrSges(6,3) * t410 + Ifges(6,4) * t420 + Ifges(6,2) * t419 + Ifges(6,6) * t480 - t441 * t421 + t481 * t423;
t389 = mrSges(5,2) * t432 - mrSges(5,3) * t413 + Ifges(5,1) * t446 + Ifges(5,4) * t445 + Ifges(5,5) * qJDD(3) - pkin(7) * t399 - qJD(3) * t438 - t487 * t400 + t490 * t401 + t458 * t437;
t388 = -mrSges(5,1) * t432 + mrSges(5,3) * t414 + Ifges(5,4) * t446 + Ifges(5,2) * t445 + Ifges(5,6) * qJDD(3) - pkin(4) * t497 + pkin(7) * t498 + qJD(3) * t439 + t490 * t400 + t487 * t401 - t459 * t437;
t377 = mrSges(4,2) * t448 - mrSges(4,3) * t435 + Ifges(4,1) * t471 + Ifges(4,4) * t472 + Ifges(4,5) * qJDD(3) - qJ(4) * t393 - qJD(3) * t461 - t484 * t388 + t486 * t389 + t460 * t505;
t376 = -mrSges(4,1) * t448 + mrSges(4,3) * t436 + Ifges(4,4) * t471 + Ifges(4,2) * t472 + Ifges(4,6) * qJDD(3) - pkin(3) * t495 + qJ(4) * t499 + qJD(3) * t462 + t486 * t388 + t484 * t389 - t460 * t506;
t375 = Ifges(3,6) * qJDD(2) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t488 * t461 + t491 * t462) * qJD(2) + t493 * Ifges(3,5) + mrSges(3,1) * t473 - Ifges(6,3) * t480 - t459 * t438 - Ifges(4,5) * t471 - Ifges(4,6) * t472 + mrSges(3,3) * t454 + t458 * t439 + t440 * t423 - t441 * t422 - Ifges(5,6) * t445 - Ifges(5,5) * t446 - mrSges(4,1) * t435 + mrSges(4,2) * t436 - Ifges(6,6) * t419 - Ifges(6,5) * t420 - mrSges(5,1) * t413 + mrSges(5,2) * t414 - mrSges(6,1) * t409 + mrSges(6,2) * t410 - pkin(4) * t399 - pkin(3) * t393 - pkin(2) * t387;
t374 = -mrSges(3,2) * t473 - mrSges(3,3) * t453 + Ifges(3,5) * qJDD(2) - t493 * Ifges(3,6) - pkin(6) * t387 - t488 * t376 + t491 * t377;
t373 = -mrSges(2,1) * t483 - mrSges(3,1) * t453 + mrSges(3,2) * t454 + mrSges(2,3) * t474 - Ifges(3,3) * qJDD(2) - pkin(1) * t381 - pkin(2) * t494 - pkin(6) * t500 - t491 * t376 - t488 * t377;
t372 = mrSges(2,2) * t483 - mrSges(2,3) * t473 - pkin(5) * t381 + t492 * t374 - t489 * t375;
t1 = [-m(1) * g(1) + t502; -m(1) * g(2) + t507; -m(1) * g(3) + m(2) * t483 + t381; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t507 + t508 * t372 - t485 * t373; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t502 + t485 * t372 + t508 * t373; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t473 - mrSges(2,2) * t474 + t489 * t374 + t492 * t375 + pkin(1) * (m(3) * t473 - t387) + pkin(5) * t501;];
tauB = t1;
