% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR3
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:23
% EndTime: 2019-12-05 16:19:27
% DurationCPUTime: 4.36s
% Computational Cost: add. (49940->258), mult. (110046->335), div. (0->0), fcn. (73622->10), ass. (0->104)
t485 = sin(pkin(8));
t487 = cos(pkin(8));
t472 = t485 * g(1) - t487 * g(2);
t473 = -t487 * g(1) - t485 * g(2);
t490 = sin(qJ(2));
t493 = cos(qJ(2));
t451 = t490 * t472 + t493 * t473;
t494 = qJD(2) ^ 2;
t447 = -t494 * pkin(2) + qJDD(2) * pkin(6) + t451;
t483 = -g(3) + qJDD(1);
t489 = sin(qJ(3));
t492 = cos(qJ(3));
t435 = -t489 * t447 + t492 * t483;
t506 = qJD(2) * qJD(3);
t504 = t492 * t506;
t470 = t489 * qJDD(2) + t504;
t430 = (-t470 + t504) * qJ(4) + (t489 * t492 * t494 + qJDD(3)) * pkin(3) + t435;
t436 = t492 * t447 + t489 * t483;
t471 = t492 * qJDD(2) - t489 * t506;
t508 = qJD(2) * t489;
t474 = qJD(3) * pkin(3) - qJ(4) * t508;
t482 = t492 ^ 2;
t431 = -t482 * t494 * pkin(3) + t471 * qJ(4) - qJD(3) * t474 + t436;
t484 = sin(pkin(9));
t486 = cos(pkin(9));
t459 = (t484 * t492 + t486 * t489) * qJD(2);
t413 = -0.2e1 * qJD(4) * t459 + t486 * t430 - t484 * t431;
t449 = t486 * t470 + t484 * t471;
t458 = (-t484 * t489 + t486 * t492) * qJD(2);
t411 = (qJD(3) * t458 - t449) * pkin(7) + (t458 * t459 + qJDD(3)) * pkin(4) + t413;
t414 = 0.2e1 * qJD(4) * t458 + t484 * t430 + t486 * t431;
t448 = -t484 * t470 + t486 * t471;
t454 = qJD(3) * pkin(4) - t459 * pkin(7);
t457 = t458 ^ 2;
t412 = -t457 * pkin(4) + t448 * pkin(7) - qJD(3) * t454 + t414;
t488 = sin(qJ(5));
t491 = cos(qJ(5));
t409 = t491 * t411 - t488 * t412;
t440 = t491 * t458 - t488 * t459;
t420 = t440 * qJD(5) + t488 * t448 + t491 * t449;
t441 = t488 * t458 + t491 * t459;
t426 = -t440 * mrSges(6,1) + t441 * mrSges(6,2);
t481 = qJD(3) + qJD(5);
t433 = -t481 * mrSges(6,2) + t440 * mrSges(6,3);
t480 = qJDD(3) + qJDD(5);
t407 = m(6) * t409 + t480 * mrSges(6,1) - t420 * mrSges(6,3) - t441 * t426 + t481 * t433;
t410 = t488 * t411 + t491 * t412;
t419 = -t441 * qJD(5) + t491 * t448 - t488 * t449;
t434 = t481 * mrSges(6,1) - t441 * mrSges(6,3);
t408 = m(6) * t410 - t480 * mrSges(6,2) + t419 * mrSges(6,3) + t440 * t426 - t481 * t434;
t399 = t491 * t407 + t488 * t408;
t443 = -t458 * mrSges(5,1) + t459 * mrSges(5,2);
t452 = -qJD(3) * mrSges(5,2) + t458 * mrSges(5,3);
t397 = m(5) * t413 + qJDD(3) * mrSges(5,1) - t449 * mrSges(5,3) + qJD(3) * t452 - t459 * t443 + t399;
t453 = qJD(3) * mrSges(5,1) - t459 * mrSges(5,3);
t499 = -t488 * t407 + t491 * t408;
t398 = m(5) * t414 - qJDD(3) * mrSges(5,2) + t448 * mrSges(5,3) - qJD(3) * t453 + t458 * t443 + t499;
t393 = t486 * t397 + t484 * t398;
t469 = (-mrSges(4,1) * t492 + mrSges(4,2) * t489) * qJD(2);
t507 = qJD(2) * t492;
t476 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t507;
t391 = m(4) * t435 + qJDD(3) * mrSges(4,1) - t470 * mrSges(4,3) + qJD(3) * t476 - t469 * t508 + t393;
t475 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t508;
t500 = -t484 * t397 + t486 * t398;
t392 = m(4) * t436 - qJDD(3) * mrSges(4,2) + t471 * mrSges(4,3) - qJD(3) * t475 + t469 * t507 + t500;
t501 = -t489 * t391 + t492 * t392;
t384 = m(3) * t451 - t494 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t501;
t450 = t493 * t472 - t490 * t473;
t497 = -qJDD(2) * pkin(2) - t450;
t446 = -t494 * pkin(6) + t497;
t432 = -t471 * pkin(3) + qJDD(4) + t474 * t508 + (-qJ(4) * t482 - pkin(6)) * t494 + t497;
t416 = -t448 * pkin(4) - t457 * pkin(7) + t459 * t454 + t432;
t498 = m(6) * t416 - t419 * mrSges(6,1) + t420 * mrSges(6,2) - t440 * t433 + t441 * t434;
t496 = m(5) * t432 - t448 * mrSges(5,1) + t449 * mrSges(5,2) - t458 * t452 + t459 * t453 + t498;
t495 = -m(4) * t446 + t471 * mrSges(4,1) - t470 * mrSges(4,2) - t475 * t508 + t476 * t507 - t496;
t403 = m(3) * t450 + qJDD(2) * mrSges(3,1) - t494 * mrSges(3,2) + t495;
t381 = t490 * t384 + t493 * t403;
t379 = m(2) * t472 + t381;
t502 = t493 * t384 - t490 * t403;
t380 = m(2) * t473 + t502;
t509 = t487 * t379 + t485 * t380;
t385 = t492 * t391 + t489 * t392;
t505 = m(3) * t483 + t385;
t503 = -t485 * t379 + t487 * t380;
t462 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t489 + Ifges(4,4) * t492) * qJD(2);
t461 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t489 + Ifges(4,2) * t492) * qJD(2);
t460 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t489 + Ifges(4,6) * t492) * qJD(2);
t439 = Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * qJD(3);
t438 = Ifges(5,4) * t459 + Ifges(5,2) * t458 + Ifges(5,6) * qJD(3);
t437 = Ifges(5,5) * t459 + Ifges(5,6) * t458 + Ifges(5,3) * qJD(3);
t423 = Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t481;
t422 = Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t481;
t421 = Ifges(6,5) * t441 + Ifges(6,6) * t440 + Ifges(6,3) * t481;
t401 = mrSges(6,2) * t416 - mrSges(6,3) * t409 + Ifges(6,1) * t420 + Ifges(6,4) * t419 + Ifges(6,5) * t480 + t440 * t421 - t481 * t422;
t400 = -mrSges(6,1) * t416 + mrSges(6,3) * t410 + Ifges(6,4) * t420 + Ifges(6,2) * t419 + Ifges(6,6) * t480 - t441 * t421 + t481 * t423;
t387 = mrSges(5,2) * t432 - mrSges(5,3) * t413 + Ifges(5,1) * t449 + Ifges(5,4) * t448 + Ifges(5,5) * qJDD(3) - pkin(7) * t399 - qJD(3) * t438 - t488 * t400 + t491 * t401 + t458 * t437;
t386 = -mrSges(5,1) * t432 + mrSges(5,3) * t414 + Ifges(5,4) * t449 + Ifges(5,2) * t448 + Ifges(5,6) * qJDD(3) - pkin(4) * t498 + pkin(7) * t499 + qJD(3) * t439 + t491 * t400 + t488 * t401 - t459 * t437;
t375 = mrSges(4,2) * t446 - mrSges(4,3) * t435 + Ifges(4,1) * t470 + Ifges(4,4) * t471 + Ifges(4,5) * qJDD(3) - qJ(4) * t393 - qJD(3) * t461 - t484 * t386 + t486 * t387 + t460 * t507;
t374 = -mrSges(4,1) * t446 + mrSges(4,3) * t436 + Ifges(4,4) * t470 + Ifges(4,2) * t471 + Ifges(4,6) * qJDD(3) - pkin(3) * t496 + qJ(4) * t500 + qJD(3) * t462 + t486 * t386 + t484 * t387 - t460 * t508;
t373 = Ifges(3,6) * qJDD(2) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t489 * t461 + t492 * t462) * qJD(2) + t494 * Ifges(3,5) - Ifges(6,3) * t480 - mrSges(3,1) * t483 - Ifges(4,5) * t470 - Ifges(4,6) * t471 + t458 * t439 - t459 * t438 - Ifges(5,6) * t448 - Ifges(5,5) * t449 + mrSges(3,3) * t451 - mrSges(4,1) * t435 + mrSges(4,2) * t436 + t440 * t423 - t441 * t422 - Ifges(6,6) * t419 - Ifges(6,5) * t420 - mrSges(5,1) * t413 + mrSges(5,2) * t414 - mrSges(6,1) * t409 + mrSges(6,2) * t410 - pkin(4) * t399 - pkin(3) * t393 - pkin(2) * t385;
t372 = mrSges(3,2) * t483 - mrSges(3,3) * t450 + Ifges(3,5) * qJDD(2) - t494 * Ifges(3,6) - pkin(6) * t385 - t489 * t374 + t492 * t375;
t371 = mrSges(2,2) * t483 - mrSges(2,3) * t472 - pkin(5) * t381 + t493 * t372 - t490 * t373;
t370 = -mrSges(2,1) * t483 + mrSges(2,3) * t473 - pkin(1) * t505 + pkin(5) * t502 + t490 * t372 + t493 * t373;
t1 = [-m(1) * g(1) + t503; -m(1) * g(2) + t509; -m(1) * g(3) + m(2) * t483 + t505; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t509 - t485 * t370 + t487 * t371; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t503 + t487 * t370 + t485 * t371; -mrSges(1,1) * g(2) + mrSges(2,1) * t472 + mrSges(3,1) * t450 + mrSges(1,2) * g(1) - mrSges(2,2) * t473 - mrSges(3,2) * t451 + Ifges(3,3) * qJDD(2) + pkin(1) * t381 + pkin(2) * t495 + pkin(6) * t501 + t492 * t374 + t489 * t375;];
tauB = t1;
