% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:37
% EndTime: 2019-12-05 16:35:49
% DurationCPUTime: 5.65s
% Computational Cost: add. (62722->269), mult. (128041->350), div. (0->0), fcn. (87146->12), ass. (0->118)
t493 = sin(pkin(9));
t496 = cos(pkin(9));
t484 = t493 * g(1) - t496 * g(2);
t485 = -t496 * g(1) - t493 * g(2);
t491 = -g(3) + qJDD(1);
t494 = sin(pkin(5));
t497 = cos(pkin(5));
t500 = sin(qJ(2));
t503 = cos(qJ(2));
t449 = -t500 * t485 + (t484 * t497 + t491 * t494) * t503;
t528 = 2 * qJD(4);
t505 = qJD(2) ^ 2;
t524 = t497 * t500;
t525 = t494 * t500;
t450 = t484 * t524 + t503 * t485 + t491 * t525;
t446 = -t505 * pkin(2) + qJDD(2) * pkin(7) + t450;
t466 = -t494 * t484 + t497 * t491;
t499 = sin(qJ(3));
t502 = cos(qJ(3));
t434 = t502 * t446 + t499 * t466;
t480 = (-pkin(3) * t502 - qJ(4) * t499) * qJD(2);
t504 = qJD(3) ^ 2;
t520 = qJD(2) * t502;
t430 = -t504 * pkin(3) + qJDD(3) * qJ(4) + t480 * t520 + t434;
t445 = -qJDD(2) * pkin(2) - t505 * pkin(7) - t449;
t519 = qJD(2) * qJD(3);
t517 = t502 * t519;
t482 = t499 * qJDD(2) + t517;
t518 = t499 * t519;
t483 = t502 * qJDD(2) - t518;
t432 = (-t482 - t517) * qJ(4) + (-t483 + t518) * pkin(3) + t445;
t492 = sin(pkin(10));
t495 = cos(pkin(10));
t521 = qJD(2) * t499;
t475 = t495 * qJD(3) - t492 * t521;
t426 = t495 * t430 + t492 * t432 + t475 * t528;
t476 = t492 * qJD(3) + t495 * t521;
t455 = -t475 * mrSges(5,1) + t476 * mrSges(5,2);
t463 = -mrSges(5,1) * t520 - t476 * mrSges(5,3);
t464 = t495 * qJDD(3) - t492 * t482;
t456 = -t475 * pkin(4) - t476 * pkin(8);
t526 = t502 ^ 2 * t505;
t424 = -pkin(4) * t526 - t483 * pkin(8) + t475 * t456 + t426;
t465 = t492 * qJDD(3) + t495 * t482;
t443 = t499 * t446;
t509 = -qJDD(3) * pkin(3) - t504 * qJ(4) + t480 * t521 + qJDD(4) + t443;
t427 = -t464 * pkin(4) - t465 * pkin(8) + (-t466 + (-pkin(4) * t476 + pkin(8) * t475) * qJD(2)) * t502 + t509;
t498 = sin(qJ(5));
t501 = cos(qJ(5));
t421 = -t498 * t424 + t501 * t427;
t457 = -t498 * t476 - t501 * t520;
t441 = t457 * qJD(5) + t501 * t465 - t498 * t483;
t458 = t501 * t476 - t498 * t520;
t442 = -t457 * mrSges(6,1) + t458 * mrSges(6,2);
t473 = qJD(5) - t475;
t447 = -t473 * mrSges(6,2) + t457 * mrSges(6,3);
t461 = qJDD(5) - t464;
t419 = m(6) * t421 + t461 * mrSges(6,1) - t441 * mrSges(6,3) - t458 * t442 + t473 * t447;
t422 = t501 * t424 + t498 * t427;
t440 = -t458 * qJD(5) - t498 * t465 - t501 * t483;
t448 = t473 * mrSges(6,1) - t458 * mrSges(6,3);
t420 = m(6) * t422 - t461 * mrSges(6,2) + t440 * mrSges(6,3) + t457 * t442 - t473 * t448;
t513 = -t498 * t419 + t501 * t420;
t413 = m(5) * t426 + t483 * mrSges(5,2) + t464 * mrSges(5,3) + t475 * t455 + t463 * t520 + t513;
t512 = t492 * t430 - t495 * t432;
t425 = -0.2e1 * qJD(4) * t476 - t512;
t462 = mrSges(5,2) * t520 + t475 * mrSges(5,3);
t423 = -pkin(8) * t526 + t483 * pkin(4) + (t528 + t456) * t476 + t512;
t508 = -m(6) * t423 + t440 * mrSges(6,1) - t441 * mrSges(6,2) + t457 * t447 - t458 * t448;
t417 = m(5) * t425 - t483 * mrSges(5,1) - t465 * mrSges(5,3) - t476 * t455 - t462 * t520 + t508;
t409 = t492 * t413 + t495 * t417;
t486 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t521;
t487 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t520;
t506 = -m(4) * t445 + t483 * mrSges(4,1) - t482 * mrSges(4,2) - t486 * t521 + t487 * t520 - t409;
t405 = m(3) * t449 + qJDD(2) * mrSges(3,1) - t505 * mrSges(3,2) + t506;
t527 = t405 * t503;
t523 = t502 * t466;
t481 = (-mrSges(4,1) * t502 + mrSges(4,2) * t499) * qJD(2);
t514 = t495 * t413 - t492 * t417;
t408 = m(4) * t434 - qJDD(3) * mrSges(4,2) + t483 * mrSges(4,3) - qJD(3) * t486 + t481 * t520 + t514;
t433 = -t443 + t523;
t414 = t501 * t419 + t498 * t420;
t429 = t509 - t523;
t507 = -m(5) * t429 + t464 * mrSges(5,1) - t465 * mrSges(5,2) + t475 * t462 - t476 * t463 - t414;
t411 = m(4) * t433 + qJDD(3) * mrSges(4,1) - t482 * mrSges(4,3) + qJD(3) * t487 - t481 * t521 + t507;
t515 = t502 * t408 - t499 * t411;
t397 = m(3) * t450 - t505 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t515;
t400 = t499 * t408 + t502 * t411;
t399 = m(3) * t466 + t400;
t387 = t397 * t524 - t494 * t399 + t497 * t527;
t385 = m(2) * t484 + t387;
t393 = t503 * t397 - t500 * t405;
t392 = m(2) * t485 + t393;
t522 = t496 * t385 + t493 * t392;
t386 = t397 * t525 + t497 * t399 + t494 * t527;
t516 = -t493 * t385 + t496 * t392;
t435 = Ifges(6,5) * t458 + Ifges(6,6) * t457 + Ifges(6,3) * t473;
t437 = Ifges(6,1) * t458 + Ifges(6,4) * t457 + Ifges(6,5) * t473;
t415 = -mrSges(6,1) * t423 + mrSges(6,3) * t422 + Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t461 - t458 * t435 + t473 * t437;
t436 = Ifges(6,4) * t458 + Ifges(6,2) * t457 + Ifges(6,6) * t473;
t416 = mrSges(6,2) * t423 - mrSges(6,3) * t421 + Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t461 + t457 * t435 - t473 * t436;
t451 = Ifges(5,5) * t476 + Ifges(5,6) * t475 - Ifges(5,3) * t520;
t452 = Ifges(5,4) * t476 + Ifges(5,2) * t475 - Ifges(5,6) * t520;
t401 = mrSges(5,2) * t429 - mrSges(5,3) * t425 + Ifges(5,1) * t465 + Ifges(5,4) * t464 - Ifges(5,5) * t483 - pkin(8) * t414 - t498 * t415 + t501 * t416 + t475 * t451 + t452 * t520;
t453 = Ifges(5,1) * t476 + Ifges(5,4) * t475 - Ifges(5,5) * t520;
t402 = -mrSges(5,1) * t429 - mrSges(6,1) * t421 + mrSges(6,2) * t422 + mrSges(5,3) * t426 + Ifges(5,4) * t465 - Ifges(6,5) * t441 + Ifges(5,2) * t464 - Ifges(5,6) * t483 - Ifges(6,6) * t440 - Ifges(6,3) * t461 - pkin(4) * t414 - t458 * t436 + t457 * t437 - t476 * t451 - t453 * t520;
t470 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t499 + Ifges(4,6) * t502) * qJD(2);
t471 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t499 + Ifges(4,2) * t502) * qJD(2);
t388 = mrSges(4,2) * t445 - mrSges(4,3) * t433 + Ifges(4,1) * t482 + Ifges(4,4) * t483 + Ifges(4,5) * qJDD(3) - qJ(4) * t409 - qJD(3) * t471 + t495 * t401 - t492 * t402 + t470 * t520;
t472 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t499 + Ifges(4,4) * t502) * qJD(2);
t389 = Ifges(4,4) * t482 + Ifges(4,6) * qJDD(3) - t470 * t521 + qJD(3) * t472 - mrSges(4,1) * t445 + mrSges(4,3) * t434 - Ifges(5,5) * t465 - Ifges(5,6) * t464 - t476 * t452 + t475 * t453 - mrSges(5,1) * t425 + mrSges(5,2) * t426 - t498 * t416 - t501 * t415 - pkin(4) * t508 - pkin(8) * t513 - pkin(3) * t409 + (Ifges(4,2) + Ifges(5,3)) * t483;
t382 = mrSges(3,2) * t466 - mrSges(3,3) * t449 + Ifges(3,5) * qJDD(2) - t505 * Ifges(3,6) - pkin(7) * t400 + t502 * t388 - t499 * t389;
t383 = Ifges(3,6) * qJDD(2) + t505 * Ifges(3,5) - mrSges(3,1) * t466 + mrSges(3,3) * t450 - Ifges(4,5) * t482 - Ifges(4,6) * t483 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t433 + mrSges(4,2) * t434 - t492 * t401 - t495 * t402 - pkin(3) * t507 - qJ(4) * t514 - pkin(2) * t400 + (-t499 * t471 + t502 * t472) * qJD(2);
t510 = pkin(6) * t393 + t382 * t500 + t383 * t503;
t381 = mrSges(3,1) * t449 - mrSges(3,2) * t450 + Ifges(3,3) * qJDD(2) + pkin(2) * t506 + pkin(7) * t515 + t499 * t388 + t502 * t389;
t380 = mrSges(2,2) * t491 - mrSges(2,3) * t484 + t503 * t382 - t500 * t383 + (-t386 * t494 - t387 * t497) * pkin(6);
t379 = -mrSges(2,1) * t491 + mrSges(2,3) * t485 - pkin(1) * t386 - t494 * t381 + t510 * t497;
t1 = [-m(1) * g(1) + t516; -m(1) * g(2) + t522; -m(1) * g(3) + m(2) * t491 + t386; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t522 - t493 * t379 + t496 * t380; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t516 + t496 * t379 + t493 * t380; -mrSges(1,1) * g(2) + mrSges(2,1) * t484 + mrSges(1,2) * g(1) - mrSges(2,2) * t485 + pkin(1) * t387 + t497 * t381 + t510 * t494;];
tauB = t1;
