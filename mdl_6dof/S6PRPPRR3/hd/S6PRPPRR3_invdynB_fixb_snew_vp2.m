% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-04 22:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:00:13
% EndTime: 2019-05-04 22:00:18
% DurationCPUTime: 4.43s
% Computational Cost: add. (60419->252), mult. (104366->314), div. (0->0), fcn. (62393->12), ass. (0->113)
t481 = sin(pkin(10));
t484 = cos(pkin(10));
t467 = g(1) * t481 - g(2) * t484;
t468 = -g(1) * t484 - g(2) * t481;
t477 = -g(3) + qJDD(1);
t482 = sin(pkin(6));
t485 = cos(pkin(6));
t488 = sin(qJ(2));
t491 = cos(qJ(2));
t436 = -t488 * t468 + (t467 * t485 + t477 * t482) * t491;
t520 = -pkin(2) - pkin(3);
t519 = -mrSges(3,1) - mrSges(4,1);
t518 = Ifges(4,4) + Ifges(3,5);
t517 = Ifges(3,6) - Ifges(4,6);
t493 = qJD(2) ^ 2;
t514 = t485 * t488;
t515 = t482 * t488;
t437 = t467 * t514 + t491 * t468 + t477 * t515;
t502 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t437;
t431 = t493 * t520 + t502;
t495 = -t493 * qJ(3) + qJDD(3) - t436;
t433 = qJDD(2) * t520 + t495;
t480 = sin(pkin(11));
t483 = cos(pkin(11));
t427 = t483 * t431 + t480 * t433;
t425 = -pkin(4) * t493 - qJDD(2) * pkin(8) + t427;
t451 = -t467 * t482 + t477 * t485;
t449 = qJDD(4) - t451;
t487 = sin(qJ(5));
t490 = cos(qJ(5));
t422 = t490 * t425 + t487 * t449;
t463 = (mrSges(6,1) * t490 - mrSges(6,2) * t487) * qJD(2);
t509 = qJD(2) * qJD(5);
t508 = t487 * t509;
t466 = -qJDD(2) * t490 + t508;
t511 = qJD(2) * t487;
t469 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t511;
t464 = (pkin(5) * t490 + pkin(9) * t487) * qJD(2);
t492 = qJD(5) ^ 2;
t510 = qJD(2) * t490;
t419 = -pkin(5) * t492 + qJDD(5) * pkin(9) - t464 * t510 + t422;
t426 = -t480 * t431 + t483 * t433;
t424 = qJDD(2) * pkin(4) - t493 * pkin(8) - t426;
t507 = t490 * t509;
t465 = -qJDD(2) * t487 - t507;
t420 = (-t465 + t507) * pkin(9) + (-t466 - t508) * pkin(5) + t424;
t486 = sin(qJ(6));
t489 = cos(qJ(6));
t416 = -t419 * t486 + t420 * t489;
t461 = qJD(5) * t489 + t486 * t511;
t444 = qJD(6) * t461 + qJDD(5) * t486 + t465 * t489;
t462 = qJD(5) * t486 - t489 * t511;
t445 = -mrSges(7,1) * t461 + mrSges(7,2) * t462;
t472 = qJD(6) + t510;
t447 = -mrSges(7,2) * t472 + mrSges(7,3) * t461;
t458 = qJDD(6) - t466;
t414 = m(7) * t416 + mrSges(7,1) * t458 - mrSges(7,3) * t444 - t445 * t462 + t447 * t472;
t417 = t419 * t489 + t420 * t486;
t443 = -qJD(6) * t462 + qJDD(5) * t489 - t465 * t486;
t448 = mrSges(7,1) * t472 - mrSges(7,3) * t462;
t415 = m(7) * t417 - mrSges(7,2) * t458 + mrSges(7,3) * t443 + t445 * t461 - t448 * t472;
t503 = -t414 * t486 + t489 * t415;
t408 = m(6) * t422 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t466 - qJD(5) * t469 - t463 * t510 + t503;
t513 = t490 * t449;
t421 = -t425 * t487 + t513;
t470 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t510;
t418 = -qJDD(5) * pkin(5) - t492 * pkin(9) - t513 + (-qJD(2) * t464 + t425) * t487;
t497 = -m(7) * t418 + t443 * mrSges(7,1) - mrSges(7,2) * t444 + t461 * t447 - t448 * t462;
t412 = m(6) * t421 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t465 + qJD(5) * t470 + t463 * t511 + t497;
t504 = t490 * t408 - t412 * t487;
t401 = m(5) * t427 - mrSges(5,1) * t493 + qJDD(2) * mrSges(5,2) + t504;
t409 = t414 * t489 + t415 * t486;
t494 = -m(6) * t424 + t466 * mrSges(6,1) - mrSges(6,2) * t465 + t469 * t511 - t470 * t510 - t409;
t406 = m(5) * t426 - qJDD(2) * mrSges(5,1) - mrSges(5,2) * t493 + t494;
t397 = t401 * t480 + t406 * t483;
t435 = -qJDD(2) * pkin(2) + t495;
t496 = -m(4) * t435 + qJDD(2) * mrSges(4,1) + t493 * mrSges(4,3) - t397;
t396 = m(3) * t436 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t493 + t496;
t516 = t396 * t491;
t434 = -pkin(2) * t493 + t502;
t505 = t483 * t401 - t480 * t406;
t500 = m(4) * t434 + qJDD(2) * mrSges(4,3) + t505;
t395 = m(3) * t437 - qJDD(2) * mrSges(3,2) + t493 * t519 + t500;
t405 = t408 * t487 + t412 * t490;
t499 = -m(5) * t449 - t405;
t404 = m(4) * t451 + t499;
t403 = m(3) * t451 + t404;
t383 = t395 * t514 - t403 * t482 + t485 * t516;
t381 = m(2) * t467 + t383;
t388 = t491 * t395 - t396 * t488;
t387 = m(2) * t468 + t388;
t512 = t484 * t381 + t481 * t387;
t382 = t395 * t515 + t485 * t403 + t482 * t516;
t506 = -t381 * t481 + t484 * t387;
t438 = Ifges(7,5) * t462 + Ifges(7,6) * t461 + Ifges(7,3) * t472;
t440 = Ifges(7,1) * t462 + Ifges(7,4) * t461 + Ifges(7,5) * t472;
t410 = -mrSges(7,1) * t418 + mrSges(7,3) * t417 + Ifges(7,4) * t444 + Ifges(7,2) * t443 + Ifges(7,6) * t458 - t438 * t462 + t440 * t472;
t439 = Ifges(7,4) * t462 + Ifges(7,2) * t461 + Ifges(7,6) * t472;
t411 = mrSges(7,2) * t418 - mrSges(7,3) * t416 + Ifges(7,1) * t444 + Ifges(7,4) * t443 + Ifges(7,5) * t458 + t438 * t461 - t439 * t472;
t453 = (Ifges(6,3) * qJD(5)) + (-Ifges(6,5) * t487 - Ifges(6,6) * t490) * qJD(2);
t454 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t487 - Ifges(6,2) * t490) * qJD(2);
t398 = mrSges(6,2) * t424 - mrSges(6,3) * t421 + Ifges(6,1) * t465 + Ifges(6,4) * t466 + Ifges(6,5) * qJDD(5) - pkin(9) * t409 - qJD(5) * t454 - t410 * t486 + t411 * t489 - t453 * t510;
t455 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t487 - Ifges(6,4) * t490) * qJD(2);
t399 = -mrSges(6,1) * t424 - mrSges(7,1) * t416 + mrSges(7,2) * t417 + mrSges(6,3) * t422 + Ifges(6,4) * t465 - Ifges(7,5) * t444 + Ifges(6,2) * t466 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t443 - Ifges(7,3) * t458 - pkin(5) * t409 + qJD(5) * t455 - t439 * t462 + t440 * t461 + t453 * t511;
t384 = mrSges(5,2) * t449 - mrSges(5,3) * t426 - Ifges(5,5) * qJDD(2) - Ifges(5,6) * t493 - pkin(8) * t405 + t398 * t490 - t399 * t487;
t389 = -Ifges(5,6) * qJDD(2) + t493 * Ifges(5,5) - mrSges(5,1) * t449 + mrSges(5,3) * t427 - Ifges(6,5) * t465 - Ifges(6,6) * t466 - Ifges(6,3) * qJDD(5) - mrSges(6,1) * t421 + mrSges(6,2) * t422 - t486 * t411 - t489 * t410 - pkin(5) * t497 - pkin(9) * t503 - pkin(4) * t405 + (t454 * t487 - t455 * t490) * qJD(2);
t377 = mrSges(4,2) * t434 + mrSges(3,3) * t437 - pkin(2) * t404 - pkin(3) * t499 - qJ(4) * t505 + qJDD(2) * t517 - t480 * t384 - t483 * t389 + t451 * t519 + t493 * t518;
t379 = mrSges(4,2) * t435 - mrSges(3,3) * t436 - qJ(3) * t404 - qJ(4) * t397 + t483 * t384 - t480 * t389 - t517 * t493 + (mrSges(3,2) - mrSges(4,3)) * t451 + t518 * qJDD(2);
t498 = pkin(7) * t388 + t377 * t491 + t379 * t488;
t378 = pkin(2) * t496 + qJ(3) * (-mrSges(4,1) * t493 + t500) + mrSges(3,1) * t436 - mrSges(3,2) * t437 - pkin(3) * t397 - mrSges(4,1) * t435 + mrSges(4,3) * t434 - pkin(8) * t504 - mrSges(5,1) * t426 + mrSges(5,2) * t427 - t487 * t398 - t490 * t399 - pkin(4) * t494 + (Ifges(3,3) + Ifges(4,2) + Ifges(5,3)) * qJDD(2);
t376 = mrSges(2,2) * t477 - mrSges(2,3) * t467 - t488 * t377 + t491 * t379 + (-t382 * t482 - t383 * t485) * pkin(7);
t375 = -mrSges(2,1) * t477 + mrSges(2,3) * t468 - pkin(1) * t382 - t482 * t378 + t485 * t498;
t1 = [-m(1) * g(1) + t506; -m(1) * g(2) + t512; -m(1) * g(3) + m(2) * t477 + t382; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t512 - t481 * t375 + t484 * t376; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t506 + t484 * t375 + t481 * t376; -mrSges(1,1) * g(2) + mrSges(2,1) * t467 + mrSges(1,2) * g(1) - mrSges(2,2) * t468 + pkin(1) * t383 + t485 * t378 + t482 * t498;];
tauB  = t1;
