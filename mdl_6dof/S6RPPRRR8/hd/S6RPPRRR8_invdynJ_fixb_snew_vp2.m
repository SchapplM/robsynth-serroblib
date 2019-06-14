% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:14:13
% EndTime: 2019-05-05 16:14:17
% DurationCPUTime: 3.09s
% Computational Cost: add. (27664->263), mult. (62761->333), div. (0->0), fcn. (44716->10), ass. (0->117)
t489 = qJD(1) ^ 2;
t483 = sin(qJ(1));
t487 = cos(qJ(1));
t507 = g(1) * t483 - t487 * g(2);
t497 = -qJ(2) * t489 + qJDD(2) - t507;
t517 = -pkin(1) - qJ(3);
t521 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t517 + t497;
t479 = cos(pkin(10));
t520 = t479 ^ 2;
t502 = -g(1) * t487 - g(2) * t483;
t519 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t502;
t518 = pkin(3) * t489;
t478 = sin(pkin(10));
t439 = t478 * g(3) + t521 * t479;
t423 = (-pkin(7) * qJDD(1) - t478 * t518) * t479 + t439;
t440 = -g(3) * t479 + t521 * t478;
t475 = t478 ^ 2;
t510 = t478 * qJDD(1);
t424 = -pkin(7) * t510 - t475 * t518 + t440;
t482 = sin(qJ(4));
t486 = cos(qJ(4));
t406 = t482 * t423 + t486 * t424;
t512 = t479 * qJD(1);
t513 = t478 * qJD(1);
t459 = -t482 * t512 - t486 * t513;
t500 = -t478 * t482 + t479 * t486;
t460 = t500 * qJD(1);
t436 = -mrSges(5,1) * t459 + mrSges(5,2) * t460;
t456 = t460 * qJD(4);
t509 = t479 * qJDD(1);
t442 = -t482 * t509 - t486 * t510 - t456;
t452 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t460;
t441 = -pkin(4) * t459 - pkin(8) * t460;
t488 = qJD(4) ^ 2;
t396 = -pkin(4) * t488 + qJDD(4) * pkin(8) + t441 * t459 + t406;
t496 = qJDD(3) + t519;
t515 = -t475 - t520;
t429 = pkin(3) * t510 + (pkin(7) * t515 + t517) * t489 + t496;
t514 = qJD(4) * t459;
t443 = qJDD(1) * t500 + t514;
t399 = (-t443 - t514) * pkin(8) + (-t442 + t456) * pkin(4) + t429;
t481 = sin(qJ(5));
t485 = cos(qJ(5));
t386 = -t396 * t481 + t485 * t399;
t446 = qJD(4) * t485 - t460 * t481;
t415 = qJD(5) * t446 + qJDD(4) * t481 + t443 * t485;
t438 = qJDD(5) - t442;
t447 = qJD(4) * t481 + t460 * t485;
t457 = qJD(5) - t459;
t383 = (t446 * t457 - t415) * pkin(9) + (t446 * t447 + t438) * pkin(5) + t386;
t387 = t485 * t396 + t481 * t399;
t414 = -qJD(5) * t447 + qJDD(4) * t485 - t443 * t481;
t427 = pkin(5) * t457 - pkin(9) * t447;
t445 = t446 ^ 2;
t384 = -pkin(5) * t445 + pkin(9) * t414 - t427 * t457 + t387;
t480 = sin(qJ(6));
t484 = cos(qJ(6));
t381 = t383 * t484 - t384 * t480;
t416 = t446 * t484 - t447 * t480;
t392 = qJD(6) * t416 + t414 * t480 + t415 * t484;
t417 = t446 * t480 + t447 * t484;
t404 = -mrSges(7,1) * t416 + mrSges(7,2) * t417;
t455 = qJD(6) + t457;
t407 = -mrSges(7,2) * t455 + mrSges(7,3) * t416;
t435 = qJDD(6) + t438;
t378 = m(7) * t381 + mrSges(7,1) * t435 - mrSges(7,3) * t392 - t404 * t417 + t407 * t455;
t382 = t383 * t480 + t384 * t484;
t391 = -qJD(6) * t417 + t414 * t484 - t415 * t480;
t408 = mrSges(7,1) * t455 - mrSges(7,3) * t417;
t379 = m(7) * t382 - mrSges(7,2) * t435 + mrSges(7,3) * t391 + t404 * t416 - t408 * t455;
t370 = t484 * t378 + t480 * t379;
t419 = -mrSges(6,1) * t446 + mrSges(6,2) * t447;
t425 = -mrSges(6,2) * t457 + mrSges(6,3) * t446;
t368 = m(6) * t386 + mrSges(6,1) * t438 - mrSges(6,3) * t415 - t419 * t447 + t425 * t457 + t370;
t426 = mrSges(6,1) * t457 - mrSges(6,3) * t447;
t503 = -t378 * t480 + t484 * t379;
t369 = m(6) * t387 - mrSges(6,2) * t438 + mrSges(6,3) * t414 + t419 * t446 - t426 * t457 + t503;
t504 = -t368 * t481 + t485 * t369;
t363 = m(5) * t406 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t442 - qJD(4) * t452 + t436 * t459 + t504;
t405 = t423 * t486 - t482 * t424;
t451 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t459;
t395 = -qJDD(4) * pkin(4) - pkin(8) * t488 + t460 * t441 - t405;
t385 = -pkin(5) * t414 - pkin(9) * t445 + t427 * t447 + t395;
t495 = m(7) * t385 - t391 * mrSges(7,1) + mrSges(7,2) * t392 - t416 * t407 + t408 * t417;
t491 = -m(6) * t395 + t414 * mrSges(6,1) - mrSges(6,2) * t415 + t446 * t425 - t426 * t447 - t495;
t374 = m(5) * t405 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t443 + qJD(4) * t451 - t436 * t460 + t491;
t516 = t482 * t363 + t486 * t374;
t364 = t485 * t368 + t481 * t369;
t506 = t515 * mrSges(4,3);
t505 = t486 * t363 - t374 * t482;
t499 = -qJDD(1) * mrSges(4,3) - t489 * (t478 * mrSges(4,1) + t479 * mrSges(4,2));
t501 = (m(4) * t439 + t479 * t499 + t516) * t479 + (m(4) * t440 + t478 * t499 + t505) * t478;
t494 = m(5) * t429 - mrSges(5,1) * t442 + t443 * mrSges(5,2) - t451 * t459 + t460 * t452 + t364;
t401 = Ifges(7,4) * t417 + Ifges(7,2) * t416 + Ifges(7,6) * t455;
t402 = Ifges(7,1) * t417 + Ifges(7,4) * t416 + Ifges(7,5) * t455;
t493 = -mrSges(7,1) * t381 + mrSges(7,2) * t382 - Ifges(7,5) * t392 - Ifges(7,6) * t391 - Ifges(7,3) * t435 - t417 * t401 + t416 * t402;
t450 = t489 * t517 + t496;
t492 = m(4) * t450 + mrSges(4,1) * t510 + mrSges(4,2) * t509 + t494;
t410 = Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t457;
t411 = Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t457;
t490 = mrSges(6,1) * t386 - mrSges(6,2) * t387 + Ifges(6,5) * t415 + Ifges(6,6) * t414 + Ifges(6,3) * t438 + pkin(5) * t370 + t447 * t410 - t446 * t411 - t493;
t462 = (Ifges(4,5) * t479 - Ifges(4,6) * t478) * qJD(1);
t458 = -qJDD(1) * pkin(1) + t497;
t454 = pkin(1) * t489 - t519;
t432 = Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * qJD(4);
t431 = Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * qJD(4);
t430 = Ifges(5,5) * t460 + Ifges(5,6) * t459 + Ifges(5,3) * qJD(4);
t409 = Ifges(6,5) * t447 + Ifges(6,6) * t446 + Ifges(6,3) * t457;
t400 = Ifges(7,5) * t417 + Ifges(7,6) * t416 + Ifges(7,3) * t455;
t372 = mrSges(7,2) * t385 - mrSges(7,3) * t381 + Ifges(7,1) * t392 + Ifges(7,4) * t391 + Ifges(7,5) * t435 + t400 * t416 - t401 * t455;
t371 = -mrSges(7,1) * t385 + mrSges(7,3) * t382 + Ifges(7,4) * t392 + Ifges(7,2) * t391 + Ifges(7,6) * t435 - t400 * t417 + t402 * t455;
t360 = mrSges(6,2) * t395 - mrSges(6,3) * t386 + Ifges(6,1) * t415 + Ifges(6,4) * t414 + Ifges(6,5) * t438 - pkin(9) * t370 - t371 * t480 + t372 * t484 + t409 * t446 - t410 * t457;
t359 = -mrSges(6,1) * t395 + mrSges(6,3) * t387 + Ifges(6,4) * t415 + Ifges(6,2) * t414 + Ifges(6,6) * t438 - pkin(5) * t495 + pkin(9) * t503 + t484 * t371 + t480 * t372 - t447 * t409 + t457 * t411;
t356 = -mrSges(5,1) * t429 + mrSges(5,3) * t406 + Ifges(5,4) * t443 + Ifges(5,2) * t442 + Ifges(5,6) * qJDD(4) - pkin(4) * t364 + qJD(4) * t432 - t460 * t430 - t490;
t355 = m(3) * t458 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t489 + t501;
t354 = mrSges(5,2) * t429 - mrSges(5,3) * t405 + Ifges(5,1) * t443 + Ifges(5,4) * t442 + Ifges(5,5) * qJDD(4) - pkin(8) * t364 - qJD(4) * t431 - t359 * t481 + t360 * t485 + t430 * t459;
t1 = [mrSges(2,1) * t507 - mrSges(2,2) * t502 + mrSges(3,2) * t458 - mrSges(3,3) * t454 + t479 * (mrSges(4,2) * t450 - mrSges(4,3) * t439 - pkin(7) * t516 + t486 * t354 - t482 * t356 - t462 * t513) - t478 * (-mrSges(4,1) * t450 + mrSges(4,3) * t440 - pkin(3) * t494 + pkin(7) * t505 + t482 * t354 + t486 * t356 - t462 * t512) - qJ(3) * t501 - pkin(1) * t355 + qJ(2) * (-m(3) * t454 + (mrSges(3,2) + t506) * t489 + t492) + (Ifges(4,1) * t520 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t479 + Ifges(4,2) * t478) * t478) * qJDD(1); t355; t489 * t506 + t492; mrSges(5,1) * t405 - mrSges(5,2) * t406 + Ifges(5,5) * t443 + Ifges(5,6) * t442 + Ifges(5,3) * qJDD(4) + pkin(4) * t491 + pkin(8) * t504 + t485 * t359 + t481 * t360 + t460 * t431 - t459 * t432; t490; -t493;];
tauJ  = t1;
