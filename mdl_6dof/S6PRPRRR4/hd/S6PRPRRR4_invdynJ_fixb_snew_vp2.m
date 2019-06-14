% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 01:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:57:59
% EndTime: 2019-05-05 00:58:03
% DurationCPUTime: 3.22s
% Computational Cost: add. (32419->265), mult. (72280->339), div. (0->0), fcn. (55327->14), ass. (0->123)
t483 = sin(pkin(11));
t486 = cos(pkin(11));
t469 = t483 * g(1) - t486 * g(2);
t481 = -g(3) + qJDD(1);
t484 = sin(pkin(6));
t487 = cos(pkin(6));
t524 = t469 * t487 + t481 * t484;
t497 = qJD(2) ^ 2;
t470 = -t486 * g(1) - t483 * g(2);
t491 = sin(qJ(2));
t495 = cos(qJ(2));
t441 = -t491 * t470 + t524 * t495;
t485 = cos(pkin(12));
t480 = t485 ^ 2;
t523 = 0.2e1 * t485;
t522 = pkin(3) * t485;
t482 = sin(pkin(12));
t521 = mrSges(4,2) * t482;
t519 = t480 * t497;
t442 = t495 * t470 + t524 * t491;
t436 = -t497 * pkin(2) + qJDD(2) * qJ(3) + t442;
t458 = -t484 * t469 + t487 * t481;
t514 = qJD(2) * qJD(3);
t516 = t485 * t458 - 0.2e1 * t482 * t514;
t413 = (-pkin(8) * qJDD(2) + t497 * t522 - t436) * t482 + t516;
t416 = t485 * t436 + t482 * t458 + t514 * t523;
t512 = t485 * qJDD(2);
t414 = -pkin(3) * t519 + pkin(8) * t512 + t416;
t490 = sin(qJ(4));
t494 = cos(qJ(4));
t400 = t490 * t413 + t494 * t414;
t462 = (-t482 * t490 + t485 * t494) * qJD(2);
t506 = t482 * t494 + t485 * t490;
t463 = t506 * qJD(2);
t445 = -t462 * mrSges(5,1) + t463 * mrSges(5,2);
t460 = t463 * qJD(4);
t513 = t482 * qJDD(2);
t449 = -t490 * t513 + t494 * t512 - t460;
t457 = qJD(4) * mrSges(5,1) - t463 * mrSges(5,3);
t448 = -t462 * pkin(4) - t463 * pkin(9);
t496 = qJD(4) ^ 2;
t393 = -t496 * pkin(4) + qJDD(4) * pkin(9) + t462 * t448 + t400;
t479 = t482 ^ 2;
t503 = qJDD(3) - t441;
t426 = (-pkin(2) - t522) * qJDD(2) + (-qJ(3) + (-t479 - t480) * pkin(8)) * t497 + t503;
t515 = t462 * qJD(4);
t450 = t506 * qJDD(2) + t515;
t404 = (-t450 - t515) * pkin(9) + (-t449 + t460) * pkin(4) + t426;
t489 = sin(qJ(5));
t493 = cos(qJ(5));
t388 = -t489 * t393 + t493 * t404;
t452 = t493 * qJD(4) - t489 * t463;
t425 = t452 * qJD(5) + t489 * qJDD(4) + t493 * t450;
t447 = qJDD(5) - t449;
t453 = t489 * qJD(4) + t493 * t463;
t461 = qJD(5) - t462;
t386 = (t452 * t461 - t425) * pkin(10) + (t452 * t453 + t447) * pkin(5) + t388;
t389 = t493 * t393 + t489 * t404;
t424 = -t453 * qJD(5) + t493 * qJDD(4) - t489 * t450;
t435 = t461 * pkin(5) - t453 * pkin(10);
t451 = t452 ^ 2;
t387 = -t451 * pkin(5) + t424 * pkin(10) - t461 * t435 + t389;
t488 = sin(qJ(6));
t492 = cos(qJ(6));
t384 = t492 * t386 - t488 * t387;
t427 = t492 * t452 - t488 * t453;
t398 = t427 * qJD(6) + t488 * t424 + t492 * t425;
t428 = t488 * t452 + t492 * t453;
t409 = -t427 * mrSges(7,1) + t428 * mrSges(7,2);
t459 = qJD(6) + t461;
t417 = -t459 * mrSges(7,2) + t427 * mrSges(7,3);
t444 = qJDD(6) + t447;
t381 = m(7) * t384 + t444 * mrSges(7,1) - t398 * mrSges(7,3) - t428 * t409 + t459 * t417;
t385 = t488 * t386 + t492 * t387;
t397 = -t428 * qJD(6) + t492 * t424 - t488 * t425;
t418 = t459 * mrSges(7,1) - t428 * mrSges(7,3);
t382 = m(7) * t385 - t444 * mrSges(7,2) + t397 * mrSges(7,3) + t427 * t409 - t459 * t418;
t373 = t492 * t381 + t488 * t382;
t429 = -t452 * mrSges(6,1) + t453 * mrSges(6,2);
t433 = -t461 * mrSges(6,2) + t452 * mrSges(6,3);
t371 = m(6) * t388 + t447 * mrSges(6,1) - t425 * mrSges(6,3) - t453 * t429 + t461 * t433 + t373;
t434 = t461 * mrSges(6,1) - t453 * mrSges(6,3);
t508 = -t488 * t381 + t492 * t382;
t372 = m(6) * t389 - t447 * mrSges(6,2) + t424 * mrSges(6,3) + t452 * t429 - t461 * t434 + t508;
t509 = -t489 * t371 + t493 * t372;
t365 = m(5) * t400 - qJDD(4) * mrSges(5,2) + t449 * mrSges(5,3) - qJD(4) * t457 + t462 * t445 + t509;
t399 = t494 * t413 - t490 * t414;
t456 = -qJD(4) * mrSges(5,2) + t462 * mrSges(5,3);
t392 = -qJDD(4) * pkin(4) - t496 * pkin(9) + t463 * t448 - t399;
t390 = -t424 * pkin(5) - t451 * pkin(10) + t453 * t435 + t392;
t504 = m(7) * t390 - t397 * mrSges(7,1) + t398 * mrSges(7,2) - t427 * t417 + t428 * t418;
t499 = -m(6) * t392 + t424 * mrSges(6,1) - t425 * mrSges(6,2) + t452 * t433 - t453 * t434 - t504;
t377 = m(5) * t399 + qJDD(4) * mrSges(5,1) - t450 * mrSges(5,3) + qJD(4) * t456 - t463 * t445 + t499;
t517 = t490 * t365 + t494 * t377;
t367 = t493 * t371 + t489 * t372;
t415 = -t482 * t436 + t516;
t505 = mrSges(4,3) * qJDD(2) + t497 * (-t485 * mrSges(4,1) + t521);
t359 = m(4) * t415 - t505 * t482 + t517;
t510 = t494 * t365 - t490 * t377;
t360 = m(4) * t416 + t505 * t485 + t510;
t511 = -t482 * t359 + t485 * t360;
t406 = Ifges(7,4) * t428 + Ifges(7,2) * t427 + Ifges(7,6) * t459;
t407 = Ifges(7,1) * t428 + Ifges(7,4) * t427 + Ifges(7,5) * t459;
t502 = -mrSges(7,1) * t384 + mrSges(7,2) * t385 - Ifges(7,5) * t398 - Ifges(7,6) * t397 - Ifges(7,3) * t444 - t428 * t406 + t427 * t407;
t501 = m(5) * t426 - t449 * mrSges(5,1) + t450 * mrSges(5,2) - t462 * t456 + t463 * t457 + t367;
t432 = -qJDD(2) * pkin(2) - t497 * qJ(3) + t503;
t500 = -m(4) * t432 + mrSges(4,1) * t512 - t501 + (t479 * t497 + t519) * mrSges(4,3);
t420 = Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t461;
t421 = Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t461;
t498 = mrSges(6,1) * t388 - mrSges(6,2) * t389 + Ifges(6,5) * t425 + Ifges(6,6) * t424 + Ifges(6,3) * t447 + pkin(5) * t373 + t453 * t420 - t452 * t421 - t502;
t439 = Ifges(5,1) * t463 + Ifges(5,4) * t462 + Ifges(5,5) * qJD(4);
t438 = Ifges(5,4) * t463 + Ifges(5,2) * t462 + Ifges(5,6) * qJD(4);
t437 = Ifges(5,5) * t463 + Ifges(5,6) * t462 + Ifges(5,3) * qJD(4);
t419 = Ifges(6,5) * t453 + Ifges(6,6) * t452 + Ifges(6,3) * t461;
t405 = Ifges(7,5) * t428 + Ifges(7,6) * t427 + Ifges(7,3) * t459;
t375 = mrSges(7,2) * t390 - mrSges(7,3) * t384 + Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t444 + t427 * t405 - t459 * t406;
t374 = -mrSges(7,1) * t390 + mrSges(7,3) * t385 + Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t444 - t428 * t405 + t459 * t407;
t366 = mrSges(4,2) * t513 - t500;
t362 = mrSges(6,2) * t392 - mrSges(6,3) * t388 + Ifges(6,1) * t425 + Ifges(6,4) * t424 + Ifges(6,5) * t447 - pkin(10) * t373 - t488 * t374 + t492 * t375 + t452 * t419 - t461 * t420;
t361 = -mrSges(6,1) * t392 + mrSges(6,3) * t389 + Ifges(6,4) * t425 + Ifges(6,2) * t424 + Ifges(6,6) * t447 - pkin(5) * t504 + pkin(10) * t508 + t492 * t374 + t488 * t375 - t453 * t419 + t461 * t421;
t357 = -mrSges(5,1) * t426 + mrSges(5,3) * t400 + Ifges(5,4) * t450 + Ifges(5,2) * t449 + Ifges(5,6) * qJDD(4) - pkin(4) * t367 + qJD(4) * t439 - t463 * t437 - t498;
t356 = mrSges(5,2) * t426 - mrSges(5,3) * t399 + Ifges(5,1) * t450 + Ifges(5,4) * t449 + Ifges(5,5) * qJDD(4) - pkin(9) * t367 - qJD(4) * t438 - t489 * t361 + t493 * t362 + t462 * t437;
t1 = [m(2) * t481 + t487 * (m(3) * t458 + t485 * t359 + t482 * t360) + (t491 * (m(3) * t442 - t497 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t511) + t495 * (-t497 * mrSges(3,2) + t500 + (mrSges(3,1) - t521) * qJDD(2) + m(3) * t441)) * t484; mrSges(3,1) * t441 - mrSges(3,2) * t442 + t482 * (mrSges(4,2) * t432 - mrSges(4,3) * t415 - pkin(8) * t517 + t494 * t356 - t490 * t357) + t485 * (-mrSges(4,1) * t432 + mrSges(4,3) * t416 - pkin(3) * t501 + pkin(8) * t510 + t490 * t356 + t494 * t357) - pkin(2) * t366 + qJ(3) * t511 + (Ifges(4,2) * t480 + Ifges(3,3) + (Ifges(4,1) * t482 + Ifges(4,4) * t523) * t482) * qJDD(2); t366; mrSges(5,1) * t399 - mrSges(5,2) * t400 + Ifges(5,5) * t450 + Ifges(5,6) * t449 + Ifges(5,3) * qJDD(4) + pkin(4) * t499 + pkin(9) * t509 + t493 * t361 + t489 * t362 + t463 * t438 - t462 * t439; t498; -t502;];
tauJ  = t1;
