% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:08:43
% EndTime: 2019-05-05 23:08:49
% DurationCPUTime: 4.08s
% Computational Cost: add. (42867->289), mult. (87068->363), div. (0->0), fcn. (58151->10), ass. (0->115)
t476 = qJD(1) ^ 2;
t493 = -pkin(1) - pkin(7);
t470 = sin(qJ(1));
t474 = cos(qJ(1));
t484 = -t474 * g(1) - t470 * g(2);
t494 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t484;
t437 = t493 * t476 - t494;
t473 = cos(qJ(3));
t491 = qJD(1) * qJD(3);
t459 = t473 * t491;
t469 = sin(qJ(3));
t453 = -t469 * qJDD(1) - t459;
t489 = t469 * t491;
t454 = t473 * qJDD(1) - t489;
t410 = (-t454 + t489) * pkin(8) + (-t453 + t459) * pkin(3) + t437;
t488 = g(1) * t470 - t474 * g(2);
t480 = -t476 * qJ(2) + qJDD(2) - t488;
t438 = t493 * qJDD(1) + t480;
t432 = -t473 * g(3) + t469 * t438;
t452 = (t469 * pkin(3) - t473 * pkin(8)) * qJD(1);
t461 = t469 * qJD(1);
t475 = qJD(3) ^ 2;
t413 = -t475 * pkin(3) + qJDD(3) * pkin(8) - t452 * t461 + t432;
t468 = sin(qJ(4));
t472 = cos(qJ(4));
t394 = t472 * t410 - t413 * t468;
t492 = t473 * qJD(1);
t449 = t472 * qJD(3) - t468 * t492;
t426 = t449 * qJD(4) + t468 * qJDD(3) + t472 * t454;
t448 = qJDD(4) - t453;
t450 = t468 * qJD(3) + t472 * t492;
t458 = t461 + qJD(4);
t384 = (t449 * t458 - t426) * qJ(5) + (t449 * t450 + t448) * pkin(4) + t394;
t395 = t468 * t410 + t472 * t413;
t425 = -t450 * qJD(4) + t472 * qJDD(3) - t468 * t454;
t434 = pkin(4) * t458 - qJ(5) * t450;
t447 = t449 ^ 2;
t386 = -pkin(4) * t447 + qJ(5) * t425 - t434 * t458 + t395;
t465 = sin(pkin(10));
t466 = cos(pkin(10));
t429 = t449 * t465 + t450 * t466;
t371 = -0.2e1 * qJD(5) * t429 + t466 * t384 - t386 * t465;
t403 = t425 * t465 + t426 * t466;
t428 = t449 * t466 - t450 * t465;
t369 = (t428 * t458 - t403) * pkin(9) + (t428 * t429 + t448) * pkin(5) + t371;
t372 = 0.2e1 * qJD(5) * t428 + t465 * t384 + t466 * t386;
t402 = t425 * t466 - t426 * t465;
t416 = pkin(5) * t458 - pkin(9) * t429;
t427 = t428 ^ 2;
t370 = -pkin(5) * t427 + pkin(9) * t402 - t416 * t458 + t372;
t467 = sin(qJ(6));
t471 = cos(qJ(6));
t367 = t369 * t471 - t467 * t370;
t405 = t428 * t471 - t467 * t429;
t380 = t405 * qJD(6) + t467 * t402 + t403 * t471;
t406 = t467 * t428 + t471 * t429;
t392 = -mrSges(7,1) * t405 + mrSges(7,2) * t406;
t457 = qJD(6) + t458;
t396 = -mrSges(7,2) * t457 + mrSges(7,3) * t405;
t446 = qJDD(6) + t448;
t362 = m(7) * t367 + mrSges(7,1) * t446 - t380 * mrSges(7,3) - t392 * t406 + t396 * t457;
t368 = t467 * t369 + t370 * t471;
t379 = -t406 * qJD(6) + t402 * t471 - t467 * t403;
t397 = mrSges(7,1) * t457 - mrSges(7,3) * t406;
t363 = m(7) * t368 - mrSges(7,2) * t446 + t379 * mrSges(7,3) + t392 * t405 - t397 * t457;
t356 = t471 * t362 + t467 * t363;
t407 = -mrSges(6,1) * t428 + mrSges(6,2) * t429;
t414 = -mrSges(6,2) * t458 + mrSges(6,3) * t428;
t354 = m(6) * t371 + mrSges(6,1) * t448 - mrSges(6,3) * t403 - t407 * t429 + t414 * t458 + t356;
t415 = mrSges(6,1) * t458 - mrSges(6,3) * t429;
t485 = -t362 * t467 + t471 * t363;
t355 = m(6) * t372 - mrSges(6,2) * t448 + mrSges(6,3) * t402 + t407 * t428 - t415 * t458 + t485;
t350 = t466 * t354 + t465 * t355;
t400 = Ifges(6,4) * t429 + Ifges(6,2) * t428 + Ifges(6,6) * t458;
t401 = Ifges(6,1) * t429 + Ifges(6,4) * t428 + Ifges(6,5) * t458;
t418 = Ifges(5,4) * t450 + Ifges(5,2) * t449 + Ifges(5,6) * t458;
t419 = Ifges(5,1) * t450 + Ifges(5,4) * t449 + Ifges(5,5) * t458;
t388 = Ifges(7,4) * t406 + Ifges(7,2) * t405 + Ifges(7,6) * t457;
t389 = Ifges(7,1) * t406 + Ifges(7,4) * t405 + Ifges(7,5) * t457;
t479 = -mrSges(7,1) * t367 + mrSges(7,2) * t368 - Ifges(7,5) * t380 - Ifges(7,6) * t379 - Ifges(7,3) * t446 - t406 * t388 + t405 * t389;
t496 = mrSges(5,1) * t394 + mrSges(6,1) * t371 - mrSges(5,2) * t395 - mrSges(6,2) * t372 + Ifges(5,5) * t426 + Ifges(6,5) * t403 + Ifges(5,6) * t425 + Ifges(6,6) * t402 + pkin(4) * t350 + pkin(5) * t356 + t429 * t400 - t428 * t401 + t450 * t418 - t449 * t419 + (Ifges(5,3) + Ifges(6,3)) * t448 - t479;
t430 = -mrSges(5,1) * t449 + mrSges(5,2) * t450;
t433 = -mrSges(5,2) * t458 + mrSges(5,3) * t449;
t348 = m(5) * t394 + mrSges(5,1) * t448 - mrSges(5,3) * t426 - t430 * t450 + t433 * t458 + t350;
t435 = mrSges(5,1) * t458 - mrSges(5,3) * t450;
t486 = -t354 * t465 + t466 * t355;
t349 = m(5) * t395 - mrSges(5,2) * t448 + mrSges(5,3) * t425 + t430 * t449 - t435 * t458 + t486;
t342 = t472 * t348 + t468 * t349;
t487 = -t468 * t348 + t472 * t349;
t431 = t469 * g(3) + t473 * t438;
t412 = -qJDD(3) * pkin(3) - t475 * pkin(8) + t452 * t492 - t431;
t393 = -t425 * pkin(4) - t447 * qJ(5) + t450 * t434 + qJDD(5) + t412;
t374 = -t402 * pkin(5) - t427 * pkin(9) + t429 * t416 + t393;
t483 = m(7) * t374 - t379 * mrSges(7,1) + t380 * mrSges(7,2) - t405 * t396 + t406 * t397;
t451 = (t469 * mrSges(4,1) + t473 * mrSges(4,2)) * qJD(1);
t455 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t461;
t456 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t492;
t365 = m(6) * t393 - t402 * mrSges(6,1) + t403 * mrSges(6,2) - t428 * t414 + t429 * t415 + t483;
t478 = -m(5) * t412 + t425 * mrSges(5,1) - t426 * mrSges(5,2) + t449 * t433 - t450 * t435 - t365;
t482 = t469 * (m(4) * t432 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t453 - qJD(3) * t456 - t451 * t461 + t487) + t473 * (m(4) * t431 + qJDD(3) * mrSges(4,1) - t454 * mrSges(4,3) + qJD(3) * t455 - t451 * t492 + t478);
t445 = (Ifges(4,5) * qJD(3)) + (t473 * Ifges(4,1) - t469 * Ifges(4,4)) * qJD(1);
t444 = (Ifges(4,6) * qJD(3)) + (t473 * Ifges(4,4) - t469 * Ifges(4,2)) * qJD(1);
t440 = -qJDD(1) * pkin(1) + t480;
t439 = t476 * pkin(1) + t494;
t417 = Ifges(5,5) * t450 + Ifges(5,6) * t449 + Ifges(5,3) * t458;
t399 = Ifges(6,5) * t429 + Ifges(6,6) * t428 + Ifges(6,3) * t458;
t387 = Ifges(7,5) * t406 + Ifges(7,6) * t405 + Ifges(7,3) * t457;
t358 = mrSges(7,2) * t374 - mrSges(7,3) * t367 + Ifges(7,1) * t380 + Ifges(7,4) * t379 + Ifges(7,5) * t446 + t387 * t405 - t388 * t457;
t357 = -mrSges(7,1) * t374 + mrSges(7,3) * t368 + Ifges(7,4) * t380 + Ifges(7,2) * t379 + Ifges(7,6) * t446 - t387 * t406 + t389 * t457;
t344 = mrSges(6,2) * t393 - mrSges(6,3) * t371 + Ifges(6,1) * t403 + Ifges(6,4) * t402 + Ifges(6,5) * t448 - pkin(9) * t356 - t467 * t357 + t471 * t358 + t428 * t399 - t458 * t400;
t343 = -mrSges(6,1) * t393 + mrSges(6,3) * t372 + Ifges(6,4) * t403 + Ifges(6,2) * t402 + Ifges(6,6) * t448 - pkin(5) * t483 + pkin(9) * t485 + t471 * t357 + t467 * t358 - t429 * t399 + t458 * t401;
t340 = m(3) * t440 + qJDD(1) * mrSges(3,2) - t476 * mrSges(3,3) + t482;
t339 = mrSges(5,2) * t412 - mrSges(5,3) * t394 + Ifges(5,1) * t426 + Ifges(5,4) * t425 + Ifges(5,5) * t448 - qJ(5) * t350 - t343 * t465 + t344 * t466 + t417 * t449 - t418 * t458;
t338 = -mrSges(5,1) * t412 + mrSges(5,3) * t395 + Ifges(5,4) * t426 + Ifges(5,2) * t425 + Ifges(5,6) * t448 - pkin(4) * t365 + qJ(5) * t486 + t466 * t343 + t465 * t344 - t450 * t417 + t458 * t419;
t1 = [mrSges(2,1) * t488 - mrSges(2,2) * t484 + mrSges(3,2) * t440 - mrSges(3,3) * t439 + t473 * (mrSges(4,2) * t437 - mrSges(4,3) * t431 + Ifges(4,1) * t454 + Ifges(4,4) * t453 + Ifges(4,5) * qJDD(3) - pkin(8) * t342 - qJD(3) * t444 - t468 * t338 + t472 * t339) - t469 * (-mrSges(4,1) * t437 + mrSges(4,3) * t432 + Ifges(4,4) * t454 + Ifges(4,2) * t453 + Ifges(4,6) * qJDD(3) - pkin(3) * t342 + qJD(3) * t445 - t496) - pkin(7) * t482 - pkin(1) * t340 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t439 + m(4) * t437 - t453 * mrSges(4,1) + t476 * mrSges(3,2) + t454 * mrSges(4,2) + t342 + qJDD(1) * mrSges(3,3) + (t455 * t469 + t456 * t473) * qJD(1)) * qJ(2); t340; Ifges(4,5) * t454 + Ifges(4,6) * t453 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t431 - mrSges(4,2) * t432 + t468 * t339 + t472 * t338 + pkin(3) * t478 + pkin(8) * t487 + (t473 * t444 + t469 * t445) * qJD(1); t496; t365; -t479;];
tauJ  = t1;
