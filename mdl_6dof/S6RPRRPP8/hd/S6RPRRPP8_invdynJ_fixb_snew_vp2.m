% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:50:58
% EndTime: 2019-05-05 21:51:01
% DurationCPUTime: 1.34s
% Computational Cost: add. (6333->244), mult. (11928->277), div. (0->0), fcn. (6695->6), ass. (0->98)
t493 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t476 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t475 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t492 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t474 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t491 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t448 = sin(qJ(1));
t450 = cos(qJ(1));
t463 = -g(1) * t450 - g(2) * t448;
t490 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t463;
t449 = cos(qJ(3));
t479 = qJD(1) * t449;
t452 = qJD(1) ^ 2;
t486 = (-pkin(1) - pkin(7));
t416 = (t486 * t452) - t490;
t447 = sin(qJ(3));
t478 = qJD(1) * qJD(3);
t467 = t449 * t478;
t435 = -qJDD(1) * t447 - t467;
t468 = t447 * t478;
t436 = qJDD(1) * t449 - t468;
t372 = (-t436 + t468) * pkin(8) + (-t435 + t467) * pkin(3) + t416;
t466 = g(1) * t448 - t450 * g(2);
t459 = -qJ(2) * t452 + qJDD(2) - t466;
t418 = t486 * qJDD(1) + t459;
t407 = -g(3) * t449 + t447 * t418;
t434 = (pkin(3) * t447 - pkin(8) * t449) * qJD(1);
t451 = qJD(3) ^ 2;
t480 = qJD(1) * t447;
t376 = -pkin(3) * t451 + qJDD(3) * pkin(8) - t434 * t480 + t407;
t446 = sin(qJ(4));
t485 = cos(qJ(4));
t369 = t485 * t372 - t446 * t376;
t431 = -t485 * qJD(3) + t446 * t479;
t432 = t446 * qJD(3) + t485 * t479;
t403 = pkin(4) * t431 - qJ(5) * t432;
t430 = qJDD(4) - t435;
t440 = qJD(4) + t480;
t439 = t440 ^ 2;
t367 = -t430 * pkin(4) - t439 * qJ(5) + t432 * t403 + qJDD(5) - t369;
t398 = -t431 * qJD(4) + t446 * qJDD(3) + t485 * t436;
t405 = -mrSges(6,2) * t431 - mrSges(6,3) * t432;
t489 = -m(6) * t367 - t398 * mrSges(6,1) - t432 * t405;
t402 = -mrSges(7,2) * t432 + mrSges(7,3) * t431;
t482 = t431 * t440;
t361 = -0.2e1 * qJD(6) * t440 + (t431 * t432 - t430) * qJ(6) + (t398 + t482) * pkin(5) + t367;
t413 = -mrSges(7,1) * t431 + mrSges(7,2) * t440;
t464 = -m(7) * t361 + t430 * mrSges(7,3) + t440 * t413;
t359 = mrSges(7,1) * t398 + t402 * t432 - t464;
t412 = mrSges(6,1) * t431 - mrSges(6,3) * t440;
t357 = mrSges(6,2) * t430 + t412 * t440 + t359 - t489;
t397 = qJD(4) * t432 - t485 * qJDD(3) + t436 * t446;
t410 = pkin(5) * t432 - qJ(6) * t440;
t429 = t431 ^ 2;
t370 = t446 * t372 + t485 * t376;
t456 = -pkin(4) * t439 + qJ(5) * t430 - t403 * t431 + t370;
t365 = -pkin(5) * t397 - qJ(6) * t429 + qJDD(6) + ((2 * qJD(5)) + t410) * t440 + t456;
t487 = -2 * qJD(5);
t366 = t440 * t487 - t456;
t414 = mrSges(6,1) * t432 + mrSges(6,2) * t440;
t411 = mrSges(7,1) * t432 - mrSges(7,3) * t440;
t472 = m(7) * t365 + t430 * mrSges(7,2) + t440 * t411;
t458 = -m(6) * t366 + t430 * mrSges(6,3) + t440 * t414 + t472;
t469 = -t476 * t431 + t493 * t432 + t475 * t440;
t470 = t492 * t431 + t476 * t432 + t474 * t440;
t481 = -t402 - t405;
t488 = -t474 * t397 + t475 * t398 + t491 * t430 + t469 * t431 + t470 * t432 + mrSges(5,1) * t369 - mrSges(5,2) * t370 + mrSges(6,2) * t367 + mrSges(7,2) * t365 - mrSges(6,3) * t366 - mrSges(7,3) * t361 - pkin(4) * t357 + qJ(5) * (t481 * t431 + (-mrSges(6,1) - mrSges(7,1)) * t397 + t458) - qJ(6) * t359;
t483 = -mrSges(7,1) - mrSges(5,3);
t404 = mrSges(5,1) * t431 + mrSges(5,2) * t432;
t408 = -mrSges(5,2) * t440 - mrSges(5,3) * t431;
t352 = m(5) * t369 + (t408 - t412) * t440 + (-t402 - t404) * t432 + (mrSges(5,1) - mrSges(6,2)) * t430 + t483 * t398 + t464 + t489;
t409 = mrSges(5,1) * t440 - mrSges(5,3) * t432;
t355 = m(5) * t370 - mrSges(5,2) * t430 - t409 * t440 + (-t404 + t481) * t431 + (-mrSges(6,1) + t483) * t397 + t458;
t349 = t485 * t352 + t446 * t355;
t471 = t474 * t431 - t475 * t432 - t491 * t440;
t465 = -t352 * t446 + t485 * t355;
t406 = g(3) * t447 + t418 * t449;
t375 = -qJDD(3) * pkin(3) - pkin(8) * t451 + t434 * t479 - t406;
t455 = (-t398 + t482) * qJ(5) + t375 + (t440 * pkin(4) + t487) * t432;
t363 = -pkin(5) * t429 + 0.2e1 * qJD(6) * t431 - t410 * t432 + (pkin(4) + qJ(6)) * t397 + t455;
t462 = m(7) * t363 - t398 * mrSges(7,2) + t397 * mrSges(7,3) - t432 * t411 + t431 * t413;
t433 = (mrSges(4,1) * t447 + mrSges(4,2) * t449) * qJD(1);
t437 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t480;
t438 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t479;
t368 = pkin(4) * t397 + t455;
t457 = -m(6) * t368 + t397 * mrSges(6,2) + t431 * t412 - t462;
t454 = -m(5) * t375 - t397 * mrSges(5,1) - t431 * t408 + (-t409 + t414) * t432 + (-mrSges(5,2) + mrSges(6,3)) * t398 + t457;
t461 = (m(4) * t407 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t435 - qJD(3) * t438 - t433 * t480 + t465) * t447 + (m(4) * t406 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t436 + qJD(3) * t437 - t433 * t479 + t454) * t449;
t425 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t449 - Ifges(4,4) * t447) * qJD(1);
t424 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t449 - Ifges(4,2) * t447) * qJD(1);
t420 = -qJDD(1) * pkin(1) + t459;
t419 = pkin(1) * t452 + t490;
t360 = -mrSges(7,1) * t397 - t402 * t431 + t472;
t356 = -mrSges(6,3) * t398 - t414 * t432 - t457;
t347 = mrSges(6,1) * t367 + mrSges(7,1) * t361 + mrSges(5,2) * t375 - mrSges(7,2) * t363 - mrSges(5,3) * t369 - mrSges(6,3) * t368 + pkin(5) * t359 - qJ(5) * t356 - t476 * t397 + t493 * t398 + t475 * t430 + t471 * t431 - t470 * t440;
t346 = m(3) * t420 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t452) + t461;
t345 = -mrSges(5,1) * t375 - mrSges(6,1) * t366 + mrSges(7,1) * t365 + mrSges(6,2) * t368 + mrSges(5,3) * t370 - mrSges(7,3) * t363 - pkin(4) * t356 + pkin(5) * t360 - qJ(6) * t462 + t492 * t397 + t476 * t398 + t474 * t430 + t471 * t432 + t469 * t440;
t1 = [mrSges(2,1) * t466 - mrSges(2,2) * t463 + mrSges(3,2) * t420 - mrSges(3,3) * t419 + t449 * (mrSges(4,2) * t416 - mrSges(4,3) * t406 + Ifges(4,1) * t436 + Ifges(4,4) * t435 + Ifges(4,5) * qJDD(3) - pkin(8) * t349 - qJD(3) * t424 - t446 * t345 + t485 * t347) - pkin(7) * t461 - pkin(1) * t346 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t419 + m(4) * t416 - mrSges(4,1) * t435 + mrSges(3,2) * t452 + mrSges(4,2) * t436 + qJDD(1) * mrSges(3,3) + t438 * t479 + t349) * qJ(2) + (qJ(2) * qJD(1) * t437 + mrSges(4,1) * t416 - mrSges(4,3) * t407 - Ifges(4,4) * t436 - Ifges(4,2) * t435 - Ifges(4,6) * qJDD(3) + pkin(3) * t349 - qJD(3) * t425 + t488) * t447; t346; Ifges(4,5) * t436 + Ifges(4,6) * t435 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t406 - mrSges(4,2) * t407 + t446 * t347 + t485 * t345 + pkin(3) * t454 + pkin(8) * t465 + (t424 * t449 + t425 * t447) * qJD(1); t488; t357; t360;];
tauJ  = t1;
