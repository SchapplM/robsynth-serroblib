% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:25:37
% EndTime: 2019-05-05 21:25:40
% DurationCPUTime: 1.51s
% Computational Cost: add. (7252->247), mult. (13639->283), div. (0->0), fcn. (7952->8), ass. (0->103)
t491 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t476 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t475 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t490 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t474 = -Ifges(5,6) + Ifges(6,5) + Ifges(7,4);
t489 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t449 = sin(qJ(1));
t451 = cos(qJ(1));
t466 = g(1) * t449 - g(2) * t451;
t429 = qJDD(1) * pkin(1) + t466;
t453 = qJD(1) ^ 2;
t462 = -g(1) * t451 - g(2) * t449;
t431 = -pkin(1) * t453 + t462;
t445 = sin(pkin(9));
t446 = cos(pkin(9));
t395 = t429 * t446 - t431 * t445;
t373 = -qJDD(1) * pkin(2) - pkin(7) * t453 - t395;
t448 = sin(qJ(3));
t450 = cos(qJ(3));
t477 = qJD(1) * qJD(3);
t467 = t450 * t477;
t433 = qJDD(1) * t448 + t467;
t468 = t448 * t477;
t434 = qJDD(1) * t450 - t468;
t362 = (-t433 - t467) * pkin(8) + (-t434 + t468) * pkin(3) + t373;
t396 = t429 * t445 + t431 * t446;
t374 = -pkin(2) * t453 + qJDD(1) * pkin(7) + t396;
t444 = -g(3) + qJDD(2);
t368 = t374 * t450 + t444 * t448;
t432 = (-pkin(3) * t450 - pkin(8) * t448) * qJD(1);
t452 = qJD(3) ^ 2;
t479 = qJD(1) * t450;
t366 = -pkin(3) * t452 + qJDD(3) * pkin(8) + t432 * t479 + t368;
t447 = sin(qJ(4));
t484 = cos(qJ(4));
t359 = t362 * t484 - t366 * t447;
t478 = t448 * qJD(1);
t427 = -qJD(3) * t484 + t447 * t478;
t428 = qJD(3) * t447 + t478 * t484;
t401 = pkin(4) * t427 - qJ(5) * t428;
t426 = qJDD(4) - t434;
t438 = -qJD(4) + t479;
t437 = t438 ^ 2;
t357 = -t426 * pkin(4) - t437 * qJ(5) + t401 * t428 + qJDD(5) - t359;
t394 = -qJD(4) * t427 + qJDD(3) * t447 + t433 * t484;
t403 = -mrSges(6,2) * t427 - mrSges(6,3) * t428;
t488 = -m(6) * t357 - mrSges(6,1) * t394 - t403 * t428;
t400 = -mrSges(7,2) * t428 + mrSges(7,3) * t427;
t481 = t427 * t438;
t485 = 2 * qJD(6);
t351 = t438 * t485 + (t427 * t428 - t426) * qJ(6) + (t394 - t481) * pkin(5) + t357;
t409 = -mrSges(7,1) * t427 - mrSges(7,2) * t438;
t463 = -m(7) * t351 + mrSges(7,3) * t426 - t409 * t438;
t349 = mrSges(7,1) * t394 + t400 * t428 - t463;
t408 = mrSges(6,1) * t427 + mrSges(6,3) * t438;
t346 = mrSges(6,2) * t426 - t408 * t438 + t349 - t488;
t393 = qJD(4) * t428 - qJDD(3) * t484 + t433 * t447;
t406 = pkin(5) * t428 + qJ(6) * t438;
t425 = t427 ^ 2;
t360 = t362 * t447 + t366 * t484;
t458 = -pkin(4) * t437 + qJ(5) * t426 - t401 * t427 + t360;
t486 = -2 * qJD(5);
t353 = -pkin(5) * t393 - qJ(6) * t425 + qJDD(6) + (t486 - t406) * t438 + t458;
t356 = 0.2e1 * qJD(5) * t438 - t458;
t410 = mrSges(6,1) * t428 - mrSges(6,2) * t438;
t407 = mrSges(7,1) * t428 + mrSges(7,3) * t438;
t472 = m(7) * t353 + mrSges(7,2) * t426 - t407 * t438;
t460 = -m(6) * t356 + mrSges(6,3) * t426 - t410 * t438 + t472;
t469 = -t427 * t476 + t428 * t491 - t438 * t475;
t470 = t427 * t490 + t428 * t476 + t438 * t474;
t480 = -t400 - t403;
t487 = t393 * t474 + t394 * t475 + t489 * t426 + t427 * t469 + t428 * t470 + mrSges(5,1) * t359 - mrSges(5,2) * t360 + mrSges(6,2) * t357 + mrSges(7,2) * t353 - mrSges(6,3) * t356 - mrSges(7,3) * t351 - pkin(4) * t346 + qJ(5) * (t480 * t427 + (-mrSges(6,1) - mrSges(7,1)) * t393 + t460) - qJ(6) * t349;
t482 = -mrSges(7,1) - mrSges(5,3);
t471 = -t427 * t474 - t428 * t475 + t438 * t489;
t430 = (-mrSges(4,1) * t450 + mrSges(4,2) * t448) * qJD(1);
t435 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t478;
t402 = mrSges(5,1) * t427 + mrSges(5,2) * t428;
t404 = mrSges(5,2) * t438 - mrSges(5,3) * t427;
t343 = m(5) * t359 + (-t404 + t408) * t438 + (-t400 - t402) * t428 + (mrSges(5,1) - mrSges(6,2)) * t426 + t482 * t394 + t463 + t488;
t405 = -mrSges(5,1) * t438 - mrSges(5,3) * t428;
t345 = m(5) * t360 - mrSges(5,2) * t426 + t405 * t438 + (-t402 + t480) * t427 + (-mrSges(6,1) + t482) * t393 + t460;
t464 = -t343 * t447 + t345 * t484;
t340 = m(4) * t368 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t434 - qJD(3) * t435 + t430 * t479 + t464;
t367 = -t374 * t448 + t444 * t450;
t436 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t479;
t365 = -qJDD(3) * pkin(3) - pkin(8) * t452 + t432 * t478 - t367;
t456 = (-t394 - t481) * qJ(5) + t365 + (-pkin(4) * t438 + t486) * t428;
t358 = pkin(4) * t393 + t456;
t355 = -pkin(5) * t425 + t427 * t485 - t406 * t428 + (pkin(4) + qJ(6)) * t393 + t456;
t461 = m(7) * t355 - t394 * mrSges(7,2) + t393 * mrSges(7,3) - t428 * t407 + t427 * t409;
t459 = -m(6) * t358 + mrSges(6,2) * t393 + t408 * t427 - t461;
t455 = -m(5) * t365 - t393 * mrSges(5,1) - t427 * t404 + (-t405 + t410) * t428 + (-mrSges(5,2) + mrSges(6,3)) * t394 + t459;
t342 = m(4) * t367 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t433 + qJD(3) * t436 - t430 * t478 + t455;
t465 = t340 * t450 - t342 * t448;
t341 = t343 * t484 + t345 * t447;
t457 = -m(4) * t373 + mrSges(4,1) * t434 - mrSges(4,2) * t433 - t435 * t478 + t436 * t479 - t341;
t419 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t448 + Ifges(4,4) * t450) * qJD(1);
t418 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t448 + Ifges(4,2) * t450) * qJD(1);
t350 = -mrSges(7,1) * t393 - t400 * t427 + t472;
t347 = -mrSges(6,3) * t394 - t410 * t428 - t459;
t338 = mrSges(6,1) * t357 + mrSges(7,1) * t351 + mrSges(5,2) * t365 - mrSges(7,2) * t355 - mrSges(5,3) * t359 - mrSges(6,3) * t358 + pkin(5) * t349 - qJ(5) * t347 - t393 * t476 + t394 * t491 + t426 * t475 + t427 * t471 + t438 * t470;
t337 = -mrSges(5,1) * t365 - mrSges(6,1) * t356 + mrSges(7,1) * t353 + mrSges(6,2) * t358 + mrSges(5,3) * t360 - mrSges(7,3) * t355 - pkin(4) * t347 + pkin(5) * t350 - qJ(6) * t461 + t393 * t490 + t394 * t476 - t426 * t474 + t428 * t471 - t438 * t469;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t466 - mrSges(2,2) * t462 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t395 - mrSges(3,2) * t396 + t448 * (mrSges(4,2) * t373 - mrSges(4,3) * t367 + Ifges(4,1) * t433 + Ifges(4,4) * t434 + Ifges(4,5) * qJDD(3) - pkin(8) * t341 - qJD(3) * t418 - t447 * t337 + t338 * t484) + pkin(2) * t457 + pkin(7) * t465 + pkin(1) * (t445 * (m(3) * t396 - mrSges(3,1) * t453 - qJDD(1) * mrSges(3,2) + t465) + t446 * (m(3) * t395 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t453 + t457)) + (-mrSges(4,1) * t373 + mrSges(4,3) * t368 + Ifges(4,4) * t433 + Ifges(4,2) * t434 + Ifges(4,6) * qJDD(3) - pkin(3) * t341 + qJD(3) * t419 - t487) * t450; m(3) * t444 + t340 * t448 + t342 * t450; Ifges(4,5) * t433 + Ifges(4,6) * t434 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t367 - mrSges(4,2) * t368 + t447 * t338 + t484 * t337 + pkin(3) * t455 + pkin(8) * t464 + (t418 * t448 - t419 * t450) * qJD(1); t487; t346; t350;];
tauJ  = t1;
