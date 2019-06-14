% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:05:07
% EndTime: 2019-05-05 18:05:10
% DurationCPUTime: 2.16s
% Computational Cost: add. (15135->264), mult. (31441->321), div. (0->0), fcn. (19847->8), ass. (0->105)
t491 = Ifges(6,1) + Ifges(7,1);
t479 = Ifges(6,4) - Ifges(7,5);
t487 = -Ifges(6,5) - Ifges(7,4);
t490 = Ifges(6,2) + Ifges(7,3);
t477 = Ifges(6,6) - Ifges(7,6);
t447 = sin(pkin(9));
t448 = cos(pkin(9));
t452 = cos(qJ(3));
t471 = qJD(1) * t452;
t433 = qJD(3) * t448 - t447 * t471;
t434 = qJD(3) * t447 + t448 * t471;
t449 = sin(qJ(5));
t481 = cos(qJ(5));
t406 = -t481 * t433 + t449 * t434;
t450 = sin(qJ(3));
t470 = qJD(1) * qJD(3);
t467 = t450 * t470;
t439 = qJDD(1) * t452 - t467;
t417 = qJDD(3) * t448 - t439 * t447;
t418 = qJDD(3) * t447 + t439 * t448;
t376 = -t406 * qJD(5) + t449 * t417 + t481 * t418;
t407 = t449 * t433 + t481 * t434;
t387 = mrSges(7,1) * t406 - mrSges(7,3) * t407;
t455 = qJD(1) ^ 2;
t482 = -pkin(1) - pkin(7);
t451 = sin(qJ(1));
t453 = cos(qJ(1));
t461 = -t453 * g(1) - t451 * g(2);
t483 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t461;
t414 = t482 * t455 - t483;
t466 = t452 * t470;
t438 = qJDD(1) * t450 + t466;
t391 = (-t439 + t467) * qJ(4) + (t438 + t466) * pkin(3) + t414;
t465 = t451 * g(1) - t453 * g(2);
t458 = -t455 * qJ(2) + qJDD(2) - t465;
t420 = t482 * qJDD(1) + t458;
t410 = -g(3) * t452 + t450 * t420;
t436 = (pkin(3) * t450 - qJ(4) * t452) * qJD(1);
t454 = qJD(3) ^ 2;
t472 = qJD(1) * t450;
t396 = -pkin(3) * t454 + qJDD(3) * qJ(4) - t436 * t472 + t410;
t369 = -0.2e1 * qJD(4) * t434 + t448 * t391 - t447 * t396;
t366 = (t433 * t472 - t418) * pkin(8) + (t433 * t434 + t438) * pkin(4) + t369;
t370 = 0.2e1 * qJD(4) * t433 + t447 * t391 + t448 * t396;
t419 = pkin(4) * t472 - pkin(8) * t434;
t432 = t433 ^ 2;
t368 = -pkin(4) * t432 + pkin(8) * t417 - t419 * t472 + t370;
t361 = t481 * t366 - t449 * t368;
t386 = pkin(5) * t406 - qJ(6) * t407;
t435 = qJDD(5) + t438;
t443 = qJD(5) + t472;
t442 = t443 ^ 2;
t360 = -t435 * pkin(5) - t442 * qJ(6) + t407 * t386 + qJDD(6) - t361;
t398 = -mrSges(7,2) * t406 + mrSges(7,3) * t443;
t462 = -m(7) * t360 + t435 * mrSges(7,1) + t443 * t398;
t356 = t376 * mrSges(7,2) + t407 * t387 - t462;
t362 = t449 * t366 + t481 * t368;
t359 = -pkin(5) * t442 + qJ(6) * t435 + 0.2e1 * qJD(6) * t443 - t386 * t406 + t362;
t375 = t407 * qJD(5) - t481 * t417 + t449 * t418;
t400 = -mrSges(7,1) * t443 + mrSges(7,2) * t407;
t468 = m(7) * t359 + t435 * mrSges(7,3) + t443 * t400;
t474 = -t479 * t406 + t491 * t407 - t487 * t443;
t475 = t490 * t406 - t479 * t407 - t477 * t443;
t486 = -Ifges(6,3) - Ifges(7,2);
t489 = -t487 * t376 - t477 * t375 - t486 * t435 + mrSges(6,1) * t361 - mrSges(7,1) * t360 - mrSges(6,2) * t362 + mrSges(7,3) * t359 - pkin(5) * t356 + qJ(6) * (-t375 * mrSges(7,2) - t406 * t387 + t468) - t475 * t407 + t474 * t406;
t480 = -mrSges(6,3) - mrSges(7,2);
t399 = mrSges(6,1) * t443 - mrSges(6,3) * t407;
t473 = -mrSges(6,1) * t406 - mrSges(6,2) * t407 - t387;
t350 = m(6) * t362 - t435 * mrSges(6,2) + t480 * t375 - t443 * t399 + t473 * t406 + t468;
t397 = -mrSges(6,2) * t443 - mrSges(6,3) * t406;
t352 = m(6) * t361 + t435 * mrSges(6,1) + t480 * t376 + t443 * t397 + t473 * t407 + t462;
t347 = t449 * t350 + t481 * t352;
t408 = -mrSges(5,1) * t433 + mrSges(5,2) * t434;
t415 = -mrSges(5,2) * t472 + mrSges(5,3) * t433;
t343 = m(5) * t369 + mrSges(5,1) * t438 - mrSges(5,3) * t418 - t408 * t434 + t415 * t472 + t347;
t416 = mrSges(5,1) * t472 - mrSges(5,3) * t434;
t463 = t481 * t350 - t352 * t449;
t344 = m(5) * t370 - mrSges(5,2) * t438 + mrSges(5,3) * t417 + t408 * t433 - t416 * t472 + t463;
t339 = t448 * t343 + t447 * t344;
t476 = t477 * t406 + t487 * t407 + t486 * t443;
t464 = -t343 * t447 + t448 * t344;
t409 = t450 * g(3) + t452 * t420;
t393 = -qJDD(3) * pkin(3) - t454 * qJ(4) + t436 * t471 + qJDD(4) - t409;
t371 = -t417 * pkin(4) - t432 * pkin(8) + t434 * t419 + t393;
t364 = -0.2e1 * qJD(6) * t407 + (t406 * t443 - t376) * qJ(6) + (t407 * t443 + t375) * pkin(5) + t371;
t357 = m(7) * t364 + t375 * mrSges(7,1) - t376 * mrSges(7,3) + t406 * t398 - t407 * t400;
t457 = m(6) * t371 + t375 * mrSges(6,1) + t376 * mrSges(6,2) + t406 * t397 + t407 * t399 + t357;
t354 = m(5) * t393 - t417 * mrSges(5,1) + t418 * mrSges(5,2) - t433 * t415 + t434 * t416 + t457;
t437 = (mrSges(4,1) * t450 + mrSges(4,2) * t452) * qJD(1);
t440 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t472;
t441 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t471;
t460 = (m(4) * t410 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t438 - qJD(3) * t441 - t437 * t472 + t464) * t450 + (m(4) * t409 + qJDD(3) * mrSges(4,1) - t439 * mrSges(4,3) + qJD(3) * t440 - t437 * t471 - t354) * t452;
t429 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t452 - Ifges(4,4) * t450) * qJD(1);
t428 = Ifges(4,6) * qJD(3) + (t452 * Ifges(4,4) - Ifges(4,2) * t450) * qJD(1);
t422 = -qJDD(1) * pkin(1) + t458;
t421 = t455 * pkin(1) + t483;
t403 = Ifges(5,1) * t434 + Ifges(5,4) * t433 + Ifges(5,5) * t472;
t402 = Ifges(5,4) * t434 + Ifges(5,2) * t433 + Ifges(5,6) * t472;
t401 = Ifges(5,5) * t434 + Ifges(5,6) * t433 + Ifges(5,3) * t472;
t346 = mrSges(6,2) * t371 + mrSges(7,2) * t360 - mrSges(6,3) * t361 - mrSges(7,3) * t364 - qJ(6) * t357 - t479 * t375 + t491 * t376 + t476 * t406 - t487 * t435 + t475 * t443;
t345 = -mrSges(6,1) * t371 - mrSges(7,1) * t364 + mrSges(7,2) * t359 + mrSges(6,3) * t362 - pkin(5) * t357 - t490 * t375 + t479 * t376 + t476 * t407 + t477 * t435 + t474 * t443;
t337 = m(3) * t422 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t455 + t460;
t336 = mrSges(5,2) * t393 - mrSges(5,3) * t369 + Ifges(5,1) * t418 + Ifges(5,4) * t417 + Ifges(5,5) * t438 - pkin(8) * t347 - t449 * t345 + t481 * t346 + t433 * t401 - t402 * t472;
t335 = -mrSges(5,1) * t393 + mrSges(5,3) * t370 + Ifges(5,4) * t418 + Ifges(5,2) * t417 + Ifges(5,6) * t438 - pkin(4) * t457 + pkin(8) * t463 + t481 * t345 + t449 * t346 - t434 * t401 + t403 * t472;
t1 = [mrSges(2,1) * t465 - mrSges(2,2) * t461 + mrSges(3,2) * t422 - mrSges(3,3) * t421 + t452 * (mrSges(4,2) * t414 - mrSges(4,3) * t409 + Ifges(4,1) * t439 - Ifges(4,4) * t438 + Ifges(4,5) * qJDD(3) - qJ(4) * t339 - qJD(3) * t428 - t335 * t447 + t336 * t448) - t450 * ((-Ifges(5,3) - Ifges(4,2)) * t438 + Ifges(4,6) * qJDD(3) + Ifges(4,4) * t439 + t433 * t403 - t434 * t402 - Ifges(5,6) * t417 - Ifges(5,5) * t418 + qJD(3) * t429 + mrSges(4,3) * t410 - mrSges(4,1) * t414 - mrSges(5,1) * t369 + mrSges(5,2) * t370 - pkin(4) * t347 - pkin(3) * t339 - t489) - pkin(7) * t460 - pkin(1) * t337 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t421 + m(4) * t414 + t438 * mrSges(4,1) + t455 * mrSges(3,2) + t439 * mrSges(4,2) + t339 + qJDD(1) * mrSges(3,3) + (t440 * t450 + t441 * t452) * qJD(1)) * qJ(2); t337; Ifges(4,5) * t439 - Ifges(4,6) * t438 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t409 - mrSges(4,2) * t410 + t447 * t336 + t448 * t335 - pkin(3) * t354 + qJ(4) * t464 + (t428 * t452 + t429 * t450) * qJD(1); t354; t489; t356;];
tauJ  = t1;
