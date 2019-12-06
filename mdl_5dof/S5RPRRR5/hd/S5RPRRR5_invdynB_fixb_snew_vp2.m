% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:21
% EndTime: 2019-12-05 18:16:24
% DurationCPUTime: 2.78s
% Computational Cost: add. (43952->228), mult. (59398->290), div. (0->0), fcn. (33553->10), ass. (0->96)
t441 = qJD(1) + qJD(3);
t447 = sin(qJ(4));
t469 = t441 * t447;
t451 = cos(qJ(4));
t468 = t441 * t451;
t449 = sin(qJ(1));
t453 = cos(qJ(1));
t429 = t453 * g(2) + t449 * g(3);
t423 = qJDD(1) * pkin(1) + t429;
t428 = t449 * g(2) - t453 * g(3);
t454 = qJD(1) ^ 2;
t424 = -t454 * pkin(1) + t428;
t444 = sin(pkin(9));
t445 = cos(pkin(9));
t406 = t445 * t423 - t444 * t424;
t404 = qJDD(1) * pkin(2) + t406;
t407 = t444 * t423 + t445 * t424;
t405 = -t454 * pkin(2) + t407;
t448 = sin(qJ(3));
t452 = cos(qJ(3));
t392 = t448 * t404 + t452 * t405;
t437 = t441 ^ 2;
t439 = qJDD(1) + qJDD(3);
t390 = -t437 * pkin(3) + t439 * pkin(7) + t392;
t443 = -g(1) + qJDD(2);
t386 = -t447 * t390 + t451 * t443;
t467 = qJD(4) * t441;
t465 = t451 * t467;
t418 = t447 * t439 + t465;
t383 = (-t418 + t465) * pkin(8) + (t437 * t447 * t451 + qJDD(4)) * pkin(4) + t386;
t387 = t451 * t390 + t447 * t443;
t419 = t451 * t439 - t447 * t467;
t427 = qJD(4) * pkin(4) - pkin(8) * t469;
t442 = t451 ^ 2;
t384 = -t442 * t437 * pkin(4) + t419 * pkin(8) - qJD(4) * t427 + t387;
t446 = sin(qJ(5));
t450 = cos(qJ(5));
t381 = t450 * t383 - t446 * t384;
t413 = (-t446 * t447 + t450 * t451) * t441;
t395 = t413 * qJD(5) + t450 * t418 + t446 * t419;
t414 = (t446 * t451 + t447 * t450) * t441;
t400 = -t413 * mrSges(6,1) + t414 * mrSges(6,2);
t440 = qJD(4) + qJD(5);
t408 = -t440 * mrSges(6,2) + t413 * mrSges(6,3);
t438 = qJDD(4) + qJDD(5);
t379 = m(6) * t381 + t438 * mrSges(6,1) - t395 * mrSges(6,3) - t414 * t400 + t440 * t408;
t382 = t446 * t383 + t450 * t384;
t394 = -t414 * qJD(5) - t446 * t418 + t450 * t419;
t409 = t440 * mrSges(6,1) - t414 * mrSges(6,3);
t380 = m(6) * t382 - t438 * mrSges(6,2) + t394 * mrSges(6,3) + t413 * t400 - t440 * t409;
t371 = t450 * t379 + t446 * t380;
t417 = (-mrSges(5,1) * t451 + mrSges(5,2) * t447) * t441;
t426 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t468;
t369 = m(5) * t386 + qJDD(4) * mrSges(5,1) - t418 * mrSges(5,3) + qJD(4) * t426 - t417 * t469 + t371;
t425 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t469;
t460 = -t446 * t379 + t450 * t380;
t370 = m(5) * t387 - qJDD(4) * mrSges(5,2) + t419 * mrSges(5,3) - qJD(4) * t425 + t417 * t468 + t460;
t461 = -t447 * t369 + t451 * t370;
t364 = m(4) * t392 - t437 * mrSges(4,1) - t439 * mrSges(4,2) + t461;
t391 = t452 * t404 - t448 * t405;
t457 = -t439 * pkin(3) - t391;
t389 = -t437 * pkin(7) + t457;
t385 = t427 * t469 - t419 * pkin(4) + (-pkin(8) * t442 - pkin(7)) * t437 + t457;
t456 = m(6) * t385 - t394 * mrSges(6,1) + t395 * mrSges(6,2) - t413 * t408 + t414 * t409;
t455 = -m(5) * t389 + t419 * mrSges(5,1) - t418 * mrSges(5,2) - t425 * t469 + t426 * t468 - t456;
t375 = m(4) * t391 + t439 * mrSges(4,1) - t437 * mrSges(4,2) + t455;
t360 = t448 * t364 + t452 * t375;
t358 = m(3) * t406 + qJDD(1) * mrSges(3,1) - t454 * mrSges(3,2) + t360;
t462 = t452 * t364 - t448 * t375;
t359 = m(3) * t407 - t454 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t462;
t352 = t445 * t358 + t444 * t359;
t365 = t451 * t369 + t447 * t370;
t466 = m(4) * t443 + t365;
t463 = -t444 * t358 + t445 * t359;
t350 = m(2) * t428 - t454 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t463;
t351 = m(2) * t429 + qJDD(1) * mrSges(2,1) - t454 * mrSges(2,2) + t352;
t464 = t453 * t350 - t449 * t351;
t459 = m(3) * t443 + t466;
t458 = -t449 * t350 - t453 * t351;
t412 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t447 + Ifges(5,4) * t451) * t441;
t411 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t447 + Ifges(5,2) * t451) * t441;
t410 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t447 + Ifges(5,6) * t451) * t441;
t398 = Ifges(6,1) * t414 + Ifges(6,4) * t413 + Ifges(6,5) * t440;
t397 = Ifges(6,4) * t414 + Ifges(6,2) * t413 + Ifges(6,6) * t440;
t396 = Ifges(6,5) * t414 + Ifges(6,6) * t413 + Ifges(6,3) * t440;
t373 = mrSges(6,2) * t385 - mrSges(6,3) * t381 + Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t438 + t413 * t396 - t440 * t397;
t372 = -mrSges(6,1) * t385 + mrSges(6,3) * t382 + Ifges(6,4) * t395 + Ifges(6,2) * t394 + Ifges(6,6) * t438 - t414 * t396 + t440 * t398;
t361 = mrSges(5,2) * t389 - mrSges(5,3) * t386 + Ifges(5,1) * t418 + Ifges(5,4) * t419 + Ifges(5,5) * qJDD(4) - pkin(8) * t371 - qJD(4) * t411 - t446 * t372 + t450 * t373 + t410 * t468;
t354 = -mrSges(5,1) * t389 + mrSges(5,3) * t387 + Ifges(5,4) * t418 + Ifges(5,2) * t419 + Ifges(5,6) * qJDD(4) - pkin(4) * t456 + pkin(8) * t460 + qJD(4) * t412 + t450 * t372 + t446 * t373 - t410 * t469;
t353 = Ifges(4,6) * t439 + t437 * Ifges(4,5) - mrSges(4,1) * t443 + mrSges(4,3) * t392 - Ifges(5,5) * t418 - Ifges(5,6) * t419 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t386 + mrSges(5,2) * t387 - Ifges(6,5) * t395 - Ifges(6,6) * t394 - Ifges(6,3) * t438 - t414 * t397 + t413 * t398 - mrSges(6,1) * t381 + mrSges(6,2) * t382 - pkin(4) * t371 - pkin(3) * t365 + (-t447 * t411 + t451 * t412) * t441;
t348 = mrSges(4,2) * t443 - mrSges(4,3) * t391 + Ifges(4,5) * t439 - t437 * Ifges(4,6) - pkin(7) * t365 - t447 * t354 + t451 * t361;
t347 = mrSges(3,2) * t443 - mrSges(3,3) * t406 + Ifges(3,5) * qJDD(1) - t454 * Ifges(3,6) - pkin(6) * t360 + t452 * t348 - t448 * t353;
t346 = -mrSges(3,1) * t443 + mrSges(3,3) * t407 + t454 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t466 + pkin(6) * t462 + t448 * t348 + t452 * t353;
t345 = -mrSges(2,2) * g(1) - mrSges(2,3) * t429 + Ifges(2,5) * qJDD(1) - t454 * Ifges(2,6) - qJ(2) * t352 - t444 * t346 + t445 * t347;
t344 = mrSges(2,1) * g(1) + mrSges(2,3) * t428 + t454 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t459 + qJ(2) * t463 + t445 * t346 + t444 * t347;
t1 = [(-m(1) - m(2)) * g(1) + t459; -m(1) * g(2) + t458; -m(1) * g(3) + t464; pkin(1) * t352 + pkin(2) * t360 + mrSges(3,1) * t406 - mrSges(3,2) * t407 + t447 * t361 + t451 * t354 + pkin(3) * t455 + pkin(7) * t461 + mrSges(4,1) * t391 - mrSges(4,2) * t392 + mrSges(2,1) * t429 - mrSges(2,2) * t428 + Ifges(4,3) * t439 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t464 - t453 * t344 - t449 * t345; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t458 - t449 * t344 + t453 * t345;];
tauB = t1;
