% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:27
% EndTime: 2019-12-31 17:24:30
% DurationCPUTime: 2.20s
% Computational Cost: add. (23613->240), mult. (51052->309), div. (0->0), fcn. (33556->8), ass. (0->95)
t437 = qJD(1) ^ 2;
t451 = pkin(2) * t437;
t432 = sin(qJ(1));
t436 = cos(qJ(1));
t422 = -t436 * g(1) - t432 * g(2);
t411 = -t437 * pkin(1) + qJDD(1) * pkin(5) + t422;
t431 = sin(qJ(2));
t450 = t431 * t411;
t435 = cos(qJ(2));
t446 = qJD(1) * qJD(2);
t416 = t431 * qJDD(1) + t435 * t446;
t383 = qJDD(2) * pkin(2) - t416 * pkin(6) - t450 + (pkin(6) * t446 + t431 * t451 - g(3)) * t435;
t399 = -t431 * g(3) + t435 * t411;
t417 = t435 * qJDD(1) - t431 * t446;
t448 = qJD(1) * t431;
t420 = qJD(2) * pkin(2) - pkin(6) * t448;
t428 = t435 ^ 2;
t384 = t417 * pkin(6) - qJD(2) * t420 - t428 * t451 + t399;
t430 = sin(qJ(3));
t434 = cos(qJ(3));
t372 = t434 * t383 - t430 * t384;
t408 = (-t430 * t431 + t434 * t435) * qJD(1);
t387 = t408 * qJD(3) + t434 * t416 + t430 * t417;
t409 = (t430 * t435 + t431 * t434) * qJD(1);
t426 = qJDD(2) + qJDD(3);
t427 = qJD(2) + qJD(3);
t364 = (t408 * t427 - t387) * pkin(7) + (t408 * t409 + t426) * pkin(3) + t372;
t373 = t430 * t383 + t434 * t384;
t386 = -t409 * qJD(3) - t430 * t416 + t434 * t417;
t402 = t427 * pkin(3) - t409 * pkin(7);
t404 = t408 ^ 2;
t365 = -t404 * pkin(3) + t386 * pkin(7) - t427 * t402 + t373;
t429 = sin(qJ(4));
t433 = cos(qJ(4));
t362 = t433 * t364 - t429 * t365;
t395 = t433 * t408 - t429 * t409;
t371 = t395 * qJD(4) + t429 * t386 + t433 * t387;
t396 = t429 * t408 + t433 * t409;
t379 = -t395 * mrSges(5,1) + t396 * mrSges(5,2);
t424 = qJD(4) + t427;
t389 = -t424 * mrSges(5,2) + t395 * mrSges(5,3);
t423 = qJDD(4) + t426;
t360 = m(5) * t362 + t423 * mrSges(5,1) - t371 * mrSges(5,3) - t396 * t379 + t424 * t389;
t363 = t429 * t364 + t433 * t365;
t370 = -t396 * qJD(4) + t433 * t386 - t429 * t387;
t390 = t424 * mrSges(5,1) - t396 * mrSges(5,3);
t361 = m(5) * t363 - t423 * mrSges(5,2) + t370 * mrSges(5,3) + t395 * t379 - t424 * t390;
t352 = t433 * t360 + t429 * t361;
t397 = -t408 * mrSges(4,1) + t409 * mrSges(4,2);
t400 = -t427 * mrSges(4,2) + t408 * mrSges(4,3);
t350 = m(4) * t372 + t426 * mrSges(4,1) - t387 * mrSges(4,3) - t409 * t397 + t427 * t400 + t352;
t401 = t427 * mrSges(4,1) - t409 * mrSges(4,3);
t442 = -t429 * t360 + t433 * t361;
t351 = m(4) * t373 - t426 * mrSges(4,2) + t386 * mrSges(4,3) + t408 * t397 - t427 * t401 + t442;
t346 = t434 * t350 + t430 * t351;
t398 = -t435 * g(3) - t450;
t415 = (-mrSges(3,1) * t435 + mrSges(3,2) * t431) * qJD(1);
t447 = qJD(1) * t435;
t419 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t447;
t344 = m(3) * t398 + qJDD(2) * mrSges(3,1) - t416 * mrSges(3,3) + qJD(2) * t419 - t415 * t448 + t346;
t418 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t448;
t443 = -t430 * t350 + t434 * t351;
t345 = m(3) * t399 - qJDD(2) * mrSges(3,2) + t417 * mrSges(3,3) - qJD(2) * t418 + t415 * t447 + t443;
t444 = -t431 * t344 + t435 * t345;
t337 = m(2) * t422 - t437 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t444;
t421 = t432 * g(1) - t436 * g(2);
t440 = -qJDD(1) * pkin(1) - t421;
t410 = -t437 * pkin(5) + t440;
t388 = -t417 * pkin(2) + t420 * t448 + (-pkin(6) * t428 - pkin(5)) * t437 + t440;
t367 = -t386 * pkin(3) - t404 * pkin(7) + t409 * t402 + t388;
t441 = m(5) * t367 - t370 * mrSges(5,1) + t371 * mrSges(5,2) - t395 * t389 + t396 * t390;
t439 = m(4) * t388 - t386 * mrSges(4,1) + t387 * mrSges(4,2) - t408 * t400 + t409 * t401 + t441;
t438 = -m(3) * t410 + t417 * mrSges(3,1) - t416 * mrSges(3,2) - t418 * t448 + t419 * t447 - t439;
t356 = m(2) * t421 + qJDD(1) * mrSges(2,1) - t437 * mrSges(2,2) + t438;
t449 = t432 * t337 + t436 * t356;
t338 = t435 * t344 + t431 * t345;
t445 = t436 * t337 - t432 * t356;
t407 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t431 + Ifges(3,4) * t435) * qJD(1);
t406 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t431 + Ifges(3,2) * t435) * qJD(1);
t405 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t431 + Ifges(3,6) * t435) * qJD(1);
t393 = Ifges(4,1) * t409 + Ifges(4,4) * t408 + Ifges(4,5) * t427;
t392 = Ifges(4,4) * t409 + Ifges(4,2) * t408 + Ifges(4,6) * t427;
t391 = Ifges(4,5) * t409 + Ifges(4,6) * t408 + Ifges(4,3) * t427;
t376 = Ifges(5,1) * t396 + Ifges(5,4) * t395 + Ifges(5,5) * t424;
t375 = Ifges(5,4) * t396 + Ifges(5,2) * t395 + Ifges(5,6) * t424;
t374 = Ifges(5,5) * t396 + Ifges(5,6) * t395 + Ifges(5,3) * t424;
t354 = mrSges(5,2) * t367 - mrSges(5,3) * t362 + Ifges(5,1) * t371 + Ifges(5,4) * t370 + Ifges(5,5) * t423 + t395 * t374 - t424 * t375;
t353 = -mrSges(5,1) * t367 + mrSges(5,3) * t363 + Ifges(5,4) * t371 + Ifges(5,2) * t370 + Ifges(5,6) * t423 - t396 * t374 + t424 * t376;
t340 = mrSges(4,2) * t388 - mrSges(4,3) * t372 + Ifges(4,1) * t387 + Ifges(4,4) * t386 + Ifges(4,5) * t426 - pkin(7) * t352 - t429 * t353 + t433 * t354 + t408 * t391 - t427 * t392;
t339 = -mrSges(4,1) * t388 + mrSges(4,3) * t373 + Ifges(4,4) * t387 + Ifges(4,2) * t386 + Ifges(4,6) * t426 - pkin(3) * t441 + pkin(7) * t442 + t433 * t353 + t429 * t354 - t409 * t391 + t427 * t393;
t334 = mrSges(3,2) * t410 - mrSges(3,3) * t398 + Ifges(3,1) * t416 + Ifges(3,4) * t417 + Ifges(3,5) * qJDD(2) - pkin(6) * t346 - qJD(2) * t406 - t430 * t339 + t434 * t340 + t405 * t447;
t333 = -mrSges(3,1) * t410 + mrSges(3,3) * t399 + Ifges(3,4) * t416 + Ifges(3,2) * t417 + Ifges(3,6) * qJDD(2) - pkin(2) * t439 + pkin(6) * t443 + qJD(2) * t407 + t434 * t339 + t430 * t340 - t405 * t448;
t332 = Ifges(2,6) * qJDD(1) + (-t431 * t406 + t435 * t407) * qJD(1) - Ifges(3,3) * qJDD(2) + t437 * Ifges(2,5) - Ifges(4,3) * t426 - Ifges(3,5) * t416 - Ifges(3,6) * t417 + mrSges(2,3) * t422 - Ifges(5,3) * t423 - t409 * t392 + t408 * t393 - mrSges(3,1) * t398 + mrSges(3,2) * t399 - Ifges(4,5) * t387 + t395 * t376 - t396 * t375 - Ifges(4,6) * t386 - Ifges(5,6) * t370 - Ifges(5,5) * t371 - mrSges(4,1) * t372 + mrSges(4,2) * t373 + mrSges(5,2) * t363 + mrSges(2,1) * g(3) - mrSges(5,1) * t362 - pkin(3) * t352 - pkin(2) * t346 - pkin(1) * t338;
t331 = -mrSges(2,2) * g(3) - mrSges(2,3) * t421 + Ifges(2,5) * qJDD(1) - t437 * Ifges(2,6) - pkin(5) * t338 - t431 * t333 + t435 * t334;
t1 = [-m(1) * g(1) + t445; -m(1) * g(2) + t449; (-m(1) - m(2)) * g(3) + t338; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t449 + t436 * t331 - t432 * t332; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t445 + t432 * t331 + t436 * t332; -mrSges(1,1) * g(2) + mrSges(2,1) * t421 + mrSges(1,2) * g(1) - mrSges(2,2) * t422 + Ifges(2,3) * qJDD(1) + pkin(1) * t438 + pkin(5) * t444 + t435 * t333 + t431 * t334;];
tauB = t1;
