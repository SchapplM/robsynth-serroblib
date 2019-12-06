% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:44
% EndTime: 2019-12-05 16:15:46
% DurationCPUTime: 2.51s
% Computational Cost: add. (33479->195), mult. (47918->250), div. (0->0), fcn. (31463->10), ass. (0->93)
t415 = qJD(2) + qJD(3);
t411 = t415 ^ 2;
t419 = cos(pkin(9));
t451 = pkin(4) * t419;
t417 = sin(pkin(9));
t450 = mrSges(5,2) * t417;
t414 = t419 ^ 2;
t449 = t411 * t414;
t412 = qJDD(2) + qJDD(3);
t448 = t412 * t419;
t434 = Ifges(5,5) * t417 + Ifges(5,6) * t419;
t447 = t411 * t434;
t418 = sin(pkin(8));
t420 = cos(pkin(8));
t400 = t418 * g(1) - t420 * g(2);
t401 = -t420 * g(1) - t418 * g(2);
t423 = sin(qJ(2));
t426 = cos(qJ(2));
t389 = t426 * t400 - t423 * t401;
t387 = qJDD(2) * pkin(2) + t389;
t390 = t423 * t400 + t426 * t401;
t427 = qJD(2) ^ 2;
t388 = -t427 * pkin(2) + t390;
t422 = sin(qJ(3));
t425 = cos(qJ(3));
t375 = t422 * t387 + t425 * t388;
t373 = -t411 * pkin(3) + t412 * qJ(4) + t375;
t416 = -g(3) + qJDD(1);
t444 = qJD(4) * t415;
t445 = t419 * t416 - 0.2e1 * t417 * t444;
t366 = (-pkin(7) * t412 + t411 * t451 - t373) * t417 + t445;
t370 = t417 * t416 + (t373 + 0.2e1 * t444) * t419;
t367 = -pkin(4) * t449 + pkin(7) * t448 + t370;
t421 = sin(qJ(5));
t424 = cos(qJ(5));
t364 = t424 * t366 - t421 * t367;
t430 = -t417 * t421 + t419 * t424;
t393 = t430 * t415;
t431 = t417 * t424 + t419 * t421;
t394 = t431 * t415;
t381 = -t393 * mrSges(6,1) + t394 * mrSges(6,2);
t383 = t393 * qJD(5) + t431 * t412;
t391 = -qJD(5) * mrSges(6,2) + t393 * mrSges(6,3);
t362 = m(6) * t364 + qJDD(5) * mrSges(6,1) - t383 * mrSges(6,3) + qJD(5) * t391 - t394 * t381;
t365 = t421 * t366 + t424 * t367;
t382 = -t394 * qJD(5) + t430 * t412;
t392 = qJD(5) * mrSges(6,1) - t394 * mrSges(6,3);
t363 = m(6) * t365 - qJDD(5) * mrSges(6,2) + t382 * mrSges(6,3) - qJD(5) * t392 + t393 * t381;
t354 = t424 * t362 + t421 * t363;
t369 = -t417 * t373 + t445;
t432 = mrSges(5,3) * t412 + (-mrSges(5,1) * t419 + t450) * t411;
t352 = m(5) * t369 - t432 * t417 + t354;
t438 = -t421 * t362 + t424 * t363;
t353 = m(5) * t370 + t432 * t419 + t438;
t439 = -t417 * t352 + t419 * t353;
t347 = m(4) * t375 - t411 * mrSges(4,1) - t412 * mrSges(4,2) + t439;
t374 = t425 * t387 - t422 * t388;
t433 = qJDD(4) - t374;
t372 = -t412 * pkin(3) - t411 * qJ(4) + t433;
t413 = t417 ^ 2;
t368 = (-pkin(3) - t451) * t412 + (-qJ(4) + (-t413 - t414) * pkin(7)) * t411 + t433;
t429 = m(6) * t368 - t382 * mrSges(6,1) + t383 * mrSges(6,2) - t393 * t391 + t394 * t392;
t428 = -m(5) * t372 + mrSges(5,1) * t448 - t429 + (t411 * t413 + t449) * mrSges(5,3);
t358 = m(4) * t374 - t411 * mrSges(4,2) + (mrSges(4,1) - t450) * t412 + t428;
t343 = t422 * t347 + t425 * t358;
t340 = m(3) * t389 + qJDD(2) * mrSges(3,1) - t427 * mrSges(3,2) + t343;
t440 = t425 * t347 - t422 * t358;
t341 = m(3) * t390 - t427 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t440;
t335 = t426 * t340 + t423 * t341;
t333 = m(2) * t400 + t335;
t441 = -t423 * t340 + t426 * t341;
t334 = m(2) * t401 + t441;
t446 = t420 * t333 + t418 * t334;
t348 = t419 * t352 + t417 * t353;
t443 = m(4) * t416 + t348;
t442 = -t418 * t333 + t420 * t334;
t437 = m(3) * t416 + t443;
t436 = Ifges(5,1) * t417 + Ifges(5,4) * t419;
t435 = Ifges(5,4) * t417 + Ifges(5,2) * t419;
t378 = Ifges(6,1) * t394 + Ifges(6,4) * t393 + Ifges(6,5) * qJD(5);
t377 = Ifges(6,4) * t394 + Ifges(6,2) * t393 + Ifges(6,6) * qJD(5);
t376 = Ifges(6,5) * t394 + Ifges(6,6) * t393 + Ifges(6,3) * qJD(5);
t356 = mrSges(6,2) * t368 - mrSges(6,3) * t364 + Ifges(6,1) * t383 + Ifges(6,4) * t382 + Ifges(6,5) * qJDD(5) - qJD(5) * t377 + t393 * t376;
t355 = -mrSges(6,1) * t368 + mrSges(6,3) * t365 + Ifges(6,4) * t383 + Ifges(6,2) * t382 + Ifges(6,6) * qJDD(5) + qJD(5) * t378 - t394 * t376;
t344 = mrSges(5,2) * t372 - mrSges(5,3) * t369 - pkin(7) * t354 - t421 * t355 + t424 * t356 + t436 * t412 + t419 * t447;
t342 = -mrSges(5,1) * t372 + mrSges(5,3) * t370 - pkin(4) * t429 + pkin(7) * t438 + t424 * t355 + t421 * t356 + t435 * t412 - t417 * t447;
t336 = -mrSges(4,1) * t416 - mrSges(5,1) * t369 - mrSges(6,1) * t364 + mrSges(5,2) * t370 + mrSges(6,2) * t365 + mrSges(4,3) * t375 - Ifges(6,5) * t383 - Ifges(6,6) * t382 - Ifges(6,3) * qJDD(5) - pkin(3) * t348 - pkin(4) * t354 - t394 * t377 + t393 * t378 + (Ifges(4,6) - t434) * t412 + (-t417 * t435 + t419 * t436 + Ifges(4,5)) * t411;
t329 = mrSges(4,2) * t416 - mrSges(4,3) * t374 + Ifges(4,5) * t412 - t411 * Ifges(4,6) - qJ(4) * t348 - t417 * t342 + t419 * t344;
t328 = mrSges(3,2) * t416 - mrSges(3,3) * t389 + Ifges(3,5) * qJDD(2) - t427 * Ifges(3,6) - pkin(6) * t343 + t425 * t329 - t422 * t336;
t327 = -mrSges(3,1) * t416 + mrSges(3,3) * t390 + t427 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t443 + pkin(6) * t440 + t422 * t329 + t425 * t336;
t326 = mrSges(2,2) * t416 - mrSges(2,3) * t400 - pkin(5) * t335 - t423 * t327 + t426 * t328;
t325 = -mrSges(2,1) * t416 + mrSges(2,3) * t401 - pkin(1) * t437 + pkin(5) * t441 + t426 * t327 + t423 * t328;
t1 = [-m(1) * g(1) + t442; -m(1) * g(2) + t446; -m(1) * g(3) + m(2) * t416 + t437; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t446 - t418 * t325 + t420 * t326; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t442 + t420 * t325 + t418 * t326; pkin(1) * t335 + mrSges(2,1) * t400 - mrSges(2,2) * t401 + pkin(2) * t343 + mrSges(3,1) * t389 - mrSges(3,2) * t390 + pkin(3) * (-t412 * t450 + t428) + qJ(4) * t439 + t417 * t344 + t419 * t342 + mrSges(4,1) * t374 - mrSges(4,2) * t375 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t412 + Ifges(3,3) * qJDD(2);];
tauB = t1;
