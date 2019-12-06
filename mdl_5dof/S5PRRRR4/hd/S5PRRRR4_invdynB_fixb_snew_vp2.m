% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:52
% EndTime: 2019-12-05 17:07:54
% DurationCPUTime: 2.53s
% Computational Cost: add. (39640->217), mult. (52923->281), div. (0->0), fcn. (33553->10), ass. (0->94)
t418 = qJD(2) + qJD(3);
t424 = sin(qJ(4));
t446 = t418 * t424;
t428 = cos(qJ(4));
t445 = t418 * t428;
t421 = sin(pkin(9));
t422 = cos(pkin(9));
t408 = t421 * g(1) - t422 * g(2);
t409 = -t422 * g(1) - t421 * g(2);
t426 = sin(qJ(2));
t430 = cos(qJ(2));
t390 = t430 * t408 - t426 * t409;
t386 = qJDD(2) * pkin(2) + t390;
t391 = t426 * t408 + t430 * t409;
t431 = qJD(2) ^ 2;
t387 = -t431 * pkin(2) + t391;
t425 = sin(qJ(3));
t429 = cos(qJ(3));
t374 = t425 * t386 + t429 * t387;
t414 = t418 ^ 2;
t416 = qJDD(2) + qJDD(3);
t372 = -t414 * pkin(3) + t416 * pkin(7) + t374;
t420 = -g(3) + qJDD(1);
t368 = -t424 * t372 + t428 * t420;
t443 = qJD(4) * t418;
t441 = t428 * t443;
t400 = t424 * t416 + t441;
t365 = (-t400 + t441) * pkin(8) + (t414 * t424 * t428 + qJDD(4)) * pkin(4) + t368;
t369 = t428 * t372 + t424 * t420;
t401 = t428 * t416 - t424 * t443;
t407 = qJD(4) * pkin(4) - pkin(8) * t446;
t419 = t428 ^ 2;
t366 = -t419 * t414 * pkin(4) + t401 * pkin(8) - qJD(4) * t407 + t369;
t423 = sin(qJ(5));
t427 = cos(qJ(5));
t363 = t427 * t365 - t423 * t366;
t395 = (-t423 * t424 + t427 * t428) * t418;
t377 = t395 * qJD(5) + t427 * t400 + t423 * t401;
t396 = (t423 * t428 + t424 * t427) * t418;
t382 = -t395 * mrSges(6,1) + t396 * mrSges(6,2);
t417 = qJD(4) + qJD(5);
t388 = -t417 * mrSges(6,2) + t395 * mrSges(6,3);
t415 = qJDD(4) + qJDD(5);
t361 = m(6) * t363 + t415 * mrSges(6,1) - t377 * mrSges(6,3) - t396 * t382 + t417 * t388;
t364 = t423 * t365 + t427 * t366;
t376 = -t396 * qJD(5) - t423 * t400 + t427 * t401;
t389 = t417 * mrSges(6,1) - t396 * mrSges(6,3);
t362 = m(6) * t364 - t415 * mrSges(6,2) + t376 * mrSges(6,3) + t395 * t382 - t417 * t389;
t353 = t427 * t361 + t423 * t362;
t399 = (-mrSges(5,1) * t428 + mrSges(5,2) * t424) * t418;
t406 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t445;
t351 = m(5) * t368 + qJDD(4) * mrSges(5,1) - t400 * mrSges(5,3) + qJD(4) * t406 - t399 * t446 + t353;
t405 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t446;
t436 = -t423 * t361 + t427 * t362;
t352 = m(5) * t369 - qJDD(4) * mrSges(5,2) + t401 * mrSges(5,3) - qJD(4) * t405 + t399 * t445 + t436;
t437 = -t424 * t351 + t428 * t352;
t346 = m(4) * t374 - t414 * mrSges(4,1) - t416 * mrSges(4,2) + t437;
t373 = t429 * t386 - t425 * t387;
t434 = -t416 * pkin(3) - t373;
t371 = -t414 * pkin(7) + t434;
t367 = t407 * t446 - t401 * pkin(4) + (-pkin(8) * t419 - pkin(7)) * t414 + t434;
t433 = m(6) * t367 - t376 * mrSges(6,1) + t377 * mrSges(6,2) - t395 * t388 + t396 * t389;
t432 = -m(5) * t371 + t401 * mrSges(5,1) - t400 * mrSges(5,2) - t405 * t446 + t406 * t445 - t433;
t357 = m(4) * t373 + t416 * mrSges(4,1) - t414 * mrSges(4,2) + t432;
t342 = t425 * t346 + t429 * t357;
t340 = m(3) * t390 + qJDD(2) * mrSges(3,1) - t431 * mrSges(3,2) + t342;
t438 = t429 * t346 - t425 * t357;
t341 = m(3) * t391 - t431 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t438;
t334 = t430 * t340 + t426 * t341;
t332 = m(2) * t408 + t334;
t439 = -t426 * t340 + t430 * t341;
t333 = m(2) * t409 + t439;
t444 = t422 * t332 + t421 * t333;
t347 = t428 * t351 + t424 * t352;
t442 = m(4) * t420 + t347;
t440 = -t421 * t332 + t422 * t333;
t435 = m(3) * t420 + t442;
t394 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t424 + Ifges(5,4) * t428) * t418;
t393 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t424 + Ifges(5,2) * t428) * t418;
t392 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t424 + Ifges(5,6) * t428) * t418;
t380 = Ifges(6,1) * t396 + Ifges(6,4) * t395 + Ifges(6,5) * t417;
t379 = Ifges(6,4) * t396 + Ifges(6,2) * t395 + Ifges(6,6) * t417;
t378 = Ifges(6,5) * t396 + Ifges(6,6) * t395 + Ifges(6,3) * t417;
t355 = mrSges(6,2) * t367 - mrSges(6,3) * t363 + Ifges(6,1) * t377 + Ifges(6,4) * t376 + Ifges(6,5) * t415 + t395 * t378 - t417 * t379;
t354 = -mrSges(6,1) * t367 + mrSges(6,3) * t364 + Ifges(6,4) * t377 + Ifges(6,2) * t376 + Ifges(6,6) * t415 - t396 * t378 + t417 * t380;
t343 = mrSges(5,2) * t371 - mrSges(5,3) * t368 + Ifges(5,1) * t400 + Ifges(5,4) * t401 + Ifges(5,5) * qJDD(4) - pkin(8) * t353 - qJD(4) * t393 - t423 * t354 + t427 * t355 + t392 * t445;
t336 = -mrSges(5,1) * t371 + mrSges(5,3) * t369 + Ifges(5,4) * t400 + Ifges(5,2) * t401 + Ifges(5,6) * qJDD(4) - pkin(4) * t433 + pkin(8) * t436 + qJD(4) * t394 + t427 * t354 + t423 * t355 - t392 * t446;
t335 = Ifges(4,6) * t416 + t414 * Ifges(4,5) - mrSges(4,1) * t420 + mrSges(4,3) * t374 - Ifges(5,5) * t400 - Ifges(5,6) * t401 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t368 + mrSges(5,2) * t369 - Ifges(6,5) * t377 - Ifges(6,6) * t376 - Ifges(6,3) * t415 - t396 * t379 + t395 * t380 - mrSges(6,1) * t363 + mrSges(6,2) * t364 - pkin(4) * t353 - pkin(3) * t347 + (-t424 * t393 + t428 * t394) * t418;
t328 = mrSges(4,2) * t420 - mrSges(4,3) * t373 + Ifges(4,5) * t416 - t414 * Ifges(4,6) - pkin(7) * t347 - t424 * t336 + t428 * t343;
t327 = mrSges(3,2) * t420 - mrSges(3,3) * t390 + Ifges(3,5) * qJDD(2) - t431 * Ifges(3,6) - pkin(6) * t342 + t429 * t328 - t425 * t335;
t326 = -mrSges(3,1) * t420 + mrSges(3,3) * t391 + t431 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t442 + pkin(6) * t438 + t425 * t328 + t429 * t335;
t325 = mrSges(2,2) * t420 - mrSges(2,3) * t408 - pkin(5) * t334 - t426 * t326 + t430 * t327;
t324 = -mrSges(2,1) * t420 + mrSges(2,3) * t409 - pkin(1) * t435 + pkin(5) * t439 + t430 * t326 + t426 * t327;
t1 = [-m(1) * g(1) + t440; -m(1) * g(2) + t444; -m(1) * g(3) + m(2) * t420 + t435; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t444 - t421 * t324 + t422 * t325; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t440 + t422 * t324 + t421 * t325; pkin(1) * t334 + mrSges(2,1) * t408 - mrSges(2,2) * t409 + pkin(2) * t342 + mrSges(3,1) * t390 - mrSges(3,2) * t391 + t424 * t343 + t428 * t336 + pkin(3) * t432 + pkin(7) * t437 + mrSges(4,1) * t373 - mrSges(4,2) * t374 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t416 + Ifges(3,3) * qJDD(2);];
tauB = t1;
