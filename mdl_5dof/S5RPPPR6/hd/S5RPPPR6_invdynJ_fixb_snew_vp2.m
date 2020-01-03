% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:44
% EndTime: 2019-12-31 17:47:45
% DurationCPUTime: 1.18s
% Computational Cost: add. (3463->195), mult. (8902->256), div. (0->0), fcn. (5074->8), ass. (0->100)
t447 = -2 * qJD(4);
t398 = qJD(1) ^ 2;
t395 = sin(qJ(1));
t397 = cos(qJ(1));
t414 = -t397 * g(1) - t395 * g(2);
t446 = -t398 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t414;
t391 = sin(pkin(7));
t430 = qJD(1) * t391;
t382 = -0.2e1 * qJD(3) * t430;
t393 = cos(pkin(7));
t421 = t395 * g(1) - t397 * g(2);
t413 = qJDD(2) - t421;
t388 = t391 ^ 2;
t389 = t393 ^ 2;
t431 = -t388 - t389;
t438 = t391 * qJ(3);
t343 = t382 + (t431 * pkin(3) - qJ(2)) * t398 + (-t438 - pkin(1) + (-pkin(2) - qJ(4)) * t393) * qJDD(1) + t413;
t429 = qJD(1) * t393;
t445 = t429 * t447 + t343;
t358 = -t391 * g(3) + t446 * t393;
t410 = -t393 * pkin(2) - t438;
t374 = t410 * qJD(1);
t349 = -t374 * t429 - t358;
t357 = -t393 * g(3) - t446 * t391;
t444 = -mrSges(4,1) - mrSges(3,3);
t390 = sin(pkin(8));
t443 = Ifges(5,1) * t390;
t392 = cos(pkin(8));
t442 = Ifges(5,4) * t392;
t441 = qJ(4) * t398;
t440 = t388 * t398;
t439 = t390 * t393;
t437 = t391 * t398;
t348 = t374 * t430 + qJDD(3) - t357;
t340 = (pkin(3) * qJDD(1) - t393 * t441) * t391 + t348;
t436 = t392 * t340;
t435 = t392 * t393;
t337 = t390 * t340 + t445 * t392;
t412 = mrSges(5,1) * t392 - mrSges(5,2) * t390;
t366 = t412 * t429;
t408 = t391 * mrSges(5,1) + mrSges(5,3) * t439;
t371 = t408 * qJD(1);
t407 = -t391 * mrSges(5,2) - mrSges(5,3) * t435;
t367 = (t392 * pkin(4) + t390 * pkin(6)) * t429;
t425 = t392 * t429;
t427 = qJDD(1) * t391;
t335 = -pkin(4) * t440 + pkin(6) * t427 - t367 * t425 + t337;
t426 = qJDD(1) * t393;
t341 = pkin(3) * t426 - t389 * t441 + qJDD(4) - t349;
t338 = ((qJDD(1) * t390 + t392 * t437) * pkin(6) + (qJDD(1) * t392 - t390 * t437) * pkin(4)) * t393 + t341;
t394 = sin(qJ(5));
t396 = cos(qJ(5));
t332 = -t394 * t335 + t396 * t338;
t404 = t391 * t396 + t394 * t439;
t364 = t404 * qJD(1);
t405 = t391 * t394 - t396 * t439;
t365 = t405 * qJD(1);
t350 = -t364 * mrSges(6,1) + t365 * mrSges(6,2);
t353 = t364 * qJD(5) + t405 * qJDD(1);
t380 = qJD(5) + t425;
t355 = -t380 * mrSges(6,2) + t364 * mrSges(6,3);
t379 = t392 * t426 + qJDD(5);
t330 = m(6) * t332 + t379 * mrSges(6,1) - t353 * mrSges(6,3) - t365 * t350 + t380 * t355;
t333 = t396 * t335 + t394 * t338;
t352 = -t365 * qJD(5) + t404 * qJDD(1);
t356 = t380 * mrSges(6,1) - t365 * mrSges(6,3);
t331 = m(6) * t333 - t379 * mrSges(6,2) + t352 * mrSges(6,3) + t364 * t350 - t380 * t356;
t417 = -t394 * t330 + t396 * t331;
t320 = m(5) * t337 + t407 * qJDD(1) + (-t366 * t435 - t371 * t391) * qJD(1) + t417;
t336 = -t445 * t390 + t436;
t372 = t407 * qJD(1);
t334 = -pkin(4) * t427 - pkin(6) * t440 - t436 + (t343 + (t447 - t367) * t429) * t390;
t400 = -m(6) * t334 + t352 * mrSges(6,1) - t353 * mrSges(6,2) + t364 * t355 - t365 * t356;
t326 = m(5) * t336 + t408 * qJDD(1) + (t366 * t439 + t372 * t391) * qJD(1) + t400;
t434 = t390 * t320 + t392 * t326;
t322 = t396 * t330 + t394 * t331;
t411 = t393 * mrSges(4,2) - t391 * mrSges(4,3);
t375 = t411 * qJD(1);
t432 = -t375 - (-t393 * mrSges(3,1) + t391 * mrSges(3,2)) * qJD(1);
t419 = t431 * mrSges(4,1);
t418 = t392 * t320 - t390 * t326;
t416 = m(4) * t348 + t434;
t415 = m(5) * t341 + t322;
t409 = -t371 * t390 + t372 * t392;
t402 = -t398 * qJ(2) + t413;
t354 = t382 + (-pkin(1) + t410) * qJDD(1) + t402;
t406 = m(4) * t354 + t418;
t403 = -Ifges(5,5) * t390 - Ifges(5,6) * t392 + Ifges(3,4) + Ifges(4,6);
t401 = (((-Ifges(5,4) * t390 - Ifges(5,2) * t392) * t390 - (-t442 - t443) * t392 - Ifges(3,6) + Ifges(4,5)) * t393 + (-Ifges(5,5) * t392 + Ifges(5,6) * t390 + Ifges(4,4) - Ifges(3,5)) * t391) * qJD(1);
t345 = Ifges(6,4) * t365 + Ifges(6,2) * t364 + Ifges(6,6) * t380;
t346 = Ifges(6,1) * t365 + Ifges(6,4) * t364 + Ifges(6,5) * t380;
t399 = mrSges(6,1) * t332 - mrSges(6,2) * t333 + Ifges(6,5) * t353 + Ifges(6,6) * t352 + Ifges(6,3) * t379 + t365 * t345 - t364 * t346;
t370 = -qJDD(1) * pkin(1) + t402;
t344 = Ifges(6,5) * t365 + Ifges(6,6) * t364 + Ifges(6,3) * t380;
t324 = mrSges(6,2) * t334 - mrSges(6,3) * t332 + Ifges(6,1) * t353 + Ifges(6,4) * t352 + Ifges(6,5) * t379 + t364 * t344 - t380 * t345;
t323 = -mrSges(6,1) * t334 + mrSges(6,3) * t333 + Ifges(6,4) * t353 + Ifges(6,2) * t352 + Ifges(6,6) * t379 - t365 * t344 + t380 * t346;
t321 = (t409 * qJD(1) + t412 * qJDD(1)) * t393 + t415;
t317 = t411 * qJDD(1) + t398 * t419 + t406;
t316 = m(3) * t370 + (t431 * mrSges(3,3) + t419) * t398 + ((-mrSges(3,1) + mrSges(4,2)) * t393 + (mrSges(3,2) - mrSges(4,3)) * t391) * qJDD(1) + t406;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t421 - mrSges(2,2) * t414 + t391 * (-qJ(3) * t317 - mrSges(3,3) * t357 + mrSges(3,2) * t370 + pkin(3) * t434 + mrSges(4,1) * t348 - mrSges(4,3) * t354 + pkin(4) * t400 + pkin(6) * t417 + mrSges(5,1) * t336 - mrSges(5,2) * t337 + t394 * t324 + t396 * t323 + (Ifges(3,1) + Ifges(4,2) + Ifges(5,3)) * t427 + (-t401 * qJD(1) + t403 * qJDD(1)) * t393) + t393 * (-mrSges(3,1) * t370 + mrSges(3,3) * t358 - mrSges(4,1) * t349 + mrSges(4,2) * t354 - t390 * (mrSges(5,2) * t341 - mrSges(5,3) * t336 - pkin(6) * t322 - t394 * t323 + t396 * t324) - t392 * (-mrSges(5,1) * t341 + mrSges(5,3) * t337 - pkin(4) * t322 - t399) + pkin(3) * t321 - qJ(4) * t418 - pkin(2) * t317 + t401 * t430 + (t391 * t403 + (t392 ^ 2 * Ifges(5,2) + Ifges(3,2) + Ifges(4,3) + (0.2e1 * t442 + t443) * t390) * t393) * qJDD(1)) - pkin(1) * t316 + qJ(2) * (-t391 * (m(3) * t357 + (t432 * qJD(1) + t444 * qJDD(1)) * t391 - t416) + (m(3) * t358 - m(4) * t349 + (t412 - t444) * t426 + (t409 - t432) * t429 + t415) * t393); t316; (qJDD(1) * mrSges(4,1) + qJD(1) * t375) * t391 + t416; t321; t399;];
tauJ = t1;
