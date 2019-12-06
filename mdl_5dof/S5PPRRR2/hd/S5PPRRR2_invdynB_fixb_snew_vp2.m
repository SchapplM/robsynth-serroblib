% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:22
% EndTime: 2019-12-05 15:14:25
% DurationCPUTime: 1.96s
% Computational Cost: add. (21984->206), mult. (39914->269), div. (0->0), fcn. (26262->10), ass. (0->89)
t387 = sin(pkin(8));
t389 = cos(pkin(8));
t378 = -t389 * g(1) - t387 * g(2);
t385 = -g(3) + qJDD(1);
t386 = sin(pkin(9));
t388 = cos(pkin(9));
t362 = -t386 * t378 + t388 * t385;
t363 = t388 * t378 + t386 * t385;
t392 = sin(qJ(3));
t395 = cos(qJ(3));
t349 = t392 * t362 + t395 * t363;
t396 = qJD(3) ^ 2;
t346 = -t396 * pkin(3) + qJDD(3) * pkin(6) + t349;
t377 = t387 * g(1) - t389 * g(2);
t376 = qJDD(2) - t377;
t391 = sin(qJ(4));
t394 = cos(qJ(4));
t342 = -t391 * t346 + t394 * t376;
t407 = qJD(3) * qJD(4);
t406 = t394 * t407;
t374 = t391 * qJDD(3) + t406;
t339 = (-t374 + t406) * pkin(7) + (t391 * t394 * t396 + qJDD(4)) * pkin(4) + t342;
t343 = t394 * t346 + t391 * t376;
t375 = t394 * qJDD(3) - t391 * t407;
t409 = qJD(3) * t391;
t381 = qJD(4) * pkin(4) - pkin(7) * t409;
t384 = t394 ^ 2;
t340 = -t384 * t396 * pkin(4) + t375 * pkin(7) - qJD(4) * t381 + t343;
t390 = sin(qJ(5));
t393 = cos(qJ(5));
t337 = t393 * t339 - t390 * t340;
t367 = (-t390 * t391 + t393 * t394) * qJD(3);
t351 = t367 * qJD(5) + t393 * t374 + t390 * t375;
t368 = (t390 * t394 + t391 * t393) * qJD(3);
t356 = -t367 * mrSges(6,1) + t368 * mrSges(6,2);
t383 = qJD(4) + qJD(5);
t360 = -t383 * mrSges(6,2) + t367 * mrSges(6,3);
t382 = qJDD(4) + qJDD(5);
t335 = m(6) * t337 + t382 * mrSges(6,1) - t351 * mrSges(6,3) - t368 * t356 + t383 * t360;
t338 = t390 * t339 + t393 * t340;
t350 = -t368 * qJD(5) - t390 * t374 + t393 * t375;
t361 = t383 * mrSges(6,1) - t368 * mrSges(6,3);
t336 = m(6) * t338 - t382 * mrSges(6,2) + t350 * mrSges(6,3) + t367 * t356 - t383 * t361;
t327 = t393 * t335 + t390 * t336;
t373 = (-mrSges(5,1) * t394 + mrSges(5,2) * t391) * qJD(3);
t408 = qJD(3) * t394;
t380 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t408;
t325 = m(5) * t342 + qJDD(4) * mrSges(5,1) - t374 * mrSges(5,3) + qJD(4) * t380 - t373 * t409 + t327;
t379 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t409;
t401 = -t390 * t335 + t393 * t336;
t326 = m(5) * t343 - qJDD(4) * mrSges(5,2) + t375 * mrSges(5,3) - qJD(4) * t379 + t373 * t408 + t401;
t402 = -t391 * t325 + t394 * t326;
t318 = m(4) * t349 - t396 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t402;
t348 = t395 * t362 - t392 * t363;
t399 = -qJDD(3) * pkin(3) - t348;
t345 = -t396 * pkin(6) + t399;
t341 = t381 * t409 - t375 * pkin(4) + (-pkin(7) * t384 - pkin(6)) * t396 + t399;
t398 = m(6) * t341 - t350 * mrSges(6,1) + t351 * mrSges(6,2) - t367 * t360 + t368 * t361;
t397 = -m(5) * t345 + t375 * mrSges(5,1) - t374 * mrSges(5,2) - t379 * t409 + t380 * t408 - t398;
t331 = m(4) * t348 + qJDD(3) * mrSges(4,1) - t396 * mrSges(4,2) + t397;
t314 = t392 * t318 + t395 * t331;
t312 = m(3) * t362 + t314;
t403 = t395 * t318 - t392 * t331;
t313 = m(3) * t363 + t403;
t404 = -t386 * t312 + t388 * t313;
t305 = m(2) * t378 + t404;
t321 = t394 * t325 + t391 * t326;
t400 = (-m(3) - m(4)) * t376 - t321;
t320 = m(2) * t377 + t400;
t410 = t387 * t305 + t389 * t320;
t306 = t388 * t312 + t386 * t313;
t405 = t389 * t305 - t387 * t320;
t366 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t391 + Ifges(5,4) * t394) * qJD(3);
t365 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t391 + Ifges(5,2) * t394) * qJD(3);
t364 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t391 + Ifges(5,6) * t394) * qJD(3);
t354 = Ifges(6,1) * t368 + Ifges(6,4) * t367 + Ifges(6,5) * t383;
t353 = Ifges(6,4) * t368 + Ifges(6,2) * t367 + Ifges(6,6) * t383;
t352 = Ifges(6,5) * t368 + Ifges(6,6) * t367 + Ifges(6,3) * t383;
t329 = mrSges(6,2) * t341 - mrSges(6,3) * t337 + Ifges(6,1) * t351 + Ifges(6,4) * t350 + Ifges(6,5) * t382 + t367 * t352 - t383 * t353;
t328 = -mrSges(6,1) * t341 + mrSges(6,3) * t338 + Ifges(6,4) * t351 + Ifges(6,2) * t350 + Ifges(6,6) * t382 - t368 * t352 + t383 * t354;
t315 = mrSges(5,2) * t345 - mrSges(5,3) * t342 + Ifges(5,1) * t374 + Ifges(5,4) * t375 + Ifges(5,5) * qJDD(4) - pkin(7) * t327 - qJD(4) * t365 - t390 * t328 + t393 * t329 + t364 * t408;
t308 = -mrSges(5,1) * t345 + mrSges(5,3) * t343 + Ifges(5,4) * t374 + Ifges(5,2) * t375 + Ifges(5,6) * qJDD(4) - pkin(4) * t398 + pkin(7) * t401 + qJD(4) * t366 + t393 * t328 + t390 * t329 - t364 * t409;
t307 = Ifges(4,6) * qJDD(3) + t396 * Ifges(4,5) - mrSges(4,1) * t376 + mrSges(4,3) * t349 - Ifges(5,5) * t374 - Ifges(5,6) * t375 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t342 + mrSges(5,2) * t343 - Ifges(6,5) * t351 - Ifges(6,6) * t350 - Ifges(6,3) * t382 - t368 * t353 + t367 * t354 - mrSges(6,1) * t337 + mrSges(6,2) * t338 - pkin(4) * t327 - pkin(3) * t321 + (-t391 * t365 + t394 * t366) * qJD(3);
t302 = mrSges(4,2) * t376 - mrSges(4,3) * t348 + Ifges(4,5) * qJDD(3) - t396 * Ifges(4,6) - pkin(6) * t321 - t391 * t308 + t394 * t315;
t301 = mrSges(3,2) * t376 - mrSges(3,3) * t362 - pkin(5) * t314 + t395 * t302 - t392 * t307;
t300 = -mrSges(2,1) * t385 - mrSges(3,1) * t362 - mrSges(4,1) * t348 + mrSges(3,2) * t363 + mrSges(4,2) * t349 + mrSges(2,3) * t378 - Ifges(4,3) * qJDD(3) - pkin(1) * t306 - pkin(2) * t314 - pkin(3) * t397 - pkin(6) * t402 - t394 * t308 - t391 * t315;
t299 = -mrSges(3,1) * t376 + mrSges(3,3) * t363 + t392 * t302 + t395 * t307 - pkin(2) * (m(4) * t376 + t321) + pkin(5) * t403;
t298 = mrSges(2,2) * t385 - mrSges(2,3) * t377 - qJ(2) * t306 - t386 * t299 + t388 * t301;
t1 = [-m(1) * g(1) + t405; -m(1) * g(2) + t410; -m(1) * g(3) + m(2) * t385 + t306; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t410 + t389 * t298 - t387 * t300; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t405 + t387 * t298 + t389 * t300; -mrSges(1,1) * g(2) + mrSges(2,1) * t377 + mrSges(1,2) * g(1) - mrSges(2,2) * t378 + pkin(1) * t400 + qJ(2) * t404 + t388 * t299 + t386 * t301;];
tauB = t1;
