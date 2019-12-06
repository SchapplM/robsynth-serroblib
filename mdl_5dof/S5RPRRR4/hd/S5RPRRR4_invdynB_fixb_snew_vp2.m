% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR4
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:29
% EndTime: 2019-12-05 18:14:31
% DurationCPUTime: 2.18s
% Computational Cost: add. (36781->185), mult. (47778->233), div. (0->0), fcn. (25560->10), ass. (0->82)
t372 = qJD(1) + qJD(3);
t363 = qJD(4) + t372;
t376 = sin(qJ(5));
t397 = t363 * t376;
t380 = cos(qJ(5));
t396 = t363 * t380;
t379 = sin(qJ(1));
t383 = cos(qJ(1));
t360 = t383 * g(2) + t379 * g(3);
t357 = qJDD(1) * pkin(1) + t360;
t359 = t379 * g(2) - t383 * g(3);
t384 = qJD(1) ^ 2;
t358 = -t384 * pkin(1) + t359;
t374 = sin(pkin(9));
t375 = cos(pkin(9));
t342 = t375 * t357 - t374 * t358;
t340 = qJDD(1) * pkin(2) + t342;
t343 = t374 * t357 + t375 * t358;
t341 = -t384 * pkin(2) + t343;
t378 = sin(qJ(3));
t382 = cos(qJ(3));
t335 = t382 * t340 - t378 * t341;
t371 = qJDD(1) + qJDD(3);
t333 = t371 * pkin(3) + t335;
t336 = t378 * t340 + t382 * t341;
t370 = t372 ^ 2;
t334 = -t370 * pkin(3) + t336;
t377 = sin(qJ(4));
t381 = cos(qJ(4));
t330 = t377 * t333 + t381 * t334;
t361 = t363 ^ 2;
t362 = qJDD(4) + t371;
t328 = -t361 * pkin(4) + t362 * pkin(8) + t330;
t373 = -g(1) + qJDD(2);
t325 = -t376 * t328 + t380 * t373;
t349 = (-mrSges(6,1) * t380 + mrSges(6,2) * t376) * t363;
t395 = qJD(5) * t363;
t350 = t376 * t362 + t380 * t395;
t356 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t396;
t323 = m(6) * t325 + qJDD(5) * mrSges(6,1) - t350 * mrSges(6,3) + qJD(5) * t356 - t349 * t397;
t326 = t380 * t328 + t376 * t373;
t351 = t380 * t362 - t376 * t395;
t355 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t397;
t324 = m(6) * t326 - qJDD(5) * mrSges(6,2) + t351 * mrSges(6,3) - qJD(5) * t355 + t349 * t396;
t389 = -t376 * t323 + t380 * t324;
t314 = m(5) * t330 - t361 * mrSges(5,1) - t362 * mrSges(5,2) + t389;
t329 = t381 * t333 - t377 * t334;
t327 = -t362 * pkin(4) - t361 * pkin(8) - t329;
t385 = -m(6) * t327 + t351 * mrSges(6,1) - t350 * mrSges(6,2) - t355 * t397 + t356 * t396;
t319 = m(5) * t329 + t362 * mrSges(5,1) - t361 * mrSges(5,2) + t385;
t311 = t377 * t314 + t381 * t319;
t308 = m(4) * t335 + t371 * mrSges(4,1) - t370 * mrSges(4,2) + t311;
t390 = t381 * t314 - t377 * t319;
t309 = m(4) * t336 - t370 * mrSges(4,1) - t371 * mrSges(4,2) + t390;
t303 = t382 * t308 + t378 * t309;
t301 = m(3) * t342 + qJDD(1) * mrSges(3,1) - t384 * mrSges(3,2) + t303;
t391 = -t378 * t308 + t382 * t309;
t302 = m(3) * t343 - t384 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t391;
t295 = t375 * t301 + t374 * t302;
t315 = t380 * t323 + t376 * t324;
t394 = m(5) * t373 + t315;
t392 = -t374 * t301 + t375 * t302;
t293 = m(2) * t359 - t384 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t392;
t294 = m(2) * t360 + qJDD(1) * mrSges(2,1) - t384 * mrSges(2,2) + t295;
t393 = t383 * t293 - t379 * t294;
t388 = m(4) * t373 + t394;
t387 = m(3) * t373 + t388;
t386 = -t379 * t293 - t383 * t294;
t346 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t376 + Ifges(6,4) * t380) * t363;
t345 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t376 + Ifges(6,2) * t380) * t363;
t344 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t376 + Ifges(6,6) * t380) * t363;
t317 = mrSges(6,2) * t327 - mrSges(6,3) * t325 + Ifges(6,1) * t350 + Ifges(6,4) * t351 + Ifges(6,5) * qJDD(5) - qJD(5) * t345 + t344 * t396;
t316 = -mrSges(6,1) * t327 + mrSges(6,3) * t326 + Ifges(6,4) * t350 + Ifges(6,2) * t351 + Ifges(6,6) * qJDD(5) + qJD(5) * t346 - t344 * t397;
t310 = -mrSges(5,1) * t373 - mrSges(6,1) * t325 + mrSges(6,2) * t326 + mrSges(5,3) * t330 + t361 * Ifges(5,5) - Ifges(6,5) * t350 + Ifges(5,6) * t362 - Ifges(6,6) * t351 - Ifges(6,3) * qJDD(5) - pkin(4) * t315 + (-t345 * t376 + t346 * t380) * t363;
t304 = mrSges(5,2) * t373 - mrSges(5,3) * t329 + Ifges(5,5) * t362 - t361 * Ifges(5,6) - pkin(8) * t315 - t376 * t316 + t380 * t317;
t297 = mrSges(4,2) * t373 - mrSges(4,3) * t335 + Ifges(4,5) * t371 - t370 * Ifges(4,6) - pkin(7) * t311 + t381 * t304 - t377 * t310;
t296 = -mrSges(4,1) * t373 + mrSges(4,3) * t336 + t370 * Ifges(4,5) + Ifges(4,6) * t371 - pkin(3) * t394 + pkin(7) * t390 + t377 * t304 + t381 * t310;
t291 = mrSges(3,2) * t373 - mrSges(3,3) * t342 + Ifges(3,5) * qJDD(1) - t384 * Ifges(3,6) - pkin(6) * t303 - t378 * t296 + t382 * t297;
t290 = -mrSges(3,1) * t373 + mrSges(3,3) * t343 + t384 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t388 + pkin(6) * t391 + t382 * t296 + t378 * t297;
t289 = -mrSges(2,2) * g(1) - mrSges(2,3) * t360 + Ifges(2,5) * qJDD(1) - t384 * Ifges(2,6) - qJ(2) * t295 - t374 * t290 + t375 * t291;
t288 = mrSges(2,1) * g(1) + mrSges(2,3) * t359 + t384 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t387 + qJ(2) * t392 + t375 * t290 + t374 * t291;
t1 = [(-m(1) - m(2)) * g(1) + t387; -m(1) * g(2) + t386; -m(1) * g(3) + t393; pkin(1) * t295 + pkin(2) * t303 + mrSges(3,1) * t342 - mrSges(3,2) * t343 + pkin(3) * t311 + mrSges(4,1) * t335 - mrSges(4,2) * t336 + pkin(8) * t389 + mrSges(5,1) * t329 - mrSges(5,2) * t330 + t376 * t317 + t380 * t316 + pkin(4) * t385 + mrSges(2,1) * t360 - mrSges(2,2) * t359 + Ifges(4,3) * t371 + Ifges(5,3) * t362 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t393 - t383 * t288 - t379 * t289; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t386 - t379 * t288 + t383 * t289;];
tauB = t1;
