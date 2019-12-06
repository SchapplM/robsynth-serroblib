% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:29
% EndTime: 2019-12-05 18:55:31
% DurationCPUTime: 1.80s
% Computational Cost: add. (16761->235), mult. (33966->302), div. (0->0), fcn. (25533->10), ass. (0->93)
t392 = cos(qJ(1));
t376 = sin(qJ(3));
t377 = sin(qJ(2));
t381 = cos(qJ(3));
t382 = cos(qJ(2));
t360 = (t382 * t376 + t377 * t381) * qJD(1);
t391 = qJD(1) * qJD(2);
t363 = t377 * qJDD(1) + t382 * t391;
t390 = t377 * t391;
t364 = t382 * qJDD(1) - t390;
t332 = -t360 * qJD(3) - t376 * t363 + t381 * t364;
t359 = (-t376 * t377 + t381 * t382) * qJD(1);
t333 = t359 * qJD(3) + t381 * t363 + t376 * t364;
t378 = sin(qJ(1));
t365 = t378 * g(1) - t392 * g(2);
t342 = -t365 + (-t364 + t390) * pkin(1);
t373 = qJD(2) + qJD(3);
t301 = (-t359 * t373 - t333) * pkin(5) + (t360 * t373 - t332) * pkin(2) + t342;
t366 = -g(1) * t392 - t378 * g(2);
t353 = -t382 * g(3) - t377 * t366;
t383 = qJD(1) ^ 2;
t345 = (t377 * t382 * t383 + qJDD(2)) * pkin(1) + t353;
t354 = -t377 * g(3) + t382 * t366;
t348 = (-t382 ^ 2 * t383 - qJD(2) ^ 2) * pkin(1) + t354;
t321 = t376 * t345 + t381 * t348;
t340 = -t359 * pkin(2) - t360 * pkin(5);
t371 = t373 ^ 2;
t372 = qJDD(2) + qJDD(3);
t309 = -t371 * pkin(2) + t372 * pkin(5) + t359 * t340 + t321;
t375 = sin(qJ(4));
t380 = cos(qJ(4));
t290 = t380 * t301 - t375 * t309;
t331 = qJDD(4) - t332;
t346 = -t375 * t360 + t380 * t373;
t347 = t380 * t360 + t375 * t373;
t288 = (t346 * t347 + t331) * pkin(3) + t290;
t291 = t375 * t301 + t380 * t309;
t355 = qJD(4) - t359;
t289 = (-t346 ^ 2 - t355 ^ 2) * pkin(3) + t291;
t374 = sin(qJ(5));
t379 = cos(qJ(5));
t286 = t379 * t288 - t374 * t289;
t311 = -t347 * qJD(4) - t375 * t333 + t380 * t372;
t312 = t346 * qJD(4) + t380 * t333 + t375 * t372;
t322 = t379 * t346 - t374 * t347;
t297 = t322 * qJD(5) + t374 * t311 + t379 * t312;
t323 = t374 * t346 + t379 * t347;
t307 = -t322 * mrSges(6,1) + t323 * mrSges(6,2);
t351 = qJD(5) + t355;
t313 = -t351 * mrSges(6,2) + t322 * mrSges(6,3);
t327 = qJDD(5) + t331;
t283 = m(6) * t286 + t327 * mrSges(6,1) - t297 * mrSges(6,3) - t323 * t307 + t351 * t313;
t287 = t374 * t288 + t379 * t289;
t296 = -t323 * qJD(5) + t379 * t311 - t374 * t312;
t314 = t351 * mrSges(6,1) - t323 * mrSges(6,3);
t284 = m(6) * t287 - t327 * mrSges(6,2) + t296 * mrSges(6,3) + t322 * t307 - t351 * t314;
t277 = t379 * t283 + t374 * t284;
t324 = -t346 * mrSges(5,1) + t347 * mrSges(5,2);
t334 = -t355 * mrSges(5,2) + t346 * mrSges(5,3);
t275 = m(5) * t290 + t331 * mrSges(5,1) - t312 * mrSges(5,3) - t347 * t324 + t355 * t334 + t277;
t335 = t355 * mrSges(5,1) - t347 * mrSges(5,3);
t276 = m(5) * t291 - t331 * mrSges(5,2) + t311 * mrSges(5,3) - t374 * t283 + t379 * t284 + t346 * t324 - t355 * t335;
t267 = t380 * t275 + t375 * t276;
t389 = -t375 * t275 + t380 * t276;
t320 = t381 * t345 - t376 * t348;
t308 = -t372 * pkin(2) - t371 * pkin(5) + t360 * t340 - t320;
t292 = (t347 * t355 - t311) * pkin(3) + t308;
t388 = m(6) * t292 - t296 * mrSges(6,1) + t297 * mrSges(6,2) - t322 * t313 + t323 * t314;
t303 = Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t351;
t304 = Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t351;
t387 = -mrSges(6,1) * t286 + mrSges(6,2) * t287 - Ifges(6,5) * t297 - Ifges(6,6) * t296 - Ifges(6,3) * t327 - t323 * t303 + t322 * t304;
t302 = Ifges(6,5) * t323 + Ifges(6,6) * t322 + Ifges(6,3) * t351;
t278 = -mrSges(6,1) * t292 + mrSges(6,3) * t287 + Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t327 - t323 * t302 + t351 * t304;
t279 = mrSges(6,2) * t292 - mrSges(6,3) * t286 + Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t327 + t322 * t302 - t351 * t303;
t315 = Ifges(5,5) * t347 + Ifges(5,6) * t346 + Ifges(5,3) * t355;
t317 = Ifges(5,1) * t347 + Ifges(5,4) * t346 + Ifges(5,5) * t355;
t269 = -mrSges(5,1) * t308 + mrSges(5,3) * t291 + Ifges(5,4) * t312 + Ifges(5,2) * t311 + Ifges(5,6) * t331 - pkin(3) * t388 + t379 * t278 + t374 * t279 - t347 * t315 + t355 * t317;
t316 = Ifges(5,4) * t347 + Ifges(5,2) * t346 + Ifges(5,6) * t355;
t271 = mrSges(5,2) * t308 - mrSges(5,3) * t290 + Ifges(5,1) * t312 + Ifges(5,4) * t311 + Ifges(5,5) * t331 - t374 * t278 + t379 * t279 + t346 * t315 - t355 * t316;
t337 = Ifges(4,4) * t360 + Ifges(4,2) * t359 + Ifges(4,6) * t373;
t338 = Ifges(4,1) * t360 + Ifges(4,4) * t359 + Ifges(4,5) * t373;
t385 = -m(5) * t308 + t311 * mrSges(5,1) - t312 * mrSges(5,2) + t346 * t334 - t347 * t335 - t388;
t386 = mrSges(4,1) * t320 - mrSges(4,2) * t321 + Ifges(4,5) * t333 + Ifges(4,6) * t332 + Ifges(4,3) * t372 + pkin(2) * t385 + pkin(5) * t389 + t380 * t269 + t375 * t271 + t360 * t337 - t359 * t338;
t384 = mrSges(5,1) * t290 - mrSges(5,2) * t291 + Ifges(5,5) * t312 + Ifges(5,6) * t311 + Ifges(5,3) * t331 + pkin(3) * t277 + t347 * t316 - t346 * t317 - t387;
t358 = Ifges(3,5) * qJD(2) + (t377 * Ifges(3,1) + t382 * Ifges(3,4)) * qJD(1);
t357 = Ifges(3,6) * qJD(2) + (t377 * Ifges(3,4) + t382 * Ifges(3,2)) * qJD(1);
t350 = t373 * mrSges(4,1) - t360 * mrSges(4,3);
t349 = -t373 * mrSges(4,2) + t359 * mrSges(4,3);
t339 = -t359 * mrSges(4,1) + t360 * mrSges(4,2);
t336 = Ifges(4,5) * t360 + Ifges(4,6) * t359 + Ifges(4,3) * t373;
t265 = -mrSges(4,1) * t342 + mrSges(4,3) * t321 + Ifges(4,4) * t333 + Ifges(4,2) * t332 + Ifges(4,6) * t372 - pkin(2) * t267 - t360 * t336 + t373 * t338 - t384;
t264 = mrSges(4,2) * t342 - mrSges(4,3) * t320 + Ifges(4,1) * t333 + Ifges(4,4) * t332 + Ifges(4,5) * t372 - pkin(5) * t267 - t375 * t269 + t380 * t271 + t359 * t336 - t373 * t337;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t365 - mrSges(2,2) * t366 + t377 * (-mrSges(3,2) * t365 - mrSges(3,3) * t353 + Ifges(3,1) * t363 + Ifges(3,4) * t364 + Ifges(3,5) * qJDD(2) - qJD(2) * t357 + t381 * t264 - t376 * t265) + t382 * (Ifges(3,4) * t363 + Ifges(3,2) * t364 + Ifges(3,6) * qJDD(2) + qJD(2) * t358 + mrSges(3,1) * t365 + mrSges(3,3) * t354 + t376 * t264 + t381 * t265 - pkin(1) * (m(4) * t342 - t332 * mrSges(4,1) + t333 * mrSges(4,2) - t359 * t349 + t360 * t350 + t267)); t386 + Ifges(3,3) * qJDD(2) + (t377 * t357 - t382 * t358) * qJD(1) + pkin(1) * (t376 * (m(4) * t321 - t372 * mrSges(4,2) + t332 * mrSges(4,3) + t359 * t339 - t373 * t350 + t389) + t381 * (m(4) * t320 + t372 * mrSges(4,1) - t333 * mrSges(4,3) - t360 * t339 + t373 * t349 + t385)) + Ifges(3,5) * t363 + Ifges(3,6) * t364 + mrSges(3,1) * t353 - mrSges(3,2) * t354; t386; t384; -t387;];
tauJ = t1;
