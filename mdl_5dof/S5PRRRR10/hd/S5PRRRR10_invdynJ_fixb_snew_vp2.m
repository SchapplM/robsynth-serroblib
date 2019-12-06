% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:24:09
% EndTime: 2019-12-05 17:24:11
% DurationCPUTime: 1.94s
% Computational Cost: add. (17138->227), mult. (35187->306), div. (0->0), fcn. (27557->14), ass. (0->106)
t373 = sin(pkin(11));
t376 = cos(pkin(11));
t365 = t373 * g(1) - t376 * g(2);
t372 = -g(3) + qJDD(1);
t375 = sin(pkin(5));
t378 = cos(pkin(5));
t406 = t365 * t378 + t372 * t375;
t374 = sin(pkin(6));
t381 = sin(qJ(3));
t385 = cos(qJ(3));
t397 = qJD(2) * qJD(3);
t359 = (-qJDD(2) * t385 + t381 * t397) * t374;
t366 = -t376 * g(1) - t373 * g(2);
t382 = sin(qJ(2));
t386 = cos(qJ(2));
t337 = -t382 * t366 + t406 * t386;
t387 = qJD(2) ^ 2;
t405 = pkin(8) * t374;
t333 = qJDD(2) * pkin(2) + t387 * t405 + t337;
t338 = t386 * t366 + t406 * t382;
t334 = -t387 * pkin(2) + qJDD(2) * t405 + t338;
t377 = cos(pkin(6));
t352 = -t375 * t365 + t378 * t372;
t403 = t352 * t374;
t309 = -t381 * t334 + (t333 * t377 + t403) * t385;
t371 = t377 * qJD(2) + qJD(3);
t398 = qJD(2) * t374;
t395 = t385 * t398;
t355 = -t371 * mrSges(4,2) + mrSges(4,3) * t395;
t356 = (-t385 * mrSges(4,1) + t381 * mrSges(4,2)) * t398;
t358 = (qJDD(2) * t381 + t385 * t397) * t374;
t370 = t377 * qJDD(2) + qJDD(3);
t400 = t377 * t381;
t310 = t333 * t400 + t385 * t334 + t381 * t403;
t357 = (-t385 * pkin(3) - t381 * pkin(9)) * t398;
t369 = t371 ^ 2;
t306 = -t369 * pkin(3) + t370 * pkin(9) + t357 * t395 + t310;
t348 = t377 * t352;
t308 = t359 * pkin(3) - t358 * pkin(9) + t348 + (-t333 + (pkin(3) * t381 - pkin(9) * t385) * t371 * qJD(2)) * t374;
t380 = sin(qJ(4));
t384 = cos(qJ(4));
t303 = t384 * t306 + t380 * t308;
t396 = t381 * t398;
t350 = t384 * t371 - t380 * t396;
t351 = t380 * t371 + t384 * t396;
t336 = -t350 * pkin(4) - t351 * pkin(10);
t353 = qJDD(4) + t359;
t364 = qJD(4) - t395;
t363 = t364 ^ 2;
t300 = -t363 * pkin(4) + t353 * pkin(10) + t350 * t336 + t303;
t305 = -t370 * pkin(3) - t369 * pkin(9) + t357 * t396 - t309;
t328 = -t351 * qJD(4) - t380 * t358 + t384 * t370;
t329 = t350 * qJD(4) + t384 * t358 + t380 * t370;
t301 = (-t350 * t364 - t329) * pkin(10) + (t351 * t364 - t328) * pkin(4) + t305;
t379 = sin(qJ(5));
t383 = cos(qJ(5));
t297 = -t379 * t300 + t383 * t301;
t339 = -t379 * t351 + t383 * t364;
t313 = t339 * qJD(5) + t383 * t329 + t379 * t353;
t340 = t383 * t351 + t379 * t364;
t318 = -t339 * mrSges(6,1) + t340 * mrSges(6,2);
t349 = qJD(5) - t350;
t320 = -t349 * mrSges(6,2) + t339 * mrSges(6,3);
t326 = qJDD(5) - t328;
t294 = m(6) * t297 + t326 * mrSges(6,1) - t313 * mrSges(6,3) - t340 * t318 + t349 * t320;
t298 = t383 * t300 + t379 * t301;
t312 = -t340 * qJD(5) - t379 * t329 + t383 * t353;
t321 = t349 * mrSges(6,1) - t340 * mrSges(6,3);
t295 = m(6) * t298 - t326 * mrSges(6,2) + t312 * mrSges(6,3) + t339 * t318 - t349 * t321;
t287 = t383 * t294 + t379 * t295;
t341 = -t364 * mrSges(5,2) + t350 * mrSges(5,3);
t342 = t364 * mrSges(5,1) - t351 * mrSges(5,3);
t390 = -m(5) * t305 + t328 * mrSges(5,1) - t329 * mrSges(5,2) + t350 * t341 - t351 * t342 - t287;
t283 = m(4) * t309 + t370 * mrSges(4,1) - t358 * mrSges(4,3) + t371 * t355 - t356 * t396 + t390;
t404 = t283 * t385;
t354 = t371 * mrSges(4,1) - mrSges(4,3) * t396;
t288 = -t379 * t294 + t383 * t295;
t335 = -t350 * mrSges(5,1) + t351 * mrSges(5,2);
t286 = m(5) * t303 - t353 * mrSges(5,2) + t328 * mrSges(5,3) + t350 * t335 - t364 * t342 + t288;
t302 = -t380 * t306 + t384 * t308;
t299 = -t353 * pkin(4) - t363 * pkin(10) + t351 * t336 - t302;
t296 = -m(6) * t299 + t312 * mrSges(6,1) - t313 * mrSges(6,2) + t339 * t320 - t340 * t321;
t292 = m(5) * t302 + t353 * mrSges(5,1) - t329 * mrSges(5,3) - t351 * t335 + t364 * t341 + t296;
t393 = t384 * t286 - t380 * t292;
t279 = m(4) * t310 - t370 * mrSges(4,2) - t359 * mrSges(4,3) - t371 * t354 + t356 * t395 + t393;
t399 = t279 * t400 + t377 * t404;
t281 = t380 * t286 + t384 * t292;
t394 = t385 * t279 - t283 * t381;
t315 = Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t349;
t316 = Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t349;
t389 = mrSges(6,1) * t297 - mrSges(6,2) * t298 + Ifges(6,5) * t313 + Ifges(6,6) * t312 + Ifges(6,3) * t326 + t340 * t315 - t339 * t316;
t314 = Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t349;
t289 = -mrSges(6,1) * t299 + mrSges(6,3) * t298 + Ifges(6,4) * t313 + Ifges(6,2) * t312 + Ifges(6,6) * t326 - t340 * t314 + t349 * t316;
t290 = mrSges(6,2) * t299 - mrSges(6,3) * t297 + Ifges(6,1) * t313 + Ifges(6,4) * t312 + Ifges(6,5) * t326 + t339 * t314 - t349 * t315;
t323 = Ifges(5,4) * t351 + Ifges(5,2) * t350 + Ifges(5,6) * t364;
t324 = Ifges(5,1) * t351 + Ifges(5,4) * t350 + Ifges(5,5) * t364;
t388 = mrSges(5,1) * t302 - mrSges(5,2) * t303 + Ifges(5,5) * t329 + Ifges(5,6) * t328 + Ifges(5,3) * t353 + pkin(4) * t296 + pkin(10) * t288 + t383 * t289 + t379 * t290 + t351 * t323 - t350 * t324;
t346 = Ifges(4,5) * t371 + (t381 * Ifges(4,1) + t385 * Ifges(4,4)) * t398;
t345 = Ifges(4,6) * t371 + (t381 * Ifges(4,4) + t385 * Ifges(4,2)) * t398;
t322 = Ifges(5,5) * t351 + Ifges(5,6) * t350 + Ifges(5,3) * t364;
t319 = -t374 * t333 + t348;
t280 = m(4) * t319 + t359 * mrSges(4,1) + t358 * mrSges(4,2) + (t354 * t381 - t355 * t385) * t398 + t281;
t276 = -mrSges(5,1) * t305 + mrSges(5,3) * t303 + Ifges(5,4) * t329 + Ifges(5,2) * t328 + Ifges(5,6) * t353 - pkin(4) * t287 - t351 * t322 + t364 * t324 - t389;
t275 = mrSges(5,2) * t305 - mrSges(5,3) * t302 + Ifges(5,1) * t329 + Ifges(5,4) * t328 + Ifges(5,5) * t353 - pkin(10) * t287 - t289 * t379 + t290 * t383 + t322 * t350 - t323 * t364;
t274 = Ifges(4,5) * t358 - Ifges(4,6) * t359 + Ifges(4,3) * t370 + mrSges(4,1) * t309 - mrSges(4,2) * t310 + t380 * t275 + t384 * t276 + pkin(3) * t390 + pkin(9) * t393 + (t345 * t381 - t346 * t385) * t398;
t1 = [m(2) * t372 + t378 * (m(3) * t352 + t377 * t280 + (t279 * t381 + t404) * t374) + (t382 * (m(3) * t338 - mrSges(3,1) * t387 - qJDD(2) * mrSges(3,2) + t394) + t386 * (m(3) * t337 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t387 - t280 * t374 + t399)) * t375; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t337 - mrSges(3,2) * t338 + t377 * t274 + pkin(2) * t399 + (t381 * (mrSges(4,2) * t319 - mrSges(4,3) * t309 + Ifges(4,1) * t358 - Ifges(4,4) * t359 + Ifges(4,5) * t370 - pkin(9) * t281 + t275 * t384 - t276 * t380 - t345 * t371) + t385 * (-mrSges(4,1) * t319 + mrSges(4,3) * t310 + Ifges(4,4) * t358 - Ifges(4,2) * t359 + Ifges(4,6) * t370 - pkin(3) * t281 + t371 * t346 - t388) - pkin(2) * t280 + pkin(8) * t394) * t374; t274; t388; t389;];
tauJ = t1;
