% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:39
% EndTime: 2019-12-05 17:19:42
% DurationCPUTime: 1.26s
% Computational Cost: add. (9643->217), mult. (18559->280), div. (0->0), fcn. (12877->12), ass. (0->95)
t357 = sin(pkin(10));
t359 = cos(pkin(10));
t348 = t357 * g(1) - t359 * g(2);
t356 = -g(3) + qJDD(1);
t358 = sin(pkin(5));
t360 = cos(pkin(5));
t386 = t348 * t360 + t356 * t358;
t349 = -t359 * g(1) - t357 * g(2);
t364 = sin(qJ(2));
t368 = cos(qJ(2));
t312 = -t364 * t349 + t386 * t368;
t313 = t368 * t349 + t386 * t364;
t370 = qJD(2) ^ 2;
t309 = -t370 * pkin(2) + qJDD(2) * pkin(7) + t313;
t328 = -t358 * t348 + t360 * t356;
t363 = sin(qJ(3));
t367 = cos(qJ(3));
t304 = t367 * t309 + t363 * t328;
t345 = (-t367 * pkin(3) - t363 * pkin(8)) * qJD(2);
t369 = qJD(3) ^ 2;
t382 = t367 * qJD(2);
t295 = -t369 * pkin(3) + qJDD(3) * pkin(8) + t345 * t382 + t304;
t308 = -qJDD(2) * pkin(2) - t370 * pkin(7) - t312;
t381 = qJD(2) * qJD(3);
t380 = t367 * t381;
t346 = t363 * qJDD(2) + t380;
t355 = t363 * t381;
t347 = t367 * qJDD(2) - t355;
t298 = (-t346 - t380) * pkin(8) + (-t347 + t355) * pkin(3) + t308;
t362 = sin(qJ(4));
t366 = cos(qJ(4));
t284 = -t362 * t295 + t366 * t298;
t383 = t363 * qJD(2);
t342 = t366 * qJD(3) - t362 * t383;
t320 = t342 * qJD(4) + t362 * qJDD(3) + t366 * t346;
t339 = qJDD(4) - t347;
t343 = t362 * qJD(3) + t366 * t383;
t354 = qJD(4) - t382;
t282 = (t342 * t354 - t320) * pkin(9) + (t342 * t343 + t339) * pkin(4) + t284;
t285 = t366 * t295 + t362 * t298;
t319 = -t343 * qJD(4) + t366 * qJDD(3) - t362 * t346;
t327 = t354 * pkin(4) - t343 * pkin(9);
t338 = t342 ^ 2;
t283 = -t338 * pkin(4) + t319 * pkin(9) - t354 * t327 + t285;
t361 = sin(qJ(5));
t365 = cos(qJ(5));
t280 = t365 * t282 - t361 * t283;
t321 = t365 * t342 - t361 * t343;
t292 = t321 * qJD(5) + t361 * t319 + t365 * t320;
t322 = t361 * t342 + t365 * t343;
t305 = -t321 * mrSges(6,1) + t322 * mrSges(6,2);
t353 = qJD(5) + t354;
t310 = -t353 * mrSges(6,2) + t321 * mrSges(6,3);
t335 = qJDD(5) + t339;
t277 = m(6) * t280 + t335 * mrSges(6,1) - t292 * mrSges(6,3) - t322 * t305 + t353 * t310;
t281 = t361 * t282 + t365 * t283;
t291 = -t322 * qJD(5) + t365 * t319 - t361 * t320;
t311 = t353 * mrSges(6,1) - t322 * mrSges(6,3);
t278 = m(6) * t281 - t335 * mrSges(6,2) + t291 * mrSges(6,3) + t321 * t305 - t353 * t311;
t270 = t365 * t277 + t361 * t278;
t344 = (-t367 * mrSges(4,1) + t363 * mrSges(4,2)) * qJD(2);
t350 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t383;
t323 = -t342 * mrSges(5,1) + t343 * mrSges(5,2);
t325 = -t354 * mrSges(5,2) + t342 * mrSges(5,3);
t268 = m(5) * t284 + t339 * mrSges(5,1) - t320 * mrSges(5,3) - t343 * t323 + t354 * t325 + t270;
t326 = t354 * mrSges(5,1) - t343 * mrSges(5,3);
t377 = -t361 * t277 + t365 * t278;
t269 = m(5) * t285 - t339 * mrSges(5,2) + t319 * mrSges(5,3) + t342 * t323 - t354 * t326 + t377;
t378 = -t362 * t268 + t366 * t269;
t265 = m(4) * t304 - qJDD(3) * mrSges(4,2) + t347 * mrSges(4,3) - qJD(3) * t350 + t344 * t382 + t378;
t303 = -t363 * t309 + t367 * t328;
t351 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t382;
t294 = -qJDD(3) * pkin(3) - t369 * pkin(8) + t345 * t383 - t303;
t286 = -t319 * pkin(4) - t338 * pkin(9) + t343 * t327 + t294;
t375 = m(6) * t286 - t291 * mrSges(6,1) + t292 * mrSges(6,2) - t321 * t310 + t322 * t311;
t372 = -m(5) * t294 + t319 * mrSges(5,1) - t320 * mrSges(5,2) + t342 * t325 - t343 * t326 - t375;
t273 = m(4) * t303 + qJDD(3) * mrSges(4,1) - t346 * mrSges(4,3) + qJD(3) * t351 - t344 * t383 + t372;
t379 = t367 * t265 - t363 * t273;
t266 = t366 * t268 + t362 * t269;
t300 = Ifges(6,4) * t322 + Ifges(6,2) * t321 + Ifges(6,6) * t353;
t301 = Ifges(6,1) * t322 + Ifges(6,4) * t321 + Ifges(6,5) * t353;
t374 = -mrSges(6,1) * t280 + mrSges(6,2) * t281 - Ifges(6,5) * t292 - Ifges(6,6) * t291 - Ifges(6,3) * t335 - t322 * t300 + t321 * t301;
t373 = -m(4) * t308 + t347 * mrSges(4,1) - t346 * mrSges(4,2) - t350 * t383 + t351 * t382 - t266;
t315 = Ifges(5,4) * t343 + Ifges(5,2) * t342 + Ifges(5,6) * t354;
t316 = Ifges(5,1) * t343 + Ifges(5,4) * t342 + Ifges(5,5) * t354;
t371 = mrSges(5,1) * t284 - mrSges(5,2) * t285 + Ifges(5,5) * t320 + Ifges(5,6) * t319 + Ifges(5,3) * t339 + pkin(4) * t270 + t343 * t315 - t342 * t316 - t374;
t334 = Ifges(4,5) * qJD(3) + (t363 * Ifges(4,1) + t367 * Ifges(4,4)) * qJD(2);
t333 = Ifges(4,6) * qJD(3) + (t363 * Ifges(4,4) + t367 * Ifges(4,2)) * qJD(2);
t314 = Ifges(5,5) * t343 + Ifges(5,6) * t342 + Ifges(5,3) * t354;
t299 = Ifges(6,5) * t322 + Ifges(6,6) * t321 + Ifges(6,3) * t353;
t272 = mrSges(6,2) * t286 - mrSges(6,3) * t280 + Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t335 + t321 * t299 - t353 * t300;
t271 = -mrSges(6,1) * t286 + mrSges(6,3) * t281 + Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t335 - t322 * t299 + t353 * t301;
t263 = mrSges(5,2) * t294 - mrSges(5,3) * t284 + Ifges(5,1) * t320 + Ifges(5,4) * t319 + Ifges(5,5) * t339 - pkin(9) * t270 - t361 * t271 + t365 * t272 + t342 * t314 - t354 * t315;
t262 = -mrSges(5,1) * t294 + mrSges(5,3) * t285 + Ifges(5,4) * t320 + Ifges(5,2) * t319 + Ifges(5,6) * t339 - pkin(4) * t375 + pkin(9) * t377 + t365 * t271 + t361 * t272 - t343 * t314 + t354 * t316;
t1 = [m(2) * t356 + t360 * (m(3) * t328 + t363 * t265 + t367 * t273) + (t364 * (m(3) * t313 - t370 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t379) + t368 * (m(3) * t312 + qJDD(2) * mrSges(3,1) - t370 * mrSges(3,2) + t373)) * t358; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t312 - mrSges(3,2) * t313 + t363 * (mrSges(4,2) * t308 - mrSges(4,3) * t303 + Ifges(4,1) * t346 + Ifges(4,4) * t347 + Ifges(4,5) * qJDD(3) - pkin(8) * t266 - qJD(3) * t333 - t362 * t262 + t366 * t263) + t367 * (-mrSges(4,1) * t308 + mrSges(4,3) * t304 + Ifges(4,4) * t346 + Ifges(4,2) * t347 + Ifges(4,6) * qJDD(3) - pkin(3) * t266 + qJD(3) * t334 - t371) + pkin(2) * t373 + pkin(7) * t379; Ifges(4,5) * t346 + Ifges(4,6) * t347 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t303 - mrSges(4,2) * t304 + t362 * t263 + t366 * t262 + pkin(3) * t372 + pkin(8) * t378 + (t363 * t333 - t367 * t334) * qJD(2); t371; -t374;];
tauJ = t1;
