% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 20:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:25:48
% EndTime: 2019-05-07 20:26:26
% DurationCPUTime: 15.07s
% Computational Cost: add. (275575->385), mult. (550015->467), div. (0->0), fcn. (391837->10), ass. (0->147)
t344 = sin(qJ(3));
t348 = cos(qJ(3));
t349 = cos(qJ(2));
t373 = qJD(1) * t349;
t345 = sin(qJ(2));
t374 = qJD(1) * t345;
t315 = -t344 * t374 + t348 * t373;
t372 = qJD(1) * qJD(2);
t324 = qJDD(1) * t345 + t349 * t372;
t325 = qJDD(1) * t349 - t345 * t372;
t287 = qJD(3) * t315 + t324 * t348 + t325 * t344;
t316 = (t344 * t349 + t345 * t348) * qJD(1);
t339 = qJD(2) + qJD(3);
t343 = sin(qJ(4));
t382 = cos(qJ(4));
t301 = t343 * t316 - t382 * t339;
t338 = qJDD(2) + qJDD(3);
t247 = -t301 * qJD(4) + t382 * t287 + t343 * t338;
t346 = sin(qJ(1));
t350 = cos(qJ(1));
t331 = -g(1) * t350 - g(2) * t346;
t351 = qJD(1) ^ 2;
t318 = -pkin(1) * t351 + qJDD(1) * pkin(7) + t331;
t378 = t345 * t318;
t381 = pkin(2) * t351;
t272 = qJDD(2) * pkin(2) - t324 * pkin(8) - t378 + (pkin(8) * t372 + t345 * t381 - g(3)) * t349;
t304 = -g(3) * t345 + t349 * t318;
t329 = qJD(2) * pkin(2) - pkin(8) * t374;
t341 = t349 ^ 2;
t273 = pkin(8) * t325 - qJD(2) * t329 - t341 * t381 + t304;
t236 = t348 * t272 - t344 * t273;
t299 = -pkin(3) * t315 - pkin(9) * t316;
t337 = t339 ^ 2;
t363 = t338 * pkin(3) + t337 * pkin(9) - t316 * t299 + t236;
t311 = qJD(4) - t315;
t379 = t301 * t311;
t386 = (-t247 + t379) * qJ(5) - t363;
t286 = -qJD(3) * t316 - t324 * t344 + t348 * t325;
t330 = t346 * g(1) - t350 * g(2);
t366 = -qJDD(1) * pkin(1) - t330;
t288 = -t325 * pkin(2) + t329 * t374 + (-pkin(8) * t341 - pkin(7)) * t351 + t366;
t224 = (-t315 * t339 - t287) * pkin(9) + (t316 * t339 - t286) * pkin(3) + t288;
t237 = t344 * t272 + t348 * t273;
t228 = -pkin(3) * t337 + pkin(9) * t338 + t299 * t315 + t237;
t213 = t382 * t224 - t343 * t228;
t214 = t343 * t224 + t382 * t228;
t302 = t382 * t316 + t343 * t339;
t246 = t302 * qJD(4) + t343 * t287 - t382 * t338;
t252 = Ifges(6,5) * t302 + Ifges(6,6) * t311 + Ifges(6,3) * t301;
t255 = Ifges(5,4) * t302 - Ifges(5,2) * t301 + Ifges(5,6) * t311;
t257 = Ifges(5,1) * t302 - Ifges(5,4) * t301 + Ifges(5,5) * t311;
t269 = mrSges(6,1) * t301 - mrSges(6,3) * t302;
t285 = qJDD(4) - t286;
t268 = pkin(4) * t301 - qJ(5) * t302;
t310 = t311 ^ 2;
t211 = -t285 * pkin(4) - t310 * qJ(5) + t302 * t268 + qJDD(5) - t213;
t203 = (-t247 - t379) * pkin(10) + (t301 * t302 - t285) * pkin(5) + t211;
t383 = 2 * qJD(5);
t209 = -pkin(4) * t310 + t285 * qJ(5) - t301 * t268 + t311 * t383 + t214;
t293 = -pkin(5) * t311 - pkin(10) * t302;
t300 = t301 ^ 2;
t204 = -pkin(5) * t300 + pkin(10) * t246 + t293 * t311 + t209;
t342 = sin(qJ(6));
t347 = cos(qJ(6));
t200 = t203 * t347 - t204 * t342;
t262 = t301 * t347 - t302 * t342;
t220 = qJD(6) * t262 + t246 * t342 + t247 * t347;
t263 = t301 * t342 + t302 * t347;
t234 = -mrSges(7,1) * t262 + mrSges(7,2) * t263;
t309 = qJD(6) - t311;
t250 = -mrSges(7,2) * t309 + mrSges(7,3) * t262;
t280 = qJDD(6) - t285;
t196 = m(7) * t200 + mrSges(7,1) * t280 - mrSges(7,3) * t220 - t234 * t263 + t250 * t309;
t201 = t203 * t342 + t204 * t347;
t219 = -qJD(6) * t263 + t246 * t347 - t247 * t342;
t251 = mrSges(7,1) * t309 - mrSges(7,3) * t263;
t197 = m(7) * t201 - mrSges(7,2) * t280 + mrSges(7,3) * t219 + t234 * t262 - t251 * t309;
t185 = t347 * t196 + t342 * t197;
t256 = Ifges(6,1) * t302 + Ifges(6,4) * t311 + Ifges(6,5) * t301;
t230 = Ifges(7,4) * t263 + Ifges(7,2) * t262 + Ifges(7,6) * t309;
t231 = Ifges(7,1) * t263 + Ifges(7,4) * t262 + Ifges(7,5) * t309;
t364 = mrSges(7,1) * t200 - mrSges(7,2) * t201 + Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t280 + t263 * t230 - t262 * t231;
t356 = mrSges(6,1) * t211 - mrSges(6,3) * t209 - Ifges(6,4) * t247 - Ifges(6,2) * t285 - Ifges(6,6) * t246 + pkin(5) * t185 - t301 * t256 + t364;
t289 = -mrSges(6,2) * t301 + mrSges(6,3) * t311;
t361 = -m(6) * t211 + t285 * mrSges(6,1) + t311 * t289 - t185;
t186 = -t342 * t196 + t347 * t197;
t292 = -mrSges(6,1) * t311 + mrSges(6,2) * t302;
t365 = m(6) * t209 + t285 * mrSges(6,3) + t311 * t292 + t186;
t385 = (t255 - t252) * t302 + mrSges(5,1) * t213 - mrSges(5,2) * t214 + Ifges(5,5) * t247 - Ifges(5,6) * t246 + Ifges(5,3) * t285 + pkin(4) * (-t247 * mrSges(6,2) - t302 * t269 + t361) + qJ(5) * (-t246 * mrSges(6,2) - t301 * t269 + t365) + t301 * t257 - t356;
t298 = -mrSges(4,1) * t315 + mrSges(4,2) * t316;
t306 = mrSges(4,1) * t339 - mrSges(4,3) * t316;
t291 = mrSges(5,1) * t311 - mrSges(5,3) * t302;
t375 = -mrSges(5,1) * t301 - mrSges(5,2) * t302 - t269;
t380 = -mrSges(5,3) - mrSges(6,2);
t180 = m(5) * t214 - t285 * mrSges(5,2) + t380 * t246 - t311 * t291 + t375 * t301 + t365;
t290 = -mrSges(5,2) * t311 - mrSges(5,3) * t301;
t182 = m(5) * t213 + t285 * mrSges(5,1) + t380 * t247 + t311 * t290 + t375 * t302 + t361;
t369 = t382 * t180 - t182 * t343;
t174 = m(4) * t237 - mrSges(4,2) * t338 + mrSges(4,3) * t286 + t298 * t315 - t306 * t339 + t369;
t305 = -mrSges(4,2) * t339 + mrSges(4,3) * t315;
t206 = -t300 * pkin(10) + (-pkin(4) - pkin(5)) * t246 + (-pkin(4) * t311 + t293 + t383) * t302 - t386;
t202 = -m(7) * t206 + t219 * mrSges(7,1) - t220 * mrSges(7,2) + t262 * t250 - t263 * t251;
t212 = -0.2e1 * qJD(5) * t302 + (t302 * t311 + t246) * pkin(4) + t386;
t194 = m(6) * t212 + mrSges(6,1) * t246 - t247 * mrSges(6,3) + t289 * t301 - t302 * t292 + t202;
t354 = m(5) * t363 - t246 * mrSges(5,1) - mrSges(5,2) * t247 - t301 * t290 - t291 * t302 - t194;
t191 = m(4) * t236 + mrSges(4,1) * t338 - mrSges(4,3) * t287 - t298 * t316 + t305 * t339 + t354;
t169 = t344 * t174 + t348 * t191;
t303 = -t349 * g(3) - t378;
t313 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t345 + Ifges(3,2) * t349) * qJD(1);
t314 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t345 + Ifges(3,4) * t349) * qJD(1);
t229 = Ifges(7,5) * t263 + Ifges(7,6) * t262 + Ifges(7,3) * t309;
t188 = -mrSges(7,1) * t206 + mrSges(7,3) * t201 + Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t280 - t229 * t263 + t231 * t309;
t189 = mrSges(7,2) * t206 - mrSges(7,3) * t200 + Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t280 + t229 * t262 - t230 * t309;
t357 = -mrSges(6,1) * t212 + mrSges(6,2) * t209 - pkin(5) * t202 - pkin(10) * t186 - t347 * t188 - t342 * t189;
t254 = Ifges(6,4) * t302 + Ifges(6,2) * t311 + Ifges(6,6) * t301;
t377 = -Ifges(5,5) * t302 + Ifges(5,6) * t301 - Ifges(5,3) * t311 - t254;
t163 = mrSges(5,1) * t363 + mrSges(5,3) * t214 - pkin(4) * t194 + (t257 + t256) * t311 + t377 * t302 + (Ifges(5,6) - Ifges(6,6)) * t285 + (Ifges(5,4) - Ifges(6,5)) * t247 + (-Ifges(5,2) - Ifges(6,3)) * t246 + t357;
t359 = mrSges(6,2) * t211 - mrSges(6,3) * t212 + Ifges(6,1) * t247 + Ifges(6,4) * t285 + Ifges(6,5) * t246 - pkin(10) * t185 - t342 * t188 + t347 * t189 + t311 * t252;
t165 = -mrSges(5,2) * t363 - mrSges(5,3) * t213 + Ifges(5,1) * t247 - Ifges(5,4) * t246 + Ifges(5,5) * t285 - qJ(5) * t194 - t311 * t255 + t377 * t301 + t359;
t295 = Ifges(4,4) * t316 + Ifges(4,2) * t315 + Ifges(4,6) * t339;
t296 = Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t339;
t358 = -mrSges(4,1) * t236 + mrSges(4,2) * t237 - Ifges(4,5) * t287 - Ifges(4,6) * t286 - Ifges(4,3) * t338 - pkin(3) * t354 - pkin(9) * t369 - t382 * t163 - t343 * t165 - t316 * t295 + t315 * t296;
t384 = mrSges(3,1) * t303 - mrSges(3,2) * t304 + Ifges(3,5) * t324 + Ifges(3,6) * t325 + Ifges(3,3) * qJDD(2) + pkin(2) * t169 + (t313 * t345 - t314 * t349) * qJD(1) - t358;
t176 = t343 * t180 + t382 * t182;
t323 = (-mrSges(3,1) * t349 + mrSges(3,2) * t345) * qJD(1);
t328 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t373;
t167 = m(3) * t303 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t324 + qJD(2) * t328 - t323 * t374 + t169;
t327 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t374;
t370 = t348 * t174 - t191 * t344;
t168 = m(3) * t304 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t325 - qJD(2) * t327 + t323 * t373 + t370;
t371 = -t167 * t345 + t349 * t168;
t294 = Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t339;
t157 = mrSges(4,2) * t288 - mrSges(4,3) * t236 + Ifges(4,1) * t287 + Ifges(4,4) * t286 + Ifges(4,5) * t338 - pkin(9) * t176 - t343 * t163 + t382 * t165 + t315 * t294 - t339 * t295;
t158 = -mrSges(4,1) * t288 + mrSges(4,3) * t237 + Ifges(4,4) * t287 + Ifges(4,2) * t286 + Ifges(4,6) * t338 - pkin(3) * t176 - t316 * t294 + t339 * t296 - t385;
t312 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t345 + Ifges(3,6) * t349) * qJD(1);
t317 = -t351 * pkin(7) + t366;
t360 = m(4) * t288 - t286 * mrSges(4,1) + mrSges(4,2) * t287 - t315 * t305 + t306 * t316 + t176;
t153 = -mrSges(3,1) * t317 + mrSges(3,3) * t304 + Ifges(3,4) * t324 + Ifges(3,2) * t325 + Ifges(3,6) * qJDD(2) - pkin(2) * t360 + pkin(8) * t370 + qJD(2) * t314 + t344 * t157 + t348 * t158 - t312 * t374;
t155 = mrSges(3,2) * t317 - mrSges(3,3) * t303 + Ifges(3,1) * t324 + Ifges(3,4) * t325 + Ifges(3,5) * qJDD(2) - pkin(8) * t169 - qJD(2) * t313 + t157 * t348 - t158 * t344 + t312 * t373;
t355 = -m(3) * t317 + t325 * mrSges(3,1) - mrSges(3,2) * t324 - t327 * t374 + t328 * t373 - t360;
t362 = mrSges(2,1) * t330 - mrSges(2,2) * t331 + Ifges(2,3) * qJDD(1) + pkin(1) * t355 + pkin(7) * t371 + t349 * t153 + t345 * t155;
t170 = m(2) * t330 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t351 + t355;
t161 = t167 * t349 + t168 * t345;
t159 = m(2) * t331 - mrSges(2,1) * t351 - qJDD(1) * mrSges(2,2) + t371;
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t331 + t351 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t161 - t384;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t330 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t351 - pkin(7) * t161 - t153 * t345 + t155 * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t350 * t151 - t346 * t156 - pkin(6) * (t159 * t346 + t170 * t350), t151, t155, t157, t165, -t254 * t301 + t359, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t346 * t151 + t350 * t156 + pkin(6) * (t159 * t350 - t170 * t346), t156, t153, t158, t163, -t302 * t252 - t356, t188; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t384, -t358, t385, Ifges(6,5) * t247 + Ifges(6,6) * t285 + Ifges(6,3) * t246 + t302 * t254 - t311 * t256 - t357, t364;];
m_new  = t1;
