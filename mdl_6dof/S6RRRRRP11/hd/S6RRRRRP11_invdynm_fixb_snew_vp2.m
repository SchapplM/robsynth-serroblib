% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:41:57
% EndTime: 2019-05-08 06:45:25
% DurationCPUTime: 101.31s
% Computational Cost: add. (1797490->411), mult. (4433537->537), div. (0->0), fcn. (3731703->14), ass. (0->172)
t347 = cos(pkin(6));
t341 = t347 * qJD(1) + qJD(2);
t344 = sin(pkin(7));
t346 = cos(pkin(7));
t345 = sin(pkin(6));
t356 = cos(qJ(2));
t380 = qJD(1) * t356;
t374 = t345 * t380;
t325 = (t341 * t344 + t346 * t374) * pkin(10);
t351 = sin(qJ(2));
t382 = qJD(1) * t345;
t397 = pkin(10) * t344;
t329 = (-pkin(2) * t356 - t351 * t397) * t382;
t379 = qJD(1) * qJD(2);
t335 = (qJDD(1) * t351 + t356 * t379) * t345;
t340 = t347 * qJDD(1) + qJDD(2);
t352 = sin(qJ(1));
t357 = cos(qJ(1));
t338 = t352 * g(1) - t357 * g(2);
t358 = qJD(1) ^ 2;
t398 = pkin(9) * t345;
t332 = qJDD(1) * pkin(1) + t358 * t398 + t338;
t339 = -t357 * g(1) - t352 * g(2);
t333 = -t358 * pkin(1) + qJDD(1) * t398 + t339;
t386 = t347 * t356;
t372 = t332 * t386 - t351 * t333;
t381 = qJD(1) * t351;
t396 = pkin(10) * t346;
t285 = -t335 * t396 + t340 * pkin(2) + t341 * t325 + (-g(3) * t356 - t329 * t381) * t345 + t372;
t375 = t345 * t381;
t328 = t341 * pkin(2) - t375 * t396;
t336 = (qJDD(1) * t356 - t351 * t379) * t345;
t369 = t336 * t346 + t340 * t344;
t387 = t347 * t351;
t383 = t332 * t387 + t356 * t333;
t286 = -t341 * t328 + (-g(3) * t351 + t329 * t380) * t345 + t369 * pkin(10) + t383;
t395 = t347 * g(3);
t291 = -t335 * t397 - t336 * pkin(2) - t395 + (-t332 + (-t325 * t356 + t328 * t351) * qJD(1)) * t345;
t350 = sin(qJ(3));
t355 = cos(qJ(3));
t246 = -t350 * t286 + (t285 * t346 + t291 * t344) * t355;
t388 = t346 * t356;
t392 = t344 * t355;
t315 = (-t350 * t351 + t355 * t388) * t382 + t341 * t392;
t303 = t315 * qJD(3) + t355 * t335 + t369 * t350;
t393 = t344 * t350;
t316 = t341 * t393 + (t350 * t388 + t351 * t355) * t382;
t326 = t346 * t341 - t344 * t374 + qJD(3);
t349 = sin(qJ(4));
t354 = cos(qJ(4));
t307 = -t349 * t316 + t354 * t326;
t317 = -t344 * t336 + t346 * t340 + qJDD(3);
t271 = t307 * qJD(4) + t354 * t303 + t349 * t317;
t308 = t354 * t316 + t349 * t326;
t314 = qJD(4) - t315;
t348 = sin(qJ(5));
t353 = cos(qJ(5));
t294 = -t348 * t308 + t353 * t314;
t302 = -t316 * qJD(3) - t350 * t335 + t369 * t355;
t301 = qJDD(4) - t302;
t245 = t294 * qJD(5) + t353 * t271 + t348 * t301;
t295 = t353 * t308 + t348 * t314;
t262 = -t294 * mrSges(7,1) + t295 * mrSges(7,2);
t389 = t346 * t350;
t247 = t285 * t389 + t355 * t286 + t291 * t393;
t305 = -t315 * pkin(3) - t316 * pkin(11);
t324 = t326 ^ 2;
t232 = -t324 * pkin(3) + t317 * pkin(11) + t315 * t305 + t247;
t260 = -t344 * t285 + t346 * t291;
t234 = (-t315 * t326 - t303) * pkin(11) + (t316 * t326 - t302) * pkin(3) + t260;
t228 = t354 * t232 + t349 * t234;
t288 = -t307 * pkin(4) - t308 * pkin(12);
t313 = t314 ^ 2;
t223 = -t313 * pkin(4) + t301 * pkin(12) + t307 * t288 + t228;
t231 = -t317 * pkin(3) - t324 * pkin(11) + t316 * t305 - t246;
t270 = -t308 * qJD(4) - t349 * t303 + t354 * t317;
t226 = (-t307 * t314 - t271) * pkin(12) + (t308 * t314 - t270) * pkin(4) + t231;
t217 = -t348 * t223 + t353 * t226;
t269 = qJDD(5) - t270;
t306 = qJD(5) - t307;
t213 = -0.2e1 * qJD(6) * t295 + (t294 * t306 - t245) * qJ(6) + (t294 * t295 + t269) * pkin(5) + t217;
t273 = -t306 * mrSges(7,2) + t294 * mrSges(7,3);
t377 = m(7) * t213 + t269 * mrSges(7,1) + t306 * t273;
t210 = -t245 * mrSges(7,3) - t295 * t262 + t377;
t218 = t353 * t223 + t348 * t226;
t244 = -t295 * qJD(5) - t348 * t271 + t353 * t301;
t255 = Ifges(6,4) * t295 + Ifges(6,2) * t294 + Ifges(6,6) * t306;
t256 = Ifges(7,1) * t295 + Ifges(7,4) * t294 + Ifges(7,5) * t306;
t257 = Ifges(6,1) * t295 + Ifges(6,4) * t294 + Ifges(6,5) * t306;
t275 = t306 * pkin(5) - t295 * qJ(6);
t293 = t294 ^ 2;
t216 = -t293 * pkin(5) + t244 * qJ(6) + 0.2e1 * qJD(6) * t294 - t306 * t275 + t218;
t254 = Ifges(7,4) * t295 + Ifges(7,2) * t294 + Ifges(7,6) * t306;
t364 = -mrSges(7,1) * t213 + mrSges(7,2) * t216 - Ifges(7,5) * t245 - Ifges(7,6) * t244 - Ifges(7,3) * t269 - t295 * t254;
t399 = mrSges(6,1) * t217 - mrSges(6,2) * t218 + Ifges(6,5) * t245 + Ifges(6,6) * t244 + Ifges(6,3) * t269 + pkin(5) * t210 + t295 * t255 - (t257 + t256) * t294 - t364;
t394 = -mrSges(6,2) - mrSges(7,2);
t391 = t345 * t351;
t390 = t345 * t356;
t263 = -t294 * mrSges(6,1) + t295 * mrSges(6,2);
t274 = -t306 * mrSges(6,2) + t294 * mrSges(6,3);
t204 = m(6) * t217 + t269 * mrSges(6,1) + t306 * t274 + (-t262 - t263) * t295 + (-mrSges(6,3) - mrSges(7,3)) * t245 + t377;
t376 = m(7) * t216 + t244 * mrSges(7,3) + t294 * t262;
t276 = t306 * mrSges(7,1) - t295 * mrSges(7,3);
t384 = -t306 * mrSges(6,1) + t295 * mrSges(6,3) - t276;
t206 = m(6) * t218 + t244 * mrSges(6,3) + t294 * t263 + t394 * t269 + t384 * t306 + t376;
t203 = -t348 * t204 + t353 * t206;
t287 = -t307 * mrSges(5,1) + t308 * mrSges(5,2);
t297 = t314 * mrSges(5,1) - t308 * mrSges(5,3);
t200 = m(5) * t228 - t301 * mrSges(5,2) + t270 * mrSges(5,3) + t307 * t287 - t314 * t297 + t203;
t227 = -t349 * t232 + t354 * t234;
t222 = -t301 * pkin(4) - t313 * pkin(12) + t308 * t288 - t227;
t220 = -t244 * pkin(5) - t293 * qJ(6) + t295 * t275 + qJDD(6) + t222;
t371 = -m(7) * t220 + t244 * mrSges(7,1) + t294 * t273;
t209 = -m(6) * t222 + t244 * mrSges(6,1) + t394 * t245 + t294 * t274 + t384 * t295 + t371;
t296 = -t314 * mrSges(5,2) + t307 * mrSges(5,3);
t208 = m(5) * t227 + t301 * mrSges(5,1) - t271 * mrSges(5,3) - t308 * t287 + t314 * t296 + t209;
t193 = t349 * t200 + t354 * t208;
t304 = -t315 * mrSges(4,1) + t316 * mrSges(4,2);
t310 = t326 * mrSges(4,1) - t316 * mrSges(4,3);
t373 = t354 * t200 - t349 * t208;
t190 = m(4) * t247 - t317 * mrSges(4,2) + t302 * mrSges(4,3) + t315 * t304 - t326 * t310 + t373;
t309 = -t326 * mrSges(4,2) + t315 * mrSges(4,3);
t192 = m(4) * t260 - t302 * mrSges(4,1) + t303 * mrSges(4,2) - t315 * t309 + t316 * t310 + t193;
t202 = t353 * t204 + t348 * t206;
t361 = -m(5) * t231 + t270 * mrSges(5,1) - t271 * mrSges(5,2) + t307 * t296 - t308 * t297 - t202;
t197 = m(4) * t246 + t317 * mrSges(4,1) - t303 * mrSges(4,3) - t316 * t304 + t326 * t309 + t361;
t179 = t190 * t393 + t346 * t192 + t197 * t392;
t180 = t346 * t355 * t197 + t190 * t389 - t344 * t192;
t311 = -g(3) * t390 + t372;
t331 = -t341 * mrSges(3,2) + mrSges(3,3) * t374;
t334 = (-mrSges(3,1) * t356 + mrSges(3,2) * t351) * t382;
t177 = m(3) * t311 + t340 * mrSges(3,1) - t335 * mrSges(3,3) + t341 * t331 - t334 * t375 + t180;
t185 = t355 * t190 - t350 * t197;
t312 = -g(3) * t391 + t383;
t330 = t341 * mrSges(3,1) - mrSges(3,3) * t375;
t184 = m(3) * t312 - t340 * mrSges(3,2) + t336 * mrSges(3,3) - t341 * t330 + t334 * t374 + t185;
t174 = -t351 * t177 + t356 * t184;
t321 = -t345 * t332 - t395;
t178 = m(3) * t321 - t336 * mrSges(3,1) + t335 * mrSges(3,2) + (t330 * t351 - t331 * t356) * t382 + t179;
t169 = t177 * t386 - t345 * t178 + t184 * t387;
t252 = Ifges(7,5) * t295 + Ifges(7,6) * t294 + Ifges(7,3) * t306;
t253 = Ifges(6,5) * t295 + Ifges(6,6) * t294 + Ifges(6,3) * t306;
t365 = -mrSges(7,1) * t220 + mrSges(7,3) * t216 + Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t269 + t306 * t256;
t194 = Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t269 + t306 * t257 - mrSges(6,1) * t222 + mrSges(6,3) * t218 - pkin(5) * (t245 * mrSges(7,2) - t371) + qJ(6) * (-t269 * mrSges(7,2) - t306 * t276 + t376) + (-pkin(5) * t276 - t252 - t253) * t295 + t365;
t363 = mrSges(7,2) * t220 - mrSges(7,3) * t213 + Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t269 + t294 * t252;
t201 = mrSges(6,2) * t222 - mrSges(6,3) * t217 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t269 - qJ(6) * t210 + t294 * t253 + (-t254 - t255) * t306 + t363;
t278 = Ifges(5,5) * t308 + Ifges(5,6) * t307 + Ifges(5,3) * t314;
t279 = Ifges(5,4) * t308 + Ifges(5,2) * t307 + Ifges(5,6) * t314;
t181 = mrSges(5,2) * t231 - mrSges(5,3) * t227 + Ifges(5,1) * t271 + Ifges(5,4) * t270 + Ifges(5,5) * t301 - pkin(12) * t202 - t348 * t194 + t353 * t201 + t307 * t278 - t314 * t279;
t280 = Ifges(5,1) * t308 + Ifges(5,4) * t307 + Ifges(5,5) * t314;
t186 = -mrSges(5,1) * t231 + mrSges(5,3) * t228 + Ifges(5,4) * t271 + Ifges(5,2) * t270 + Ifges(5,6) * t301 - pkin(4) * t202 - t308 * t278 + t314 * t280 - t399;
t298 = Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t326;
t299 = Ifges(4,4) * t316 + Ifges(4,2) * t315 + Ifges(4,6) * t326;
t171 = mrSges(4,2) * t260 - mrSges(4,3) * t246 + Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t317 - pkin(11) * t193 + t354 * t181 - t349 * t186 + t315 * t298 - t326 * t299;
t300 = Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t326;
t359 = mrSges(5,1) * t227 - mrSges(5,2) * t228 + Ifges(5,5) * t271 + Ifges(5,6) * t270 + Ifges(5,3) * t301 + pkin(4) * t209 + pkin(12) * t203 + t353 * t194 + t348 * t201 + t308 * t279 - t307 * t280;
t175 = -mrSges(4,1) * t260 + mrSges(4,3) * t247 + Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t317 - pkin(3) * t193 - t316 * t298 + t326 * t300 - t359;
t366 = pkin(10) * t185 + t171 * t350 + t175 * t355;
t170 = mrSges(4,1) * t246 - mrSges(4,2) * t247 + Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * t317 + pkin(3) * t361 + pkin(11) * t373 + t349 * t181 + t354 * t186 + t316 * t299 - t315 * t300;
t319 = Ifges(3,6) * t341 + (Ifges(3,4) * t351 + Ifges(3,2) * t356) * t382;
t320 = Ifges(3,5) * t341 + (Ifges(3,1) * t351 + Ifges(3,4) * t356) * t382;
t161 = mrSges(3,1) * t311 - mrSges(3,2) * t312 + Ifges(3,5) * t335 + Ifges(3,6) * t336 + Ifges(3,3) * t340 + pkin(2) * t180 + t346 * t170 + (t319 * t351 - t320 * t356) * t382 + t366 * t344;
t318 = Ifges(3,3) * t341 + (Ifges(3,5) * t351 + Ifges(3,6) * t356) * t382;
t163 = -mrSges(3,1) * t321 + mrSges(3,3) * t312 + Ifges(3,4) * t335 + Ifges(3,2) * t336 + Ifges(3,6) * t340 - pkin(2) * t179 - t344 * t170 - t318 * t375 + t341 * t320 + t366 * t346;
t165 = t318 * t374 + mrSges(3,2) * t321 - mrSges(3,3) * t311 + Ifges(3,1) * t335 + Ifges(3,4) * t336 + Ifges(3,5) * t340 + t355 * t171 - t350 * t175 - t341 * t319 + (-t179 * t344 - t180 * t346) * pkin(10);
t362 = mrSges(2,1) * t338 - mrSges(2,2) * t339 + Ifges(2,3) * qJDD(1) + pkin(1) * t169 + t347 * t161 + t163 * t390 + t165 * t391 + t174 * t398;
t172 = m(2) * t339 - t358 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t174;
t168 = t347 * t178 + (t177 * t356 + t184 * t351) * t345;
t166 = m(2) * t338 + qJDD(1) * mrSges(2,1) - t358 * mrSges(2,2) + t169;
t159 = -mrSges(2,2) * g(3) - mrSges(2,3) * t338 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) - t351 * t163 + t356 * t165 + (-t168 * t345 - t169 * t347) * pkin(9);
t158 = mrSges(2,1) * g(3) + mrSges(2,3) * t339 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t168 - t345 * t161 + (pkin(9) * t174 + t163 * t356 + t165 * t351) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t159 - t352 * t158 - pkin(8) * (t357 * t166 + t352 * t172), t159, t165, t171, t181, t201, -t306 * t254 + t363; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t159 + t357 * t158 + pkin(8) * (-t352 * t166 + t357 * t172), t158, t163, t175, t186, t194, -t295 * t252 + t365; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t161, t170, t359, t399, -t294 * t256 - t364;];
m_new  = t1;
