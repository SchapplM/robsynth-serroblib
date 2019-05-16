% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:01:46
% EndTime: 2019-05-06 02:03:29
% DurationCPUTime: 83.90s
% Computational Cost: add. (1328252->397), mult. (4125125->532), div. (0->0), fcn. (3514147->14), ass. (0->175)
t344 = sin(pkin(12));
t346 = sin(pkin(6));
t347 = cos(pkin(12));
t349 = cos(pkin(6));
t352 = sin(qJ(3));
t348 = cos(pkin(7));
t356 = cos(qJ(3));
t394 = t348 * t356;
t345 = sin(pkin(7));
t399 = t345 * t356;
t362 = t346 * (-t344 * t352 + t347 * t394) + t349 * t399;
t317 = t362 * qJD(1);
t395 = t348 * t352;
t400 = t345 * t352;
t364 = t349 * t400 + (t344 * t356 + t347 * t395) * t346;
t318 = t364 * qJD(1);
t306 = -t318 * qJD(3) + qJDD(1) * t362;
t397 = t346 * t348;
t329 = (t345 * t349 + t347 * t397) * qJD(1) * pkin(9);
t353 = sin(qJ(1));
t357 = cos(qJ(1));
t341 = -g(1) * t357 - g(2) * t353;
t358 = qJD(1) ^ 2;
t403 = qJ(2) * t346;
t333 = -pkin(1) * t358 + qJDD(1) * t403 + t341;
t407 = pkin(9) * t344;
t377 = -pkin(2) * t347 - t345 * t407;
t391 = qJD(1) * t346;
t404 = pkin(9) * qJDD(1);
t372 = qJD(1) * t377 * t391 + t348 * t404;
t340 = t353 * g(1) - g(2) * t357;
t332 = qJDD(1) * pkin(1) + t358 * t403 + t340;
t385 = qJD(2) * t391;
t396 = t347 * t349;
t398 = t346 * t347;
t378 = -g(3) * t398 + t332 * t396 - 0.2e1 * t344 * t385;
t287 = (pkin(2) * qJDD(1) + qJD(1) * t329) * t349 + (-t346 * t372 - t333) * t344 + t378;
t334 = (pkin(2) * t349 - t397 * t407) * qJD(1);
t401 = t344 * t349;
t386 = t332 * t401 + (t333 + 0.2e1 * t385) * t347;
t288 = (-qJD(1) * t334 + t345 * t404) * t349 + (-g(3) * t344 + t347 * t372) * t346 + t386;
t384 = -t349 * g(3) + qJDD(2);
t297 = (-t332 + t377 * qJDD(1) + (-t329 * t347 + t334 * t344) * qJD(1)) * t346 + t384;
t237 = -t352 * t288 + (t287 * t348 + t297 * t345) * t356;
t307 = t317 * qJD(3) + qJDD(1) * t364;
t373 = -t345 * t398 + t348 * t349;
t330 = qJD(1) * t373 + qJD(3);
t351 = sin(qJ(4));
t355 = cos(qJ(4));
t311 = -t318 * t351 + t330 * t355;
t327 = qJDD(1) * t373 + qJDD(3);
t283 = qJD(4) * t311 + t307 * t355 + t327 * t351;
t312 = t318 * t355 + t330 * t351;
t316 = qJD(4) - t317;
t350 = sin(qJ(5));
t354 = cos(qJ(5));
t295 = -t312 * t350 + t316 * t354;
t303 = qJDD(4) - t306;
t251 = qJD(5) * t295 + t283 * t354 + t303 * t350;
t296 = t312 * t354 + t316 * t350;
t264 = -mrSges(7,1) * t295 + mrSges(7,2) * t296;
t238 = t287 * t395 + t356 * t288 + t297 * t400;
t305 = -pkin(3) * t317 - pkin(10) * t318;
t326 = t330 ^ 2;
t234 = -pkin(3) * t326 + pkin(10) * t327 + t305 * t317 + t238;
t262 = -t345 * t287 + t348 * t297;
t236 = (-t317 * t330 - t307) * pkin(10) + (t318 * t330 - t306) * pkin(3) + t262;
t230 = t355 * t234 + t351 * t236;
t290 = -pkin(4) * t311 - pkin(11) * t312;
t315 = t316 ^ 2;
t225 = -pkin(4) * t315 + pkin(11) * t303 + t290 * t311 + t230;
t233 = -t327 * pkin(3) - t326 * pkin(10) + t318 * t305 - t237;
t282 = -qJD(4) * t312 - t307 * t351 + t327 * t355;
t228 = (-t311 * t316 - t283) * pkin(11) + (t312 * t316 - t282) * pkin(4) + t233;
t219 = -t350 * t225 + t354 * t228;
t280 = qJDD(5) - t282;
t308 = qJD(5) - t311;
t215 = -0.2e1 * qJD(6) * t296 + (t295 * t308 - t251) * qJ(6) + (t295 * t296 + t280) * pkin(5) + t219;
t267 = -mrSges(7,2) * t308 + mrSges(7,3) * t295;
t388 = m(7) * t215 + t280 * mrSges(7,1) + t308 * t267;
t212 = -t251 * mrSges(7,3) - t296 * t264 + t388;
t220 = t354 * t225 + t350 * t228;
t250 = -qJD(5) * t296 - t283 * t350 + t303 * t354;
t257 = Ifges(6,4) * t296 + Ifges(6,2) * t295 + Ifges(6,6) * t308;
t258 = Ifges(7,1) * t296 + Ifges(7,4) * t295 + Ifges(7,5) * t308;
t259 = Ifges(6,1) * t296 + Ifges(6,4) * t295 + Ifges(6,5) * t308;
t269 = pkin(5) * t308 - qJ(6) * t296;
t294 = t295 ^ 2;
t218 = -pkin(5) * t294 + qJ(6) * t250 + 0.2e1 * qJD(6) * t295 - t269 * t308 + t220;
t256 = Ifges(7,4) * t296 + Ifges(7,2) * t295 + Ifges(7,6) * t308;
t369 = -mrSges(7,1) * t215 + mrSges(7,2) * t218 - Ifges(7,5) * t251 - Ifges(7,6) * t250 - Ifges(7,3) * t280 - t296 * t256;
t408 = mrSges(6,1) * t219 - mrSges(6,2) * t220 + Ifges(6,5) * t251 + Ifges(6,6) * t250 + Ifges(6,3) * t280 + pkin(5) * t212 + t296 * t257 - (t259 + t258) * t295 - t369;
t406 = -mrSges(6,2) - mrSges(7,2);
t405 = Ifges(3,3) * t349;
t402 = t344 * t346;
t265 = -mrSges(6,1) * t295 + mrSges(6,2) * t296;
t268 = -mrSges(6,2) * t308 + mrSges(6,3) * t295;
t206 = m(6) * t219 + t280 * mrSges(6,1) + t308 * t268 + (-t264 - t265) * t296 + (-mrSges(6,3) - mrSges(7,3)) * t251 + t388;
t387 = m(7) * t218 + t250 * mrSges(7,3) + t295 * t264;
t270 = mrSges(7,1) * t308 - mrSges(7,3) * t296;
t392 = -mrSges(6,1) * t308 + mrSges(6,3) * t296 - t270;
t208 = m(6) * t220 + t250 * mrSges(6,3) + t295 * t265 + t280 * t406 + t308 * t392 + t387;
t205 = -t206 * t350 + t354 * t208;
t289 = -mrSges(5,1) * t311 + mrSges(5,2) * t312;
t299 = mrSges(5,1) * t316 - mrSges(5,3) * t312;
t202 = m(5) * t230 - mrSges(5,2) * t303 + mrSges(5,3) * t282 + t289 * t311 - t299 * t316 + t205;
t229 = -t351 * t234 + t236 * t355;
t224 = -pkin(4) * t303 - pkin(11) * t315 + t312 * t290 - t229;
t222 = -pkin(5) * t250 - qJ(6) * t294 + t269 * t296 + qJDD(6) + t224;
t382 = -m(7) * t222 + t250 * mrSges(7,1) + t295 * t267;
t211 = -m(6) * t224 + t250 * mrSges(6,1) + t251 * t406 + t295 * t268 + t296 * t392 + t382;
t298 = -mrSges(5,2) * t316 + mrSges(5,3) * t311;
t210 = m(5) * t229 + t303 * mrSges(5,1) - t283 * mrSges(5,3) - t312 * t289 + t316 * t298 + t211;
t195 = t351 * t202 + t355 * t210;
t304 = -mrSges(4,1) * t317 + mrSges(4,2) * t318;
t314 = mrSges(4,1) * t330 - mrSges(4,3) * t318;
t383 = t355 * t202 - t210 * t351;
t192 = m(4) * t238 - mrSges(4,2) * t327 + mrSges(4,3) * t306 + t304 * t317 - t314 * t330 + t383;
t313 = -mrSges(4,2) * t330 + mrSges(4,3) * t317;
t194 = m(4) * t262 - mrSges(4,1) * t306 + mrSges(4,2) * t307 - t313 * t317 + t314 * t318 + t195;
t204 = t206 * t354 + t208 * t350;
t361 = -m(5) * t233 + t282 * mrSges(5,1) - mrSges(5,2) * t283 + t311 * t298 - t299 * t312 - t204;
t199 = m(4) * t237 + mrSges(4,1) * t327 - mrSges(4,3) * t307 - t304 * t318 + t313 * t330 + t361;
t181 = t192 * t400 + t348 * t194 + t199 * t399;
t182 = t192 * t395 - t345 * t194 + t199 * t394;
t309 = -t344 * t333 + t378;
t381 = -mrSges(3,1) * t347 + mrSges(3,2) * t344;
t331 = t381 * t391;
t375 = -mrSges(3,2) * t349 + mrSges(3,3) * t398;
t336 = t375 * qJD(1);
t376 = mrSges(3,1) * t349 - mrSges(3,3) * t402;
t179 = m(3) * t309 + t376 * qJDD(1) + (-t331 * t402 + t336 * t349) * qJD(1) + t182;
t187 = t356 * t192 - t352 * t199;
t310 = -g(3) * t402 + t386;
t335 = t376 * qJD(1);
t186 = m(3) * t310 + t375 * qJDD(1) + (t331 * t398 - t335 * t349) * qJD(1) + t187;
t176 = -t179 * t344 + t347 * t186;
t319 = -t346 * t332 + t384;
t180 = m(3) * t319 + (t381 * qJDD(1) + (t335 * t344 - t336 * t347) * qJD(1)) * t346 + t181;
t171 = t179 * t396 - t180 * t346 + t186 * t401;
t380 = Ifges(3,5) * t344 + Ifges(3,6) * t347;
t254 = Ifges(7,5) * t296 + Ifges(7,6) * t295 + Ifges(7,3) * t308;
t255 = Ifges(6,5) * t296 + Ifges(6,6) * t295 + Ifges(6,3) * t308;
t370 = -mrSges(7,1) * t222 + mrSges(7,3) * t218 + Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t280 + t308 * t258;
t196 = Ifges(6,4) * t251 + Ifges(6,2) * t250 + Ifges(6,6) * t280 + t308 * t259 - mrSges(6,1) * t224 + mrSges(6,3) * t220 - pkin(5) * (t251 * mrSges(7,2) - t382) + qJ(6) * (-t280 * mrSges(7,2) - t308 * t270 + t387) + (-pkin(5) * t270 - t254 - t255) * t296 + t370;
t368 = mrSges(7,2) * t222 - mrSges(7,3) * t215 + Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t280 + t295 * t254;
t203 = mrSges(6,2) * t224 - mrSges(6,3) * t219 + Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t280 - qJ(6) * t212 + t295 * t255 + (-t256 - t257) * t308 + t368;
t276 = Ifges(5,5) * t312 + Ifges(5,6) * t311 + Ifges(5,3) * t316;
t277 = Ifges(5,4) * t312 + Ifges(5,2) * t311 + Ifges(5,6) * t316;
t183 = mrSges(5,2) * t233 - mrSges(5,3) * t229 + Ifges(5,1) * t283 + Ifges(5,4) * t282 + Ifges(5,5) * t303 - pkin(11) * t204 - t196 * t350 + t203 * t354 + t276 * t311 - t277 * t316;
t278 = Ifges(5,1) * t312 + Ifges(5,4) * t311 + Ifges(5,5) * t316;
t188 = -mrSges(5,1) * t233 + mrSges(5,3) * t230 + Ifges(5,4) * t283 + Ifges(5,2) * t282 + Ifges(5,6) * t303 - pkin(4) * t204 - t312 * t276 + t316 * t278 - t408;
t300 = Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t330;
t301 = Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t330;
t173 = mrSges(4,2) * t262 - mrSges(4,3) * t237 + Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t327 - pkin(10) * t195 + t183 * t355 - t188 * t351 + t300 * t317 - t301 * t330;
t302 = Ifges(4,1) * t318 + Ifges(4,4) * t317 + Ifges(4,5) * t330;
t359 = mrSges(5,1) * t229 - mrSges(5,2) * t230 + Ifges(5,5) * t283 + Ifges(5,6) * t282 + Ifges(5,3) * t303 + pkin(4) * t211 + pkin(11) * t205 + t354 * t196 + t350 * t203 + t312 * t277 - t311 * t278;
t177 = -mrSges(4,1) * t262 + mrSges(4,3) * t238 + Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t327 - pkin(3) * t195 - t318 * t300 + t330 * t302 - t359;
t371 = pkin(9) * t187 + t173 * t352 + t177 * t356;
t367 = Ifges(3,5) * t349 + (Ifges(3,1) * t344 + Ifges(3,4) * t347) * t346;
t366 = Ifges(3,6) * t349 + (Ifges(3,4) * t344 + Ifges(3,2) * t347) * t346;
t172 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t327 + pkin(3) * t361 + pkin(10) * t383 + t351 * t183 + t355 * t188 + t318 * t301 - t317 * t302;
t323 = t366 * qJD(1);
t324 = t367 * qJD(1);
t163 = qJDD(1) * t405 + mrSges(3,1) * t309 - mrSges(3,2) * t310 + pkin(2) * t182 + t348 * t172 + t371 * t345 + (t380 * qJDD(1) + (t323 * t344 - t324 * t347) * qJD(1)) * t346;
t322 = (t346 * t380 + t405) * qJD(1);
t165 = -mrSges(3,1) * t319 + mrSges(3,3) * t310 - pkin(2) * t181 - t345 * t172 + (-t322 * t402 + t324 * t349) * qJD(1) + t371 * t348 + t366 * qJDD(1);
t167 = mrSges(3,2) * t319 - mrSges(3,3) * t309 + t356 * t173 - t352 * t177 + (t322 * t398 - t323 * t349) * qJD(1) + (-t181 * t345 - t182 * t348) * pkin(9) + t367 * qJDD(1);
t365 = mrSges(2,1) * t340 - mrSges(2,2) * t341 + Ifges(2,3) * qJDD(1) + pkin(1) * t171 + t349 * t163 + t165 * t398 + t167 * t402 + t176 * t403;
t174 = m(2) * t341 - mrSges(2,1) * t358 - qJDD(1) * mrSges(2,2) + t176;
t170 = t349 * t180 + (t179 * t347 + t186 * t344) * t346;
t168 = m(2) * t340 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t358 + t171;
t161 = -mrSges(2,2) * g(3) - mrSges(2,3) * t340 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) - t344 * t165 + t347 * t167 + (-t170 * t346 - t171 * t349) * qJ(2);
t160 = mrSges(2,1) * g(3) + mrSges(2,3) * t341 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t170 - t346 * t163 + (qJ(2) * t176 + t165 * t347 + t167 * t344) * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t161 - t353 * t160 - pkin(8) * (t168 * t357 + t174 * t353), t161, t167, t173, t183, t203, -t256 * t308 + t368; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t353 * t161 + t357 * t160 + pkin(8) * (-t168 * t353 + t174 * t357), t160, t165, t177, t188, t196, -t296 * t254 + t370; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t365, t365, t163, t172, t359, t408, -t295 * t258 - t369;];
m_new  = t1;
