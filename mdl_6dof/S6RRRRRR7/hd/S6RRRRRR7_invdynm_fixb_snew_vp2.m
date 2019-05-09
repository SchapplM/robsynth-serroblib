% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 12:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 12:17:01
% EndTime: 2019-05-08 12:20:33
% DurationCPUTime: 103.73s
% Computational Cost: add. (1930903->399), mult. (4114122->516), div. (0->0), fcn. (3352932->14), ass. (0->167)
t340 = sin(pkin(6));
t346 = sin(qJ(2));
t352 = cos(qJ(2));
t369 = qJD(1) * qJD(2);
t327 = (-qJDD(1) * t352 + t346 * t369) * t340;
t378 = pkin(8) * t340;
t341 = cos(pkin(6));
t377 = t341 * g(3);
t376 = t340 * t346;
t375 = t340 * t352;
t374 = t341 * t346;
t373 = t341 * t352;
t371 = qJD(1) * t340;
t325 = (-pkin(2) * t352 - pkin(9) * t346) * t371;
t336 = t341 * qJD(1) + qJD(2);
t334 = t336 ^ 2;
t335 = t341 * qJDD(1) + qJDD(2);
t370 = qJD(1) * t352;
t347 = sin(qJ(1));
t353 = cos(qJ(1));
t332 = t347 * g(1) - t353 * g(2);
t354 = qJD(1) ^ 2;
t322 = qJDD(1) * pkin(1) + t354 * t378 + t332;
t333 = -t353 * g(1) - t347 * g(2);
t323 = -t354 * pkin(1) + qJDD(1) * t378 + t333;
t372 = t322 * t374 + t352 * t323;
t273 = -t334 * pkin(2) + t335 * pkin(9) + (-g(3) * t346 + t325 * t370) * t340 + t372;
t326 = (qJDD(1) * t346 + t352 * t369) * t340;
t274 = t327 * pkin(2) - t326 * pkin(9) - t377 + (-t322 + (pkin(2) * t346 - pkin(9) * t352) * t336 * qJD(1)) * t340;
t345 = sin(qJ(3));
t351 = cos(qJ(3));
t250 = t351 * t273 + t345 * t274;
t368 = t346 * t371;
t314 = t351 * t336 - t345 * t368;
t315 = t345 * t336 + t351 * t368;
t298 = -t314 * pkin(3) - t315 * pkin(10);
t319 = qJDD(3) + t327;
t367 = t340 * t370;
t331 = qJD(3) - t367;
t329 = t331 ^ 2;
t242 = -t329 * pkin(3) + t319 * pkin(10) + t314 * t298 + t250;
t295 = -g(3) * t375 + t322 * t373 - t346 * t323;
t272 = -t335 * pkin(2) - t334 * pkin(9) + t325 * t368 - t295;
t293 = -t315 * qJD(3) - t345 * t326 + t351 * t335;
t294 = t314 * qJD(3) + t351 * t326 + t345 * t335;
t246 = (-t314 * t331 - t294) * pkin(10) + (t315 * t331 - t293) * pkin(3) + t272;
t344 = sin(qJ(4));
t350 = cos(qJ(4));
t225 = -t344 * t242 + t350 * t246;
t300 = -t344 * t315 + t350 * t331;
t260 = t300 * qJD(4) + t350 * t294 + t344 * t319;
t291 = qJDD(4) - t293;
t301 = t350 * t315 + t344 * t331;
t313 = qJD(4) - t314;
t216 = (t300 * t313 - t260) * pkin(11) + (t300 * t301 + t291) * pkin(4) + t225;
t226 = t350 * t242 + t344 * t246;
t259 = -t301 * qJD(4) - t344 * t294 + t350 * t319;
t282 = t313 * pkin(4) - t301 * pkin(11);
t299 = t300 ^ 2;
t224 = -t299 * pkin(4) + t259 * pkin(11) - t313 * t282 + t226;
t343 = sin(qJ(5));
t349 = cos(qJ(5));
t210 = t349 * t216 - t343 * t224;
t276 = t349 * t300 - t343 * t301;
t239 = t276 * qJD(5) + t343 * t259 + t349 * t260;
t277 = t343 * t300 + t349 * t301;
t286 = qJDD(5) + t291;
t311 = qJD(5) + t313;
t207 = (t276 * t311 - t239) * pkin(12) + (t276 * t277 + t286) * pkin(5) + t210;
t211 = t343 * t216 + t349 * t224;
t238 = -t277 * qJD(5) + t349 * t259 - t343 * t260;
t263 = t311 * pkin(5) - t277 * pkin(12);
t275 = t276 ^ 2;
t208 = -t275 * pkin(5) + t238 * pkin(12) - t311 * t263 + t211;
t342 = sin(qJ(6));
t348 = cos(qJ(6));
t205 = t348 * t207 - t342 * t208;
t255 = t348 * t276 - t342 * t277;
t222 = t255 * qJD(6) + t342 * t238 + t348 * t239;
t256 = t342 * t276 + t348 * t277;
t234 = -t255 * mrSges(7,1) + t256 * mrSges(7,2);
t308 = qJD(6) + t311;
t247 = -t308 * mrSges(7,2) + t255 * mrSges(7,3);
t285 = qJDD(6) + t286;
t201 = m(7) * t205 + t285 * mrSges(7,1) - t222 * mrSges(7,3) - t256 * t234 + t308 * t247;
t206 = t342 * t207 + t348 * t208;
t221 = -t256 * qJD(6) + t348 * t238 - t342 * t239;
t248 = t308 * mrSges(7,1) - t256 * mrSges(7,3);
t202 = m(7) * t206 - t285 * mrSges(7,2) + t221 * mrSges(7,3) + t255 * t234 - t308 * t248;
t193 = t348 * t201 + t342 * t202;
t257 = -t276 * mrSges(6,1) + t277 * mrSges(6,2);
t261 = -t311 * mrSges(6,2) + t276 * mrSges(6,3);
t190 = m(6) * t210 + t286 * mrSges(6,1) - t239 * mrSges(6,3) - t277 * t257 + t311 * t261 + t193;
t262 = t311 * mrSges(6,1) - t277 * mrSges(6,3);
t364 = -t342 * t201 + t348 * t202;
t191 = m(6) * t211 - t286 * mrSges(6,2) + t238 * mrSges(6,3) + t276 * t257 - t311 * t262 + t364;
t186 = t349 * t190 + t343 * t191;
t278 = -t300 * mrSges(5,1) + t301 * mrSges(5,2);
t280 = -t313 * mrSges(5,2) + t300 * mrSges(5,3);
t184 = m(5) * t225 + t291 * mrSges(5,1) - t260 * mrSges(5,3) - t301 * t278 + t313 * t280 + t186;
t281 = t313 * mrSges(5,1) - t301 * mrSges(5,3);
t365 = -t343 * t190 + t349 * t191;
t185 = m(5) * t226 - t291 * mrSges(5,2) + t259 * mrSges(5,3) + t300 * t278 - t313 * t281 + t365;
t180 = -t344 * t184 + t350 * t185;
t297 = -t314 * mrSges(4,1) + t315 * mrSges(4,2);
t303 = t331 * mrSges(4,1) - t315 * mrSges(4,3);
t178 = m(4) * t250 - t319 * mrSges(4,2) + t293 * mrSges(4,3) + t314 * t297 - t331 * t303 + t180;
t249 = -t345 * t273 + t351 * t274;
t241 = -t319 * pkin(3) - t329 * pkin(10) + t315 * t298 - t249;
t227 = -t259 * pkin(4) - t299 * pkin(11) + t301 * t282 + t241;
t213 = -t238 * pkin(5) - t275 * pkin(12) + t277 * t263 + t227;
t363 = m(7) * t213 - t221 * mrSges(7,1) + t222 * mrSges(7,2) - t255 * t247 + t256 * t248;
t359 = m(6) * t227 - t238 * mrSges(6,1) + t239 * mrSges(6,2) - t276 * t261 + t277 * t262 + t363;
t203 = -m(5) * t241 + t259 * mrSges(5,1) - t260 * mrSges(5,2) + t300 * t280 - t301 * t281 - t359;
t302 = -t331 * mrSges(4,2) + t314 * mrSges(4,3);
t197 = m(4) * t249 + t319 * mrSges(4,1) - t294 * mrSges(4,3) - t315 * t297 + t331 * t302 + t203;
t173 = t345 * t178 + t351 * t197;
t296 = -g(3) * t376 + t372;
t320 = t336 * mrSges(3,1) - mrSges(3,3) * t368;
t324 = (-mrSges(3,1) * t352 + mrSges(3,2) * t346) * t371;
t366 = t351 * t178 - t345 * t197;
t171 = m(3) * t296 - t335 * mrSges(3,2) - t327 * mrSges(3,3) - t336 * t320 + t324 * t367 + t366;
t321 = -t336 * mrSges(3,2) + mrSges(3,3) * t367;
t179 = t350 * t184 + t344 * t185;
t358 = -m(4) * t272 + t293 * mrSges(4,1) - t294 * mrSges(4,2) + t314 * t302 - t315 * t303 - t179;
t175 = m(3) * t295 + t335 * mrSges(3,1) - t326 * mrSges(3,3) + t336 * t321 - t324 * t368 + t358;
t165 = t352 * t171 - t346 * t175;
t307 = -t340 * t322 - t377;
t172 = m(3) * t307 + t327 * mrSges(3,1) + t326 * mrSges(3,2) + (t320 * t346 - t321 * t352) * t371 + t173;
t162 = t171 * t374 - t340 * t172 + t175 * t373;
t229 = Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * t308;
t231 = Ifges(7,1) * t256 + Ifges(7,4) * t255 + Ifges(7,5) * t308;
t194 = -mrSges(7,1) * t213 + mrSges(7,3) * t206 + Ifges(7,4) * t222 + Ifges(7,2) * t221 + Ifges(7,6) * t285 - t256 * t229 + t308 * t231;
t230 = Ifges(7,4) * t256 + Ifges(7,2) * t255 + Ifges(7,6) * t308;
t195 = mrSges(7,2) * t213 - mrSges(7,3) * t205 + Ifges(7,1) * t222 + Ifges(7,4) * t221 + Ifges(7,5) * t285 + t255 * t229 - t308 * t230;
t251 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t311;
t253 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t311;
t181 = -mrSges(6,1) * t227 + mrSges(6,3) * t211 + Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t286 - pkin(5) * t363 + pkin(12) * t364 + t348 * t194 + t342 * t195 - t277 * t251 + t311 * t253;
t252 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t311;
t182 = mrSges(6,2) * t227 - mrSges(6,3) * t210 + Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t286 - pkin(12) * t193 - t342 * t194 + t348 * t195 + t276 * t251 - t311 * t252;
t264 = Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t313;
t266 = Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t313;
t167 = -mrSges(5,1) * t241 + mrSges(5,3) * t226 + Ifges(5,4) * t260 + Ifges(5,2) * t259 + Ifges(5,6) * t291 - pkin(4) * t359 + pkin(11) * t365 + t349 * t181 + t343 * t182 - t301 * t264 + t313 * t266;
t265 = Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * t313;
t168 = mrSges(5,2) * t241 - mrSges(5,3) * t225 + Ifges(5,1) * t260 + Ifges(5,4) * t259 + Ifges(5,5) * t291 - pkin(11) * t186 - t343 * t181 + t349 * t182 + t300 * t264 - t313 * t265;
t287 = Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t331;
t288 = Ifges(4,4) * t315 + Ifges(4,2) * t314 + Ifges(4,6) * t331;
t158 = mrSges(4,2) * t272 - mrSges(4,3) * t249 + Ifges(4,1) * t294 + Ifges(4,4) * t293 + Ifges(4,5) * t319 - pkin(10) * t179 - t344 * t167 + t350 * t168 + t314 * t287 - t331 * t288;
t289 = Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t331;
t360 = -mrSges(7,1) * t205 + mrSges(7,2) * t206 - Ifges(7,5) * t222 - Ifges(7,6) * t221 - Ifges(7,3) * t285 - t256 * t230 + t255 * t231;
t357 = -mrSges(6,1) * t210 + mrSges(6,2) * t211 - Ifges(6,5) * t239 - Ifges(6,6) * t238 - Ifges(6,3) * t286 - pkin(5) * t193 - t277 * t252 + t276 * t253 + t360;
t355 = mrSges(5,1) * t225 - mrSges(5,2) * t226 + Ifges(5,5) * t260 + Ifges(5,6) * t259 + Ifges(5,3) * t291 + pkin(4) * t186 + t301 * t265 - t300 * t266 - t357;
t166 = -mrSges(4,1) * t272 + mrSges(4,3) * t250 + Ifges(4,4) * t294 + Ifges(4,2) * t293 + Ifges(4,6) * t319 - pkin(3) * t179 - t315 * t287 + t331 * t289 - t355;
t305 = Ifges(3,6) * t336 + (Ifges(3,4) * t346 + Ifges(3,2) * t352) * t371;
t306 = Ifges(3,5) * t336 + (Ifges(3,1) * t346 + Ifges(3,4) * t352) * t371;
t153 = Ifges(3,5) * t326 - Ifges(3,6) * t327 + Ifges(3,3) * t335 + mrSges(3,1) * t295 - mrSges(3,2) * t296 + t345 * t158 + t351 * t166 + pkin(2) * t358 + pkin(9) * t366 + (t305 * t346 - t306 * t352) * t371;
t304 = Ifges(3,3) * t336 + (Ifges(3,5) * t346 + Ifges(3,6) * t352) * t371;
t155 = mrSges(3,2) * t307 - mrSges(3,3) * t295 + Ifges(3,1) * t326 - Ifges(3,4) * t327 + Ifges(3,5) * t335 - pkin(9) * t173 + t351 * t158 - t345 * t166 + t304 * t367 - t336 * t305;
t356 = mrSges(4,1) * t249 - mrSges(4,2) * t250 + Ifges(4,5) * t294 + Ifges(4,6) * t293 + Ifges(4,3) * t319 + pkin(3) * t203 + pkin(10) * t180 + t350 * t167 + t344 * t168 + t315 * t288 - t314 * t289;
t157 = -mrSges(3,1) * t307 + mrSges(3,3) * t296 + Ifges(3,4) * t326 - Ifges(3,2) * t327 + Ifges(3,6) * t335 - pkin(2) * t173 - t304 * t368 + t336 * t306 - t356;
t361 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t162 + t341 * t153 + t155 * t376 + t157 * t375 + t165 * t378;
t163 = m(2) * t333 - t354 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t165;
t161 = t341 * t172 + (t171 * t346 + t175 * t352) * t340;
t159 = m(2) * t332 + qJDD(1) * mrSges(2,1) - t354 * mrSges(2,2) + t162;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t354 * Ifges(2,6) + t352 * t155 - t346 * t157 + (-t161 * t340 - t162 * t341) * pkin(8);
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t354 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t161 - t340 * t153 + (pkin(8) * t165 + t155 * t346 + t157 * t352) * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t353 * t151 - t347 * t150 - pkin(7) * (t353 * t159 + t347 * t163), t151, t155, t158, t168, t182, t195; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t347 * t151 + t353 * t150 + pkin(7) * (-t347 * t159 + t353 * t163), t150, t157, t166, t167, t181, t194; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, t153, t356, t355, -t357, -t360;];
m_new  = t1;
