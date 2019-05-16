% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 22:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:12:04
% EndTime: 2019-05-07 22:14:53
% DurationCPUTime: 86.26s
% Computational Cost: add. (1623131->399), mult. (3473841->517), div. (0->0), fcn. (2839989->14), ass. (0->165)
t381 = cos(qJ(4));
t343 = sin(pkin(6));
t380 = pkin(8) * t343;
t345 = cos(pkin(6));
t379 = t345 * g(3);
t349 = sin(qJ(2));
t378 = t343 * t349;
t353 = cos(qJ(2));
t377 = t343 * t353;
t376 = t345 * t349;
t375 = t345 * t353;
t373 = qJD(1) * t343;
t326 = (-pkin(2) * t353 - pkin(9) * t349) * t373;
t338 = qJD(1) * t345 + qJD(2);
t336 = t338 ^ 2;
t337 = qJDD(1) * t345 + qJDD(2);
t372 = qJD(1) * t353;
t350 = sin(qJ(1));
t354 = cos(qJ(1));
t333 = t350 * g(1) - t354 * g(2);
t355 = qJD(1) ^ 2;
t323 = qJDD(1) * pkin(1) + t355 * t380 + t333;
t334 = -t354 * g(1) - t350 * g(2);
t371 = qJDD(1) * t343;
t324 = -t355 * pkin(1) + pkin(8) * t371 + t334;
t374 = t323 * t376 + t353 * t324;
t279 = -t336 * pkin(2) + t337 * pkin(9) + (-g(3) * t349 + t326 * t372) * t343 + t374;
t327 = (qJD(2) * t372 + qJDD(1) * t349) * t343;
t370 = t349 * t373;
t328 = -qJD(2) * t370 + t353 * t371;
t280 = -t328 * pkin(2) - t327 * pkin(9) - t379 + (-t323 + (pkin(2) * t349 - pkin(9) * t353) * t338 * qJD(1)) * t343;
t348 = sin(qJ(3));
t352 = cos(qJ(3));
t246 = -t279 * t348 + t352 * t280;
t315 = t338 * t352 - t348 * t370;
t296 = t315 * qJD(3) + t327 * t352 + t348 * t337;
t316 = t348 * t338 + t352 * t370;
t320 = qJDD(3) - t328;
t369 = t343 * t372;
t332 = qJD(3) - t369;
t233 = (t315 * t332 - t296) * pkin(10) + (t315 * t316 + t320) * pkin(3) + t246;
t247 = t352 * t279 + t348 * t280;
t295 = -t316 * qJD(3) - t348 * t327 + t337 * t352;
t305 = pkin(3) * t332 - pkin(10) * t316;
t314 = t315 ^ 2;
t240 = -pkin(3) * t314 + pkin(10) * t295 - t305 * t332 + t247;
t347 = sin(qJ(4));
t225 = t347 * t233 + t240 * t381;
t301 = t347 * t315 + t316 * t381;
t259 = qJD(4) * t301 - t295 * t381 + t296 * t347;
t300 = -t315 * t381 + t316 * t347;
t274 = mrSges(5,1) * t300 + mrSges(5,2) * t301;
t330 = qJD(4) + t332;
t287 = mrSges(5,1) * t330 - mrSges(5,3) * t301;
t319 = qJDD(4) + t320;
t273 = pkin(4) * t300 - qJ(5) * t301;
t329 = t330 ^ 2;
t219 = -pkin(4) * t329 + qJ(5) * t319 - t273 * t300 + t225;
t297 = -g(3) * t377 + t323 * t375 - t349 * t324;
t278 = -t337 * pkin(2) - t336 * pkin(9) + t326 * t370 - t297;
t242 = -t295 * pkin(3) - t314 * pkin(10) + t316 * t305 + t278;
t260 = -t300 * qJD(4) + t347 * t295 + t296 * t381;
t222 = (t300 * t330 - t260) * qJ(5) + (t301 * t330 + t259) * pkin(4) + t242;
t342 = sin(pkin(12));
t344 = cos(pkin(12));
t285 = t301 * t344 + t330 * t342;
t214 = -0.2e1 * qJD(5) * t285 - t219 * t342 + t344 * t222;
t249 = t260 * t344 + t319 * t342;
t284 = -t301 * t342 + t330 * t344;
t212 = (t284 * t300 - t249) * pkin(11) + (t284 * t285 + t259) * pkin(5) + t214;
t215 = 0.2e1 * qJD(5) * t284 + t344 * t219 + t342 * t222;
t248 = -t260 * t342 + t319 * t344;
t268 = pkin(5) * t300 - pkin(11) * t285;
t283 = t284 ^ 2;
t213 = -pkin(5) * t283 + pkin(11) * t248 - t268 * t300 + t215;
t346 = sin(qJ(6));
t351 = cos(qJ(6));
t210 = t212 * t351 - t213 * t346;
t263 = t284 * t351 - t285 * t346;
t230 = qJD(6) * t263 + t248 * t346 + t249 * t351;
t264 = t284 * t346 + t285 * t351;
t241 = -mrSges(7,1) * t263 + mrSges(7,2) * t264;
t299 = qJD(6) + t300;
t244 = -mrSges(7,2) * t299 + mrSges(7,3) * t263;
t258 = qJDD(6) + t259;
t205 = m(7) * t210 + mrSges(7,1) * t258 - mrSges(7,3) * t230 - t241 * t264 + t244 * t299;
t211 = t212 * t346 + t213 * t351;
t229 = -qJD(6) * t264 + t248 * t351 - t249 * t346;
t245 = mrSges(7,1) * t299 - mrSges(7,3) * t264;
t206 = m(7) * t211 - mrSges(7,2) * t258 + mrSges(7,3) * t229 + t241 * t263 - t245 * t299;
t197 = t351 * t205 + t346 * t206;
t265 = -mrSges(6,1) * t284 + mrSges(6,2) * t285;
t266 = -mrSges(6,2) * t300 + mrSges(6,3) * t284;
t195 = m(6) * t214 + mrSges(6,1) * t259 - mrSges(6,3) * t249 - t265 * t285 + t266 * t300 + t197;
t267 = mrSges(6,1) * t300 - mrSges(6,3) * t285;
t365 = -t205 * t346 + t351 * t206;
t196 = m(6) * t215 - mrSges(6,2) * t259 + mrSges(6,3) * t248 + t265 * t284 - t267 * t300 + t365;
t366 = -t195 * t342 + t344 * t196;
t188 = m(5) * t225 - mrSges(5,2) * t319 - mrSges(5,3) * t259 - t274 * t300 - t287 * t330 + t366;
t224 = t233 * t381 - t347 * t240;
t286 = -mrSges(5,2) * t330 - mrSges(5,3) * t300;
t218 = -t319 * pkin(4) - t329 * qJ(5) + t301 * t273 + qJDD(5) - t224;
t216 = -t248 * pkin(5) - t283 * pkin(11) + t285 * t268 + t218;
t364 = m(7) * t216 - t229 * mrSges(7,1) + mrSges(7,2) * t230 - t263 * t244 + t245 * t264;
t359 = -m(6) * t218 + t248 * mrSges(6,1) - mrSges(6,2) * t249 + t284 * t266 - t267 * t285 - t364;
t201 = m(5) * t224 + mrSges(5,1) * t319 - mrSges(5,3) * t260 - t274 * t301 + t286 * t330 + t359;
t179 = t347 * t188 + t201 * t381;
t302 = -mrSges(4,1) * t315 + mrSges(4,2) * t316;
t303 = -mrSges(4,2) * t332 + mrSges(4,3) * t315;
t177 = m(4) * t246 + mrSges(4,1) * t320 - mrSges(4,3) * t296 - t302 * t316 + t303 * t332 + t179;
t304 = mrSges(4,1) * t332 - mrSges(4,3) * t316;
t367 = t188 * t381 - t201 * t347;
t178 = m(4) * t247 - mrSges(4,2) * t320 + mrSges(4,3) * t295 + t302 * t315 - t304 * t332 + t367;
t172 = t352 * t177 + t348 * t178;
t190 = t344 * t195 + t342 * t196;
t298 = -g(3) * t378 + t374;
t321 = mrSges(3,1) * t338 - mrSges(3,3) * t370;
t325 = (-mrSges(3,1) * t353 + mrSges(3,2) * t349) * t373;
t368 = -t348 * t177 + t352 * t178;
t170 = m(3) * t298 - t337 * mrSges(3,2) + t328 * mrSges(3,3) - t338 * t321 + t325 * t369 + t368;
t322 = -t338 * mrSges(3,2) + mrSges(3,3) * t369;
t361 = m(5) * t242 + t259 * mrSges(5,1) + mrSges(5,2) * t260 + t300 * t286 + t287 * t301 + t190;
t358 = -m(4) * t278 + t295 * mrSges(4,1) - mrSges(4,2) * t296 + t315 * t303 - t304 * t316 - t361;
t185 = m(3) * t297 + mrSges(3,1) * t337 - mrSges(3,3) * t327 + t322 * t338 - t325 * t370 + t358;
t166 = t353 * t170 - t185 * t349;
t309 = -t343 * t323 - t379;
t171 = m(3) * t309 - t328 * mrSges(3,1) + t327 * mrSges(3,2) + (t321 * t349 - t322 * t353) * t373 + t172;
t163 = t170 * t376 - t171 * t343 + t185 * t375;
t234 = Ifges(7,5) * t264 + Ifges(7,6) * t263 + Ifges(7,3) * t299;
t236 = Ifges(7,1) * t264 + Ifges(7,4) * t263 + Ifges(7,5) * t299;
t198 = -mrSges(7,1) * t216 + mrSges(7,3) * t211 + Ifges(7,4) * t230 + Ifges(7,2) * t229 + Ifges(7,6) * t258 - t234 * t264 + t236 * t299;
t235 = Ifges(7,4) * t264 + Ifges(7,2) * t263 + Ifges(7,6) * t299;
t199 = mrSges(7,2) * t216 - mrSges(7,3) * t210 + Ifges(7,1) * t230 + Ifges(7,4) * t229 + Ifges(7,5) * t258 + t234 * t263 - t235 * t299;
t250 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t300;
t252 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t300;
t181 = -mrSges(6,1) * t218 + mrSges(6,3) * t215 + Ifges(6,4) * t249 + Ifges(6,2) * t248 + Ifges(6,6) * t259 - pkin(5) * t364 + pkin(11) * t365 + t351 * t198 + t346 * t199 - t285 * t250 + t300 * t252;
t251 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t300;
t183 = mrSges(6,2) * t218 - mrSges(6,3) * t214 + Ifges(6,1) * t249 + Ifges(6,4) * t248 + Ifges(6,5) * t259 - pkin(11) * t197 - t198 * t346 + t199 * t351 + t250 * t284 - t251 * t300;
t269 = Ifges(5,5) * t301 - Ifges(5,6) * t300 + Ifges(5,3) * t330;
t270 = Ifges(5,4) * t301 - Ifges(5,2) * t300 + Ifges(5,6) * t330;
t167 = mrSges(5,2) * t242 - mrSges(5,3) * t224 + Ifges(5,1) * t260 - Ifges(5,4) * t259 + Ifges(5,5) * t319 - qJ(5) * t190 - t181 * t342 + t183 * t344 - t269 * t300 - t270 * t330;
t271 = Ifges(5,1) * t301 - Ifges(5,4) * t300 + Ifges(5,5) * t330;
t362 = -mrSges(7,1) * t210 + mrSges(7,2) * t211 - Ifges(7,5) * t230 - Ifges(7,6) * t229 - Ifges(7,3) * t258 - t264 * t235 + t263 * t236;
t357 = -mrSges(6,1) * t214 + mrSges(6,2) * t215 - Ifges(6,5) * t249 - Ifges(6,6) * t248 - pkin(5) * t197 - t285 * t251 + t284 * t252 + t362;
t173 = (-Ifges(5,2) - Ifges(6,3)) * t259 + t330 * t271 + Ifges(5,6) * t319 - t301 * t269 + Ifges(5,4) * t260 - mrSges(5,1) * t242 + mrSges(5,3) * t225 - pkin(4) * t190 + t357;
t289 = Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t332;
t291 = Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t332;
t156 = -mrSges(4,1) * t278 + mrSges(4,3) * t247 + Ifges(4,4) * t296 + Ifges(4,2) * t295 + Ifges(4,6) * t320 - pkin(3) * t361 + pkin(10) * t367 + t347 * t167 + t173 * t381 - t316 * t289 + t332 * t291;
t290 = Ifges(4,4) * t316 + Ifges(4,2) * t315 + Ifges(4,6) * t332;
t159 = mrSges(4,2) * t278 - mrSges(4,3) * t246 + Ifges(4,1) * t296 + Ifges(4,4) * t295 + Ifges(4,5) * t320 - pkin(10) * t179 + t167 * t381 - t347 * t173 + t315 * t289 - t332 * t290;
t307 = Ifges(3,6) * t338 + (Ifges(3,4) * t349 + Ifges(3,2) * t353) * t373;
t308 = Ifges(3,5) * t338 + (Ifges(3,1) * t349 + Ifges(3,4) * t353) * t373;
t153 = Ifges(3,5) * t327 + Ifges(3,6) * t328 + Ifges(3,3) * t337 + mrSges(3,1) * t297 - mrSges(3,2) * t298 + t348 * t159 + t352 * t156 + pkin(2) * t358 + pkin(9) * t368 + (t307 * t349 - t308 * t353) * t373;
t306 = Ifges(3,3) * t338 + (Ifges(3,5) * t349 + Ifges(3,6) * t353) * t373;
t155 = mrSges(3,2) * t309 - mrSges(3,3) * t297 + Ifges(3,1) * t327 + Ifges(3,4) * t328 + Ifges(3,5) * t337 - pkin(9) * t172 - t348 * t156 + t352 * t159 + t306 * t369 - t338 * t307;
t360 = -mrSges(5,1) * t224 + mrSges(5,2) * t225 - Ifges(5,5) * t260 + Ifges(5,6) * t259 - Ifges(5,3) * t319 - pkin(4) * t359 - qJ(5) * t366 - t344 * t181 - t342 * t183 - t301 * t270 - t300 * t271;
t356 = mrSges(4,1) * t246 - mrSges(4,2) * t247 + Ifges(4,5) * t296 + Ifges(4,6) * t295 + Ifges(4,3) * t320 + pkin(3) * t179 + t316 * t290 - t315 * t291 - t360;
t158 = -mrSges(3,1) * t309 + mrSges(3,3) * t298 + Ifges(3,4) * t327 + Ifges(3,2) * t328 + Ifges(3,6) * t337 - pkin(2) * t172 - t306 * t370 + t338 * t308 - t356;
t363 = mrSges(2,1) * t333 - mrSges(2,2) * t334 + Ifges(2,3) * qJDD(1) + pkin(1) * t163 + t345 * t153 + t155 * t378 + t158 * t377 + t166 * t380;
t164 = m(2) * t334 - t355 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t166;
t162 = t345 * t171 + (t170 * t349 + t185 * t353) * t343;
t160 = m(2) * t333 + qJDD(1) * mrSges(2,1) - t355 * mrSges(2,2) + t163;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t333 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) + t353 * t155 - t349 * t158 + (-t162 * t343 - t163 * t345) * pkin(8);
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t334 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t343 * t153 + (pkin(8) * t166 + t155 * t349 + t158 * t353) * t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t151 - t350 * t150 - pkin(7) * (t354 * t160 + t350 * t164), t151, t155, t159, t167, t183, t199; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t350 * t151 + t354 * t150 + pkin(7) * (-t350 * t160 + t354 * t164), t150, t158, t156, t173, t181, t198; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t363, t363, t153, t356, -t360, Ifges(6,3) * t259 - t357, -t362;];
m_new  = t1;
