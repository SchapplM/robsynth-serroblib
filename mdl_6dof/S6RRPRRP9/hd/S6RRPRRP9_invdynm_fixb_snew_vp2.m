% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:27:41
% EndTime: 2019-05-06 18:28:51
% DurationCPUTime: 37.70s
% Computational Cost: add. (654218->394), mult. (1483194->499), div. (0->0), fcn. (1190545->12), ass. (0->155)
t341 = cos(pkin(6));
t334 = qJD(1) * t341 + qJD(2);
t338 = sin(pkin(11));
t340 = cos(pkin(11));
t344 = sin(qJ(2));
t339 = sin(pkin(6));
t371 = qJD(1) * t339;
t366 = t344 * t371;
t312 = t334 * t340 - t338 * t366;
t313 = t334 * t338 + t340 * t366;
t343 = sin(qJ(4));
t347 = cos(qJ(4));
t294 = t312 * t347 - t313 * t343;
t348 = cos(qJ(2));
t370 = qJD(1) * t348;
t324 = (qJD(2) * t370 + qJDD(1) * t344) * t339;
t333 = qJDD(1) * t341 + qJDD(2);
t300 = -t324 * t338 + t333 * t340;
t301 = t324 * t340 + t333 * t338;
t263 = qJD(4) * t294 + t300 * t343 + t301 * t347;
t295 = t312 * t343 + t313 * t347;
t365 = t339 * t370;
t328 = qJD(4) - t365;
t342 = sin(qJ(5));
t346 = cos(qJ(5));
t282 = -t295 * t342 + t328 * t346;
t369 = qJDD(1) * t339;
t325 = -qJD(2) * t366 + t348 * t369;
t317 = qJDD(4) - t325;
t238 = qJD(5) * t282 + t263 * t346 + t317 * t342;
t283 = t295 * t346 + t328 * t342;
t252 = -mrSges(7,1) * t282 + mrSges(7,2) * t283;
t322 = (-pkin(2) * t348 - qJ(3) * t344) * t371;
t332 = t334 ^ 2;
t345 = sin(qJ(1));
t349 = cos(qJ(1));
t329 = g(1) * t345 - g(2) * t349;
t350 = qJD(1) ^ 2;
t381 = pkin(8) * t339;
t320 = qJDD(1) * pkin(1) + t350 * t381 + t329;
t330 = -g(1) * t349 - g(2) * t345;
t321 = -pkin(1) * t350 + pkin(8) * t369 + t330;
t376 = t341 * t344;
t372 = t320 * t376 + t321 * t348;
t279 = -t332 * pkin(2) + t333 * qJ(3) + (-g(3) * t344 + t322 * t370) * t339 + t372;
t380 = t341 * g(3);
t280 = -t325 * pkin(2) - t380 - t324 * qJ(3) + (-t320 + (pkin(2) * t344 - qJ(3) * t348) * t334 * qJD(1)) * t339;
t226 = -0.2e1 * qJD(3) * t313 - t338 * t279 + t280 * t340;
t219 = (-t312 * t365 - t301) * pkin(9) + (t312 * t313 - t325) * pkin(3) + t226;
t227 = 0.2e1 * qJD(3) * t312 + t279 * t340 + t280 * t338;
t302 = -pkin(3) * t365 - pkin(9) * t313;
t311 = t312 ^ 2;
t222 = -pkin(3) * t311 + pkin(9) * t300 + t302 * t365 + t227;
t214 = t219 * t343 + t222 * t347;
t274 = -pkin(4) * t294 - pkin(10) * t295;
t327 = t328 ^ 2;
t211 = -pkin(4) * t327 + pkin(10) * t317 + t274 * t294 + t214;
t375 = t341 * t348;
t377 = t339 * t348;
t291 = -g(3) * t377 + t320 * t375 - t321 * t344;
t278 = -pkin(2) * t333 - qJ(3) * t332 + t322 * t366 + qJDD(3) - t291;
t239 = -pkin(3) * t300 - pkin(9) * t311 + t302 * t313 + t278;
t262 = -qJD(4) * t295 + t300 * t347 - t301 * t343;
t217 = (-t294 * t328 - t263) * pkin(10) + (t295 * t328 - t262) * pkin(4) + t239;
t205 = -t342 * t211 + t217 * t346;
t261 = qJDD(5) - t262;
t293 = qJD(5) - t294;
t201 = -0.2e1 * qJD(6) * t283 + (t282 * t293 - t238) * qJ(6) + (t282 * t283 + t261) * pkin(5) + t205;
t264 = -mrSges(7,2) * t293 + mrSges(7,3) * t282;
t368 = m(7) * t201 + mrSges(7,1) * t261 + t264 * t293;
t198 = -t238 * mrSges(7,3) - t283 * t252 + t368;
t206 = t211 * t346 + t217 * t342;
t237 = -qJD(5) * t283 - t263 * t342 + t317 * t346;
t244 = Ifges(6,4) * t283 + Ifges(6,2) * t282 + Ifges(6,6) * t293;
t245 = Ifges(7,1) * t283 + Ifges(7,4) * t282 + Ifges(7,5) * t293;
t246 = Ifges(6,1) * t283 + Ifges(6,4) * t282 + Ifges(6,5) * t293;
t266 = pkin(5) * t293 - qJ(6) * t283;
t281 = t282 ^ 2;
t204 = -pkin(5) * t281 + qJ(6) * t237 + 0.2e1 * qJD(6) * t282 - t266 * t293 + t206;
t243 = Ifges(7,4) * t283 + Ifges(7,2) * t282 + Ifges(7,6) * t293;
t359 = -mrSges(7,1) * t201 + mrSges(7,2) * t204 - Ifges(7,5) * t238 - Ifges(7,6) * t237 - Ifges(7,3) * t261 - t243 * t283;
t382 = mrSges(6,1) * t205 - mrSges(6,2) * t206 + Ifges(6,5) * t238 + Ifges(6,6) * t237 + Ifges(6,3) * t261 + pkin(5) * t198 + t283 * t244 - (t246 + t245) * t282 - t359;
t379 = -mrSges(6,2) - mrSges(7,2);
t378 = t339 * t344;
t273 = -mrSges(5,1) * t294 + mrSges(5,2) * t295;
t285 = mrSges(5,1) * t328 - mrSges(5,3) * t295;
t253 = -mrSges(6,1) * t282 + mrSges(6,2) * t283;
t265 = -mrSges(6,2) * t293 + mrSges(6,3) * t282;
t190 = m(6) * t205 + t261 * mrSges(6,1) + t293 * t265 + (-t252 - t253) * t283 + (-mrSges(6,3) - mrSges(7,3)) * t238 + t368;
t367 = m(7) * t204 + mrSges(7,3) * t237 + t252 * t282;
t267 = mrSges(7,1) * t293 - mrSges(7,3) * t283;
t373 = -mrSges(6,1) * t293 + mrSges(6,3) * t283 - t267;
t195 = m(6) * t206 + mrSges(6,3) * t237 + t253 * t282 + t261 * t379 + t293 * t373 + t367;
t362 = -t190 * t342 + t195 * t346;
t183 = m(5) * t214 - mrSges(5,2) * t317 + mrSges(5,3) * t262 + t273 * t294 - t285 * t328 + t362;
t213 = t219 * t347 - t222 * t343;
t284 = -mrSges(5,2) * t328 + mrSges(5,3) * t294;
t210 = -pkin(4) * t317 - pkin(10) * t327 + t274 * t295 - t213;
t208 = -pkin(5) * t237 - qJ(6) * t281 + t266 * t283 + qJDD(6) + t210;
t361 = -m(7) * t208 + mrSges(7,1) * t237 + t264 * t282;
t354 = -m(6) * t210 + mrSges(6,1) * t237 + t238 * t379 + t265 * t282 + t283 * t373 + t361;
t192 = m(5) * t213 + mrSges(5,1) * t317 - mrSges(5,3) * t263 - t273 * t295 + t284 * t328 + t354;
t176 = t183 * t343 + t192 * t347;
t296 = -mrSges(4,1) * t312 + mrSges(4,2) * t313;
t298 = mrSges(4,2) * t365 + mrSges(4,3) * t312;
t174 = m(4) * t226 - mrSges(4,1) * t325 - mrSges(4,3) * t301 - t296 * t313 - t298 * t365 + t176;
t299 = -mrSges(4,1) * t365 - mrSges(4,3) * t313;
t363 = t183 * t347 - t192 * t343;
t175 = m(4) * t227 + mrSges(4,2) * t325 + mrSges(4,3) * t300 + t296 * t312 + t299 * t365 + t363;
t169 = t174 * t340 + t175 * t338;
t187 = t190 * t346 + t195 * t342;
t292 = -g(3) * t378 + t372;
t318 = mrSges(3,1) * t334 - mrSges(3,3) * t366;
t323 = (-mrSges(3,1) * t348 + mrSges(3,2) * t344) * t371;
t364 = -t174 * t338 + t175 * t340;
t167 = m(3) * t292 - mrSges(3,2) * t333 + mrSges(3,3) * t325 - t318 * t334 + t323 * t365 + t364;
t319 = -mrSges(3,2) * t334 + mrSges(3,3) * t365;
t356 = m(5) * t239 - mrSges(5,1) * t262 + mrSges(5,2) * t263 - t284 * t294 + t285 * t295 + t187;
t353 = -m(4) * t278 + mrSges(4,1) * t300 - mrSges(4,2) * t301 + t298 * t312 - t299 * t313 - t356;
t180 = m(3) * t291 + mrSges(3,1) * t333 - mrSges(3,3) * t324 + t319 * t334 - t323 * t366 + t353;
t163 = t167 * t348 - t180 * t344;
t306 = -t339 * t320 - t380;
t168 = m(3) * t306 - t325 * mrSges(3,1) + t324 * mrSges(3,2) + (t318 * t344 - t319 * t348) * t371 + t169;
t159 = t167 * t376 - t168 * t339 + t180 * t375;
t360 = -mrSges(7,1) * t208 + mrSges(7,3) * t204 + Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t261 + t245 * t293;
t241 = Ifges(7,5) * t283 + Ifges(7,6) * t282 + Ifges(7,3) * t293;
t358 = mrSges(7,2) * t208 - mrSges(7,3) * t201 + Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t261 + t241 * t282;
t242 = Ifges(6,5) * t283 + Ifges(6,6) * t282 + Ifges(6,3) * t293;
t178 = Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t261 + t293 * t246 - mrSges(6,1) * t210 + mrSges(6,3) * t206 - pkin(5) * (t238 * mrSges(7,2) - t361) + qJ(6) * (-mrSges(7,2) * t261 - t267 * t293 + t367) + (-pkin(5) * t267 - t241 - t242) * t283 + t360;
t185 = mrSges(6,2) * t210 - mrSges(6,3) * t205 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t261 - qJ(6) * t198 + t242 * t282 + (-t243 - t244) * t293 + t358;
t269 = Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t328;
t270 = Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t328;
t164 = mrSges(5,2) * t239 - mrSges(5,3) * t213 + Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t317 - pkin(10) * t187 - t178 * t342 + t185 * t346 + t269 * t294 - t270 * t328;
t271 = Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t328;
t170 = -mrSges(5,1) * t239 + mrSges(5,3) * t214 + Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * t317 - pkin(4) * t187 - t295 * t269 + t328 * t271 - t382;
t286 = Ifges(4,5) * t313 + Ifges(4,6) * t312 - Ifges(4,3) * t365;
t288 = Ifges(4,1) * t313 + Ifges(4,4) * t312 - Ifges(4,5) * t365;
t155 = -mrSges(4,1) * t278 + mrSges(4,3) * t227 + Ifges(4,4) * t301 + Ifges(4,2) * t300 - Ifges(4,6) * t325 - pkin(3) * t356 + pkin(9) * t363 + t343 * t164 + t347 * t170 - t313 * t286 - t288 * t365;
t287 = Ifges(4,4) * t313 + Ifges(4,2) * t312 - Ifges(4,6) * t365;
t160 = mrSges(4,2) * t278 - mrSges(4,3) * t226 + Ifges(4,1) * t301 + Ifges(4,4) * t300 - Ifges(4,5) * t325 - pkin(9) * t176 + t164 * t347 - t170 * t343 + t286 * t312 + t287 * t365;
t304 = Ifges(3,6) * t334 + (Ifges(3,4) * t344 + Ifges(3,2) * t348) * t371;
t305 = Ifges(3,5) * t334 + (Ifges(3,1) * t344 + Ifges(3,4) * t348) * t371;
t150 = Ifges(3,5) * t324 + Ifges(3,6) * t325 + Ifges(3,3) * t333 + mrSges(3,1) * t291 - mrSges(3,2) * t292 + t338 * t160 + t340 * t155 + pkin(2) * t353 + qJ(3) * t364 + (t304 * t344 - t305 * t348) * t371;
t303 = Ifges(3,3) * t334 + (Ifges(3,5) * t344 + Ifges(3,6) * t348) * t371;
t152 = mrSges(3,2) * t306 - mrSges(3,3) * t291 + Ifges(3,1) * t324 + Ifges(3,4) * t325 + Ifges(3,5) * t333 - qJ(3) * t169 - t155 * t338 + t160 * t340 + t303 * t365 - t304 * t334;
t355 = -mrSges(5,1) * t213 + mrSges(5,2) * t214 - Ifges(5,5) * t263 - Ifges(5,6) * t262 - Ifges(5,3) * t317 - pkin(4) * t354 - pkin(10) * t362 - t178 * t346 - t185 * t342 - t270 * t295 + t294 * t271;
t351 = -mrSges(4,1) * t226 + mrSges(4,2) * t227 - Ifges(4,5) * t301 - Ifges(4,6) * t300 - pkin(3) * t176 - t313 * t287 + t312 * t288 + t355;
t154 = t351 - t303 * t366 + Ifges(3,4) * t324 + Ifges(3,6) * t333 + t334 * t305 - pkin(2) * t169 - mrSges(3,1) * t306 + mrSges(3,3) * t292 + (Ifges(3,2) + Ifges(4,3)) * t325;
t357 = mrSges(2,1) * t329 - mrSges(2,2) * t330 + Ifges(2,3) * qJDD(1) + pkin(1) * t159 + t150 * t341 + t152 * t378 + t154 * t377 + t163 * t381;
t161 = m(2) * t330 - mrSges(2,1) * t350 - qJDD(1) * mrSges(2,2) + t163;
t158 = t341 * t168 + (t167 * t344 + t180 * t348) * t339;
t156 = m(2) * t329 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t350 + t159;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t329 + Ifges(2,5) * qJDD(1) - t350 * Ifges(2,6) + t348 * t152 - t344 * t154 + (-t158 * t339 - t159 * t341) * pkin(8);
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t330 + t350 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t158 - t339 * t150 + (pkin(8) * t163 + t152 * t344 + t154 * t348) * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t349 * t148 - t345 * t147 - pkin(7) * (t156 * t349 + t161 * t345), t148, t152, t160, t164, t185, -t243 * t293 + t358; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t345 * t148 + t349 * t147 + pkin(7) * (-t156 * t345 + t161 * t349), t147, t154, t155, t170, t178, -t283 * t241 + t360; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t357, t357, t150, -Ifges(4,3) * t325 - t351, -t355, t382, -t282 * t245 - t359;];
m_new  = t1;
