% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 09:22:16
% EndTime: 2019-05-08 09:24:53
% DurationCPUTime: 64.02s
% Computational Cost: add. (1170246->385), mult. (2429895->480), div. (0->0), fcn. (1825294->12), ass. (0->154)
t333 = sin(qJ(1));
t339 = cos(qJ(1));
t320 = t333 * g(1) - t339 * g(2);
t341 = qJD(1) ^ 2;
t302 = -qJDD(1) * pkin(1) - t341 * pkin(7) - t320;
t332 = sin(qJ(2));
t338 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t357 = t338 * t358;
t313 = t332 * qJDD(1) + t357;
t324 = t332 * t358;
t314 = t338 * qJDD(1) - t324;
t267 = (-t313 - t357) * pkin(8) + (-t314 + t324) * pkin(2) + t302;
t321 = -t339 * g(1) - t333 * g(2);
t303 = -t341 * pkin(1) + qJDD(1) * pkin(7) + t321;
t289 = -t332 * g(3) + t338 * t303;
t312 = (-pkin(2) * t338 - pkin(8) * t332) * qJD(1);
t340 = qJD(2) ^ 2;
t359 = t338 * qJD(1);
t270 = -t340 * pkin(2) + qJDD(2) * pkin(8) + t312 * t359 + t289;
t331 = sin(qJ(3));
t337 = cos(qJ(3));
t251 = t337 * t267 - t331 * t270;
t360 = qJD(1) * t332;
t309 = t337 * qJD(2) - t331 * t360;
t281 = t309 * qJD(3) + t331 * qJDD(2) + t337 * t313;
t308 = qJDD(3) - t314;
t310 = t331 * qJD(2) + t337 * t360;
t323 = qJD(3) - t359;
t231 = (t309 * t323 - t281) * pkin(9) + (t309 * t310 + t308) * pkin(3) + t251;
t252 = t331 * t267 + t337 * t270;
t280 = -t310 * qJD(3) + t337 * qJDD(2) - t331 * t313;
t290 = t323 * pkin(3) - t310 * pkin(9);
t307 = t309 ^ 2;
t233 = -t307 * pkin(3) + t280 * pkin(9) - t323 * t290 + t252;
t330 = sin(qJ(4));
t336 = cos(qJ(4));
t212 = t336 * t231 - t330 * t233;
t283 = t336 * t309 - t330 * t310;
t250 = t283 * qJD(4) + t330 * t280 + t336 * t281;
t284 = t330 * t309 + t336 * t310;
t304 = qJDD(4) + t308;
t322 = qJD(4) + t323;
t207 = (t283 * t322 - t250) * pkin(10) + (t283 * t284 + t304) * pkin(4) + t212;
t213 = t330 * t231 + t336 * t233;
t249 = -t284 * qJD(4) + t336 * t280 - t330 * t281;
t273 = t322 * pkin(4) - t284 * pkin(10);
t282 = t283 ^ 2;
t209 = -t282 * pkin(4) + t249 * pkin(10) - t322 * t273 + t213;
t329 = sin(qJ(5));
t335 = cos(qJ(5));
t195 = t335 * t207 - t329 * t209;
t262 = t335 * t283 - t329 * t284;
t227 = t262 * qJD(5) + t329 * t249 + t335 * t250;
t263 = t329 * t283 + t335 * t284;
t298 = qJDD(5) + t304;
t317 = qJD(5) + t322;
t192 = (t262 * t317 - t227) * pkin(11) + (t262 * t263 + t298) * pkin(5) + t195;
t196 = t329 * t207 + t335 * t209;
t226 = -t263 * qJD(5) + t335 * t249 - t329 * t250;
t255 = t317 * pkin(5) - t263 * pkin(11);
t261 = t262 ^ 2;
t193 = -t261 * pkin(5) + t226 * pkin(11) - t317 * t255 + t196;
t328 = sin(qJ(6));
t334 = cos(qJ(6));
t191 = t328 * t192 + t334 * t193;
t288 = -t338 * g(3) - t332 * t303;
t269 = -qJDD(2) * pkin(2) - t340 * pkin(8) + t312 * t360 - t288;
t244 = -t280 * pkin(3) - t307 * pkin(9) + t310 * t290 + t269;
t215 = -t249 * pkin(4) - t282 * pkin(10) + t284 * t273 + t244;
t198 = -t226 * pkin(5) - t261 * pkin(11) + t263 * t255 + t215;
t242 = t328 * t262 + t334 * t263;
t203 = -t242 * qJD(6) + t334 * t226 - t328 * t227;
t241 = t334 * t262 - t328 * t263;
t204 = t241 * qJD(6) + t328 * t226 + t334 * t227;
t315 = qJD(6) + t317;
t216 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t315;
t218 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t315;
t295 = qJDD(6) + t298;
t180 = -mrSges(7,1) * t198 + mrSges(7,3) * t191 + Ifges(7,4) * t204 + Ifges(7,2) * t203 + Ifges(7,6) * t295 - t242 * t216 + t315 * t218;
t190 = t334 * t192 - t328 * t193;
t217 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t315;
t181 = mrSges(7,2) * t198 - mrSges(7,3) * t190 + Ifges(7,1) * t204 + Ifges(7,4) * t203 + Ifges(7,5) * t295 + t241 * t216 - t315 * t217;
t236 = Ifges(6,5) * t263 + Ifges(6,6) * t262 + Ifges(6,3) * t317;
t238 = Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t317;
t234 = -t315 * mrSges(7,2) + t241 * mrSges(7,3);
t235 = t315 * mrSges(7,1) - t242 * mrSges(7,3);
t352 = m(7) * t198 - t203 * mrSges(7,1) + t204 * mrSges(7,2) - t241 * t234 + t242 * t235;
t221 = -t241 * mrSges(7,1) + t242 * mrSges(7,2);
t187 = m(7) * t190 + t295 * mrSges(7,1) - t204 * mrSges(7,3) - t242 * t221 + t315 * t234;
t188 = m(7) * t191 - t295 * mrSges(7,2) + t203 * mrSges(7,3) + t241 * t221 - t315 * t235;
t353 = -t328 * t187 + t334 * t188;
t164 = -mrSges(6,1) * t215 + mrSges(6,3) * t196 + Ifges(6,4) * t227 + Ifges(6,2) * t226 + Ifges(6,6) * t298 - pkin(5) * t352 + pkin(11) * t353 + t334 * t180 + t328 * t181 - t263 * t236 + t317 * t238;
t179 = t334 * t187 + t328 * t188;
t237 = Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t317;
t165 = mrSges(6,2) * t215 - mrSges(6,3) * t195 + Ifges(6,1) * t227 + Ifges(6,4) * t226 + Ifges(6,5) * t298 - pkin(11) * t179 - t328 * t180 + t334 * t181 + t262 * t236 - t317 * t237;
t256 = Ifges(5,5) * t284 + Ifges(5,6) * t283 + Ifges(5,3) * t322;
t258 = Ifges(5,1) * t284 + Ifges(5,4) * t283 + Ifges(5,5) * t322;
t253 = -t317 * mrSges(6,2) + t262 * mrSges(6,3);
t254 = t317 * mrSges(6,1) - t263 * mrSges(6,3);
t349 = m(6) * t215 - t226 * mrSges(6,1) + t227 * mrSges(6,2) - t262 * t253 + t263 * t254 + t352;
t243 = -t262 * mrSges(6,1) + t263 * mrSges(6,2);
t176 = m(6) * t195 + t298 * mrSges(6,1) - t227 * mrSges(6,3) - t263 * t243 + t317 * t253 + t179;
t177 = m(6) * t196 - t298 * mrSges(6,2) + t226 * mrSges(6,3) + t262 * t243 - t317 * t254 + t353;
t354 = -t329 * t176 + t335 * t177;
t158 = -mrSges(5,1) * t244 + mrSges(5,3) * t213 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * t304 - pkin(4) * t349 + pkin(10) * t354 + t335 * t164 + t329 * t165 - t284 * t256 + t322 * t258;
t172 = t335 * t176 + t329 * t177;
t257 = Ifges(5,4) * t284 + Ifges(5,2) * t283 + Ifges(5,6) * t322;
t159 = mrSges(5,2) * t244 - mrSges(5,3) * t212 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * t304 - pkin(10) * t172 - t329 * t164 + t335 * t165 + t283 * t256 - t322 * t257;
t274 = Ifges(4,5) * t310 + Ifges(4,6) * t309 + Ifges(4,3) * t323;
t276 = Ifges(4,1) * t310 + Ifges(4,4) * t309 + Ifges(4,5) * t323;
t271 = -t322 * mrSges(5,2) + t283 * mrSges(5,3);
t272 = t322 * mrSges(5,1) - t284 * mrSges(5,3);
t346 = m(5) * t244 - t249 * mrSges(5,1) + t250 * mrSges(5,2) - t283 * t271 + t284 * t272 + t349;
t264 = -t283 * mrSges(5,1) + t284 * mrSges(5,2);
t169 = m(5) * t212 + t304 * mrSges(5,1) - t250 * mrSges(5,3) - t284 * t264 + t322 * t271 + t172;
t170 = m(5) * t213 - t304 * mrSges(5,2) + t249 * mrSges(5,3) + t283 * t264 - t322 * t272 + t354;
t355 = -t330 * t169 + t336 * t170;
t145 = -mrSges(4,1) * t269 + mrSges(4,3) * t252 + Ifges(4,4) * t281 + Ifges(4,2) * t280 + Ifges(4,6) * t308 - pkin(3) * t346 + pkin(9) * t355 + t336 * t158 + t330 * t159 - t310 * t274 + t323 * t276;
t163 = t336 * t169 + t330 * t170;
t275 = Ifges(4,4) * t310 + Ifges(4,2) * t309 + Ifges(4,6) * t323;
t146 = mrSges(4,2) * t269 - mrSges(4,3) * t251 + Ifges(4,1) * t281 + Ifges(4,4) * t280 + Ifges(4,5) * t308 - pkin(9) * t163 - t330 * t158 + t336 * t159 + t309 * t274 - t323 * t275;
t285 = -t309 * mrSges(4,1) + t310 * mrSges(4,2);
t286 = -t323 * mrSges(4,2) + t309 * mrSges(4,3);
t161 = m(4) * t251 + t308 * mrSges(4,1) - t281 * mrSges(4,3) - t310 * t285 + t323 * t286 + t163;
t287 = t323 * mrSges(4,1) - t310 * mrSges(4,3);
t162 = m(4) * t252 - t308 * mrSges(4,2) + t280 * mrSges(4,3) + t309 * t285 - t323 * t287 + t355;
t157 = -t331 * t161 + t337 * t162;
t183 = -m(4) * t269 + t280 * mrSges(4,1) - t281 * mrSges(4,2) + t309 * t286 - t310 * t287 - t346;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t332 + Ifges(3,2) * t338) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t332 + Ifges(3,4) * t338) * qJD(1);
t361 = mrSges(3,1) * t288 - mrSges(3,2) * t289 + Ifges(3,5) * t313 + Ifges(3,6) * t314 + Ifges(3,3) * qJDD(2) + pkin(2) * t183 + pkin(8) * t157 + t337 * t145 + t331 * t146 + (t332 * t300 - t338 * t301) * qJD(1);
t311 = (-mrSges(3,1) * t338 + mrSges(3,2) * t332) * qJD(1);
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t360;
t155 = m(3) * t289 - qJDD(2) * mrSges(3,2) + t314 * mrSges(3,3) - qJD(2) * t318 + t311 * t359 + t157;
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t182 = m(3) * t288 + qJDD(2) * mrSges(3,1) - t313 * mrSges(3,3) + qJD(2) * t319 - t311 * t360 + t183;
t356 = t338 * t155 - t332 * t182;
t156 = t337 * t161 + t331 * t162;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t332 + Ifges(3,6) * t338) * qJD(1);
t144 = mrSges(3,2) * t302 - mrSges(3,3) * t288 + Ifges(3,1) * t313 + Ifges(3,4) * t314 + Ifges(3,5) * qJDD(2) - pkin(8) * t156 - qJD(2) * t300 - t331 * t145 + t337 * t146 + t299 * t359;
t348 = -mrSges(7,1) * t190 + mrSges(7,2) * t191 - Ifges(7,5) * t204 - Ifges(7,6) * t203 - Ifges(7,3) * t295 - t242 * t217 + t241 * t218;
t345 = -mrSges(6,1) * t195 + mrSges(6,2) * t196 - Ifges(6,5) * t227 - Ifges(6,6) * t226 - Ifges(6,3) * t298 - pkin(5) * t179 - t263 * t237 + t262 * t238 + t348;
t343 = -mrSges(5,1) * t212 + mrSges(5,2) * t213 - Ifges(5,5) * t250 - Ifges(5,6) * t249 - Ifges(5,3) * t304 - pkin(4) * t172 - t284 * t257 + t283 * t258 + t345;
t342 = mrSges(4,1) * t251 - mrSges(4,2) * t252 + Ifges(4,5) * t281 + Ifges(4,6) * t280 + Ifges(4,3) * t308 + pkin(3) * t163 + t310 * t275 - t309 * t276 - t343;
t148 = -mrSges(3,1) * t302 + mrSges(3,3) * t289 + Ifges(3,4) * t313 + Ifges(3,2) * t314 + Ifges(3,6) * qJDD(2) - pkin(2) * t156 + qJD(2) * t301 - t299 * t360 - t342;
t347 = -m(3) * t302 + t314 * mrSges(3,1) - t313 * mrSges(3,2) - t318 * t360 + t319 * t359 - t156;
t350 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t347 + pkin(7) * t356 + t332 * t144 + t338 * t148;
t152 = m(2) * t320 + qJDD(1) * mrSges(2,1) - t341 * mrSges(2,2) + t347;
t151 = t332 * t155 + t338 * t182;
t149 = m(2) * t321 - t341 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t356;
t142 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t341 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t361;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - t341 * Ifges(2,6) - pkin(7) * t151 + t338 * t144 - t332 * t148;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t141 - t333 * t142 - pkin(6) * (t333 * t149 + t339 * t152), t141, t144, t146, t159, t165, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t333 * t141 + t339 * t142 + pkin(6) * (t339 * t149 - t333 * t152), t142, t148, t145, t158, t164, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t361, t342, -t343, -t345, -t348;];
m_new  = t1;
