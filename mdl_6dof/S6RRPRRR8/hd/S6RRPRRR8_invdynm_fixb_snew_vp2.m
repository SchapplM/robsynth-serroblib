% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 22:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:37:27
% EndTime: 2019-05-06 22:38:50
% DurationCPUTime: 57.20s
% Computational Cost: add. (1028928->385), mult. (2268015->482), div. (0->0), fcn. (1693660->12), ass. (0->152)
t334 = sin(qJ(1));
t339 = cos(qJ(1));
t320 = t334 * g(1) - t339 * g(2);
t341 = qJD(1) ^ 2;
t302 = -qJDD(1) * pkin(1) - t341 * pkin(7) - t320;
t333 = sin(qJ(2));
t338 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t357 = t338 * t358;
t314 = qJDD(1) * t333 + t357;
t324 = t333 * t358;
t315 = qJDD(1) * t338 - t324;
t267 = (-t314 - t357) * qJ(3) + (-t315 + t324) * pkin(2) + t302;
t321 = -g(1) * t339 - g(2) * t334;
t303 = -pkin(1) * t341 + qJDD(1) * pkin(7) + t321;
t284 = -g(3) * t333 + t338 * t303;
t312 = (-pkin(2) * t338 - qJ(3) * t333) * qJD(1);
t340 = qJD(2) ^ 2;
t359 = qJD(1) * t338;
t270 = -pkin(2) * t340 + qJDD(2) * qJ(3) + t312 * t359 + t284;
t328 = sin(pkin(11));
t329 = cos(pkin(11));
t360 = qJD(1) * t333;
t309 = qJD(2) * t328 + t329 * t360;
t244 = -0.2e1 * qJD(3) * t309 + t329 * t267 - t328 * t270;
t289 = qJDD(2) * t328 + t314 * t329;
t308 = qJD(2) * t329 - t328 * t360;
t231 = (-t308 * t359 - t289) * pkin(8) + (t308 * t309 - t315) * pkin(3) + t244;
t245 = 0.2e1 * qJD(3) * t308 + t328 * t267 + t329 * t270;
t288 = qJDD(2) * t329 - t314 * t328;
t290 = -pkin(3) * t359 - pkin(8) * t309;
t307 = t308 ^ 2;
t233 = -pkin(3) * t307 + pkin(8) * t288 + t290 * t359 + t245;
t332 = sin(qJ(4));
t337 = cos(qJ(4));
t212 = t337 * t231 - t332 * t233;
t280 = t308 * t337 - t309 * t332;
t255 = qJD(4) * t280 + t288 * t332 + t289 * t337;
t281 = t308 * t332 + t309 * t337;
t311 = qJDD(4) - t315;
t323 = qJD(4) - t359;
t207 = (t280 * t323 - t255) * pkin(9) + (t280 * t281 + t311) * pkin(4) + t212;
t213 = t332 * t231 + t337 * t233;
t254 = -qJD(4) * t281 + t288 * t337 - t289 * t332;
t273 = pkin(4) * t323 - pkin(9) * t281;
t279 = t280 ^ 2;
t209 = -pkin(4) * t279 + pkin(9) * t254 - t273 * t323 + t213;
t331 = sin(qJ(5));
t336 = cos(qJ(5));
t195 = t336 * t207 - t331 * t209;
t262 = t280 * t336 - t281 * t331;
t227 = qJD(5) * t262 + t254 * t331 + t255 * t336;
t263 = t280 * t331 + t281 * t336;
t305 = qJDD(5) + t311;
t322 = qJD(5) + t323;
t192 = (t262 * t322 - t227) * pkin(10) + (t262 * t263 + t305) * pkin(5) + t195;
t196 = t331 * t207 + t336 * t209;
t226 = -qJD(5) * t263 + t254 * t336 - t255 * t331;
t253 = pkin(5) * t322 - pkin(10) * t263;
t261 = t262 ^ 2;
t193 = -pkin(5) * t261 + pkin(10) * t226 - t253 * t322 + t196;
t330 = sin(qJ(6));
t335 = cos(qJ(6));
t191 = t192 * t330 + t193 * t335;
t283 = -t338 * g(3) - t333 * t303;
t269 = -qJDD(2) * pkin(2) - qJ(3) * t340 + t312 * t360 + qJDD(3) - t283;
t246 = -pkin(3) * t288 - pkin(8) * t307 + t309 * t290 + t269;
t220 = -pkin(4) * t254 - pkin(9) * t279 + t281 * t273 + t246;
t198 = -pkin(5) * t226 - pkin(10) * t261 + t253 * t263 + t220;
t242 = t262 * t330 + t263 * t335;
t203 = -qJD(6) * t242 + t226 * t335 - t227 * t330;
t241 = t262 * t335 - t263 * t330;
t204 = qJD(6) * t241 + t226 * t330 + t227 * t335;
t317 = qJD(6) + t322;
t214 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t317;
t216 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t317;
t296 = qJDD(6) + t305;
t180 = -mrSges(7,1) * t198 + mrSges(7,3) * t191 + Ifges(7,4) * t204 + Ifges(7,2) * t203 + Ifges(7,6) * t296 - t214 * t242 + t216 * t317;
t190 = t192 * t335 - t193 * t330;
t215 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t317;
t181 = mrSges(7,2) * t198 - mrSges(7,3) * t190 + Ifges(7,1) * t204 + Ifges(7,4) * t203 + Ifges(7,5) * t296 + t214 * t241 - t215 * t317;
t236 = Ifges(6,5) * t263 + Ifges(6,6) * t262 + Ifges(6,3) * t322;
t238 = Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t322;
t234 = -mrSges(7,2) * t317 + mrSges(7,3) * t241;
t235 = mrSges(7,1) * t317 - mrSges(7,3) * t242;
t352 = m(7) * t198 - t203 * mrSges(7,1) + t204 * mrSges(7,2) - t241 * t234 + t242 * t235;
t221 = -mrSges(7,1) * t241 + mrSges(7,2) * t242;
t185 = m(7) * t190 + mrSges(7,1) * t296 - mrSges(7,3) * t204 - t221 * t242 + t234 * t317;
t186 = m(7) * t191 - mrSges(7,2) * t296 + mrSges(7,3) * t203 + t221 * t241 - t235 * t317;
t353 = -t185 * t330 + t335 * t186;
t164 = -mrSges(6,1) * t220 + mrSges(6,3) * t196 + Ifges(6,4) * t227 + Ifges(6,2) * t226 + Ifges(6,6) * t305 - pkin(5) * t352 + pkin(10) * t353 + t335 * t180 + t330 * t181 - t263 * t236 + t322 * t238;
t179 = t335 * t185 + t330 * t186;
t237 = Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t322;
t165 = mrSges(6,2) * t220 - mrSges(6,3) * t195 + Ifges(6,1) * t227 + Ifges(6,4) * t226 + Ifges(6,5) * t305 - pkin(10) * t179 - t180 * t330 + t181 * t335 + t236 * t262 - t237 * t322;
t256 = Ifges(5,5) * t281 + Ifges(5,6) * t280 + Ifges(5,3) * t323;
t258 = Ifges(5,1) * t281 + Ifges(5,4) * t280 + Ifges(5,5) * t323;
t248 = -mrSges(6,2) * t322 + mrSges(6,3) * t262;
t249 = mrSges(6,1) * t322 - mrSges(6,3) * t263;
t349 = m(6) * t220 - t226 * mrSges(6,1) + t227 * mrSges(6,2) - t262 * t248 + t263 * t249 + t352;
t243 = -mrSges(6,1) * t262 + mrSges(6,2) * t263;
t176 = m(6) * t195 + mrSges(6,1) * t305 - mrSges(6,3) * t227 - t243 * t263 + t248 * t322 + t179;
t177 = m(6) * t196 - mrSges(6,2) * t305 + mrSges(6,3) * t226 + t243 * t262 - t249 * t322 + t353;
t354 = -t176 * t331 + t336 * t177;
t158 = -mrSges(5,1) * t246 + mrSges(5,3) * t213 + Ifges(5,4) * t255 + Ifges(5,2) * t254 + Ifges(5,6) * t311 - pkin(4) * t349 + pkin(9) * t354 + t336 * t164 + t331 * t165 - t281 * t256 + t323 * t258;
t172 = t336 * t176 + t331 * t177;
t257 = Ifges(5,4) * t281 + Ifges(5,2) * t280 + Ifges(5,6) * t323;
t159 = mrSges(5,2) * t246 - mrSges(5,3) * t212 + Ifges(5,1) * t255 + Ifges(5,4) * t254 + Ifges(5,5) * t311 - pkin(9) * t172 - t164 * t331 + t165 * t336 + t256 * t280 - t257 * t323;
t274 = Ifges(4,5) * t309 + Ifges(4,6) * t308 - Ifges(4,3) * t359;
t276 = Ifges(4,1) * t309 + Ifges(4,4) * t308 - Ifges(4,5) * t359;
t271 = -mrSges(5,2) * t323 + mrSges(5,3) * t280;
t272 = mrSges(5,1) * t323 - mrSges(5,3) * t281;
t346 = m(5) * t246 - t254 * mrSges(5,1) + t255 * mrSges(5,2) - t280 * t271 + t281 * t272 + t349;
t264 = -mrSges(5,1) * t280 + mrSges(5,2) * t281;
t169 = m(5) * t212 + mrSges(5,1) * t311 - mrSges(5,3) * t255 - t264 * t281 + t271 * t323 + t172;
t170 = m(5) * t213 - mrSges(5,2) * t311 + mrSges(5,3) * t254 + t264 * t280 - t272 * t323 + t354;
t355 = -t169 * t332 + t337 * t170;
t145 = -mrSges(4,1) * t269 + mrSges(4,3) * t245 + Ifges(4,4) * t289 + Ifges(4,2) * t288 - Ifges(4,6) * t315 - pkin(3) * t346 + pkin(8) * t355 + t337 * t158 + t332 * t159 - t309 * t274 - t276 * t359;
t163 = t337 * t169 + t332 * t170;
t275 = Ifges(4,4) * t309 + Ifges(4,2) * t308 - Ifges(4,6) * t359;
t146 = mrSges(4,2) * t269 - mrSges(4,3) * t244 + Ifges(4,1) * t289 + Ifges(4,4) * t288 - Ifges(4,5) * t315 - pkin(8) * t163 - t158 * t332 + t159 * t337 + t274 * t308 + t275 * t359;
t282 = -mrSges(4,1) * t308 + mrSges(4,2) * t309;
t286 = mrSges(4,2) * t359 + mrSges(4,3) * t308;
t161 = m(4) * t244 - mrSges(4,1) * t315 - mrSges(4,3) * t289 - t282 * t309 - t286 * t359 + t163;
t287 = -mrSges(4,1) * t359 - mrSges(4,3) * t309;
t162 = m(4) * t245 + mrSges(4,2) * t315 + mrSges(4,3) * t288 + t282 * t308 + t287 * t359 + t355;
t157 = -t161 * t328 + t329 * t162;
t188 = -m(4) * t269 + t288 * mrSges(4,1) - t289 * mrSges(4,2) + t308 * t286 - t309 * t287 - t346;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t333 + Ifges(3,2) * t338) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t333 + Ifges(3,4) * t338) * qJD(1);
t361 = mrSges(3,1) * t283 - mrSges(3,2) * t284 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t188 + qJ(3) * t157 + t329 * t145 + t328 * t146 + (t300 * t333 - t301 * t338) * qJD(1);
t313 = (-mrSges(3,1) * t338 + mrSges(3,2) * t333) * qJD(1);
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t360;
t155 = m(3) * t284 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t318 + t313 * t359 + t157;
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t187 = m(3) * t283 + qJDD(2) * mrSges(3,1) - t314 * mrSges(3,3) + qJD(2) * t319 - t313 * t360 + t188;
t356 = t338 * t155 - t187 * t333;
t156 = t161 * t329 + t162 * t328;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t333 + Ifges(3,6) * t338) * qJD(1);
t144 = mrSges(3,2) * t302 - mrSges(3,3) * t283 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - qJ(3) * t156 - qJD(2) * t300 - t145 * t328 + t146 * t329 + t299 * t359;
t348 = -mrSges(7,1) * t190 + mrSges(7,2) * t191 - Ifges(7,5) * t204 - Ifges(7,6) * t203 - Ifges(7,3) * t296 - t242 * t215 + t241 * t216;
t345 = -mrSges(6,1) * t195 + mrSges(6,2) * t196 - Ifges(6,5) * t227 - Ifges(6,6) * t226 - Ifges(6,3) * t305 - pkin(5) * t179 - t263 * t237 + t262 * t238 + t348;
t343 = -mrSges(5,1) * t212 + mrSges(5,2) * t213 - Ifges(5,5) * t255 - Ifges(5,6) * t254 - Ifges(5,3) * t311 - pkin(4) * t172 - t281 * t257 + t280 * t258 + t345;
t342 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t289 + Ifges(4,6) * t288 + pkin(3) * t163 + t309 * t275 - t308 * t276 - t343;
t148 = -t299 * t360 - t342 + (Ifges(3,2) + Ifges(4,3)) * t315 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t314 + qJD(2) * t301 - mrSges(3,1) * t302 + mrSges(3,3) * t284 - pkin(2) * t156;
t347 = -m(3) * t302 + t315 * mrSges(3,1) - mrSges(3,2) * t314 - t318 * t360 + t319 * t359 - t156;
t350 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t347 + pkin(7) * t356 + t333 * t144 + t338 * t148;
t152 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t341 + t347;
t151 = t155 * t333 + t187 * t338;
t149 = m(2) * t321 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t356;
t142 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t341 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t361;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t341 - pkin(7) * t151 + t144 * t338 - t148 * t333;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t141 - t334 * t142 - pkin(6) * (t149 * t334 + t152 * t339), t141, t144, t146, t159, t165, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t141 + t339 * t142 + pkin(6) * (t149 * t339 - t152 * t334), t142, t148, t145, t158, t164, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t361, -Ifges(4,3) * t315 + t342, -t343, -t345, -t348;];
m_new  = t1;
