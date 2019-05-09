% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:56:46
% EndTime: 2019-05-07 11:00:07
% DurationCPUTime: 60.50s
% Computational Cost: add. (1084947->384), mult. (2302590->482), div. (0->0), fcn. (1706118->12), ass. (0->152)
t334 = sin(qJ(1));
t339 = cos(qJ(1));
t320 = t334 * g(1) - t339 * g(2);
t341 = qJD(1) ^ 2;
t302 = -qJDD(1) * pkin(1) - t341 * pkin(7) - t320;
t333 = sin(qJ(2));
t338 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t357 = t338 * t358;
t314 = t333 * qJDD(1) + t357;
t324 = t333 * t358;
t315 = t338 * qJDD(1) - t324;
t267 = (-t314 - t357) * pkin(8) + (-t315 + t324) * pkin(2) + t302;
t321 = -t339 * g(1) - t334 * g(2);
t303 = -t341 * pkin(1) + qJDD(1) * pkin(7) + t321;
t292 = -t333 * g(3) + t338 * t303;
t313 = (-pkin(2) * t338 - pkin(8) * t333) * qJD(1);
t340 = qJD(2) ^ 2;
t359 = t338 * qJD(1);
t270 = -t340 * pkin(2) + qJDD(2) * pkin(8) + t313 * t359 + t292;
t332 = sin(qJ(3));
t337 = cos(qJ(3));
t246 = t337 * t267 - t332 * t270;
t360 = qJD(1) * t333;
t310 = t337 * qJD(2) - t332 * t360;
t283 = t310 * qJD(3) + t332 * qJDD(2) + t337 * t314;
t309 = qJDD(3) - t315;
t311 = t332 * qJD(2) + t337 * t360;
t323 = qJD(3) - t359;
t231 = (t310 * t323 - t283) * qJ(4) + (t310 * t311 + t309) * pkin(3) + t246;
t247 = t332 * t267 + t337 * t270;
t282 = -t311 * qJD(3) + t337 * qJDD(2) - t332 * t314;
t289 = t323 * pkin(3) - t311 * qJ(4);
t308 = t310 ^ 2;
t233 = -t308 * pkin(3) + t282 * qJ(4) - t323 * t289 + t247;
t328 = sin(pkin(11));
t329 = cos(pkin(11));
t286 = t328 * t310 + t329 * t311;
t212 = -0.2e1 * qJD(4) * t286 + t329 * t231 - t328 * t233;
t258 = t328 * t282 + t329 * t283;
t285 = t329 * t310 - t328 * t311;
t201 = (t285 * t323 - t258) * pkin(9) + (t285 * t286 + t309) * pkin(4) + t212;
t213 = 0.2e1 * qJD(4) * t285 + t328 * t231 + t329 * t233;
t257 = t329 * t282 - t328 * t283;
t273 = t323 * pkin(4) - t286 * pkin(9);
t284 = t285 ^ 2;
t209 = -t284 * pkin(4) + t257 * pkin(9) - t323 * t273 + t213;
t331 = sin(qJ(5));
t336 = cos(qJ(5));
t195 = t336 * t201 - t331 * t209;
t262 = t336 * t285 - t331 * t286;
t227 = t262 * qJD(5) + t331 * t257 + t336 * t258;
t263 = t331 * t285 + t336 * t286;
t305 = qJDD(5) + t309;
t322 = qJD(5) + t323;
t192 = (t262 * t322 - t227) * pkin(10) + (t262 * t263 + t305) * pkin(5) + t195;
t196 = t331 * t201 + t336 * t209;
t226 = -t263 * qJD(5) + t336 * t257 - t331 * t258;
t250 = t322 * pkin(5) - t263 * pkin(10);
t261 = t262 ^ 2;
t193 = -t261 * pkin(5) + t226 * pkin(10) - t322 * t250 + t196;
t330 = sin(qJ(6));
t335 = cos(qJ(6));
t191 = t330 * t192 + t335 * t193;
t291 = -t338 * g(3) - t333 * t303;
t269 = -qJDD(2) * pkin(2) - t340 * pkin(8) + t313 * t360 - t291;
t244 = -t282 * pkin(3) - t308 * qJ(4) + t311 * t289 + qJDD(4) + t269;
t215 = -t257 * pkin(4) - t284 * pkin(9) + t286 * t273 + t244;
t198 = -t226 * pkin(5) - t261 * pkin(10) + t263 * t250 + t215;
t242 = t330 * t262 + t335 * t263;
t206 = -t242 * qJD(6) + t335 * t226 - t330 * t227;
t241 = t335 * t262 - t330 * t263;
t207 = t241 * qJD(6) + t330 * t226 + t335 * t227;
t317 = qJD(6) + t322;
t216 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t317;
t218 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t317;
t298 = qJDD(6) + t305;
t180 = -mrSges(7,1) * t198 + mrSges(7,3) * t191 + Ifges(7,4) * t207 + Ifges(7,2) * t206 + Ifges(7,6) * t298 - t242 * t216 + t317 * t218;
t190 = t335 * t192 - t330 * t193;
t217 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t317;
t181 = mrSges(7,2) * t198 - mrSges(7,3) * t190 + Ifges(7,1) * t207 + Ifges(7,4) * t206 + Ifges(7,5) * t298 + t241 * t216 - t317 * t217;
t236 = Ifges(6,5) * t263 + Ifges(6,6) * t262 + Ifges(6,3) * t322;
t238 = Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t322;
t234 = -t317 * mrSges(7,2) + t241 * mrSges(7,3);
t235 = t317 * mrSges(7,1) - t242 * mrSges(7,3);
t352 = m(7) * t198 - t206 * mrSges(7,1) + t207 * mrSges(7,2) - t241 * t234 + t242 * t235;
t221 = -t241 * mrSges(7,1) + t242 * mrSges(7,2);
t185 = m(7) * t190 + t298 * mrSges(7,1) - t207 * mrSges(7,3) - t242 * t221 + t317 * t234;
t186 = m(7) * t191 - t298 * mrSges(7,2) + t206 * mrSges(7,3) + t241 * t221 - t317 * t235;
t353 = -t330 * t185 + t335 * t186;
t164 = -mrSges(6,1) * t215 + mrSges(6,3) * t196 + Ifges(6,4) * t227 + Ifges(6,2) * t226 + Ifges(6,6) * t305 - pkin(5) * t352 + pkin(10) * t353 + t335 * t180 + t330 * t181 - t263 * t236 + t322 * t238;
t179 = t335 * t185 + t330 * t186;
t237 = Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t322;
t165 = mrSges(6,2) * t215 - mrSges(6,3) * t195 + Ifges(6,1) * t227 + Ifges(6,4) * t226 + Ifges(6,5) * t305 - pkin(10) * t179 - t330 * t180 + t335 * t181 + t262 * t236 - t322 * t237;
t254 = Ifges(5,5) * t286 + Ifges(5,6) * t285 + Ifges(5,3) * t323;
t256 = Ifges(5,1) * t286 + Ifges(5,4) * t285 + Ifges(5,5) * t323;
t248 = -t322 * mrSges(6,2) + t262 * mrSges(6,3);
t249 = t322 * mrSges(6,1) - t263 * mrSges(6,3);
t349 = m(6) * t215 - t226 * mrSges(6,1) + t227 * mrSges(6,2) - t262 * t248 + t263 * t249 + t352;
t243 = -t262 * mrSges(6,1) + t263 * mrSges(6,2);
t176 = m(6) * t195 + t305 * mrSges(6,1) - t227 * mrSges(6,3) - t263 * t243 + t322 * t248 + t179;
t177 = m(6) * t196 - t305 * mrSges(6,2) + t226 * mrSges(6,3) + t262 * t243 - t322 * t249 + t353;
t354 = -t331 * t176 + t336 * t177;
t158 = -mrSges(5,1) * t244 + mrSges(5,3) * t213 + Ifges(5,4) * t258 + Ifges(5,2) * t257 + Ifges(5,6) * t309 - pkin(4) * t349 + pkin(9) * t354 + t336 * t164 + t331 * t165 - t286 * t254 + t323 * t256;
t172 = t336 * t176 + t331 * t177;
t255 = Ifges(5,4) * t286 + Ifges(5,2) * t285 + Ifges(5,6) * t323;
t159 = mrSges(5,2) * t244 - mrSges(5,3) * t212 + Ifges(5,1) * t258 + Ifges(5,4) * t257 + Ifges(5,5) * t309 - pkin(9) * t172 - t331 * t164 + t336 * t165 + t285 * t254 - t323 * t255;
t274 = Ifges(4,5) * t311 + Ifges(4,6) * t310 + Ifges(4,3) * t323;
t276 = Ifges(4,1) * t311 + Ifges(4,4) * t310 + Ifges(4,5) * t323;
t271 = -t323 * mrSges(5,2) + t285 * mrSges(5,3);
t272 = t323 * mrSges(5,1) - t286 * mrSges(5,3);
t346 = m(5) * t244 - t257 * mrSges(5,1) + t258 * mrSges(5,2) - t285 * t271 + t286 * t272 + t349;
t264 = -t285 * mrSges(5,1) + t286 * mrSges(5,2);
t169 = m(5) * t212 + t309 * mrSges(5,1) - t258 * mrSges(5,3) - t286 * t264 + t323 * t271 + t172;
t170 = m(5) * t213 - t309 * mrSges(5,2) + t257 * mrSges(5,3) + t285 * t264 - t323 * t272 + t354;
t355 = -t328 * t169 + t329 * t170;
t145 = -mrSges(4,1) * t269 + mrSges(4,3) * t247 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * t309 - pkin(3) * t346 + qJ(4) * t355 + t329 * t158 + t328 * t159 - t311 * t274 + t323 * t276;
t163 = t329 * t169 + t328 * t170;
t275 = Ifges(4,4) * t311 + Ifges(4,2) * t310 + Ifges(4,6) * t323;
t146 = mrSges(4,2) * t269 - mrSges(4,3) * t246 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * t309 - qJ(4) * t163 - t328 * t158 + t329 * t159 + t310 * t274 - t323 * t275;
t287 = -t310 * mrSges(4,1) + t311 * mrSges(4,2);
t288 = -t323 * mrSges(4,2) + t310 * mrSges(4,3);
t161 = m(4) * t246 + t309 * mrSges(4,1) - t283 * mrSges(4,3) - t311 * t287 + t323 * t288 + t163;
t290 = t323 * mrSges(4,1) - t311 * mrSges(4,3);
t162 = m(4) * t247 - t309 * mrSges(4,2) + t282 * mrSges(4,3) + t310 * t287 - t323 * t290 + t355;
t157 = -t332 * t161 + t337 * t162;
t188 = -m(4) * t269 + t282 * mrSges(4,1) - t283 * mrSges(4,2) + t310 * t288 - t311 * t290 - t346;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t333 + Ifges(3,2) * t338) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t333 + Ifges(3,4) * t338) * qJD(1);
t361 = mrSges(3,1) * t291 - mrSges(3,2) * t292 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t188 + pkin(8) * t157 + t337 * t145 + t332 * t146 + (t333 * t300 - t338 * t301) * qJD(1);
t312 = (-mrSges(3,1) * t338 + mrSges(3,2) * t333) * qJD(1);
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t360;
t155 = m(3) * t292 - qJDD(2) * mrSges(3,2) + t315 * mrSges(3,3) - qJD(2) * t318 + t312 * t359 + t157;
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t187 = m(3) * t291 + qJDD(2) * mrSges(3,1) - t314 * mrSges(3,3) + qJD(2) * t319 - t312 * t360 + t188;
t356 = t338 * t155 - t333 * t187;
t156 = t337 * t161 + t332 * t162;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t333 + Ifges(3,6) * t338) * qJD(1);
t144 = mrSges(3,2) * t302 - mrSges(3,3) * t291 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t156 - qJD(2) * t300 - t332 * t145 + t337 * t146 + t299 * t359;
t348 = -mrSges(7,1) * t190 + mrSges(7,2) * t191 - Ifges(7,5) * t207 - Ifges(7,6) * t206 - Ifges(7,3) * t298 - t242 * t217 + t241 * t218;
t345 = -mrSges(6,1) * t195 + mrSges(6,2) * t196 - Ifges(6,5) * t227 - Ifges(6,6) * t226 - Ifges(6,3) * t305 - pkin(5) * t179 - t263 * t237 + t262 * t238 + t348;
t343 = -mrSges(5,1) * t212 + mrSges(5,2) * t213 - Ifges(5,5) * t258 - Ifges(5,6) * t257 - Ifges(5,3) * t309 - pkin(4) * t172 - t286 * t255 + t285 * t256 + t345;
t342 = mrSges(4,1) * t246 - mrSges(4,2) * t247 + Ifges(4,5) * t283 + Ifges(4,6) * t282 + Ifges(4,3) * t309 + pkin(3) * t163 + t311 * t275 - t310 * t276 - t343;
t148 = -mrSges(3,1) * t302 + mrSges(3,3) * t292 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t156 + qJD(2) * t301 - t299 * t360 - t342;
t347 = -m(3) * t302 + t315 * mrSges(3,1) - t314 * mrSges(3,2) - t318 * t360 + t319 * t359 - t156;
t350 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t347 + pkin(7) * t356 + t333 * t144 + t338 * t148;
t152 = m(2) * t320 + qJDD(1) * mrSges(2,1) - t341 * mrSges(2,2) + t347;
t151 = t333 * t155 + t338 * t187;
t149 = m(2) * t321 - t341 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t356;
t142 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t341 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t361;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - t341 * Ifges(2,6) - pkin(7) * t151 + t338 * t144 - t333 * t148;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t141 - t334 * t142 - pkin(6) * (t334 * t149 + t339 * t152), t141, t144, t146, t159, t165, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t141 + t339 * t142 + pkin(6) * (t339 * t149 - t334 * t152), t142, t148, t145, t158, t164, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t361, t342, -t343, -t345, -t348;];
m_new  = t1;
