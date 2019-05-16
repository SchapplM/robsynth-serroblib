% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:45:25
% EndTime: 2019-05-07 20:47:37
% DurationCPUTime: 61.92s
% Computational Cost: add. (1127355->384), mult. (2358622->482), div. (0->0), fcn. (1754680->12), ass. (0->152)
t334 = sin(qJ(1));
t339 = cos(qJ(1));
t320 = t334 * g(1) - t339 * g(2);
t341 = qJD(1) ^ 2;
t303 = -qJDD(1) * pkin(1) - t341 * pkin(7) - t320;
t333 = sin(qJ(2));
t338 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t357 = t338 * t358;
t314 = qJDD(1) * t333 + t357;
t324 = t333 * t358;
t315 = qJDD(1) * t338 - t324;
t269 = (-t314 - t357) * pkin(8) + (-t315 + t324) * pkin(2) + t303;
t321 = -g(1) * t339 - g(2) * t334;
t304 = -pkin(1) * t341 + qJDD(1) * pkin(7) + t321;
t291 = -g(3) * t333 + t338 * t304;
t313 = (-pkin(2) * t338 - pkin(8) * t333) * qJD(1);
t340 = qJD(2) ^ 2;
t359 = qJD(1) * t338;
t272 = -pkin(2) * t340 + qJDD(2) * pkin(8) + t313 * t359 + t291;
t332 = sin(qJ(3));
t337 = cos(qJ(3));
t251 = t337 * t269 - t332 * t272;
t360 = qJD(1) * t333;
t310 = qJD(2) * t337 - t332 * t360;
t283 = qJD(3) * t310 + qJDD(2) * t332 + t314 * t337;
t309 = qJDD(3) - t315;
t311 = qJD(2) * t332 + t337 * t360;
t323 = qJD(3) - t359;
t231 = (t310 * t323 - t283) * pkin(9) + (t310 * t311 + t309) * pkin(3) + t251;
t252 = t332 * t269 + t337 * t272;
t282 = -qJD(3) * t311 + qJDD(2) * t337 - t314 * t332;
t292 = pkin(3) * t323 - pkin(9) * t311;
t308 = t310 ^ 2;
t233 = -pkin(3) * t308 + pkin(9) * t282 - t292 * t323 + t252;
t331 = sin(qJ(4));
t336 = cos(qJ(4));
t212 = t336 * t231 - t331 * t233;
t285 = t310 * t336 - t311 * t331;
t250 = qJD(4) * t285 + t282 * t331 + t283 * t336;
t286 = t310 * t331 + t311 * t336;
t305 = qJDD(4) + t309;
t322 = qJD(4) + t323;
t201 = (t285 * t322 - t250) * qJ(5) + (t285 * t286 + t305) * pkin(4) + t212;
t213 = t331 * t231 + t336 * t233;
t249 = -qJD(4) * t286 + t282 * t336 - t283 * t331;
t274 = pkin(4) * t322 - qJ(5) * t286;
t284 = t285 ^ 2;
t209 = -pkin(4) * t284 + qJ(5) * t249 - t274 * t322 + t213;
t328 = sin(pkin(11));
t329 = cos(pkin(11));
t265 = t285 * t328 + t286 * t329;
t195 = -0.2e1 * qJD(5) * t265 + t329 * t201 - t328 * t209;
t228 = t249 * t328 + t250 * t329;
t264 = t285 * t329 - t286 * t328;
t192 = (t264 * t322 - t228) * pkin(10) + (t264 * t265 + t305) * pkin(5) + t195;
t196 = 0.2e1 * qJD(5) * t264 + t328 * t201 + t329 * t209;
t227 = t249 * t329 - t250 * t328;
t255 = pkin(5) * t322 - pkin(10) * t265;
t263 = t264 ^ 2;
t193 = -pkin(5) * t263 + pkin(10) * t227 - t255 * t322 + t196;
t330 = sin(qJ(6));
t335 = cos(qJ(6));
t191 = t192 * t330 + t193 * t335;
t290 = -t338 * g(3) - t333 * t304;
t271 = -qJDD(2) * pkin(2) - pkin(8) * t340 + t313 * t360 - t290;
t244 = -pkin(3) * t282 - pkin(9) * t308 + t311 * t292 + t271;
t215 = -pkin(4) * t249 - qJ(5) * t284 + t286 * t274 + qJDD(5) + t244;
t198 = -pkin(5) * t227 - pkin(10) * t263 + t255 * t265 + t215;
t242 = t264 * t330 + t265 * t335;
t206 = -qJD(6) * t242 + t227 * t335 - t228 * t330;
t241 = t264 * t335 - t265 * t330;
t207 = qJD(6) * t241 + t227 * t330 + t228 * t335;
t317 = qJD(6) + t322;
t216 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t317;
t218 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t317;
t299 = qJDD(6) + t305;
t180 = -mrSges(7,1) * t198 + mrSges(7,3) * t191 + Ifges(7,4) * t207 + Ifges(7,2) * t206 + Ifges(7,6) * t299 - t216 * t242 + t218 * t317;
t190 = t192 * t335 - t193 * t330;
t217 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t317;
t181 = mrSges(7,2) * t198 - mrSges(7,3) * t190 + Ifges(7,1) * t207 + Ifges(7,4) * t206 + Ifges(7,5) * t299 + t216 * t241 - t217 * t317;
t236 = Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t322;
t238 = Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t322;
t234 = -mrSges(7,2) * t317 + mrSges(7,3) * t241;
t235 = mrSges(7,1) * t317 - mrSges(7,3) * t242;
t352 = m(7) * t198 - t206 * mrSges(7,1) + t207 * mrSges(7,2) - t241 * t234 + t242 * t235;
t221 = -mrSges(7,1) * t241 + mrSges(7,2) * t242;
t185 = m(7) * t190 + mrSges(7,1) * t299 - mrSges(7,3) * t207 - t221 * t242 + t234 * t317;
t186 = m(7) * t191 - mrSges(7,2) * t299 + mrSges(7,3) * t206 + t221 * t241 - t235 * t317;
t353 = -t185 * t330 + t335 * t186;
t164 = -mrSges(6,1) * t215 + mrSges(6,3) * t196 + Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t305 - pkin(5) * t352 + pkin(10) * t353 + t335 * t180 + t330 * t181 - t265 * t236 + t322 * t238;
t179 = t335 * t185 + t330 * t186;
t237 = Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t322;
t165 = mrSges(6,2) * t215 - mrSges(6,3) * t195 + Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t305 - pkin(10) * t179 - t180 * t330 + t181 * t335 + t236 * t264 - t237 * t322;
t256 = Ifges(5,5) * t286 + Ifges(5,6) * t285 + Ifges(5,3) * t322;
t258 = Ifges(5,1) * t286 + Ifges(5,4) * t285 + Ifges(5,5) * t322;
t253 = -mrSges(6,2) * t322 + mrSges(6,3) * t264;
t254 = mrSges(6,1) * t322 - mrSges(6,3) * t265;
t349 = m(6) * t215 - t227 * mrSges(6,1) + t228 * mrSges(6,2) - t264 * t253 + t265 * t254 + t352;
t243 = -mrSges(6,1) * t264 + mrSges(6,2) * t265;
t176 = m(6) * t195 + mrSges(6,1) * t305 - mrSges(6,3) * t228 - t243 * t265 + t253 * t322 + t179;
t177 = m(6) * t196 - mrSges(6,2) * t305 + mrSges(6,3) * t227 + t243 * t264 - t254 * t322 + t353;
t354 = -t176 * t328 + t329 * t177;
t158 = -mrSges(5,1) * t244 + mrSges(5,3) * t213 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * t305 - pkin(4) * t349 + qJ(5) * t354 + t329 * t164 + t328 * t165 - t286 * t256 + t322 * t258;
t172 = t329 * t176 + t328 * t177;
t257 = Ifges(5,4) * t286 + Ifges(5,2) * t285 + Ifges(5,6) * t322;
t159 = mrSges(5,2) * t244 - mrSges(5,3) * t212 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * t305 - qJ(5) * t172 - t164 * t328 + t165 * t329 + t256 * t285 - t257 * t322;
t276 = Ifges(4,5) * t311 + Ifges(4,6) * t310 + Ifges(4,3) * t323;
t278 = Ifges(4,1) * t311 + Ifges(4,4) * t310 + Ifges(4,5) * t323;
t273 = -mrSges(5,2) * t322 + mrSges(5,3) * t285;
t275 = mrSges(5,1) * t322 - mrSges(5,3) * t286;
t346 = m(5) * t244 - t249 * mrSges(5,1) + t250 * mrSges(5,2) - t285 * t273 + t286 * t275 + t349;
t266 = -mrSges(5,1) * t285 + mrSges(5,2) * t286;
t169 = m(5) * t212 + mrSges(5,1) * t305 - mrSges(5,3) * t250 - t266 * t286 + t273 * t322 + t172;
t170 = m(5) * t213 - mrSges(5,2) * t305 + mrSges(5,3) * t249 + t266 * t285 - t275 * t322 + t354;
t355 = -t169 * t331 + t336 * t170;
t145 = -mrSges(4,1) * t271 + mrSges(4,3) * t252 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * t309 - pkin(3) * t346 + pkin(9) * t355 + t336 * t158 + t331 * t159 - t311 * t276 + t323 * t278;
t163 = t336 * t169 + t331 * t170;
t277 = Ifges(4,4) * t311 + Ifges(4,2) * t310 + Ifges(4,6) * t323;
t146 = mrSges(4,2) * t271 - mrSges(4,3) * t251 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * t309 - pkin(9) * t163 - t158 * t331 + t159 * t336 + t276 * t310 - t277 * t323;
t287 = -mrSges(4,1) * t310 + mrSges(4,2) * t311;
t288 = -mrSges(4,2) * t323 + mrSges(4,3) * t310;
t161 = m(4) * t251 + mrSges(4,1) * t309 - mrSges(4,3) * t283 - t287 * t311 + t288 * t323 + t163;
t289 = mrSges(4,1) * t323 - mrSges(4,3) * t311;
t162 = m(4) * t252 - mrSges(4,2) * t309 + mrSges(4,3) * t282 + t287 * t310 - t289 * t323 + t355;
t157 = -t161 * t332 + t337 * t162;
t188 = -m(4) * t271 + t282 * mrSges(4,1) - t283 * mrSges(4,2) + t310 * t288 - t311 * t289 - t346;
t301 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t333 + Ifges(3,2) * t338) * qJD(1);
t302 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t333 + Ifges(3,4) * t338) * qJD(1);
t361 = mrSges(3,1) * t290 - mrSges(3,2) * t291 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t188 + pkin(8) * t157 + t337 * t145 + t332 * t146 + (t301 * t333 - t302 * t338) * qJD(1);
t312 = (-mrSges(3,1) * t338 + mrSges(3,2) * t333) * qJD(1);
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t360;
t155 = m(3) * t291 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t318 + t312 * t359 + t157;
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t187 = m(3) * t290 + qJDD(2) * mrSges(3,1) - t314 * mrSges(3,3) + qJD(2) * t319 - t312 * t360 + t188;
t356 = t338 * t155 - t187 * t333;
t156 = t161 * t337 + t162 * t332;
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t333 + Ifges(3,6) * t338) * qJD(1);
t144 = mrSges(3,2) * t303 - mrSges(3,3) * t290 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t156 - qJD(2) * t301 - t145 * t332 + t146 * t337 + t300 * t359;
t348 = -mrSges(7,1) * t190 + mrSges(7,2) * t191 - Ifges(7,5) * t207 - Ifges(7,6) * t206 - Ifges(7,3) * t299 - t242 * t217 + t241 * t218;
t345 = -mrSges(6,1) * t195 + mrSges(6,2) * t196 - Ifges(6,5) * t228 - Ifges(6,6) * t227 - Ifges(6,3) * t305 - pkin(5) * t179 - t265 * t237 + t264 * t238 + t348;
t343 = -mrSges(5,1) * t212 + mrSges(5,2) * t213 - Ifges(5,5) * t250 - Ifges(5,6) * t249 - Ifges(5,3) * t305 - pkin(4) * t172 - t286 * t257 + t285 * t258 + t345;
t342 = mrSges(4,1) * t251 - mrSges(4,2) * t252 + Ifges(4,5) * t283 + Ifges(4,6) * t282 + Ifges(4,3) * t309 + pkin(3) * t163 + t311 * t277 - t310 * t278 - t343;
t148 = -mrSges(3,1) * t303 + mrSges(3,3) * t291 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t156 + qJD(2) * t302 - t300 * t360 - t342;
t347 = -m(3) * t303 + t315 * mrSges(3,1) - mrSges(3,2) * t314 - t318 * t360 + t319 * t359 - t156;
t350 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t347 + pkin(7) * t356 + t333 * t144 + t338 * t148;
t152 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t341 + t347;
t151 = t155 * t333 + t187 * t338;
t149 = m(2) * t321 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t356;
t142 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t341 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t361;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t341 - pkin(7) * t151 + t144 * t338 - t148 * t333;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t141 - t334 * t142 - pkin(6) * (t149 * t334 + t152 * t339), t141, t144, t146, t159, t165, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t141 + t339 * t142 + pkin(6) * (t149 * t339 - t152 * t334), t142, t148, t145, t158, t164, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t361, t342, -t343, -t345, -t348;];
m_new  = t1;
