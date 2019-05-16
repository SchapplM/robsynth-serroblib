% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:40:24
% EndTime: 2019-05-07 18:41:43
% DurationCPUTime: 39.54s
% Computational Cost: add. (730545->395), mult. (1568332->500), div. (0->0), fcn. (1248567->12), ass. (0->156)
t330 = sin(pkin(6));
t334 = sin(qJ(2));
t338 = cos(qJ(2));
t358 = qJD(1) * qJD(2);
t317 = (-qJDD(1) * t338 + t334 * t358) * t330;
t372 = -2 * qJD(5);
t371 = pkin(8) * t330;
t331 = cos(pkin(6));
t370 = t331 * g(3);
t369 = -mrSges(6,3) - mrSges(7,2);
t368 = cos(pkin(11));
t367 = t330 * t334;
t366 = t330 * t338;
t365 = t331 * t334;
t364 = t331 * t338;
t360 = qJD(1) * t330;
t315 = (-pkin(2) * t338 - pkin(9) * t334) * t360;
t325 = t331 * qJD(1) + qJD(2);
t323 = t325 ^ 2;
t324 = t331 * qJDD(1) + qJDD(2);
t359 = qJD(1) * t338;
t335 = sin(qJ(1));
t339 = cos(qJ(1));
t321 = t335 * g(1) - t339 * g(2);
t340 = qJD(1) ^ 2;
t312 = qJDD(1) * pkin(1) + t340 * t371 + t321;
t322 = -t339 * g(1) - t335 * g(2);
t313 = -t340 * pkin(1) + qJDD(1) * t371 + t322;
t361 = t312 * t365 + t338 * t313;
t264 = -t323 * pkin(2) + t324 * pkin(9) + (-g(3) * t334 + t315 * t359) * t330 + t361;
t316 = (qJDD(1) * t334 + t338 * t358) * t330;
t265 = t317 * pkin(2) - t316 * pkin(9) - t370 + (-t312 + (pkin(2) * t334 - pkin(9) * t338) * t325 * qJD(1)) * t330;
t333 = sin(qJ(3));
t337 = cos(qJ(3));
t232 = t337 * t264 + t333 * t265;
t356 = t334 * t360;
t305 = t337 * t325 - t333 * t356;
t306 = t333 * t325 + t337 * t356;
t290 = -t305 * pkin(3) - t306 * pkin(10);
t309 = qJDD(3) + t317;
t355 = t330 * t359;
t320 = qJD(3) - t355;
t319 = t320 ^ 2;
t213 = -t319 * pkin(3) + t309 * pkin(10) + t305 * t290 + t232;
t287 = -g(3) * t366 + t312 * t364 - t334 * t313;
t263 = -t324 * pkin(2) - t323 * pkin(9) + t315 * t356 - t287;
t285 = -t306 * qJD(3) - t333 * t316 + t337 * t324;
t286 = t305 * qJD(3) + t337 * t316 + t333 * t324;
t218 = (-t305 * t320 - t286) * pkin(10) + (t306 * t320 - t285) * pkin(3) + t263;
t332 = sin(qJ(4));
t336 = cos(qJ(4));
t208 = -t332 * t213 + t336 * t218;
t293 = -t332 * t306 + t336 * t320;
t250 = t293 * qJD(4) + t336 * t286 + t332 * t309;
t283 = qJDD(4) - t285;
t294 = t336 * t306 + t332 * t320;
t304 = qJD(4) - t305;
t205 = (t293 * t304 - t250) * qJ(5) + (t293 * t294 + t283) * pkin(4) + t208;
t209 = t336 * t213 + t332 * t218;
t249 = -t294 * qJD(4) - t332 * t286 + t336 * t309;
t272 = t304 * pkin(4) - t294 * qJ(5);
t292 = t293 ^ 2;
t207 = -t292 * pkin(4) + t249 * qJ(5) - t304 * t272 + t209;
t329 = sin(pkin(11));
t267 = -t293 * t368 + t329 * t294;
t201 = t329 * t205 + t368 * t207 + t267 * t372;
t228 = -t249 * t368 + t329 * t250;
t268 = t329 * t293 + t294 * t368;
t253 = t304 * mrSges(6,1) - t268 * mrSges(6,3);
t242 = pkin(5) * t267 - qJ(6) * t268;
t303 = t304 ^ 2;
t196 = -t303 * pkin(5) + t283 * qJ(6) + 0.2e1 * qJD(6) * t304 - t267 * t242 + t201;
t254 = -t304 * mrSges(7,1) + t268 * mrSges(7,2);
t357 = m(7) * t196 + t283 * mrSges(7,3) + t304 * t254;
t243 = mrSges(7,1) * t267 - mrSges(7,3) * t268;
t362 = -mrSges(6,1) * t267 - mrSges(6,2) * t268 - t243;
t183 = m(6) * t201 - t283 * mrSges(6,2) + t228 * t369 - t304 * t253 + t267 * t362 + t357;
t350 = t205 * t368 - t329 * t207;
t200 = t268 * t372 + t350;
t229 = t329 * t249 + t250 * t368;
t252 = -t304 * mrSges(6,2) - t267 * mrSges(6,3);
t198 = -t283 * pkin(5) - t303 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t242) * t268 - t350;
t251 = -t267 * mrSges(7,2) + t304 * mrSges(7,3);
t352 = -m(7) * t198 + t283 * mrSges(7,1) + t304 * t251;
t185 = m(6) * t200 + t283 * mrSges(6,1) + t229 * t369 + t304 * t252 + t268 * t362 + t352;
t178 = t329 * t183 + t368 * t185;
t269 = -t293 * mrSges(5,1) + t294 * mrSges(5,2);
t271 = -t304 * mrSges(5,2) + t293 * mrSges(5,3);
t176 = m(5) * t208 + t283 * mrSges(5,1) - t250 * mrSges(5,3) - t294 * t269 + t304 * t271 + t178;
t273 = t304 * mrSges(5,1) - t294 * mrSges(5,3);
t353 = t368 * t183 - t329 * t185;
t177 = m(5) * t209 - t283 * mrSges(5,2) + t249 * mrSges(5,3) + t293 * t269 - t304 * t273 + t353;
t174 = -t332 * t176 + t336 * t177;
t289 = -t305 * mrSges(4,1) + t306 * mrSges(4,2);
t296 = t320 * mrSges(4,1) - t306 * mrSges(4,3);
t172 = m(4) * t232 - t309 * mrSges(4,2) + t285 * mrSges(4,3) + t305 * t289 - t320 * t296 + t174;
t231 = -t333 * t264 + t337 * t265;
t212 = -t309 * pkin(3) - t319 * pkin(10) + t306 * t290 - t231;
t210 = -t249 * pkin(4) - t292 * qJ(5) + t294 * t272 + qJDD(5) + t212;
t203 = -0.2e1 * qJD(6) * t268 + (t267 * t304 - t229) * qJ(6) + (t268 * t304 + t228) * pkin(5) + t210;
t193 = m(7) * t203 + t228 * mrSges(7,1) - t229 * mrSges(7,3) + t267 * t251 - t268 * t254;
t345 = m(6) * t210 + t228 * mrSges(6,1) + t229 * mrSges(6,2) + t267 * t252 + t268 * t253 + t193;
t188 = -m(5) * t212 + t249 * mrSges(5,1) - t250 * mrSges(5,2) + t293 * t271 - t294 * t273 - t345;
t295 = -t320 * mrSges(4,2) + t305 * mrSges(4,3);
t187 = m(4) * t231 + t309 * mrSges(4,1) - t286 * mrSges(4,3) - t306 * t289 + t320 * t295 + t188;
t167 = t333 * t172 + t337 * t187;
t235 = Ifges(7,4) * t268 + Ifges(7,2) * t304 + Ifges(7,6) * t267;
t363 = -Ifges(6,5) * t268 + Ifges(6,6) * t267 - Ifges(6,3) * t304 - t235;
t288 = -g(3) * t367 + t361;
t310 = t325 * mrSges(3,1) - mrSges(3,3) * t356;
t314 = (-mrSges(3,1) * t338 + mrSges(3,2) * t334) * t360;
t354 = t337 * t172 - t333 * t187;
t165 = m(3) * t288 - t324 * mrSges(3,2) - t317 * mrSges(3,3) - t325 * t310 + t314 * t355 + t354;
t311 = -t325 * mrSges(3,2) + mrSges(3,3) * t355;
t173 = t336 * t176 + t332 * t177;
t344 = -m(4) * t263 + t285 * mrSges(4,1) - t286 * mrSges(4,2) + t305 * t295 - t306 * t296 - t173;
t169 = m(3) * t287 + t324 * mrSges(3,1) - t316 * mrSges(3,3) + t325 * t311 - t314 * t356 + t344;
t159 = t338 * t165 - t334 * t169;
t300 = -t330 * t312 - t370;
t166 = m(3) * t300 + t317 * mrSges(3,1) + t316 * mrSges(3,2) + (t310 * t334 - t311 * t338) * t360 + t167;
t156 = t165 * t365 - t330 * t166 + t169 * t364;
t351 = -mrSges(7,1) * t203 + mrSges(7,2) * t196;
t233 = Ifges(7,5) * t268 + Ifges(7,6) * t304 + Ifges(7,3) * t267;
t348 = mrSges(7,2) * t198 - mrSges(7,3) * t203 + Ifges(7,1) * t229 + Ifges(7,4) * t283 + Ifges(7,5) * t228 + t304 * t233;
t237 = Ifges(7,1) * t268 + Ifges(7,4) * t304 + Ifges(7,5) * t267;
t238 = Ifges(6,1) * t268 - Ifges(6,4) * t267 + Ifges(6,5) * t304;
t179 = -mrSges(6,1) * t210 + mrSges(6,3) * t201 - pkin(5) * t193 + (t237 + t238) * t304 + (Ifges(6,6) - Ifges(7,6)) * t283 + t363 * t268 + (Ifges(6,4) - Ifges(7,5)) * t229 + (-Ifges(6,2) - Ifges(7,3)) * t228 + t351;
t236 = Ifges(6,4) * t268 - Ifges(6,2) * t267 + Ifges(6,6) * t304;
t180 = mrSges(6,2) * t210 - mrSges(6,3) * t200 + Ifges(6,1) * t229 - Ifges(6,4) * t228 + Ifges(6,5) * t283 - qJ(6) * t193 - t304 * t236 + t267 * t363 + t348;
t255 = Ifges(5,5) * t294 + Ifges(5,6) * t293 + Ifges(5,3) * t304;
t257 = Ifges(5,1) * t294 + Ifges(5,4) * t293 + Ifges(5,5) * t304;
t161 = -mrSges(5,1) * t212 + mrSges(5,3) * t209 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * t283 - pkin(4) * t345 + qJ(5) * t353 + t179 * t368 + t329 * t180 - t294 * t255 + t304 * t257;
t256 = Ifges(5,4) * t294 + Ifges(5,2) * t293 + Ifges(5,6) * t304;
t162 = mrSges(5,2) * t212 - mrSges(5,3) * t208 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * t283 - qJ(5) * t178 - t329 * t179 + t180 * t368 + t293 * t255 - t304 * t256;
t279 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * t320;
t280 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * t320;
t152 = mrSges(4,2) * t263 - mrSges(4,3) * t231 + Ifges(4,1) * t286 + Ifges(4,4) * t285 + Ifges(4,5) * t309 - pkin(10) * t173 - t332 * t161 + t336 * t162 + t305 * t279 - t320 * t280;
t281 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * t320;
t346 = mrSges(7,1) * t198 - mrSges(7,3) * t196 - Ifges(7,4) * t229 - Ifges(7,2) * t283 - Ifges(7,6) * t228 + t268 * t233 - t267 * t237;
t343 = mrSges(6,2) * t201 - t267 * t238 - qJ(6) * (-t228 * mrSges(7,2) - t267 * t243 + t357) - pkin(5) * (-t229 * mrSges(7,2) - t268 * t243 + t352) - mrSges(6,1) * t200 - t268 * t236 + Ifges(6,6) * t228 - Ifges(6,5) * t229 - Ifges(6,3) * t283 + t346;
t341 = mrSges(5,1) * t208 - mrSges(5,2) * t209 + Ifges(5,5) * t250 + Ifges(5,6) * t249 + Ifges(5,3) * t283 + pkin(4) * t178 + t294 * t256 - t293 * t257 - t343;
t160 = -mrSges(4,1) * t263 + mrSges(4,3) * t232 + Ifges(4,4) * t286 + Ifges(4,2) * t285 + Ifges(4,6) * t309 - pkin(3) * t173 - t306 * t279 + t320 * t281 - t341;
t298 = Ifges(3,6) * t325 + (Ifges(3,4) * t334 + Ifges(3,2) * t338) * t360;
t299 = Ifges(3,5) * t325 + (Ifges(3,1) * t334 + Ifges(3,4) * t338) * t360;
t147 = Ifges(3,5) * t316 - Ifges(3,6) * t317 + Ifges(3,3) * t324 + mrSges(3,1) * t287 - mrSges(3,2) * t288 + t333 * t152 + t337 * t160 + pkin(2) * t344 + pkin(9) * t354 + (t298 * t334 - t299 * t338) * t360;
t297 = Ifges(3,3) * t325 + (Ifges(3,5) * t334 + Ifges(3,6) * t338) * t360;
t149 = mrSges(3,2) * t300 - mrSges(3,3) * t287 + Ifges(3,1) * t316 - Ifges(3,4) * t317 + Ifges(3,5) * t324 - pkin(9) * t167 + t337 * t152 - t333 * t160 + t297 * t355 - t325 * t298;
t342 = mrSges(4,1) * t231 - mrSges(4,2) * t232 + Ifges(4,5) * t286 + Ifges(4,6) * t285 + Ifges(4,3) * t309 + pkin(3) * t188 + pkin(10) * t174 + t336 * t161 + t332 * t162 + t306 * t280 - t305 * t281;
t151 = -mrSges(3,1) * t300 + mrSges(3,3) * t288 + Ifges(3,4) * t316 - Ifges(3,2) * t317 + Ifges(3,6) * t324 - pkin(2) * t167 - t297 * t356 + t325 * t299 - t342;
t347 = mrSges(2,1) * t321 - mrSges(2,2) * t322 + Ifges(2,3) * qJDD(1) + pkin(1) * t156 + t331 * t147 + t149 * t367 + t151 * t366 + t159 * t371;
t157 = m(2) * t322 - t340 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t159;
t155 = t331 * t166 + (t165 * t334 + t169 * t338) * t330;
t153 = m(2) * t321 + qJDD(1) * mrSges(2,1) - t340 * mrSges(2,2) + t156;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t321 + Ifges(2,5) * qJDD(1) - t340 * Ifges(2,6) + t338 * t149 - t334 * t151 + (-t155 * t330 - t156 * t331) * pkin(8);
t144 = mrSges(2,1) * g(3) + mrSges(2,3) * t322 + t340 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t155 - t330 * t147 + (pkin(8) * t159 + t149 * t334 + t151 * t338) * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t145 - t335 * t144 - pkin(7) * (t339 * t153 + t335 * t157), t145, t149, t152, t162, t180, -t267 * t235 + t348; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t145 + t339 * t144 + pkin(7) * (-t335 * t153 + t339 * t157), t144, t151, t160, t161, t179, -t346; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t347, t347, t147, t342, t341, -t343, Ifges(7,5) * t229 + Ifges(7,6) * t283 + Ifges(7,3) * t228 + t268 * t235 - t304 * t237 - t351;];
m_new  = t1;
