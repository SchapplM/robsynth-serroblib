% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:47:07
% EndTime: 2019-05-08 04:47:58
% DurationCPUTime: 20.17s
% Computational Cost: add. (383315->383), mult. (762051->467), div. (0->0), fcn. (553095->10), ass. (0->144)
t333 = sin(qJ(2));
t337 = cos(qJ(2));
t360 = qJD(1) * qJD(2);
t314 = qJDD(1) * t333 + t337 * t360;
t334 = sin(qJ(1));
t338 = cos(qJ(1));
t321 = -g(1) * t338 - g(2) * t334;
t339 = qJD(1) ^ 2;
t308 = -pkin(1) * t339 + qJDD(1) * pkin(7) + t321;
t365 = t333 * t308;
t367 = pkin(2) * t339;
t264 = qJDD(2) * pkin(2) - t314 * pkin(8) - t365 + (pkin(8) * t360 + t333 * t367 - g(3)) * t337;
t294 = -g(3) * t333 + t308 * t337;
t315 = qJDD(1) * t337 - t333 * t360;
t362 = qJD(1) * t333;
t319 = qJD(2) * pkin(2) - pkin(8) * t362;
t329 = t337 ^ 2;
t265 = pkin(8) * t315 - qJD(2) * t319 - t329 * t367 + t294;
t332 = sin(qJ(3));
t336 = cos(qJ(3));
t241 = t264 * t332 + t265 * t336;
t306 = (t332 * t337 + t333 * t336) * qJD(1);
t277 = -qJD(3) * t306 - t314 * t332 + t315 * t336;
t361 = qJD(1) * t337;
t305 = -t332 * t362 + t336 * t361;
t287 = -mrSges(4,1) * t305 + mrSges(4,2) * t306;
t327 = qJD(2) + qJD(3);
t296 = mrSges(4,1) * t327 - mrSges(4,3) * t306;
t326 = qJDD(2) + qJDD(3);
t278 = qJD(3) * t305 + t314 * t336 + t315 * t332;
t320 = t334 * g(1) - g(2) * t338;
t351 = -qJDD(1) * pkin(1) - t320;
t279 = -t315 * pkin(2) + t319 * t362 + (-pkin(8) * t329 - pkin(7)) * t339 + t351;
t222 = (-t305 * t327 - t278) * pkin(9) + (t306 * t327 - t277) * pkin(3) + t279;
t288 = -pkin(3) * t305 - pkin(9) * t306;
t325 = t327 ^ 2;
t225 = -pkin(3) * t325 + pkin(9) * t326 + t288 * t305 + t241;
t331 = sin(qJ(4));
t335 = cos(qJ(4));
t204 = t222 * t335 - t331 * t225;
t291 = -t306 * t331 + t327 * t335;
t247 = qJD(4) * t291 + t278 * t335 + t326 * t331;
t276 = qJDD(4) - t277;
t292 = t306 * t335 + t327 * t331;
t301 = qJD(4) - t305;
t201 = (t291 * t301 - t247) * pkin(10) + (t291 * t292 + t276) * pkin(4) + t204;
t205 = t222 * t331 + t225 * t335;
t246 = -qJD(4) * t292 - t278 * t331 + t326 * t335;
t282 = pkin(4) * t301 - pkin(10) * t292;
t290 = t291 ^ 2;
t203 = -pkin(4) * t290 + pkin(10) * t246 - t282 * t301 + t205;
t330 = sin(qJ(5));
t368 = cos(qJ(5));
t197 = t201 * t330 + t203 * t368;
t258 = t291 * t330 + t292 * t368;
t216 = qJD(5) * t258 - t246 * t368 + t247 * t330;
t299 = qJD(5) + t301;
t250 = mrSges(6,1) * t299 - mrSges(6,3) * t258;
t257 = -t291 * t368 + t292 * t330;
t271 = qJDD(5) + t276;
t236 = pkin(5) * t257 - qJ(6) * t258;
t297 = t299 ^ 2;
t192 = -pkin(5) * t297 + qJ(6) * t271 + 0.2e1 * qJD(6) * t299 - t236 * t257 + t197;
t251 = -mrSges(7,1) * t299 + mrSges(7,2) * t258;
t359 = m(7) * t192 + mrSges(7,3) * t271 + t251 * t299;
t237 = mrSges(7,1) * t257 - mrSges(7,3) * t258;
t363 = -mrSges(6,1) * t257 - mrSges(6,2) * t258 - t237;
t366 = -mrSges(6,3) - mrSges(7,2);
t178 = m(6) * t197 - t271 * mrSges(6,2) + t216 * t366 - t299 * t250 + t257 * t363 + t359;
t196 = t201 * t368 - t203 * t330;
t217 = -qJD(5) * t257 + t246 * t330 + t247 * t368;
t249 = -mrSges(6,2) * t299 - mrSges(6,3) * t257;
t194 = -pkin(5) * t271 - qJ(6) * t297 + t236 * t258 + qJDD(6) - t196;
t248 = -mrSges(7,2) * t257 + mrSges(7,3) * t299;
t354 = -m(7) * t194 + mrSges(7,1) * t271 + t248 * t299;
t180 = m(6) * t196 + t271 * mrSges(6,1) + t217 * t366 + t299 * t249 + t258 * t363 + t354;
t175 = t178 * t330 + t180 * t368;
t262 = -mrSges(5,1) * t291 + mrSges(5,2) * t292;
t280 = -mrSges(5,2) * t301 + mrSges(5,3) * t291;
t171 = m(5) * t204 + mrSges(5,1) * t276 - mrSges(5,3) * t247 - t262 * t292 + t280 * t301 + t175;
t281 = mrSges(5,1) * t301 - mrSges(5,3) * t292;
t355 = t178 * t368 - t180 * t330;
t172 = m(5) * t205 - mrSges(5,2) * t276 + mrSges(5,3) * t246 + t262 * t291 - t281 * t301 + t355;
t356 = -t171 * t331 + t172 * t335;
t164 = m(4) * t241 - mrSges(4,2) * t326 + mrSges(4,3) * t277 + t287 * t305 - t296 * t327 + t356;
t240 = t336 * t264 - t265 * t332;
t295 = -mrSges(4,2) * t327 + mrSges(4,3) * t305;
t224 = -t326 * pkin(3) - t325 * pkin(9) + t288 * t306 - t240;
t206 = -t246 * pkin(4) - t290 * pkin(10) + t282 * t292 + t224;
t199 = -0.2e1 * qJD(6) * t258 + (t257 * t299 - t217) * qJ(6) + (t258 * t299 + t216) * pkin(5) + t206;
t189 = m(7) * t199 + mrSges(7,1) * t216 - mrSges(7,3) * t217 + t248 * t257 - t251 * t258;
t345 = m(6) * t206 + mrSges(6,1) * t216 + mrSges(6,2) * t217 + t249 * t257 + t250 * t258 + t189;
t342 = -m(5) * t224 + mrSges(5,1) * t246 - mrSges(5,2) * t247 + t280 * t291 - t281 * t292 - t345;
t182 = m(4) * t240 + mrSges(4,1) * t326 - mrSges(4,3) * t278 - t287 * t306 + t295 * t327 + t342;
t159 = t164 * t332 + t182 * t336;
t293 = -t337 * g(3) - t365;
t303 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t333 + Ifges(3,2) * t337) * qJD(1);
t304 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t333 + Ifges(3,4) * t337) * qJD(1);
t231 = Ifges(7,1) * t258 + Ifges(7,4) * t299 + Ifges(7,5) * t257;
t232 = Ifges(6,1) * t258 - Ifges(6,4) * t257 + Ifges(6,5) * t299;
t353 = -mrSges(7,1) * t199 + mrSges(7,2) * t192;
t229 = Ifges(7,4) * t258 + Ifges(7,2) * t299 + Ifges(7,6) * t257;
t364 = -Ifges(6,5) * t258 + Ifges(6,6) * t257 - Ifges(6,3) * t299 - t229;
t173 = -mrSges(6,1) * t206 + mrSges(6,3) * t197 - pkin(5) * t189 + (t231 + t232) * t299 + (Ifges(6,6) - Ifges(7,6)) * t271 + t364 * t258 + (Ifges(6,4) - Ifges(7,5)) * t217 + (-Ifges(6,2) - Ifges(7,3)) * t216 + t353;
t230 = Ifges(6,4) * t258 - Ifges(6,2) * t257 + Ifges(6,6) * t299;
t227 = Ifges(7,5) * t258 + Ifges(7,6) * t299 + Ifges(7,3) * t257;
t350 = mrSges(7,2) * t194 - mrSges(7,3) * t199 + Ifges(7,1) * t217 + Ifges(7,4) * t271 + Ifges(7,5) * t216 + t227 * t299;
t174 = mrSges(6,2) * t206 - mrSges(6,3) * t196 + Ifges(6,1) * t217 - Ifges(6,4) * t216 + Ifges(6,5) * t271 - qJ(6) * t189 - t299 * t230 + t257 * t364 + t350;
t252 = Ifges(5,5) * t292 + Ifges(5,6) * t291 + Ifges(5,3) * t301;
t254 = Ifges(5,1) * t292 + Ifges(5,4) * t291 + Ifges(5,5) * t301;
t153 = -mrSges(5,1) * t224 + mrSges(5,3) * t205 + Ifges(5,4) * t247 + Ifges(5,2) * t246 + Ifges(5,6) * t276 - pkin(4) * t345 + pkin(10) * t355 + t173 * t368 + t330 * t174 - t292 * t252 + t301 * t254;
t253 = Ifges(5,4) * t292 + Ifges(5,2) * t291 + Ifges(5,6) * t301;
t155 = mrSges(5,2) * t224 - mrSges(5,3) * t204 + Ifges(5,1) * t247 + Ifges(5,4) * t246 + Ifges(5,5) * t276 - pkin(10) * t175 - t173 * t330 + t174 * t368 + t252 * t291 - t253 * t301;
t284 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * t327;
t285 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * t327;
t346 = -mrSges(4,1) * t240 + mrSges(4,2) * t241 - Ifges(4,5) * t278 - Ifges(4,6) * t277 - Ifges(4,3) * t326 - pkin(3) * t342 - pkin(9) * t356 - t153 * t335 - t155 * t331 - t284 * t306 + t305 * t285;
t369 = mrSges(3,1) * t293 - mrSges(3,2) * t294 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t159 + (t303 * t333 - t304 * t337) * qJD(1) - t346;
t166 = t171 * t335 + t172 * t331;
t313 = (-mrSges(3,1) * t337 + mrSges(3,2) * t333) * qJD(1);
t318 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t361;
t157 = m(3) * t293 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t314 + qJD(2) * t318 - t313 * t362 + t159;
t317 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t362;
t357 = t164 * t336 - t182 * t332;
t158 = m(3) * t294 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t317 + t313 * t361 + t357;
t358 = -t157 * t333 + t158 * t337;
t283 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * t327;
t147 = mrSges(4,2) * t279 - mrSges(4,3) * t240 + Ifges(4,1) * t278 + Ifges(4,4) * t277 + Ifges(4,5) * t326 - pkin(9) * t166 - t153 * t331 + t155 * t335 + t283 * t305 - t284 * t327;
t348 = mrSges(7,1) * t194 - mrSges(7,3) * t192 - Ifges(7,4) * t217 - Ifges(7,2) * t271 - Ifges(7,6) * t216 + t227 * t258 - t231 * t257;
t343 = mrSges(6,2) * t197 - t257 * t232 - qJ(6) * (-t216 * mrSges(7,2) - t257 * t237 + t359) - pkin(5) * (-t217 * mrSges(7,2) - t258 * t237 + t354) - mrSges(6,1) * t196 + Ifges(6,6) * t216 - Ifges(6,5) * t217 - t258 * t230 - Ifges(6,3) * t271 + t348;
t340 = mrSges(5,1) * t204 - mrSges(5,2) * t205 + Ifges(5,5) * t247 + Ifges(5,6) * t246 + Ifges(5,3) * t276 + pkin(4) * t175 + t253 * t292 - t254 * t291 - t343;
t151 = -mrSges(4,1) * t279 + mrSges(4,3) * t241 + Ifges(4,4) * t278 + Ifges(4,2) * t277 + Ifges(4,6) * t326 - pkin(3) * t166 - t283 * t306 + t285 * t327 - t340;
t302 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t333 + Ifges(3,6) * t337) * qJD(1);
t307 = -t339 * pkin(7) + t351;
t347 = m(4) * t279 - mrSges(4,1) * t277 + mrSges(4,2) * t278 - t295 * t305 + t296 * t306 + t166;
t143 = -mrSges(3,1) * t307 + mrSges(3,3) * t294 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t347 + pkin(8) * t357 + qJD(2) * t304 + t332 * t147 + t336 * t151 - t302 * t362;
t146 = mrSges(3,2) * t307 - mrSges(3,3) * t293 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t159 - qJD(2) * t303 + t147 * t336 - t151 * t332 + t302 * t361;
t344 = -m(3) * t307 + mrSges(3,1) * t315 - mrSges(3,2) * t314 - t317 * t362 + t318 * t361 - t347;
t349 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t344 + pkin(7) * t358 + t143 * t337 + t146 * t333;
t160 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t339 + t344;
t150 = t157 * t337 + t158 * t333;
t148 = m(2) * t321 - mrSges(2,1) * t339 - qJDD(1) * mrSges(2,2) + t358;
t144 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t339 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t150 - t369;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t339 - pkin(7) * t150 - t143 * t333 + t146 * t337;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t338 * t141 - t334 * t144 - pkin(6) * (t148 * t334 + t160 * t338), t141, t146, t147, t155, t174, -t229 * t257 + t350; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t141 + t338 * t144 + pkin(6) * (t148 * t338 - t160 * t334), t144, t143, t151, t153, t173, -t348; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t349, t349, t369, -t346, t340, -t343, Ifges(7,5) * t217 + Ifges(7,6) * t271 + Ifges(7,3) * t216 + t258 * t229 - t299 * t231 - t353;];
m_new  = t1;
