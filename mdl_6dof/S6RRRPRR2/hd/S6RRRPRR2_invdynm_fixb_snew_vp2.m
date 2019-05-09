% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR2
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
% Datum: 2019-05-07 10:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:56:36
% EndTime: 2019-05-07 09:57:37
% DurationCPUTime: 45.23s
% Computational Cost: add. (803658->386), mult. (1805995->488), div. (0->0), fcn. (1357987->12), ass. (0->153)
t339 = sin(qJ(2));
t344 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t318 = t339 * qJDD(1) + t344 * t366;
t340 = sin(qJ(1));
t345 = cos(qJ(1));
t325 = -t345 * g(1) - t340 * g(2);
t346 = qJD(1) ^ 2;
t313 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t325;
t370 = t339 * t313;
t371 = pkin(2) * t346;
t273 = qJDD(2) * pkin(2) - t318 * pkin(8) - t370 + (pkin(8) * t366 + t339 * t371 - g(3)) * t344;
t300 = -t339 * g(3) + t344 * t313;
t319 = t344 * qJDD(1) - t339 * t366;
t369 = qJD(1) * t339;
t323 = qJD(2) * pkin(2) - pkin(8) * t369;
t333 = t344 ^ 2;
t274 = t319 * pkin(8) - qJD(2) * t323 - t333 * t371 + t300;
t338 = sin(qJ(3));
t343 = cos(qJ(3));
t245 = t343 * t273 - t338 * t274;
t310 = (-t338 * t339 + t343 * t344) * qJD(1);
t283 = t310 * qJD(3) + t343 * t318 + t338 * t319;
t311 = (t338 * t344 + t339 * t343) * qJD(1);
t330 = qJDD(2) + qJDD(3);
t331 = qJD(2) + qJD(3);
t226 = (t310 * t331 - t283) * qJ(4) + (t310 * t311 + t330) * pkin(3) + t245;
t246 = t338 * t273 + t343 * t274;
t282 = -t311 * qJD(3) - t338 * t318 + t343 * t319;
t302 = t331 * pkin(3) - t311 * qJ(4);
t306 = t310 ^ 2;
t230 = -t306 * pkin(3) + t282 * qJ(4) - t331 * t302 + t246;
t334 = sin(pkin(11));
t335 = cos(pkin(11));
t297 = t334 * t310 + t335 * t311;
t210 = -0.2e1 * qJD(4) * t297 + t335 * t226 - t334 * t230;
t296 = t335 * t310 - t334 * t311;
t211 = 0.2e1 * qJD(4) * t296 + t334 * t226 + t335 * t230;
t256 = t335 * t282 - t334 * t283;
t267 = -t296 * mrSges(5,1) + t297 * mrSges(5,2);
t286 = t331 * mrSges(5,1) - t297 * mrSges(5,3);
t268 = -t296 * pkin(4) - t297 * pkin(9);
t329 = t331 ^ 2;
t208 = -t329 * pkin(4) + t330 * pkin(9) + t296 * t268 + t211;
t324 = t340 * g(1) - t345 * g(2);
t358 = -qJDD(1) * pkin(1) - t324;
t284 = -t319 * pkin(2) + t323 * t369 + (-pkin(8) * t333 - pkin(7)) * t346 + t358;
t235 = -t282 * pkin(3) - t306 * qJ(4) + t311 * t302 + qJDD(4) + t284;
t257 = t334 * t282 + t335 * t283;
t214 = (-t296 * t331 - t257) * pkin(9) + (t297 * t331 - t256) * pkin(4) + t235;
t337 = sin(qJ(5));
t342 = cos(qJ(5));
t203 = -t337 * t208 + t342 * t214;
t280 = -t337 * t297 + t342 * t331;
t233 = t280 * qJD(5) + t342 * t257 + t337 * t330;
t255 = qJDD(5) - t256;
t281 = t342 * t297 + t337 * t331;
t290 = qJD(5) - t296;
t201 = (t280 * t290 - t233) * pkin(10) + (t280 * t281 + t255) * pkin(5) + t203;
t204 = t342 * t208 + t337 * t214;
t232 = -t281 * qJD(5) - t337 * t257 + t342 * t330;
t261 = t290 * pkin(5) - t281 * pkin(10);
t276 = t280 ^ 2;
t202 = -t276 * pkin(5) + t232 * pkin(10) - t290 * t261 + t204;
t336 = sin(qJ(6));
t341 = cos(qJ(6));
t199 = t341 * t201 - t336 * t202;
t249 = t341 * t280 - t336 * t281;
t219 = t249 * qJD(6) + t336 * t232 + t341 * t233;
t250 = t336 * t280 + t341 * t281;
t227 = -t249 * mrSges(7,1) + t250 * mrSges(7,2);
t287 = qJD(6) + t290;
t236 = -t287 * mrSges(7,2) + t249 * mrSges(7,3);
t248 = qJDD(6) + t255;
t194 = m(7) * t199 + t248 * mrSges(7,1) - t219 * mrSges(7,3) - t250 * t227 + t287 * t236;
t200 = t336 * t201 + t341 * t202;
t218 = -t250 * qJD(6) + t341 * t232 - t336 * t233;
t237 = t287 * mrSges(7,1) - t250 * mrSges(7,3);
t195 = m(7) * t200 - t248 * mrSges(7,2) + t218 * mrSges(7,3) + t249 * t227 - t287 * t237;
t186 = t341 * t194 + t336 * t195;
t258 = -t280 * mrSges(6,1) + t281 * mrSges(6,2);
t259 = -t290 * mrSges(6,2) + t280 * mrSges(6,3);
t184 = m(6) * t203 + t255 * mrSges(6,1) - t233 * mrSges(6,3) - t281 * t258 + t290 * t259 + t186;
t260 = t290 * mrSges(6,1) - t281 * mrSges(6,3);
t361 = -t336 * t194 + t341 * t195;
t185 = m(6) * t204 - t255 * mrSges(6,2) + t232 * mrSges(6,3) + t280 * t258 - t290 * t260 + t361;
t362 = -t337 * t184 + t342 * t185;
t177 = m(5) * t211 - t330 * mrSges(5,2) + t256 * mrSges(5,3) + t296 * t267 - t331 * t286 + t362;
t285 = -t331 * mrSges(5,2) + t296 * mrSges(5,3);
t207 = -t330 * pkin(4) - t329 * pkin(9) + t297 * t268 - t210;
t205 = -t232 * pkin(5) - t276 * pkin(10) + t281 * t261 + t207;
t355 = m(7) * t205 - t218 * mrSges(7,1) + t219 * mrSges(7,2) - t249 * t236 + t250 * t237;
t351 = -m(6) * t207 + t232 * mrSges(6,1) - t233 * mrSges(6,2) + t280 * t259 - t281 * t260 - t355;
t190 = m(5) * t210 + t330 * mrSges(5,1) - t257 * mrSges(5,3) - t297 * t267 + t331 * t285 + t351;
t168 = t334 * t177 + t335 * t190;
t298 = -t310 * mrSges(4,1) + t311 * mrSges(4,2);
t301 = -t331 * mrSges(4,2) + t310 * mrSges(4,3);
t165 = m(4) * t245 + t330 * mrSges(4,1) - t283 * mrSges(4,3) - t311 * t298 + t331 * t301 + t168;
t303 = t331 * mrSges(4,1) - t311 * mrSges(4,3);
t363 = t335 * t177 - t334 * t190;
t166 = m(4) * t246 - t330 * mrSges(4,2) + t282 * mrSges(4,3) + t310 * t298 - t331 * t303 + t363;
t160 = t343 * t165 + t338 * t166;
t299 = -t344 * g(3) - t370;
t308 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t344) * qJD(1);
t309 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t344) * qJD(1);
t292 = Ifges(4,4) * t311 + Ifges(4,2) * t310 + Ifges(4,6) * t331;
t293 = Ifges(4,1) * t311 + Ifges(4,4) * t310 + Ifges(4,5) * t331;
t221 = Ifges(7,5) * t250 + Ifges(7,6) * t249 + Ifges(7,3) * t287;
t223 = Ifges(7,1) * t250 + Ifges(7,4) * t249 + Ifges(7,5) * t287;
t187 = -mrSges(7,1) * t205 + mrSges(7,3) * t200 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t248 - t250 * t221 + t287 * t223;
t222 = Ifges(7,4) * t250 + Ifges(7,2) * t249 + Ifges(7,6) * t287;
t188 = mrSges(7,2) * t205 - mrSges(7,3) * t199 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t248 + t249 * t221 - t287 * t222;
t238 = Ifges(6,5) * t281 + Ifges(6,6) * t280 + Ifges(6,3) * t290;
t240 = Ifges(6,1) * t281 + Ifges(6,4) * t280 + Ifges(6,5) * t290;
t170 = -mrSges(6,1) * t207 + mrSges(6,3) * t204 + Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t255 - pkin(5) * t355 + pkin(10) * t361 + t341 * t187 + t336 * t188 - t281 * t238 + t290 * t240;
t239 = Ifges(6,4) * t281 + Ifges(6,2) * t280 + Ifges(6,6) * t290;
t172 = mrSges(6,2) * t207 - mrSges(6,3) * t203 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t255 - pkin(10) * t186 - t336 * t187 + t341 * t188 + t280 * t238 - t290 * t239;
t263 = Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t331;
t264 = Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t331;
t353 = -mrSges(5,1) * t210 + mrSges(5,2) * t211 - Ifges(5,5) * t257 - Ifges(5,6) * t256 - Ifges(5,3) * t330 - pkin(4) * t351 - pkin(9) * t362 - t342 * t170 - t337 * t172 - t297 * t263 + t296 * t264;
t350 = -mrSges(4,1) * t245 + mrSges(4,2) * t246 - Ifges(4,5) * t283 - Ifges(4,6) * t282 - Ifges(4,3) * t330 - pkin(3) * t168 - t311 * t292 + t310 * t293 + t353;
t372 = mrSges(3,1) * t299 - mrSges(3,2) * t300 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t160 + (t339 * t308 - t344 * t309) * qJD(1) - t350;
t179 = t342 * t184 + t337 * t185;
t368 = qJD(1) * t344;
t317 = (-mrSges(3,1) * t344 + mrSges(3,2) * t339) * qJD(1);
t322 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t368;
t158 = m(3) * t299 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,3) + qJD(2) * t322 - t317 * t369 + t160;
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t369;
t364 = -t338 * t165 + t343 * t166;
t159 = m(3) * t300 - qJDD(2) * mrSges(3,2) + t319 * mrSges(3,3) - qJD(2) * t321 + t317 * t368 + t364;
t365 = -t339 * t158 + t344 * t159;
t357 = m(5) * t235 - t256 * mrSges(5,1) + t257 * mrSges(5,2) - t296 * t285 + t297 * t286 + t179;
t262 = Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t331;
t156 = mrSges(5,2) * t235 - mrSges(5,3) * t210 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * t330 - pkin(9) * t179 - t337 * t170 + t342 * t172 + t296 * t262 - t331 * t263;
t354 = -mrSges(7,1) * t199 + mrSges(7,2) * t200 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t248 - t250 * t222 + t249 * t223;
t348 = mrSges(6,1) * t203 - mrSges(6,2) * t204 + Ifges(6,5) * t233 + Ifges(6,6) * t232 + Ifges(6,3) * t255 + pkin(5) * t186 + t281 * t239 - t280 * t240 - t354;
t161 = -mrSges(5,1) * t235 + mrSges(5,3) * t211 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * t330 - pkin(4) * t179 - t297 * t262 + t331 * t264 - t348;
t291 = Ifges(4,5) * t311 + Ifges(4,6) * t310 + Ifges(4,3) * t331;
t151 = -mrSges(4,1) * t284 + mrSges(4,3) * t246 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * t330 - pkin(3) * t357 + qJ(4) * t363 + t334 * t156 + t335 * t161 - t311 * t291 + t331 * t293;
t152 = mrSges(4,2) * t284 - mrSges(4,3) * t245 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * t330 - qJ(4) * t168 + t335 * t156 - t334 * t161 + t310 * t291 - t331 * t292;
t307 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t344) * qJD(1);
t312 = -t346 * pkin(7) + t358;
t352 = m(4) * t284 - t282 * mrSges(4,1) + t283 * mrSges(4,2) - t310 * t301 + t311 * t303 + t357;
t147 = -mrSges(3,1) * t312 + mrSges(3,3) * t300 + Ifges(3,4) * t318 + Ifges(3,2) * t319 + Ifges(3,6) * qJDD(2) - pkin(2) * t352 + pkin(8) * t364 + qJD(2) * t309 + t343 * t151 + t338 * t152 - t307 * t369;
t149 = mrSges(3,2) * t312 - mrSges(3,3) * t299 + Ifges(3,1) * t318 + Ifges(3,4) * t319 + Ifges(3,5) * qJDD(2) - pkin(8) * t160 - qJD(2) * t308 - t338 * t151 + t343 * t152 + t307 * t368;
t349 = -m(3) * t312 + t319 * mrSges(3,1) - t318 * mrSges(3,2) - t321 * t369 + t322 * t368 - t352;
t356 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t349 + pkin(7) * t365 + t344 * t147 + t339 * t149;
t173 = m(2) * t324 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t349;
t155 = t344 * t158 + t339 * t159;
t153 = m(2) * t325 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t365;
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t155 - t372;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t155 - t339 * t147 + t344 * t149;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t145 - t340 * t150 - pkin(6) * (t340 * t153 + t345 * t173), t145, t149, t152, t156, t172, t188; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t145 + t345 * t150 + pkin(6) * (t345 * t153 - t340 * t173), t150, t147, t151, t161, t170, t187; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t372, -t350, -t353, t348, -t354;];
m_new  = t1;
