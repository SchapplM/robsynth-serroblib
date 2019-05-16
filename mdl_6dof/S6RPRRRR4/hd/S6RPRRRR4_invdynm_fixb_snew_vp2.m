% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:12:44
% EndTime: 2019-05-06 03:13:42
% DurationCPUTime: 48.72s
% Computational Cost: add. (802195->364), mult. (2006195->455), div. (0->0), fcn. (1626457->12), ass. (0->154)
t338 = qJD(1) ^ 2;
t332 = sin(qJ(1));
t337 = cos(qJ(1));
t308 = -t337 * g(1) - t332 * g(2);
t301 = -t338 * pkin(1) + qJDD(1) * qJ(2) + t308;
t326 = sin(pkin(11));
t327 = cos(pkin(11));
t365 = qJD(1) * qJD(2);
t363 = -t327 * g(3) - 0.2e1 * t326 * t365;
t370 = pkin(2) * t327;
t274 = (-pkin(7) * qJDD(1) + t338 * t370 - t301) * t326 + t363;
t292 = -t326 * g(3) + (t301 + 0.2e1 * t365) * t327;
t364 = qJDD(1) * t327;
t322 = t327 ^ 2;
t368 = t322 * t338;
t275 = -pkin(2) * t368 + pkin(7) * t364 + t292;
t331 = sin(qJ(3));
t336 = cos(qJ(3));
t253 = t336 * t274 - t331 * t275;
t352 = t326 * t336 + t327 * t331;
t351 = -t326 * t331 + t327 * t336;
t299 = t351 * qJD(1);
t366 = t299 * qJD(3);
t290 = t352 * qJDD(1) + t366;
t300 = t352 * qJD(1);
t231 = (-t290 + t366) * pkin(8) + (t299 * t300 + qJDD(3)) * pkin(3) + t253;
t254 = t331 * t274 + t336 * t275;
t289 = -t300 * qJD(3) + t351 * qJDD(1);
t295 = qJD(3) * pkin(3) - t300 * pkin(8);
t298 = t299 ^ 2;
t240 = -t298 * pkin(3) + t289 * pkin(8) - qJD(3) * t295 + t254;
t330 = sin(qJ(4));
t335 = cos(qJ(4));
t212 = t335 * t231 - t330 * t240;
t280 = t335 * t299 - t330 * t300;
t249 = t280 * qJD(4) + t330 * t289 + t335 * t290;
t281 = t330 * t299 + t335 * t300;
t320 = qJDD(3) + qJDD(4);
t323 = qJD(3) + qJD(4);
t203 = (t280 * t323 - t249) * pkin(9) + (t280 * t281 + t320) * pkin(4) + t212;
t213 = t330 * t231 + t335 * t240;
t248 = -t281 * qJD(4) + t335 * t289 - t330 * t290;
t272 = t323 * pkin(4) - t281 * pkin(9);
t276 = t280 ^ 2;
t205 = -t276 * pkin(4) + t248 * pkin(9) - t323 * t272 + t213;
t329 = sin(qJ(5));
t334 = cos(qJ(5));
t201 = t329 * t203 + t334 * t205;
t265 = t329 * t280 + t334 * t281;
t219 = -t265 * qJD(5) + t334 * t248 - t329 * t249;
t264 = t334 * t280 - t329 * t281;
t237 = -t264 * mrSges(6,1) + t265 * mrSges(6,2);
t318 = qJD(5) + t323;
t256 = t318 * mrSges(6,1) - t265 * mrSges(6,3);
t317 = qJDD(5) + t320;
t239 = -t264 * pkin(5) - t265 * pkin(10);
t316 = t318 ^ 2;
t197 = -t316 * pkin(5) + t317 * pkin(10) + t264 * t239 + t201;
t321 = t326 ^ 2;
t307 = t332 * g(1) - t337 * g(2);
t357 = qJDD(2) - t307;
t288 = (-pkin(1) - t370) * qJDD(1) + (-qJ(2) + (-t321 - t322) * pkin(7)) * t338 + t357;
t243 = -t289 * pkin(3) - t298 * pkin(8) + t300 * t295 + t288;
t210 = -t248 * pkin(4) - t276 * pkin(9) + t281 * t272 + t243;
t220 = t264 * qJD(5) + t329 * t248 + t334 * t249;
t198 = t210 + (-t264 * t318 - t220) * pkin(10) + (t265 * t318 - t219) * pkin(5);
t328 = sin(qJ(6));
t333 = cos(qJ(6));
t194 = -t328 * t197 + t333 * t198;
t250 = -t328 * t265 + t333 * t318;
t208 = t250 * qJD(6) + t333 * t220 + t328 * t317;
t218 = qJDD(6) - t219;
t251 = t333 * t265 + t328 * t318;
t226 = -t250 * mrSges(7,1) + t251 * mrSges(7,2);
t260 = qJD(6) - t264;
t229 = -t260 * mrSges(7,2) + t250 * mrSges(7,3);
t190 = m(7) * t194 + t218 * mrSges(7,1) - t208 * mrSges(7,3) - t251 * t226 + t260 * t229;
t195 = t333 * t197 + t328 * t198;
t207 = -t251 * qJD(6) - t328 * t220 + t333 * t317;
t230 = t260 * mrSges(7,1) - t251 * mrSges(7,3);
t191 = m(7) * t195 - t218 * mrSges(7,2) + t207 * mrSges(7,3) + t250 * t226 - t260 * t230;
t358 = -t328 * t190 + t333 * t191;
t177 = m(6) * t201 - t317 * mrSges(6,2) + t219 * mrSges(6,3) + t264 * t237 - t318 * t256 + t358;
t200 = t334 * t203 - t329 * t205;
t255 = -t318 * mrSges(6,2) + t264 * mrSges(6,3);
t196 = -t317 * pkin(5) - t316 * pkin(10) + t265 * t239 - t200;
t347 = -m(7) * t196 + t207 * mrSges(7,1) - t208 * mrSges(7,2) + t250 * t229 - t251 * t230;
t186 = m(6) * t200 + t317 * mrSges(6,1) - t220 * mrSges(6,3) - t265 * t237 + t318 * t255 + t347;
t172 = t329 * t177 + t334 * t186;
t266 = -t280 * mrSges(5,1) + t281 * mrSges(5,2);
t270 = -t323 * mrSges(5,2) + t280 * mrSges(5,3);
t169 = m(5) * t212 + t320 * mrSges(5,1) - t249 * mrSges(5,3) - t281 * t266 + t323 * t270 + t172;
t271 = t323 * mrSges(5,1) - t281 * mrSges(5,3);
t359 = t334 * t177 - t329 * t186;
t170 = m(5) * t213 - t320 * mrSges(5,2) + t248 * mrSges(5,3) + t280 * t266 - t323 * t271 + t359;
t163 = t335 * t169 + t330 * t170;
t284 = -t299 * mrSges(4,1) + t300 * mrSges(4,2);
t293 = -qJD(3) * mrSges(4,2) + t299 * mrSges(4,3);
t160 = m(4) * t253 + qJDD(3) * mrSges(4,1) - t290 * mrSges(4,3) + qJD(3) * t293 - t300 * t284 + t163;
t294 = qJD(3) * mrSges(4,1) - t300 * mrSges(4,3);
t360 = -t330 * t169 + t335 * t170;
t161 = m(4) * t254 - qJDD(3) * mrSges(4,2) + t289 * mrSges(4,3) - qJD(3) * t294 + t299 * t284 + t360;
t155 = t336 * t160 + t331 * t161;
t291 = -t326 * t301 + t363;
t278 = Ifges(4,4) * t300 + Ifges(4,2) * t299 + Ifges(4,6) * qJD(3);
t279 = Ifges(4,1) * t300 + Ifges(4,4) * t299 + Ifges(4,5) * qJD(3);
t258 = Ifges(5,4) * t281 + Ifges(5,2) * t280 + Ifges(5,6) * t323;
t259 = Ifges(5,1) * t281 + Ifges(5,4) * t280 + Ifges(5,5) * t323;
t221 = Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t260;
t223 = Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t260;
t183 = -mrSges(7,1) * t196 + mrSges(7,3) * t195 + Ifges(7,4) * t208 + Ifges(7,2) * t207 + Ifges(7,6) * t218 - t251 * t221 + t260 * t223;
t222 = Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t260;
t184 = mrSges(7,2) * t196 - mrSges(7,3) * t194 + Ifges(7,1) * t208 + Ifges(7,4) * t207 + Ifges(7,5) * t218 + t250 * t221 - t260 * t222;
t233 = Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t318;
t234 = Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t318;
t345 = -mrSges(6,1) * t200 + mrSges(6,2) * t201 - Ifges(6,5) * t220 - Ifges(6,6) * t219 - Ifges(6,3) * t317 - pkin(5) * t347 - pkin(10) * t358 - t333 * t183 - t328 * t184 - t265 * t233 + t264 * t234;
t342 = -mrSges(5,1) * t212 + mrSges(5,2) * t213 - Ifges(5,5) * t249 - Ifges(5,6) * t248 - Ifges(5,3) * t320 - pkin(4) * t172 - t281 * t258 + t280 * t259 + t345;
t340 = mrSges(4,1) * t253 - mrSges(4,2) * t254 + Ifges(4,5) * t290 + Ifges(4,6) * t289 + Ifges(4,3) * qJDD(3) + pkin(3) * t163 + t300 * t278 - t299 * t279 - t342;
t355 = Ifges(3,4) * t326 + Ifges(3,2) * t327;
t356 = Ifges(3,1) * t326 + Ifges(3,4) * t327;
t371 = -mrSges(3,1) * t291 + mrSges(3,2) * t292 - pkin(2) * t155 - (t326 * t355 - t327 * t356) * t338 - t340;
t369 = t326 * mrSges(3,2);
t179 = t333 * t190 + t328 * t191;
t354 = Ifges(3,5) * t326 + Ifges(3,6) * t327;
t367 = t338 * t354;
t350 = mrSges(3,3) * qJDD(1) + t338 * (-mrSges(3,1) * t327 + t369);
t153 = m(3) * t291 - t350 * t326 + t155;
t361 = -t331 * t160 + t336 * t161;
t154 = m(3) * t292 + t350 * t327 + t361;
t362 = -t326 * t153 + t327 * t154;
t349 = m(6) * t210 - t219 * mrSges(6,1) + t220 * mrSges(6,2) - t264 * t255 + t265 * t256 + t179;
t232 = Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t318;
t164 = mrSges(6,2) * t210 - mrSges(6,3) * t200 + Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t317 - pkin(10) * t179 - t328 * t183 + t333 * t184 + t264 * t232 - t318 * t233;
t344 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t208 + Ifges(7,6) * t207 + Ifges(7,3) * t218 + t251 * t222 - t250 * t223;
t165 = -mrSges(6,1) * t210 + mrSges(6,3) * t201 + Ifges(6,4) * t220 + Ifges(6,2) * t219 + Ifges(6,6) * t317 - pkin(5) * t179 - t265 * t232 + t318 * t234 - t344;
t257 = Ifges(5,5) * t281 + Ifges(5,6) * t280 + Ifges(5,3) * t323;
t151 = -mrSges(5,1) * t243 + mrSges(5,3) * t213 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t320 - pkin(4) * t349 + pkin(9) * t359 + t329 * t164 + t334 * t165 - t281 * t257 + t323 * t259;
t156 = mrSges(5,2) * t243 - mrSges(5,3) * t212 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t320 - pkin(9) * t172 + t334 * t164 - t329 * t165 + t280 * t257 - t323 * t258;
t277 = Ifges(4,5) * t300 + Ifges(4,6) * t299 + Ifges(4,3) * qJD(3);
t346 = m(5) * t243 - t248 * mrSges(5,1) + t249 * mrSges(5,2) - t280 * t270 + t281 * t271 + t349;
t146 = -mrSges(4,1) * t288 + mrSges(4,3) * t254 + Ifges(4,4) * t290 + Ifges(4,2) * t289 + Ifges(4,6) * qJDD(3) - pkin(3) * t346 + pkin(8) * t360 + qJD(3) * t279 + t335 * t151 + t330 * t156 - t300 * t277;
t147 = mrSges(4,2) * t288 - mrSges(4,3) * t253 + Ifges(4,1) * t290 + Ifges(4,4) * t289 + Ifges(4,5) * qJDD(3) - pkin(8) * t163 - qJD(3) * t278 - t330 * t151 + t335 * t156 + t299 * t277;
t297 = -qJDD(1) * pkin(1) - t338 * qJ(2) + t357;
t343 = m(4) * t288 - t289 * mrSges(4,1) + t290 * mrSges(4,2) - t299 * t293 + t300 * t294 + t346;
t142 = -mrSges(3,1) * t297 + mrSges(3,3) * t292 - pkin(2) * t343 + pkin(7) * t361 + t355 * qJDD(1) + t336 * t146 + t331 * t147 - t326 * t367;
t144 = mrSges(3,2) * t297 - mrSges(3,3) * t291 - pkin(7) * t155 + t356 * qJDD(1) - t331 * t146 + t336 * t147 + t327 * t367;
t341 = -m(3) * t297 + mrSges(3,1) * t364 - t343 + (t321 * t338 + t368) * mrSges(3,3);
t348 = -mrSges(2,2) * t308 + qJ(2) * t362 + t327 * t142 + t326 * t144 + pkin(1) * (-qJDD(1) * t369 + t341) + mrSges(2,1) * t307 + Ifges(2,3) * qJDD(1);
t173 = t341 + (mrSges(2,1) - t369) * qJDD(1) - t338 * mrSges(2,2) + m(2) * t307;
t150 = t327 * t153 + t326 * t154;
t148 = m(2) * t308 - t338 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t362;
t145 = mrSges(2,1) * g(3) - pkin(1) * t150 + (Ifges(2,6) - t354) * qJDD(1) + t338 * Ifges(2,5) + mrSges(2,3) * t308 + t371;
t140 = -mrSges(2,2) * g(3) - mrSges(2,3) * t307 + Ifges(2,5) * qJDD(1) - t338 * Ifges(2,6) - qJ(2) * t150 - t326 * t142 + t327 * t144;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t140 - t332 * t145 - pkin(6) * (t332 * t148 + t337 * t173), t140, t144, t147, t156, t164, t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t332 * t140 + t337 * t145 + pkin(6) * (t337 * t148 - t332 * t173), t145, t142, t146, t151, t165, t183; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t348, t348, t354 * qJDD(1) - t371, t340, -t342, -t345, t344;];
m_new  = t1;
