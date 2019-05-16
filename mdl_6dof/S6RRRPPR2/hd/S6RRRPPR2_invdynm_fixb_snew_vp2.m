% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:19:59
% EndTime: 2019-05-07 04:20:29
% DurationCPUTime: 17.62s
% Computational Cost: add. (283991->389), mult. (647485->470), div. (0->0), fcn. (469668->10), ass. (0->149)
t340 = sin(qJ(2));
t344 = cos(qJ(2));
t369 = qJD(1) * qJD(2);
t319 = qJDD(1) * t340 + t344 * t369;
t341 = sin(qJ(1));
t345 = cos(qJ(1));
t326 = -g(1) * t345 - g(2) * t341;
t346 = qJD(1) ^ 2;
t314 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t326;
t376 = t314 * t340;
t380 = pkin(2) * t346;
t269 = qJDD(2) * pkin(2) - pkin(8) * t319 - t376 + (pkin(8) * t369 + t340 * t380 - g(3)) * t344;
t301 = -g(3) * t340 + t344 * t314;
t320 = qJDD(1) * t344 - t340 * t369;
t372 = qJD(1) * t340;
t324 = qJD(2) * pkin(2) - pkin(8) * t372;
t336 = t344 ^ 2;
t270 = pkin(8) * t320 - qJD(2) * t324 - t336 * t380 + t301;
t339 = sin(qJ(3));
t343 = cos(qJ(3));
t234 = t343 * t269 - t270 * t339;
t311 = (-t339 * t340 + t343 * t344) * qJD(1);
t279 = qJD(3) * t311 + t319 * t343 + t320 * t339;
t312 = (t339 * t344 + t340 * t343) * qJD(1);
t333 = qJDD(2) + qJDD(3);
t334 = qJD(2) + qJD(3);
t213 = (t311 * t334 - t279) * qJ(4) + (t311 * t312 + t333) * pkin(3) + t234;
t235 = t339 * t269 + t343 * t270;
t278 = -qJD(3) * t312 - t319 * t339 + t320 * t343;
t303 = pkin(3) * t334 - qJ(4) * t312;
t307 = t311 ^ 2;
t216 = -pkin(3) * t307 + qJ(4) * t278 - t303 * t334 + t235;
t337 = sin(pkin(10));
t378 = cos(pkin(10));
t298 = t337 * t311 + t312 * t378;
t384 = -2 * qJD(4);
t210 = t213 * t378 - t337 * t216 + t298 * t384;
t297 = -t311 * t378 + t312 * t337;
t288 = t297 * t384;
t375 = t337 * t213 + t378 * t216;
t211 = t288 + t375;
t247 = -t278 * t378 + t279 * t337;
t248 = t337 * t278 + t279 * t378;
t255 = Ifges(5,4) * t298 - Ifges(5,2) * t297 + Ifges(5,6) * t334;
t264 = -mrSges(6,2) * t297 - mrSges(6,3) * t298;
t283 = mrSges(6,1) * t297 - mrSges(6,3) * t334;
t262 = pkin(4) * t297 - qJ(5) * t298;
t332 = t334 ^ 2;
t206 = -t333 * pkin(4) - t332 * qJ(5) + t298 * t262 + qJDD(5) - t210;
t377 = t297 * t334;
t200 = (t297 * t298 - t333) * pkin(9) + (t248 + t377) * pkin(5) + t206;
t285 = pkin(5) * t298 - pkin(9) * t334;
t294 = t297 ^ 2;
t325 = t341 * g(1) - t345 * g(2);
t362 = -qJDD(1) * pkin(1) - t325;
t280 = -t320 * pkin(2) + t324 * t372 + (-pkin(8) * t336 - pkin(7)) * t346 + t362;
t224 = -t278 * pkin(3) - t307 * qJ(4) + t312 * t303 + qJDD(4) + t280;
t381 = -2 * qJD(5);
t350 = (-t248 + t377) * qJ(5) + t224 + (t334 * pkin(4) + t381) * t298;
t203 = (pkin(4) + pkin(9)) * t247 + t350 - t298 * t285 - t294 * pkin(5);
t338 = sin(qJ(6));
t342 = cos(qJ(6));
t197 = t200 * t342 - t203 * t338;
t276 = t297 * t342 - t334 * t338;
t222 = qJD(6) * t276 + t247 * t338 + t333 * t342;
t246 = qJDD(6) + t248;
t277 = t297 * t338 + t334 * t342;
t249 = -mrSges(7,1) * t276 + mrSges(7,2) * t277;
t290 = qJD(6) + t298;
t250 = -mrSges(7,2) * t290 + mrSges(7,3) * t276;
t194 = m(7) * t197 + mrSges(7,1) * t246 - mrSges(7,3) * t222 - t249 * t277 + t250 * t290;
t198 = t200 * t338 + t203 * t342;
t221 = -qJD(6) * t277 + t247 * t342 - t333 * t338;
t251 = mrSges(7,1) * t290 - mrSges(7,3) * t277;
t195 = m(7) * t198 - mrSges(7,2) * t246 + mrSges(7,3) * t221 + t249 * t276 - t251 * t290;
t182 = t194 * t342 + t195 * t338;
t361 = pkin(4) * t332 - qJ(5) * t333 - t375;
t202 = -pkin(5) * t247 - pkin(9) * t294 - t262 * t297 + t288 + ((2 * qJD(5)) + t285) * t334 - t361;
t225 = Ifges(7,5) * t277 + Ifges(7,6) * t276 + Ifges(7,3) * t290;
t227 = Ifges(7,1) * t277 + Ifges(7,4) * t276 + Ifges(7,5) * t290;
t185 = -mrSges(7,1) * t202 + mrSges(7,3) * t198 + Ifges(7,4) * t222 + Ifges(7,2) * t221 + Ifges(7,6) * t246 - t225 * t277 + t227 * t290;
t226 = Ifges(7,4) * t277 + Ifges(7,2) * t276 + Ifges(7,6) * t290;
t186 = mrSges(7,2) * t202 - mrSges(7,3) * t197 + Ifges(7,1) * t222 + Ifges(7,4) * t221 + Ifges(7,5) * t246 + t225 * t276 - t226 * t290;
t204 = t334 * t381 + ((2 * qJD(4)) + t262) * t297 + t361;
t252 = Ifges(6,5) * t334 - Ifges(6,6) * t298 + Ifges(6,3) * t297;
t355 = -mrSges(6,2) * t206 + mrSges(6,3) * t204 - Ifges(6,1) * t333 + Ifges(6,4) * t248 - Ifges(6,5) * t247 + pkin(9) * t182 + t338 * t185 - t342 * t186 + t298 * t252;
t199 = -m(7) * t202 + mrSges(7,1) * t221 - t222 * mrSges(7,2) + t250 * t276 - t277 * t251;
t284 = mrSges(6,1) * t298 + mrSges(6,2) * t334;
t356 = -m(6) * t204 + t333 * mrSges(6,3) + t334 * t284 - t199;
t359 = -m(6) * t206 - t248 * mrSges(6,1) - t298 * t264 - t182;
t254 = Ifges(6,4) * t334 - Ifges(6,2) * t298 + Ifges(6,6) * t297;
t373 = Ifges(5,1) * t298 - Ifges(5,4) * t297 + Ifges(5,5) * t334 - t254;
t386 = -mrSges(5,2) * t211 + pkin(4) * (-mrSges(6,2) * t333 - t283 * t334 + t359) + qJ(5) * (-mrSges(6,1) * t247 - t264 * t297 + t356) + mrSges(5,1) * t210 + t298 * t255 - Ifges(5,6) * t247 + Ifges(5,5) * t248 + Ifges(5,3) * t333 - t355 + t373 * t297;
t263 = mrSges(5,1) * t297 + mrSges(5,2) * t298;
t281 = -mrSges(5,2) * t334 - mrSges(5,3) * t297;
t178 = m(5) * t210 - mrSges(5,3) * t248 - t263 * t298 + (t281 - t283) * t334 + (mrSges(5,1) - mrSges(6,2)) * t333 + t359;
t282 = mrSges(5,1) * t334 - mrSges(5,3) * t298;
t189 = m(5) * t211 - mrSges(5,2) * t333 - t282 * t334 + (-t263 - t264) * t297 + (-mrSges(5,3) - mrSges(6,1)) * t247 + t356;
t174 = t378 * t178 + t337 * t189;
t292 = Ifges(4,4) * t312 + Ifges(4,2) * t311 + Ifges(4,6) * t334;
t293 = Ifges(4,1) * t312 + Ifges(4,4) * t311 + Ifges(4,5) * t334;
t385 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t279 + Ifges(4,6) * t278 + Ifges(4,3) * t333 + pkin(3) * t174 + t312 * t292 - t311 * t293 + t386;
t299 = -mrSges(4,1) * t311 + mrSges(4,2) * t312;
t302 = -mrSges(4,2) * t334 + mrSges(4,3) * t311;
t171 = m(4) * t234 + mrSges(4,1) * t333 - mrSges(4,3) * t279 - t299 * t312 + t302 * t334 + t174;
t304 = mrSges(4,1) * t334 - mrSges(4,3) * t312;
t365 = -t178 * t337 + t378 * t189;
t172 = m(4) * t235 - mrSges(4,2) * t333 + mrSges(4,3) * t278 + t299 * t311 - t304 * t334 + t365;
t166 = t343 * t171 + t339 * t172;
t300 = -g(3) * t344 - t376;
t309 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t340 + Ifges(3,2) * t344) * qJD(1);
t310 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t340 + Ifges(3,4) * t344) * qJD(1);
t383 = mrSges(3,1) * t300 - mrSges(3,2) * t301 + Ifges(3,5) * t319 + Ifges(3,6) * t320 + Ifges(3,3) * qJDD(2) + pkin(2) * t166 + (t340 * t309 - t344 * t310) * qJD(1) + t385;
t379 = Ifges(5,4) + Ifges(6,6);
t183 = -t338 * t194 + t342 * t195;
t256 = Ifges(6,1) * t334 - Ifges(6,4) * t298 + Ifges(6,5) * t297;
t374 = -Ifges(5,5) * t298 + Ifges(5,6) * t297 - Ifges(5,3) * t334 - t256;
t371 = qJD(1) * t344;
t318 = (-mrSges(3,1) * t344 + mrSges(3,2) * t340) * qJD(1);
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t371;
t164 = m(3) * t300 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t319 + qJD(2) * t323 - t318 * t372 + t166;
t322 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t372;
t366 = -t339 * t171 + t343 * t172;
t165 = m(3) * t301 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t320 - qJD(2) * t322 + t318 * t371 + t366;
t367 = -t164 * t340 + t344 * t165;
t208 = t247 * pkin(4) + t350;
t181 = m(6) * t208 - t247 * mrSges(6,2) - t248 * mrSges(6,3) - t297 * t283 - t298 * t284 + t183;
t354 = -mrSges(6,1) * t204 + mrSges(6,2) * t208 - pkin(5) * t199 - pkin(9) * t183 - t342 * t185 - t338 * t186;
t162 = -mrSges(5,1) * t224 + mrSges(5,3) * t211 - pkin(4) * t181 + t373 * t334 + (Ifges(5,6) - Ifges(6,5)) * t333 + t374 * t298 + t379 * t248 + (-Ifges(5,2) - Ifges(6,3)) * t247 + t354;
t358 = mrSges(7,1) * t197 - mrSges(7,2) * t198 + Ifges(7,5) * t222 + Ifges(7,6) * t221 + Ifges(7,3) * t246 + t277 * t226 - t276 * t227;
t353 = mrSges(6,1) * t206 - mrSges(6,3) * t208 + pkin(5) * t182 + t358;
t167 = t353 + (-t255 + t252) * t334 + (Ifges(5,5) - Ifges(6,4)) * t333 + t374 * t297 + (Ifges(5,1) + Ifges(6,2)) * t248 - t379 * t247 + mrSges(5,2) * t224 - mrSges(5,3) * t210 - qJ(5) * t181;
t291 = Ifges(4,5) * t312 + Ifges(4,6) * t311 + Ifges(4,3) * t334;
t357 = m(5) * t224 + t247 * mrSges(5,1) + t248 * mrSges(5,2) + t297 * t281 + t298 * t282 + t181;
t157 = -mrSges(4,1) * t280 + mrSges(4,3) * t235 + Ifges(4,4) * t279 + Ifges(4,2) * t278 + Ifges(4,6) * t333 - pkin(3) * t357 + qJ(4) * t365 + t162 * t378 + t337 * t167 - t312 * t291 + t334 * t293;
t158 = mrSges(4,2) * t280 - mrSges(4,3) * t234 + Ifges(4,1) * t279 + Ifges(4,4) * t278 + Ifges(4,5) * t333 - qJ(4) * t174 - t337 * t162 + t167 * t378 + t311 * t291 - t334 * t292;
t308 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t340 + Ifges(3,6) * t344) * qJD(1);
t313 = -pkin(7) * t346 + t362;
t352 = m(4) * t280 - t278 * mrSges(4,1) + t279 * mrSges(4,2) - t311 * t302 + t312 * t304 + t357;
t153 = -mrSges(3,1) * t313 + mrSges(3,3) * t301 + Ifges(3,4) * t319 + Ifges(3,2) * t320 + Ifges(3,6) * qJDD(2) - pkin(2) * t352 + pkin(8) * t366 + qJD(2) * t310 + t343 * t157 + t339 * t158 - t308 * t372;
t155 = mrSges(3,2) * t313 - mrSges(3,3) * t300 + Ifges(3,1) * t319 + Ifges(3,4) * t320 + Ifges(3,5) * qJDD(2) - pkin(8) * t166 - qJD(2) * t309 - t157 * t339 + t158 * t343 + t308 * t371;
t349 = -m(3) * t313 + t320 * mrSges(3,1) - t319 * mrSges(3,2) - t322 * t372 + t323 * t371 - t352;
t360 = mrSges(2,1) * t325 - mrSges(2,2) * t326 + Ifges(2,3) * qJDD(1) + pkin(1) * t349 + pkin(7) * t367 + t344 * t153 + t340 * t155;
t175 = m(2) * t325 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t349;
t161 = t164 * t344 + t165 * t340;
t159 = m(2) * t326 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t367;
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t326 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t161 - t383;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t325 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t161 - t153 * t340 + t155 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t151 - t341 * t156 - pkin(6) * (t159 * t341 + t175 * t345), t151, t155, t158, t167, -t297 * t254 - t355, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t341 * t151 + t345 * t156 + pkin(6) * (t159 * t345 - t341 * t175), t156, t153, t157, t162, Ifges(6,4) * t333 - Ifges(6,2) * t248 + Ifges(6,6) * t247 - t334 * t252 + t297 * t256 - t353, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t360, t360, t383, t385, t386, Ifges(6,5) * t333 - Ifges(6,6) * t248 + Ifges(6,3) * t247 + t334 * t254 + t298 * t256 - t354, t358;];
m_new  = t1;
