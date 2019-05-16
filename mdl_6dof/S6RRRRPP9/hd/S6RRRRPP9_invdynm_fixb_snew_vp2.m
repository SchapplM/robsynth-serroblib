% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:12:04
% EndTime: 2019-05-07 19:12:39
% DurationCPUTime: 15.56s
% Computational Cost: add. (271519->402), mult. (575528->479), div. (0->0), fcn. (445123->10), ass. (0->154)
t330 = sin(pkin(6));
t334 = sin(qJ(2));
t337 = cos(qJ(2));
t360 = qJD(1) * qJD(2);
t316 = (-qJDD(1) * t337 + t334 * t360) * t330;
t362 = qJD(1) * t330;
t314 = (-pkin(2) * t337 - pkin(9) * t334) * t362;
t331 = cos(pkin(6));
t326 = qJD(1) * t331 + qJD(2);
t324 = t326 ^ 2;
t325 = qJDD(1) * t331 + qJDD(2);
t361 = qJD(1) * t337;
t335 = sin(qJ(1));
t338 = cos(qJ(1));
t322 = t335 * g(1) - g(2) * t338;
t339 = qJD(1) ^ 2;
t377 = pkin(8) * t330;
t311 = qJDD(1) * pkin(1) + t339 * t377 + t322;
t323 = -g(1) * t338 - g(2) * t335;
t312 = -pkin(1) * t339 + qJDD(1) * t377 + t323;
t369 = t331 * t334;
t363 = t311 * t369 + t337 * t312;
t252 = -t324 * pkin(2) + t325 * pkin(9) + (-g(3) * t334 + t314 * t361) * t330 + t363;
t315 = (qJDD(1) * t334 + t337 * t360) * t330;
t376 = t331 * g(3);
t253 = t316 * pkin(2) - t315 * pkin(9) - t376 + (-t311 + (pkin(2) * t334 - pkin(9) * t337) * t326 * qJD(1)) * t330;
t333 = sin(qJ(3));
t336 = cos(qJ(3));
t213 = t336 * t252 + t333 * t253;
t358 = t334 * t362;
t304 = t336 * t326 - t333 * t358;
t305 = t326 * t333 + t336 * t358;
t288 = -pkin(3) * t304 - pkin(10) * t305;
t308 = qJDD(3) + t316;
t357 = t330 * t361;
t321 = qJD(3) - t357;
t320 = t321 ^ 2;
t209 = -pkin(3) * t320 + pkin(10) * t308 + t288 * t304 + t213;
t368 = t331 * t337;
t370 = t330 * t337;
t285 = -g(3) * t370 + t311 * t368 - t334 * t312;
t251 = -t325 * pkin(2) - t324 * pkin(9) + t314 * t358 - t285;
t283 = -qJD(3) * t305 - t315 * t333 + t325 * t336;
t284 = qJD(3) * t304 + t315 * t336 + t325 * t333;
t211 = (-t304 * t321 - t284) * pkin(10) + (t305 * t321 - t283) * pkin(3) + t251;
t332 = sin(qJ(4));
t379 = cos(qJ(4));
t204 = -t332 * t209 + t379 * t211;
t291 = t305 * t332 - t379 * t321;
t292 = t379 * t305 + t332 * t321;
t258 = pkin(4) * t291 - qJ(5) * t292;
t281 = qJDD(4) - t283;
t302 = qJD(4) - t304;
t301 = t302 ^ 2;
t202 = -t281 * pkin(4) - t301 * qJ(5) + t292 * t258 + qJDD(5) - t204;
t232 = -t291 * qJD(4) + t379 * t284 + t332 * t308;
t260 = -mrSges(6,2) * t291 - mrSges(6,3) * t292;
t382 = -m(6) * t202 - t232 * mrSges(6,1) - t292 * t260;
t257 = -mrSges(7,2) * t292 + mrSges(7,3) * t291;
t373 = t291 * t302;
t193 = -0.2e1 * qJD(6) * t302 + (t291 * t292 - t281) * qJ(6) + (t232 + t373) * pkin(5) + t202;
t265 = -mrSges(7,1) * t291 + mrSges(7,2) * t302;
t355 = -m(7) * t193 + t281 * mrSges(7,3) + t302 * t265;
t189 = t232 * mrSges(7,1) + t292 * t257 - t355;
t205 = t379 * t209 + t332 * t211;
t231 = qJD(4) * t292 + t284 * t332 - t379 * t308;
t238 = Ifges(7,4) * t302 + Ifges(7,2) * t291 + Ifges(7,6) * t292;
t240 = Ifges(5,4) * t292 - Ifges(5,2) * t291 + Ifges(5,6) * t302;
t264 = mrSges(6,1) * t291 - mrSges(6,3) * t302;
t348 = -t301 * pkin(4) + t281 * qJ(5) - t291 * t258 + t205;
t380 = -2 * qJD(5);
t200 = t302 * t380 - t348;
t236 = Ifges(6,5) * t302 - Ifges(6,6) * t292 + Ifges(6,3) * t291;
t262 = pkin(5) * t292 - qJ(6) * t302;
t290 = t291 ^ 2;
t196 = -t231 * pkin(5) - t290 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t262) * t302 + t348;
t235 = Ifges(7,5) * t302 + Ifges(7,6) * t291 + Ifges(7,3) * t292;
t350 = -mrSges(7,2) * t196 + mrSges(7,3) * t193 - Ifges(7,1) * t281 - Ifges(7,4) * t231 - Ifges(7,5) * t232 - t291 * t235;
t342 = -mrSges(6,2) * t202 + mrSges(6,3) * t200 - Ifges(6,1) * t281 + Ifges(6,4) * t232 - Ifges(6,5) * t231 + qJ(6) * t189 + t292 * t236 + t350;
t266 = mrSges(6,1) * t292 + mrSges(6,2) * t302;
t263 = mrSges(7,1) * t292 - mrSges(7,3) * t302;
t359 = m(7) * t196 + t281 * mrSges(7,2) + t302 * t263;
t353 = -m(6) * t200 + t281 * mrSges(6,3) + t302 * t266 + t359;
t364 = -t257 - t260;
t239 = Ifges(6,4) * t302 - Ifges(6,2) * t292 + Ifges(6,6) * t291;
t365 = Ifges(5,1) * t292 - Ifges(5,4) * t291 + Ifges(5,5) * t302 - t239;
t381 = t365 * t291 + (t240 - t238) * t292 + mrSges(5,1) * t204 - mrSges(5,2) * t205 + Ifges(5,5) * t232 - Ifges(5,6) * t231 + Ifges(5,3) * t281 + pkin(4) * (-t281 * mrSges(6,2) - t302 * t264 - t189 + t382) + qJ(5) * (t364 * t291 + (-mrSges(6,1) - mrSges(7,1)) * t231 + t353) - t342;
t375 = -mrSges(7,1) - mrSges(5,3);
t374 = Ifges(5,4) + Ifges(6,6);
t372 = t292 * t238;
t371 = t330 * t334;
t259 = mrSges(5,1) * t291 + mrSges(5,2) * t292;
t267 = -mrSges(5,2) * t302 - mrSges(5,3) * t291;
t182 = m(5) * t204 + (-t264 + t267) * t302 + (-t257 - t259) * t292 + (mrSges(5,1) - mrSges(6,2)) * t281 + t375 * t232 + t355 + t382;
t268 = mrSges(5,1) * t302 - mrSges(5,3) * t292;
t184 = m(5) * t205 - t281 * mrSges(5,2) - t302 * t268 + (-t259 + t364) * t291 + (-mrSges(6,1) + t375) * t231 + t353;
t179 = -t182 * t332 + t379 * t184;
t287 = -mrSges(4,1) * t304 + mrSges(4,2) * t305;
t294 = mrSges(4,1) * t321 - mrSges(4,3) * t305;
t177 = m(4) * t213 - mrSges(4,2) * t308 + mrSges(4,3) * t283 + t287 * t304 - t294 * t321 + t179;
t212 = -t333 * t252 + t336 * t253;
t208 = -t308 * pkin(3) - t320 * pkin(10) + t305 * t288 - t212;
t344 = (-t232 + t373) * qJ(5) + t208 + (t302 * pkin(4) + t380) * t292;
t199 = -t290 * pkin(5) + 0.2e1 * qJD(6) * t291 - t292 * t262 + (pkin(4) + qJ(6)) * t231 + t344;
t190 = m(7) * t199 - t232 * mrSges(7,2) + t231 * mrSges(7,3) - t292 * t263 + t291 * t265;
t203 = t231 * pkin(4) + t344;
t349 = -m(6) * t203 + t231 * mrSges(6,2) + t291 * t264 - t190;
t185 = -m(5) * t208 - t231 * mrSges(5,1) - t291 * t267 + (t266 - t268) * t292 + (-mrSges(5,2) + mrSges(6,3)) * t232 + t349;
t293 = -mrSges(4,2) * t321 + mrSges(4,3) * t304;
t181 = m(4) * t212 + t308 * mrSges(4,1) - t284 * mrSges(4,3) - t305 * t287 + t321 * t293 + t185;
t171 = t333 * t177 + t336 * t181;
t241 = Ifges(7,1) * t302 + Ifges(7,4) * t291 + Ifges(7,5) * t292;
t242 = Ifges(6,1) * t302 - Ifges(6,4) * t292 + Ifges(6,5) * t291;
t366 = t241 + t242;
t286 = -g(3) * t371 + t363;
t309 = mrSges(3,1) * t326 - mrSges(3,3) * t358;
t313 = (-mrSges(3,1) * t337 + mrSges(3,2) * t334) * t362;
t356 = t336 * t177 - t181 * t333;
t169 = m(3) * t286 - mrSges(3,2) * t325 - mrSges(3,3) * t316 - t309 * t326 + t313 * t357 + t356;
t310 = -mrSges(3,2) * t326 + mrSges(3,3) * t357;
t178 = t379 * t182 + t332 * t184;
t345 = -m(4) * t251 + t283 * mrSges(4,1) - t284 * mrSges(4,2) + t304 * t293 - t305 * t294 - t178;
t174 = m(3) * t285 + t325 * mrSges(3,1) - t315 * mrSges(3,3) + t326 * t310 - t313 * t358 + t345;
t165 = t337 * t169 - t174 * t334;
t298 = -t330 * t311 - t376;
t170 = m(3) * t298 + t316 * mrSges(3,1) + t315 * mrSges(3,2) + (t309 * t334 - t310 * t337) * t362 + t171;
t160 = t169 * t369 - t170 * t330 + t174 * t368;
t352 = mrSges(7,1) * t196 - mrSges(7,3) * t199 - Ifges(7,4) * t281 - Ifges(7,2) * t231 - Ifges(7,6) * t232 - t292 * t241;
t351 = -mrSges(7,1) * t193 + mrSges(7,2) * t199 - Ifges(7,5) * t281 - Ifges(7,6) * t231 - Ifges(7,3) * t232 - t302 * t238;
t188 = -t232 * mrSges(6,3) - t292 * t266 - t349;
t237 = Ifges(5,5) * t292 - Ifges(5,6) * t291 + Ifges(5,3) * t302;
t343 = mrSges(6,1) * t200 - mrSges(6,2) * t203 + pkin(5) * (t231 * mrSges(7,1) + t291 * t257 - t359) + qJ(6) * t190 - t352;
t166 = -t343 + (t235 + t365) * t302 + mrSges(5,3) * t205 - mrSges(5,1) * t208 + (-Ifges(5,2) - Ifges(6,3)) * t231 + t374 * t232 - pkin(4) * t188 + (-t237 - t242) * t292 + (Ifges(5,6) - Ifges(6,5)) * t281;
t346 = -mrSges(6,1) * t202 + mrSges(6,3) * t203 - pkin(5) * t189 + t351;
t172 = -t346 + (-t240 + t236) * t302 + (-t237 - t366) * t291 + (Ifges(5,5) - Ifges(6,4)) * t281 + (Ifges(5,1) + Ifges(6,2)) * t232 - t374 * t231 - mrSges(5,3) * t204 + mrSges(5,2) * t208 - qJ(5) * t188;
t277 = Ifges(4,5) * t305 + Ifges(4,6) * t304 + Ifges(4,3) * t321;
t278 = Ifges(4,4) * t305 + Ifges(4,2) * t304 + Ifges(4,6) * t321;
t161 = mrSges(4,2) * t251 - mrSges(4,3) * t212 + Ifges(4,1) * t284 + Ifges(4,4) * t283 + Ifges(4,5) * t308 - pkin(10) * t178 - t332 * t166 + t379 * t172 + t304 * t277 - t321 * t278;
t279 = Ifges(4,1) * t305 + Ifges(4,4) * t304 + Ifges(4,5) * t321;
t162 = -mrSges(4,1) * t251 + mrSges(4,3) * t213 + Ifges(4,4) * t284 + Ifges(4,2) * t283 + Ifges(4,6) * t308 - pkin(3) * t178 - t305 * t277 + t321 * t279 - t381;
t296 = Ifges(3,6) * t326 + (Ifges(3,4) * t334 + Ifges(3,2) * t337) * t362;
t297 = Ifges(3,5) * t326 + (Ifges(3,1) * t334 + Ifges(3,4) * t337) * t362;
t152 = Ifges(3,5) * t315 - Ifges(3,6) * t316 + Ifges(3,3) * t325 + mrSges(3,1) * t285 - mrSges(3,2) * t286 + t333 * t161 + t336 * t162 + pkin(2) * t345 + pkin(9) * t356 + (t296 * t334 - t297 * t337) * t362;
t295 = Ifges(3,3) * t326 + (Ifges(3,5) * t334 + Ifges(3,6) * t337) * t362;
t154 = mrSges(3,2) * t298 - mrSges(3,3) * t285 + Ifges(3,1) * t315 - Ifges(3,4) * t316 + Ifges(3,5) * t325 - pkin(9) * t171 + t161 * t336 - t162 * t333 + t295 * t357 - t296 * t326;
t341 = mrSges(4,1) * t212 - mrSges(4,2) * t213 + Ifges(4,5) * t284 + Ifges(4,6) * t283 + Ifges(4,3) * t308 + pkin(3) * t185 + pkin(10) * t179 + t379 * t166 + t332 * t172 + t305 * t278 - t304 * t279;
t156 = -mrSges(3,1) * t298 + mrSges(3,3) * t286 + Ifges(3,4) * t315 - Ifges(3,2) * t316 + Ifges(3,6) * t325 - pkin(2) * t171 - t295 * t358 + t326 * t297 - t341;
t347 = mrSges(2,1) * t322 - mrSges(2,2) * t323 + Ifges(2,3) * qJDD(1) + pkin(1) * t160 + t331 * t152 + t154 * t371 + t156 * t370 + t165 * t377;
t163 = m(2) * t323 - mrSges(2,1) * t339 - qJDD(1) * mrSges(2,2) + t165;
t159 = t331 * t170 + (t169 * t334 + t174 * t337) * t330;
t157 = m(2) * t322 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t339 + t160;
t150 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - t339 * Ifges(2,6) + t337 * t154 - t334 * t156 + (-t159 * t330 - t160 * t331) * pkin(8);
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t323 + t339 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t159 - t330 * t152 + (pkin(8) * t165 + t154 * t334 + t156 * t337) * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t338 * t150 - t335 * t149 - pkin(7) * (t157 * t338 + t163 * t335), t150, t154, t161, t172, -t291 * t239 - t342 - t372, -t350 - t372; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t150 + t338 * t149 + pkin(7) * (-t157 * t335 + t163 * t338), t149, t156, t162, t166, Ifges(6,4) * t281 - Ifges(6,2) * t232 + Ifges(6,6) * t231 - t302 * t236 + t366 * t291 + t346, -t302 * t235 - t352; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t347, t347, t152, t341, t381, t343 + (-t235 + t239) * t302 + t292 * t242 + Ifges(6,5) * t281 - Ifges(6,6) * t232 + Ifges(6,3) * t231, -t291 * t241 - t351;];
m_new  = t1;
