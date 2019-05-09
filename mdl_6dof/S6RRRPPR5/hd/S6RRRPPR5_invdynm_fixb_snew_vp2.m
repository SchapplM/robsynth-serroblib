% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 05:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:00:07
% EndTime: 2019-05-07 05:02:00
% DurationCPUTime: 80.14s
% Computational Cost: add. (1465710->398), mult. (3265808->519), div. (0->0), fcn. (2643381->14), ass. (0->164)
t388 = -2 * qJD(4);
t347 = sin(pkin(6));
t352 = sin(qJ(2));
t356 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t333 = (-qJDD(1) * t356 + t352 * t376) * t347;
t379 = qJD(1) * t347;
t331 = (-pkin(2) * t356 - pkin(9) * t352) * t379;
t349 = cos(pkin(6));
t341 = qJD(1) * t349 + qJD(2);
t339 = t341 ^ 2;
t340 = qJDD(1) * t349 + qJDD(2);
t378 = qJD(1) * t356;
t353 = sin(qJ(1));
t357 = cos(qJ(1));
t337 = t353 * g(1) - g(2) * t357;
t358 = qJD(1) ^ 2;
t387 = pkin(8) * t347;
t328 = qJDD(1) * pkin(1) + t358 * t387 + t337;
t338 = -g(1) * t357 - g(2) * t353;
t329 = -pkin(1) * t358 + qJDD(1) * t387 + t338;
t382 = t349 * t352;
t380 = t328 * t382 + t356 * t329;
t283 = -t339 * pkin(2) + t340 * pkin(9) + (-g(3) * t352 + t331 * t378) * t347 + t380;
t332 = (qJDD(1) * t352 + t356 * t376) * t347;
t386 = t349 * g(3);
t284 = t333 * pkin(2) - t332 * pkin(9) - t386 + (-t328 + (pkin(2) * t352 - pkin(9) * t356) * t341 * qJD(1)) * t347;
t351 = sin(qJ(3));
t355 = cos(qJ(3));
t249 = -t351 * t283 + t355 * t284;
t375 = t352 * t379;
t320 = t341 * t355 - t351 * t375;
t300 = qJD(3) * t320 + t332 * t355 + t340 * t351;
t321 = t341 * t351 + t355 * t375;
t325 = qJDD(3) + t333;
t374 = t347 * t378;
t336 = qJD(3) - t374;
t237 = (t320 * t336 - t300) * qJ(4) + (t320 * t321 + t325) * pkin(3) + t249;
t250 = t355 * t283 + t351 * t284;
t299 = -qJD(3) * t321 - t332 * t351 + t340 * t355;
t310 = pkin(3) * t336 - qJ(4) * t321;
t319 = t320 ^ 2;
t244 = -pkin(3) * t319 + qJ(4) * t299 - t310 * t336 + t250;
t346 = sin(pkin(11));
t385 = cos(pkin(11));
t307 = t346 * t320 + t385 * t321;
t225 = t385 * t237 - t346 * t244 + t307 * t388;
t384 = t347 * t352;
t383 = t347 * t356;
t381 = t349 * t356;
t306 = -t385 * t320 + t321 * t346;
t226 = t346 * t237 + t385 * t244 + t306 * t388;
t271 = -t385 * t299 + t300 * t346;
t278 = mrSges(5,1) * t306 + mrSges(5,2) * t307;
t291 = mrSges(5,1) * t336 - mrSges(5,3) * t307;
t277 = pkin(4) * t306 - qJ(5) * t307;
t335 = t336 ^ 2;
t223 = -pkin(4) * t335 + qJ(5) * t325 - t277 * t306 + t226;
t301 = -g(3) * t383 + t328 * t381 - t352 * t329;
t282 = -t340 * pkin(2) - t339 * pkin(9) + t331 * t375 - t301;
t246 = -t299 * pkin(3) - t319 * qJ(4) + t321 * t310 + qJDD(4) + t282;
t272 = t346 * t299 + t385 * t300;
t229 = (t306 * t336 - t272) * qJ(5) + (t307 * t336 + t271) * pkin(4) + t246;
t345 = sin(pkin(12));
t348 = cos(pkin(12));
t289 = t307 * t348 + t336 * t345;
t218 = -0.2e1 * qJD(5) * t289 - t345 * t223 + t348 * t229;
t259 = t272 * t348 + t325 * t345;
t288 = -t307 * t345 + t336 * t348;
t216 = (t288 * t306 - t259) * pkin(10) + (t288 * t289 + t271) * pkin(5) + t218;
t219 = 0.2e1 * qJD(5) * t288 + t348 * t223 + t345 * t229;
t258 = -t272 * t345 + t325 * t348;
t265 = pkin(5) * t306 - pkin(10) * t289;
t287 = t288 ^ 2;
t217 = -pkin(5) * t287 + pkin(10) * t258 - t265 * t306 + t219;
t350 = sin(qJ(6));
t354 = cos(qJ(6));
t214 = t216 * t354 - t217 * t350;
t260 = t288 * t354 - t289 * t350;
t235 = qJD(6) * t260 + t258 * t350 + t259 * t354;
t261 = t288 * t350 + t289 * t354;
t245 = -mrSges(7,1) * t260 + mrSges(7,2) * t261;
t305 = qJD(6) + t306;
t247 = -mrSges(7,2) * t305 + mrSges(7,3) * t260;
t270 = qJDD(6) + t271;
t209 = m(7) * t214 + mrSges(7,1) * t270 - mrSges(7,3) * t235 - t245 * t261 + t247 * t305;
t215 = t216 * t350 + t217 * t354;
t234 = -qJD(6) * t261 + t258 * t354 - t259 * t350;
t248 = mrSges(7,1) * t305 - mrSges(7,3) * t261;
t210 = m(7) * t215 - mrSges(7,2) * t270 + mrSges(7,3) * t234 + t245 * t260 - t248 * t305;
t201 = t354 * t209 + t350 * t210;
t262 = -mrSges(6,1) * t288 + mrSges(6,2) * t289;
t263 = -mrSges(6,2) * t306 + mrSges(6,3) * t288;
t199 = m(6) * t218 + mrSges(6,1) * t271 - mrSges(6,3) * t259 - t262 * t289 + t263 * t306 + t201;
t264 = mrSges(6,1) * t306 - mrSges(6,3) * t289;
t370 = -t209 * t350 + t354 * t210;
t200 = m(6) * t219 - mrSges(6,2) * t271 + mrSges(6,3) * t258 + t262 * t288 - t264 * t306 + t370;
t371 = -t199 * t345 + t348 * t200;
t192 = m(5) * t226 - mrSges(5,2) * t325 - mrSges(5,3) * t271 - t278 * t306 - t291 * t336 + t371;
t290 = -mrSges(5,2) * t336 - mrSges(5,3) * t306;
t222 = -t325 * pkin(4) - t335 * qJ(5) + t307 * t277 + qJDD(5) - t225;
t220 = -t258 * pkin(5) - t287 * pkin(10) + t289 * t265 + t222;
t367 = m(7) * t220 - t234 * mrSges(7,1) + mrSges(7,2) * t235 - t260 * t247 + t248 * t261;
t362 = -m(6) * t222 + t258 * mrSges(6,1) - mrSges(6,2) * t259 + t288 * t263 - t264 * t289 - t367;
t205 = m(5) * t225 + mrSges(5,1) * t325 - mrSges(5,3) * t272 - t278 * t307 + t290 * t336 + t362;
t183 = t346 * t192 + t385 * t205;
t308 = -mrSges(4,1) * t320 + mrSges(4,2) * t321;
t309 = -mrSges(4,2) * t336 + mrSges(4,3) * t320;
t181 = m(4) * t249 + mrSges(4,1) * t325 - mrSges(4,3) * t300 - t308 * t321 + t309 * t336 + t183;
t311 = mrSges(4,1) * t336 - mrSges(4,3) * t321;
t372 = t385 * t192 - t205 * t346;
t182 = m(4) * t250 - mrSges(4,2) * t325 + mrSges(4,3) * t299 + t308 * t320 - t311 * t336 + t372;
t176 = t355 * t181 + t351 * t182;
t194 = t348 * t199 + t345 * t200;
t302 = -g(3) * t384 + t380;
t326 = mrSges(3,1) * t341 - mrSges(3,3) * t375;
t330 = (-mrSges(3,1) * t356 + mrSges(3,2) * t352) * t379;
t373 = -t181 * t351 + t355 * t182;
t174 = m(3) * t302 - mrSges(3,2) * t340 - mrSges(3,3) * t333 - t326 * t341 + t330 * t374 + t373;
t327 = -mrSges(3,2) * t341 + mrSges(3,3) * t374;
t364 = m(5) * t246 + t271 * mrSges(5,1) + mrSges(5,2) * t272 + t306 * t290 + t291 * t307 + t194;
t361 = -m(4) * t282 + t299 * mrSges(4,1) - mrSges(4,2) * t300 + t320 * t309 - t311 * t321 - t364;
t189 = m(3) * t301 + mrSges(3,1) * t340 - mrSges(3,3) * t332 + t327 * t341 - t330 * t375 + t361;
t170 = t356 * t174 - t189 * t352;
t315 = -t347 * t328 - t386;
t175 = m(3) * t315 + t333 * mrSges(3,1) + t332 * mrSges(3,2) + (t326 * t352 - t327 * t356) * t379 + t176;
t167 = t174 * t382 - t175 * t347 + t189 * t381;
t238 = Ifges(7,5) * t261 + Ifges(7,6) * t260 + Ifges(7,3) * t305;
t240 = Ifges(7,1) * t261 + Ifges(7,4) * t260 + Ifges(7,5) * t305;
t202 = -mrSges(7,1) * t220 + mrSges(7,3) * t215 + Ifges(7,4) * t235 + Ifges(7,2) * t234 + Ifges(7,6) * t270 - t238 * t261 + t240 * t305;
t239 = Ifges(7,4) * t261 + Ifges(7,2) * t260 + Ifges(7,6) * t305;
t203 = mrSges(7,2) * t220 - mrSges(7,3) * t214 + Ifges(7,1) * t235 + Ifges(7,4) * t234 + Ifges(7,5) * t270 + t238 * t260 - t239 * t305;
t251 = Ifges(6,5) * t289 + Ifges(6,6) * t288 + Ifges(6,3) * t306;
t253 = Ifges(6,1) * t289 + Ifges(6,4) * t288 + Ifges(6,5) * t306;
t185 = -mrSges(6,1) * t222 + mrSges(6,3) * t219 + Ifges(6,4) * t259 + Ifges(6,2) * t258 + Ifges(6,6) * t271 - pkin(5) * t367 + pkin(10) * t370 + t354 * t202 + t350 * t203 - t289 * t251 + t306 * t253;
t252 = Ifges(6,4) * t289 + Ifges(6,2) * t288 + Ifges(6,6) * t306;
t187 = mrSges(6,2) * t222 - mrSges(6,3) * t218 + Ifges(6,1) * t259 + Ifges(6,4) * t258 + Ifges(6,5) * t271 - pkin(10) * t201 - t202 * t350 + t203 * t354 + t251 * t288 - t252 * t306;
t273 = Ifges(5,5) * t307 - Ifges(5,6) * t306 + Ifges(5,3) * t336;
t274 = Ifges(5,4) * t307 - Ifges(5,2) * t306 + Ifges(5,6) * t336;
t171 = mrSges(5,2) * t246 - mrSges(5,3) * t225 + Ifges(5,1) * t272 - Ifges(5,4) * t271 + Ifges(5,5) * t325 - qJ(5) * t194 - t185 * t345 + t187 * t348 - t273 * t306 - t274 * t336;
t275 = Ifges(5,1) * t307 - Ifges(5,4) * t306 + Ifges(5,5) * t336;
t365 = -mrSges(7,1) * t214 + mrSges(7,2) * t215 - Ifges(7,5) * t235 - Ifges(7,6) * t234 - Ifges(7,3) * t270 - t261 * t239 + t260 * t240;
t360 = -mrSges(6,1) * t218 + mrSges(6,2) * t219 - Ifges(6,5) * t259 - Ifges(6,6) * t258 - pkin(5) * t201 - t289 * t252 + t288 * t253 + t365;
t177 = -pkin(4) * t194 + mrSges(5,3) * t226 + t360 + (-Ifges(5,2) - Ifges(6,3)) * t271 - mrSges(5,1) * t246 + Ifges(5,4) * t272 - t307 * t273 + Ifges(5,6) * t325 + t336 * t275;
t293 = Ifges(4,5) * t321 + Ifges(4,6) * t320 + Ifges(4,3) * t336;
t295 = Ifges(4,1) * t321 + Ifges(4,4) * t320 + Ifges(4,5) * t336;
t160 = -mrSges(4,1) * t282 + mrSges(4,3) * t250 + Ifges(4,4) * t300 + Ifges(4,2) * t299 + Ifges(4,6) * t325 - pkin(3) * t364 + qJ(4) * t372 + t346 * t171 + t385 * t177 - t321 * t293 + t336 * t295;
t294 = Ifges(4,4) * t321 + Ifges(4,2) * t320 + Ifges(4,6) * t336;
t163 = mrSges(4,2) * t282 - mrSges(4,3) * t249 + Ifges(4,1) * t300 + Ifges(4,4) * t299 + Ifges(4,5) * t325 - qJ(4) * t183 + t385 * t171 - t346 * t177 + t320 * t293 - t336 * t294;
t313 = Ifges(3,6) * t341 + (Ifges(3,4) * t352 + Ifges(3,2) * t356) * t379;
t314 = Ifges(3,5) * t341 + (Ifges(3,1) * t352 + Ifges(3,4) * t356) * t379;
t157 = Ifges(3,5) * t332 - Ifges(3,6) * t333 + Ifges(3,3) * t340 + mrSges(3,1) * t301 - mrSges(3,2) * t302 + t351 * t163 + t355 * t160 + pkin(2) * t361 + pkin(9) * t373 + (t313 * t352 - t314 * t356) * t379;
t312 = Ifges(3,3) * t341 + (Ifges(3,5) * t352 + Ifges(3,6) * t356) * t379;
t159 = mrSges(3,2) * t315 - mrSges(3,3) * t301 + Ifges(3,1) * t332 - Ifges(3,4) * t333 + Ifges(3,5) * t340 - pkin(9) * t176 - t160 * t351 + t163 * t355 + t312 * t374 - t313 * t341;
t363 = -mrSges(5,1) * t225 + mrSges(5,2) * t226 - Ifges(5,5) * t272 + Ifges(5,6) * t271 - Ifges(5,3) * t325 - pkin(4) * t362 - qJ(5) * t371 - t348 * t185 - t345 * t187 - t307 * t274 - t306 * t275;
t359 = mrSges(4,1) * t249 - mrSges(4,2) * t250 + Ifges(4,5) * t300 + Ifges(4,6) * t299 + Ifges(4,3) * t325 + pkin(3) * t183 + t321 * t294 - t320 * t295 - t363;
t162 = -mrSges(3,1) * t315 + mrSges(3,3) * t302 + Ifges(3,4) * t332 - Ifges(3,2) * t333 + Ifges(3,6) * t340 - pkin(2) * t176 - t312 * t375 + t341 * t314 - t359;
t366 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t167 + t349 * t157 + t159 * t384 + t162 * t383 + t170 * t387;
t168 = m(2) * t338 - mrSges(2,1) * t358 - qJDD(1) * mrSges(2,2) + t170;
t166 = t349 * t175 + (t174 * t352 + t189 * t356) * t347;
t164 = m(2) * t337 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t358 + t167;
t155 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - t358 * Ifges(2,6) + t356 * t159 - t352 * t162 + (-t166 * t347 - t167 * t349) * pkin(8);
t154 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t347 * t157 + (pkin(8) * t170 + t159 * t352 + t162 * t356) * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t357 * t155 - t353 * t154 - pkin(7) * (t164 * t357 + t168 * t353), t155, t159, t163, t171, t187, t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t353 * t155 + t357 * t154 + pkin(7) * (-t164 * t353 + t168 * t357), t154, t162, t160, t177, t185, t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t157, t359, -t363, Ifges(6,3) * t271 - t360, -t365;];
m_new  = t1;
