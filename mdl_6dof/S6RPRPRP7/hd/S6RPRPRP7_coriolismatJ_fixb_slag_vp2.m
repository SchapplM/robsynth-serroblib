% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:16
% EndTime: 2019-03-09 03:21:26
% DurationCPUTime: 5.01s
% Computational Cost: add. (10550->463), mult. (19521->628), div. (0->0), fcn. (20265->6), ass. (0->238)
t266 = sin(qJ(5));
t263 = t266 ^ 2;
t268 = cos(qJ(5));
t264 = t268 ^ 2;
t337 = t263 + t264;
t435 = (mrSges(6,3) + mrSges(7,3)) * t337;
t265 = sin(pkin(9));
t247 = pkin(3) * t265 + pkin(8);
t344 = qJ(6) + t247;
t216 = t344 * t266;
t217 = t344 * t268;
t425 = t216 * t266 + t217 * t268;
t434 = m(7) * t425;
t267 = sin(qJ(3));
t363 = cos(pkin(9));
t389 = cos(qJ(3));
t230 = t265 * t267 - t363 * t389;
t231 = t265 * t389 + t267 * t363;
t433 = t231 * mrSges(5,1) - t230 * mrSges(5,2);
t249 = t267 * pkin(3) + qJ(2);
t432 = m(5) * t249;
t380 = Ifges(6,6) + Ifges(7,6);
t430 = t380 * t268;
t429 = mrSges(7,3) / 0.2e1;
t353 = t230 * t266;
t151 = -t231 * mrSges(7,2) + mrSges(7,3) * t353;
t428 = -t151 / 0.2e1;
t248 = -pkin(3) * t363 - pkin(4);
t381 = t268 * pkin(5);
t237 = t248 - t381;
t383 = m(7) * t237;
t427 = Ifges(6,5) + Ifges(7,5);
t426 = Ifges(7,3) + Ifges(6,3);
t408 = m(7) * pkin(5);
t326 = -mrSges(7,1) - t408;
t269 = -pkin(1) - pkin(7);
t346 = t267 * t269;
t235 = -qJ(4) * t267 + t346;
t323 = t389 * t269;
t236 = -qJ(4) * t389 + t323;
t424 = t363 * t235 + t265 * t236;
t256 = t266 * mrSges(7,1);
t257 = t268 * mrSges(7,2);
t338 = t257 + t256;
t260 = Ifges(7,4) * t268;
t242 = t266 * Ifges(7,1) + t260;
t261 = Ifges(6,4) * t268;
t243 = t266 * Ifges(6,1) + t261;
t422 = t230 * t265 + t363 * t231;
t333 = t389 * pkin(3);
t146 = -t230 * pkin(4) + t231 * pkin(8) + t333;
t162 = t235 * t265 - t363 * t236;
t61 = t268 * t146 + t162 * t266;
t62 = t266 * t146 - t162 * t268;
t298 = -t266 * t61 + t268 * t62;
t350 = t231 * t268;
t40 = -pkin(5) * t230 + qJ(6) * t350 + t61;
t351 = t231 * t266;
t50 = qJ(6) * t351 + t62;
t420 = -t266 * t40 + t268 * t50;
t369 = t268 * mrSges(6,1);
t375 = t266 * mrSges(6,2);
t239 = -t369 + t375;
t419 = m(6) * t248 + t239;
t379 = Ifges(6,4) * t266;
t304 = -Ifges(6,1) * t268 + t379;
t101 = t231 * Ifges(6,5) + t230 * t304;
t328 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t378 = Ifges(7,4) * t266;
t303 = -Ifges(7,1) * t268 + t378;
t99 = t231 * Ifges(7,5) + t230 * t303;
t418 = -t328 * t231 - t101 / 0.2e1 - t99 / 0.2e1;
t177 = t230 ^ 2;
t417 = t231 ^ 2;
t416 = m(5) / 0.2e1;
t414 = m(6) / 0.2e1;
t413 = m(6) / 0.4e1;
t412 = m(7) / 0.2e1;
t411 = m(7) / 0.4e1;
t410 = -pkin(5) / 0.2e1;
t409 = m(5) * pkin(3);
t407 = mrSges(6,2) / 0.2e1;
t406 = mrSges(7,2) / 0.2e1;
t405 = -t40 / 0.2e1;
t352 = t230 * t268;
t145 = pkin(4) * t231 + pkin(8) * t230 + t249;
t59 = t268 * t145 - t266 * t424;
t48 = qJ(6) * t352 + t59;
t404 = t48 / 0.2e1;
t209 = mrSges(7,2) * t353;
t135 = -mrSges(7,1) * t352 + t209;
t402 = t135 / 0.2e1;
t136 = t230 * t239;
t401 = t136 / 0.2e1;
t139 = t338 * t230;
t400 = -t139 / 0.2e1;
t153 = -mrSges(7,1) * t230 + mrSges(7,3) * t350;
t399 = -t153 / 0.2e1;
t398 = t162 / 0.2e1;
t397 = -t230 / 0.2e1;
t395 = -t231 / 0.2e1;
t368 = t268 * mrSges(7,1);
t374 = t266 * mrSges(7,2);
t238 = -t368 + t374;
t394 = -t238 / 0.2e1;
t393 = -t247 / 0.2e1;
t392 = -t266 / 0.2e1;
t391 = t266 / 0.2e1;
t390 = t268 / 0.2e1;
t334 = pkin(5) * t353;
t106 = t162 - t334;
t388 = m(7) * t106;
t107 = -pkin(5) * t351 + t424;
t387 = m(7) * t107;
t385 = m(7) * t230;
t384 = m(7) * t231;
t382 = pkin(5) * t266;
t306 = t266 * mrSges(6,1) + t268 * mrSges(6,2);
t140 = t230 * t306;
t60 = t145 * t266 + t268 * t424;
t299 = t266 * t59 - t268 * t60;
t39 = t231 * pkin(5) + t48;
t49 = qJ(6) * t353 + t60;
t300 = t266 * t39 - t268 * t49;
t155 = t231 * mrSges(7,1) + mrSges(7,3) * t352;
t156 = t231 * mrSges(6,1) + mrSges(6,3) * t352;
t341 = t155 + t156;
t152 = -t231 * mrSges(6,2) + mrSges(6,3) * t353;
t342 = t151 + t152;
t358 = t162 * t230;
t9 = (mrSges(5,3) * t230 + t139 + t140) * t230 + (mrSges(5,3) * t231 + t266 * t341 - t268 * t342) * t231 + m(5) * (-t231 * t424 - t358) + m(7) * (-t106 * t230 + t231 * t300) + m(6) * (t231 * t299 - t358);
t377 = qJD(1) * t9;
t100 = -t230 * Ifges(6,5) + t231 * t304;
t137 = t338 * t231;
t138 = t306 * t231;
t149 = mrSges(7,2) * t230 + mrSges(7,3) * t351;
t150 = mrSges(6,2) * t230 + mrSges(6,3) * t351;
t154 = -mrSges(6,1) * t230 + mrSges(6,3) * t350;
t213 = t231 * mrSges(5,2);
t327 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t301 = Ifges(7,2) * t266 - t260;
t95 = t231 * Ifges(7,6) + t230 * t301;
t302 = Ifges(6,2) * t266 - t261;
t97 = t231 * Ifges(6,6) + t230 * t302;
t286 = t95 / 0.2e1 + t97 / 0.2e1 + t327 * t231;
t94 = -t230 * Ifges(7,6) + t231 * t301;
t96 = -t230 * Ifges(6,6) + t231 * t302;
t98 = -t230 * Ifges(7,5) + t231 * t303;
t1 = -t249 * t213 + t61 * t156 - t162 * t138 - t424 * t140 + t49 * t149 + t60 * t150 + t50 * t151 + t62 * t152 + t39 * t153 + t59 * t154 + t40 * t155 - t106 * t137 - t107 * t139 + m(6) * (t162 * t424 + t59 * t61 + t60 * t62) + m(7) * (t106 * t107 + t39 * t40 + t49 * t50) + (-t249 * mrSges(5,1) - Ifges(5,4) * t230 + (-t98 / 0.2e1 - t100 / 0.2e1 + t328 * t230) * t268 + (t94 / 0.2e1 + t96 / 0.2e1 - t327 * t230) * t266) * t230 + (Ifges(5,4) * t231 + t418 * t268 + t286 * t266 + (Ifges(5,1) - Ifges(5,2) - t426) * t230) * t231 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t267) * t267 + (qJ(2) * mrSges(4,1) + (-Ifges(4,1) + Ifges(4,2)) * t267 - Ifges(4,4) * t389) * t389 + (t432 + t433) * t333;
t376 = t1 * qJD(1);
t240 = t268 * Ifges(7,2) + t378;
t141 = t230 * t240;
t241 = t268 * Ifges(6,2) + t379;
t142 = t230 * t241;
t143 = t230 * t242;
t144 = t230 * t243;
t325 = m(7) * (-t39 + t48);
t4 = t106 * t135 + t162 * t136 + t48 * t151 + t59 * t152 - t60 * t156 + (t325 - t155) * t49 + ((t142 / 0.2e1 + t141 / 0.2e1 - t39 * mrSges(7,3) - t59 * mrSges(6,3) - t418) * t266 + (-t143 / 0.2e1 - t144 / 0.2e1 + t49 * mrSges(7,3) + t60 * mrSges(6,3) + (t139 - t388) * pkin(5) + t286) * t268) * t230;
t365 = t4 * qJD(1);
t313 = t264 / 0.2e1 + t263 / 0.2e1;
t307 = t313 * mrSges(6,3);
t281 = mrSges(7,3) * t313 + t307;
t316 = t155 / 0.2e1 + t156 / 0.2e1;
t317 = t428 - t152 / 0.2e1;
t318 = t402 + t401;
t329 = t407 + t406;
t330 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t7 = (t231 * t317 + t329) * t266 + (t231 * t281 + t318) * t230 + (-t316 * t231 + (t177 * t410 + t231 * t404 + t39 * t395 + t410) * m(7) - t330) * t268;
t364 = t7 * qJD(1);
t310 = t337 * t247;
t272 = (-t230 * t248 - t231 * t310) * t414 + (-t230 * t237 - t231 * t425) * t412 + t230 * t394 + t239 * t397 + (t230 * t363 - t231 * t265) * t409 / 0.2e1 + t395 * t435;
t275 = (t266 * t62 + t268 * t61) * t414 + (t266 * t50 + t268 * t40) * t412 + t333 * t416 + (t149 + t150) * t391 + (t153 + t154) * t390;
t10 = t230 * mrSges(5,1) + t213 + t272 - t275;
t362 = qJD(1) * t10;
t284 = t267 * mrSges(4,1) + mrSges(4,2) * t389 + t433;
t290 = m(7) * (t266 * t49 + t39 * t268);
t15 = mrSges(3,3) + t341 * t268 + t342 * t266 + (m(4) + m(3)) * qJ(2) + t432 + t290 + m(6) * (t266 * t60 + t268 * t59) + t284;
t361 = qJD(1) * t15;
t22 = (t266 * t151 + t268 * t155 + t290) * t230;
t360 = qJD(1) * t22;
t349 = t247 * t268;
t312 = m(7) * t337;
t282 = -m(5) / 0.2e1 - m(6) * t337 / 0.2e1 - t312 / 0.2e1;
t335 = t411 + t413;
t283 = (-t177 - t417) * t416 + 0.2e1 * t335 * (-t337 * t417 - t177);
t26 = t282 + t283;
t348 = t26 * qJD(1);
t314 = -t263 / 0.4e1 - t264 / 0.4e1;
t68 = 0.2e1 * (0.1e1 / 0.4e1 - t314) * t385;
t345 = t68 * qJD(1);
t343 = t151 * t390 + t155 * t392;
t336 = t408 / 0.2e1;
t331 = -t385 / 0.2e1;
t220 = -m(7) * t382 - t338;
t322 = t351 / 0.2e1;
t321 = -t350 / 0.2e1;
t258 = Ifges(7,5) * t268;
t259 = Ifges(6,5) * t268;
t315 = t259 / 0.4e1 + t258 / 0.4e1;
t309 = t325 / 0.2e1;
t308 = t312 / 0.2e1;
t27 = 0.4e1 * t335 * (0.1e1 - t337) * t231 * t230;
t271 = (-t138 / 0.2e1 - t137 / 0.2e1 + t317 * t268 + t316 * t266) * t230 + (-t140 / 0.2e1 + t400 + (t149 / 0.2e1 + t150 / 0.2e1) * t268 + (t399 - t154 / 0.2e1) * t266) * t231 + 0.2e1 * ((t106 + t420) * t411 + (t162 + t298) * t413) * t231 + 0.2e1 * ((t107 + t300) * t411 + (t424 + t299) * t413) * t230;
t296 = -t216 * t268 + t217 * t266;
t289 = m(7) * t296;
t6 = -t289 / 0.2e1 + t271;
t297 = t6 * qJD(1) + t27 * qJD(2);
t74 = t326 * t352 + t209;
t294 = qJD(1) * t74 - qJD(3) * t220;
t293 = t336 + t330;
t292 = -t155 / 0.2e1 + t309;
t287 = -t156 / 0.2e1 + t292;
t12 = (t231 * t329 + t317) * t268 + (t231 * t293 - t287) * t266;
t288 = t12 * qJD(1);
t102 = t337 * mrSges(7,3) + t434;
t277 = (t230 * t296 - t300) * t412 + t343;
t20 = (t257 / 0.2e1 + t256 / 0.2e1) * t231 - t387 / 0.2e1 + t277;
t66 = 0.2e1 * (0.1e1 / 0.4e1 + t314) * t384;
t285 = qJD(1) * t20 - qJD(2) * t66 + qJD(3) * t102;
t280 = mrSges(7,1) * t405 + t50 * t406 - t61 * mrSges(6,1) / 0.2e1 + t62 * t407;
t23 = -t248 * t306 - t238 * t382 - t237 * t338 + (-t242 / 0.2e1 - t243 / 0.2e1 - t261 / 0.2e1 - t260 / 0.2e1) * t268 + (-pkin(5) * t383 + t241 / 0.2e1 + t240 / 0.2e1 + (Ifges(7,4) / 0.2e1 + Ifges(6,4) / 0.2e1) * t266 + (-Ifges(7,1) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t268) * t266;
t273 = t247 * t307 + (t216 * t429 + t243 / 0.4e1 + t242 / 0.4e1 + t261 / 0.4e1 + t260 / 0.4e1 + (-Ifges(6,2) / 0.4e1 - Ifges(7,2) / 0.4e1) * t266) * t266 + (t217 * t429 + t241 / 0.4e1 + t240 / 0.4e1 + (-Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1) * t268 + (Ifges(6,4) / 0.4e1 + Ifges(7,4) / 0.4e1) * t266 + (t394 - t383 / 0.2e1) * pkin(5)) * t268;
t274 = t292 * t217 + (t156 * t393 + t142 / 0.4e1 + t141 / 0.4e1 + mrSges(6,2) * t398 + t101 / 0.4e1 + t99 / 0.4e1 + (t404 - t39 / 0.2e1) * mrSges(7,3)) * t268 + t106 * t338 / 0.2e1 + t216 * t428 + t237 * t402 + t248 * t401;
t276 = t152 * t393 + mrSges(6,1) * t398 + t144 / 0.4e1 + t143 / 0.4e1 - t97 / 0.4e1 - t95 / 0.4e1 + (t400 + t388 / 0.2e1) * pkin(5);
t3 = (m(7) * t405 + t399) * pkin(5) + t276 * t266 + (t328 * t268 + (-0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t266 + t315) * t231 + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + t273) * t230 + t274 + t280;
t279 = -t3 * qJD(1) + t23 * qJD(3);
t278 = t266 * t293 + t268 * t329;
t69 = t231 * t308 + t384 / 0.2e1;
t67 = t230 * t308 + t331;
t29 = t230 * t278 + t334 * t412 + (t306 + t338) * t230 / 0.2e1;
t25 = -t282 + t283;
t19 = t387 / 0.2e1 + mrSges(7,2) * t321 - mrSges(7,1) * t351 / 0.2e1 + t277;
t13 = t152 * t390 + t156 * t392 + t231 * t278 + t266 * t309 + t343;
t11 = t272 + t275;
t8 = -t375 / 0.2e1 + t369 / 0.2e1 - t374 / 0.2e1 + t368 / 0.2e1 + t268 * t336 + (t331 * t381 + t318) * t230 + (t281 * t230 + t317 * t266 + t268 * t287) * t231;
t5 = t289 / 0.2e1 + t271;
t2 = t315 * t231 + t273 * t230 + t274 + pkin(5) * t153 / 0.2e1 + t40 * t336 + ((-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t231 + t276) * t266 - t280 + t426 * t397 + t380 * t322 + t427 * t321;
t14 = [qJD(2) * t15 + qJD(3) * t1 + qJD(4) * t9 + qJD(5) * t4 + qJD(6) * t22, qJD(3) * t5 + qJD(4) * t25 + qJD(5) * t8 + t361, t376 + t5 * qJD(2) + (m(7) * (t107 * t237 - t216 * t40 + t217 * t50) - Ifges(4,5) * t267 - t248 * t138 - Ifges(5,5) * t231 - t237 * t137 + t107 * t238 + Ifges(5,6) * t230 - t216 * t153 + t217 * t149 + t150 * t349 - Ifges(4,6) * t389 - mrSges(4,1) * t346 - mrSges(4,2) * t323 + t422 * mrSges(5,3) * pkin(3) + (t266 * t427 + t430) * t397 + (t100 + t98) * t391 + (t96 + t94) * t390 + (t241 + t240) * t322 + (t243 + t242) * t321 + (m(6) * t298 - t266 * t154) * t247 - (t265 * t409 - mrSges(5,2)) * t162 + (-t363 * t409 - mrSges(5,1) + t419) * t424 + t420 * mrSges(7,3) + t298 * mrSges(6,3)) * qJD(3) + t11 * qJD(4) + t2 * qJD(5) + t19 * qJD(6), qJD(2) * t25 + qJD(3) * t11 + qJD(5) * t13 + qJD(6) * t67 + t377, t8 * qJD(2) + t2 * qJD(3) + t13 * qJD(4) + t365 + (-mrSges(6,1) * t60 - mrSges(6,2) * t59 - mrSges(7,2) * t48 + (t430 + (-mrSges(7,3) * pkin(5) + t427) * t266) * t230 + t326 * t49) * qJD(5), qJD(3) * t19 + qJD(4) * t67 + t360; qJD(3) * t6 + qJD(4) * t26 + qJD(5) * t7 - t361, t27 * qJD(3), t29 * qJD(5) + t69 * qJD(6) + t297 + (-t422 * t409 - t284 + (t238 + t419 + t383) * t231 + (-m(6) * t310 - t434 - t435) * t230) * qJD(3), t348, t364 + t29 * qJD(3) + ((mrSges(6,2) + mrSges(7,2)) * t266 + (-mrSges(6,1) + t326) * t268) * qJD(5) * t231, t69 * qJD(3); -qJD(2) * t6 + qJD(4) * t10 + qJD(5) * t3 + qJD(6) * t20 - t376, -qJD(6) * t66 - t297, -qJD(5) * t23 + qJD(6) * t102, t362, -t279 + (-mrSges(6,1) * t349 + mrSges(7,2) * t216 - mrSges(7,3) * t381 + t258 + t259 + (mrSges(6,2) * t247 - t380) * t266 + t326 * t217) * qJD(5), t285; -qJD(2) * t26 - qJD(3) * t10 - qJD(5) * t12 + qJD(6) * t68 - t377, -t348, -t362, 0 (-t306 + t220) * qJD(5) - t288, t345; -qJD(2) * t7 - qJD(3) * t3 + qJD(4) * t12 - qJD(6) * t74 - t365, -t364, t220 * qJD(6) + t279, t288, 0, -t294; -qJD(3) * t20 - qJD(4) * t68 + qJD(5) * t74 - t360, t66 * qJD(3), -qJD(5) * t220 - t285, -t345, t294, 0;];
Cq  = t14;
