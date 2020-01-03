% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:58
% EndTime: 2019-12-31 20:17:08
% DurationCPUTime: 5.74s
% Computational Cost: add. (14141->384), mult. (27700->520), div. (0->0), fcn. (32103->8), ass. (0->223)
t216 = sin(qJ(5));
t218 = cos(qJ(5));
t263 = mrSges(6,1) * t218 - mrSges(6,2) * t216;
t346 = sin(qJ(2));
t284 = t346 * pkin(6);
t246 = -qJ(3) * t346 - t284;
t348 = cos(qJ(2));
t285 = t348 * pkin(6);
t247 = qJ(3) * t348 + t285;
t301 = sin(pkin(9));
t302 = cos(pkin(9));
t178 = t246 * t301 + t247 * t302;
t236 = t301 * t346 - t302 * t348;
t146 = -pkin(7) * t236 + t178;
t217 = sin(qJ(4));
t347 = cos(qJ(4));
t197 = -t301 * t348 - t302 * t346;
t378 = t246 * t302 - t247 * t301;
t397 = t197 * pkin(7) + t378;
t406 = t146 * t347 + t217 * t397;
t412 = t406 * t263;
t413 = t406 * mrSges(5,1);
t92 = t146 * t217 - t347 * t397;
t419 = t92 * mrSges(5,2);
t421 = -t412 - t413 + t419;
t231 = -t197 * t347 - t217 * t236;
t393 = -t231 / 0.2e1;
t420 = 0.2e1 * t393;
t208 = pkin(2) * t302 + pkin(3);
t273 = pkin(2) * t301;
t187 = t208 * t217 + t273 * t347;
t418 = t187 * t92;
t417 = t216 * t92;
t416 = t218 * t92;
t340 = t406 * t92;
t414 = -t412 / 0.2e1 - t413 / 0.2e1;
t172 = t197 * t217 - t236 * t347;
t315 = t218 * mrSges(6,2);
t320 = t216 * mrSges(6,1);
t203 = t315 + t320;
t398 = t203 * t172;
t411 = m(6) * t406 + t398;
t410 = -t197 * mrSges(4,1) - t236 * mrSges(4,2);
t115 = t203 * t231;
t409 = t406 * t115 + t92 * t398;
t213 = t346 * pkin(2);
t408 = m(4) * t213;
t209 = -pkin(2) * t348 - pkin(1);
t182 = pkin(3) * t236 + t209;
t349 = t218 / 0.2e1;
t352 = -t216 / 0.2e1;
t210 = Ifges(6,5) * t218;
t330 = Ifges(6,6) * t216;
t384 = t210 - t330;
t392 = t231 / 0.2e1;
t334 = Ifges(6,4) * t216;
t262 = Ifges(6,1) * t218 - t334;
t395 = Ifges(6,5) * t231 + t172 * t262;
t211 = Ifges(6,4) * t218;
t261 = -Ifges(6,2) * t216 + t211;
t396 = Ifges(6,6) * t231 + t172 * t261;
t400 = -t172 / 0.2e1;
t407 = Ifges(5,4) * t420 + t395 * t349 + t396 * t352 + t182 * mrSges(5,1) + t384 * t392 + (Ifges(5,2) + Ifges(6,3)) * t400;
t363 = -t172 / 0.4e1;
t399 = t172 / 0.2e1;
t360 = t172 / 0.4e1;
t350 = -t218 / 0.2e1;
t344 = pkin(4) * t231;
t124 = -pkin(8) * t172 + t344;
t390 = mrSges(6,1) * t231;
t389 = mrSges(6,2) * t231;
t204 = Ifges(6,5) * t216 + Ifges(6,6) * t218;
t296 = t231 * t204;
t205 = Ifges(6,2) * t218 + t334;
t292 = t216 * t205;
t193 = -t292 / 0.2e1;
t206 = Ifges(6,1) * t216 + t211;
t289 = t218 * t206;
t385 = t289 / 0.2e1 + t193;
t214 = t216 ^ 2;
t215 = t218 ^ 2;
t288 = t214 + t215;
t70 = -Ifges(6,5) * t172 + t231 * t262;
t310 = t218 * t70;
t67 = -Ifges(6,6) * t172 + t231 * t261;
t317 = t216 * t67;
t383 = t310 / 0.4e1 - t317 / 0.4e1;
t46 = t124 * t218 + t417;
t47 = t124 * t216 - t416;
t382 = -t216 * t46 + t218 * t47;
t381 = t310 / 0.2e1 - t317 / 0.2e1 + Ifges(5,1) * t392 + Ifges(5,4) * t172;
t380 = t261 * t350 + t262 * t352;
t379 = (-t210 / 0.2e1 + t330 / 0.2e1) * t172;
t377 = 0.2e1 * m(6);
t376 = m(5) / 0.2e1;
t375 = m(6) / 0.2e1;
t374 = pkin(4) / 0.2e1;
t373 = -pkin(8) / 0.2e1;
t372 = pkin(8) / 0.2e1;
t371 = m(4) * pkin(2);
t370 = -t46 / 0.2e1;
t369 = t47 / 0.2e1;
t368 = t92 / 0.2e1;
t186 = t208 * t347 - t217 * t273;
t184 = -pkin(4) - t186;
t357 = t184 / 0.2e1;
t185 = pkin(8) + t187;
t356 = -t185 / 0.2e1;
t355 = t204 / 0.2e1;
t354 = -t206 / 0.4e1;
t353 = t210 / 0.4e1;
t351 = t216 / 0.2e1;
t345 = pkin(4) * t398;
t343 = pkin(4) * t203;
t333 = Ifges(5,5) * t172;
t331 = Ifges(5,6) * t231;
t293 = t216 * t172;
t117 = -mrSges(6,3) * t293 - t389;
t319 = t216 * mrSges(6,3);
t283 = t231 * t319;
t118 = mrSges(6,2) * t172 - t283;
t290 = t218 * t172;
t120 = -mrSges(6,3) * t290 + t390;
t314 = t218 * mrSges(6,3);
t121 = -mrSges(6,1) * t172 - t231 * t314;
t168 = t172 * mrSges(5,2);
t183 = -pkin(3) * t197 + t213;
t234 = pkin(2) * t236;
t235 = t236 ^ 2;
t282 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t90 = -pkin(4) * t172 - pkin(8) * t231 + t182;
t38 = -t216 * t406 + t218 * t90;
t39 = t216 * t90 + t218 * t406;
t95 = t124 + t183;
t40 = t218 * t95 + t417;
t41 = t216 * t95 - t416;
t1 = t235 * Ifges(4,4) + m(6) * (t38 * t40 + t39 * t41 + t340) + t39 * t117 + t41 * t118 + t38 * t120 + t40 * t121 + (-mrSges(4,2) * t213 - Ifges(4,4) * t197 + (-Ifges(4,2) + Ifges(4,1)) * t236) * t197 + (t183 * mrSges(5,2) + Ifges(5,1) * t399 + t407) * t231 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t348) * t348 + (-pkin(1) * mrSges(3,1) + mrSges(4,1) * t234 + (Ifges(3,1) - Ifges(3,2)) * t348 - Ifges(3,4) * t346) * t346 + (-t183 * mrSges(5,1) - t282 * t231 + t400 * t384 + t381) * t172 + (t408 + t410) * t209 + t409 + (m(5) * t183 + t168) * t182;
t329 = t1 * qJD(1);
t328 = t231 * mrSges(5,1);
t327 = t231 * mrSges(5,3);
t326 = t172 * mrSges(5,3);
t325 = t172 * mrSges(6,3);
t324 = t185 * t46;
t323 = t185 * t47;
t322 = t186 * t38;
t321 = t186 * t39;
t318 = t216 * t40;
t316 = t216 * t395;
t313 = t218 * t39;
t312 = t218 * t41;
t311 = t218 * t396;
t116 = -t172 * t319 - t389;
t119 = -t172 * t314 + t390;
t4 = m(6) * (t38 * t46 + t39 * t47 + t340) + t47 * t118 + t39 * t116 + t46 * t121 + t38 * t119 + t407 * t231 - (-t182 * mrSges(5,2) - t379 + (-Ifges(5,1) / 0.2e1 + t282) * t231 - t381) * t172 + t409;
t309 = t4 * qJD(1);
t112 = t263 * t231;
t5 = -t92 * t112 + t39 * t121 + (t70 * t351 + t67 * t349 - t172 * t355 + mrSges(6,3) * t313 + (-t206 * t350 + t193) * t231) * t231 + (-t118 - t283) * t38;
t306 = t5 * qJD(1);
t304 = t92 * t231;
t10 = t121 * t293 - t172 * t326 - t118 * t290 - m(6) * (t304 + (-t216 * t38 + t313) * t172) - m(5) * (t172 * t406 + t304) - m(4) * (-t178 * t236 + t197 * t378) + (-t197 ^ 2 - t235) * mrSges(4,3) + (-t327 - t115) * t231;
t300 = qJD(1) * t10;
t268 = t288 * t172;
t272 = t263 * t392;
t219 = (t172 * t187 - t186 * t231) * t376 + (t184 * t231 + t185 * t268) * t375 - t272 + (t197 * t302 - t236 * t301) * t371 / 0.2e1 + t288 * t325 / 0.2e1;
t225 = t183 * t376 + (t216 * t41 + t218 * t40) * t375 + t117 * t351 + t120 * t349 + t408 / 0.2e1;
t11 = -t168 + t219 - t225 - t328 - t410;
t299 = t11 * qJD(1);
t270 = t215 / 0.2e1 + t214 / 0.2e1;
t226 = t270 * t325 - t168 / 0.2e1 + (pkin(8) * t268 - t344) * t375 - t272;
t230 = -m(6) * (t216 * t47 + t218 * t46) / 0.2e1 + mrSges(5,2) * t400 + t116 * t352 + t119 * t350;
t14 = mrSges(5,1) * t420 + t226 + t230;
t298 = t14 * qJD(1);
t253 = -t315 / 0.2e1 - t320 / 0.2e1;
t245 = t253 * t172;
t250 = t118 * t350 + t121 * t351;
t17 = t245 + t250;
t297 = t17 * qJD(1);
t295 = t184 * t203;
t294 = t216 * t120;
t291 = t218 * t117;
t287 = mrSges(6,3) * t318;
t286 = mrSges(6,3) * t312;
t281 = Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1;
t279 = t210 / 0.2e1;
t274 = t203 * t368;
t271 = t360 + t363;
t267 = t288 * t186;
t266 = -t380 + t385;
t264 = mrSges(6,3) * t270;
t259 = t312 - t318;
t233 = -t186 * mrSges(5,2) + mrSges(6,3) * t267 + (-mrSges(5,1) - t263) * t187;
t27 = m(6) * (t184 * t187 + t185 * t267) + t233;
t229 = (t184 * t406 + t418) * t375 + t398 * t357 + t187 * t115 / 0.2e1 + t414;
t241 = t395 / 0.4e1 + t119 * t356 - t186 * t121 / 0.2e1;
t242 = t396 / 0.4e1 + t185 * t116 / 0.2e1 + t186 * t118 / 0.2e1;
t3 = t345 / 0.2e1 + (t392 + t393) * Ifges(5,6) + (t399 + t400) * Ifges(5,5) + (t368 - t92 / 0.2e1) * mrSges(5,2) + (m(6) * t374 + mrSges(5,1) / 0.2e1 + t263 / 0.2e1) * t406 + (t117 * t373 - t396 / 0.4e1 - t271 * t206 + (t369 - t41 / 0.2e1) * mrSges(6,3) + (-pkin(8) * t41 / 0.4e1 + t323 / 0.4e1 + t321 / 0.4e1) * t377 + t242) * t218 + (t120 * t372 - t395 / 0.4e1 + t271 * t205 + (t370 + t40 / 0.2e1) * mrSges(6,3) + (pkin(8) * t40 / 0.4e1 - t324 / 0.4e1 - t322 / 0.4e1) * t377 + t241) * t216 + t229;
t257 = t3 * qJD(1) + t27 * qJD(2);
t224 = (t354 - t211 / 0.4e1 + t281 * t216) * t216 + (-0.3e1 / 0.4e1 * t334 - t205 / 0.4e1 - t281 * t218) * t218;
t220 = (-t185 * t264 + t224) * t231 + t112 * t357 + t274;
t239 = Ifges(6,3) * t393 - t40 * mrSges(6,1) / 0.2e1 + t41 * mrSges(6,2) / 0.2e1;
t7 = -t172 * t353 + (t118 * t356 - t67 / 0.4e1 + (t360 + t399) * Ifges(6,6)) * t216 + (t121 * t356 + t70 / 0.4e1 + Ifges(6,5) * t400) * t218 + t220 + t239;
t97 = t266 + t295;
t256 = qJD(1) * t7 + qJD(2) * t97;
t254 = mrSges(6,1) * t370 + mrSges(6,2) * t369;
t252 = t384 * t363;
t251 = t118 * t352 + t121 * t350;
t249 = t292 / 0.2e1 - t289 / 0.2e1;
t244 = t253 * t186;
t228 = t249 + t380;
t126 = t228 + t343;
t52 = (t374 - t184 / 0.2e1) * t203 + t244 + t228;
t221 = -pkin(8) * t264 + t224;
t227 = t251 * pkin(8) - pkin(4) * t112 / 0.2e1 + t274 + t383;
t9 = -(-0.3e1 / 0.4e1 * t330 + t353 + t279) * t172 + (-Ifges(6,3) / 0.2e1 + t221) * t231 + t227 + t254;
t240 = t9 * qJD(1) - t52 * qJD(2) - t126 * qJD(4);
t53 = -t343 / 0.2e1 + t295 / 0.2e1 + t244 + t266;
t18 = t245 - t250;
t15 = mrSges(5,1) * t392 - t328 / 0.2e1 + t226 - t230;
t13 = t219 + t225;
t8 = Ifges(6,3) * t392 + t221 * t231 + t227 + t252 - t254 - t379;
t6 = t252 + t172 * t279 - Ifges(6,6) * t293 / 0.2e1 + t251 * t185 + t220 - t239 + t383;
t2 = t414 + ((-t322 - t324) * t375 + mrSges(6,3) * t370 + t205 * t363 + t241) * t216 + (-Ifges(5,6) / 0.2e1 + t204 / 0.4e1) * t231 + t316 / 0.4e1 + t292 * t363 - t345 / 0.2e1 + t286 / 0.2e1 + t296 / 0.4e1 + t333 / 0.2e1 - t331 / 0.2e1 + t419 + t311 / 0.4e1 + t229 - t287 / 0.2e1 + (-pkin(4) * t406 + t259 * pkin(8)) * t375 + Ifges(5,5) * t399 + t289 * t360 + t291 * t372 + t294 * t373 + (-t172 * t354 + mrSges(6,3) * t369 + (t321 + t323) * t375 + t242) * t218;
t12 = [qJD(2) * t1 - qJD(3) * t10 + qJD(4) * t4 - qJD(5) * t5, t329 + (-Ifges(4,5) * t236 + (t197 * t273 + t234 * t302) * mrSges(4,3) - Ifges(3,6) * t346 + Ifges(3,5) * t348 + t316 / 0.2e1 + t286 + t296 / 0.2e1 + t333 - t331 + m(5) * (-t186 * t406 - t418) + t311 / 0.2e1 + (m(6) * t259 + t291 - t294) * t185 - t186 * t326 - t187 * t327 + (-t178 * t302 + t301 * t378) * t371 - t178 * mrSges(4,1) + t411 * t184 + t421 - t287 + t385 * t172 + mrSges(3,2) * t284 + Ifges(4,6) * t197 - t378 * mrSges(4,2) - mrSges(3,1) * t285) * qJD(2) + t13 * qJD(3) + t2 * qJD(4) + t6 * qJD(5), qJD(2) * t13 + qJD(4) * t15 + qJD(5) * t18 - t300, t2 * qJD(2) + t15 * qJD(3) + t8 * qJD(5) + t309 + (t395 * t351 + t396 * t349 - (-Ifges(5,5) + t249) * t172 + (t355 - Ifges(5,6)) * t231 - t411 * pkin(4) + (m(6) * t382 + t218 * t116 - t216 * t119) * pkin(8) + t382 * mrSges(6,3) + t421) * qJD(4), -t306 + t6 * qJD(2) + t18 * qJD(3) + t8 * qJD(4) + (-t39 * mrSges(6,1) - t38 * mrSges(6,2) - t296) * qJD(5); qJD(3) * t11 + qJD(4) * t3 + qJD(5) * t7 - t329, qJD(4) * t27 + qJD(5) * t97, t299, (m(6) * (-pkin(4) * t187 + pkin(8) * t267) + t233) * qJD(4) + t53 * qJD(5) + t257, t53 * qJD(4) + (-t185 * t263 + t384) * qJD(5) + t256; -qJD(2) * t11 - qJD(4) * t14 - qJD(5) * t17 + t300, -t299, 0, -t298, -qJD(5) * t203 - t297; -qJD(2) * t3 + qJD(3) * t14 + qJD(5) * t9 - t309, -qJD(5) * t52 - t257, t298, -t126 * qJD(5), (-pkin(8) * t263 + t384) * qJD(5) + t240; -qJD(2) * t7 + qJD(3) * t17 - qJD(4) * t9 + t306, qJD(4) * t52 - t256, t297, -t240, 0;];
Cq = t12;
