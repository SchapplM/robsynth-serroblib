% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:01
% EndTime: 2019-03-09 06:57:18
% DurationCPUTime: 7.94s
% Computational Cost: add. (12857->579), mult. (29512->791), div. (0->0), fcn. (20504->10), ass. (0->280)
t264 = sin(qJ(5));
t365 = -t264 / 0.2e1;
t268 = cos(qJ(5));
t265 = sin(qJ(4));
t254 = pkin(3) * t265 + pkin(9);
t360 = -pkin(10) - t254;
t302 = qJD(5) * t360;
t269 = cos(qJ(4));
t353 = pkin(3) * qJD(4);
t312 = t269 * t353;
t270 = cos(qJ(3));
t319 = qJD(1) * t270;
t266 = sin(qJ(3));
t320 = qJD(1) * t266;
t228 = -t265 * t320 + t269 * t319;
t327 = t228 * t264;
t314 = pkin(10) * t327;
t251 = sin(pkin(11)) * pkin(1) + pkin(7);
t244 = t251 * qJD(1);
t301 = pkin(8) * qJD(1) + t244;
t315 = t266 * qJD(2);
t210 = t270 * t301 + t315;
t200 = t265 * t210;
t258 = t270 * qJD(2);
t285 = t301 * t266;
t209 = t258 - t285;
t146 = t209 * t269 - t200;
t238 = t265 * t270 + t269 * t266;
t229 = t238 * qJD(1);
t190 = pkin(4) * t229 - pkin(9) * t228;
t166 = pkin(3) * t320 + t190;
t83 = t268 * t146 + t264 * t166;
t425 = t264 * t302 + t268 * t312 + t314 - t83;
t326 = t228 * t268;
t297 = t229 * pkin(5) - pkin(10) * t326;
t82 = -t146 * t264 + t268 * t166;
t424 = -t264 * t312 + t268 * t302 - t297 - t82;
t260 = qJD(3) + qJD(4);
t423 = t260 * Ifges(5,6) / 0.2e1;
t383 = -pkin(10) - pkin(9);
t309 = qJD(5) * t383;
t202 = qJD(3) * pkin(3) + t209;
t140 = t202 * t269 - t200;
t85 = t268 * t140 + t264 * t190;
t422 = t264 * t309 + t314 - t85;
t84 = -t140 * t264 + t268 * t190;
t421 = t268 * t309 - t297 - t84;
t203 = -t229 * t264 + t260 * t268;
t225 = qJD(5) - t228;
t204 = t229 * t268 + t260 * t264;
t357 = Ifges(6,4) * t204;
t116 = Ifges(6,2) * t203 + Ifges(6,6) * t225 + t357;
t127 = -pkin(4) * t260 - t140;
t295 = mrSges(6,1) * t264 + mrSges(6,2) * t268;
t420 = t116 * t365 + t127 * t295;
t310 = -cos(pkin(11)) * pkin(1) - pkin(2);
t242 = -pkin(3) * t270 + t310;
t230 = qJD(1) * t242;
t337 = t260 * Ifges(5,5);
t419 = t230 * mrSges(5,2) + t337 / 0.2e1;
t340 = t229 * Ifges(5,4);
t418 = t423 + t340 / 0.2e1 + t228 * Ifges(5,2) / 0.2e1;
t263 = sin(qJ(6));
t267 = cos(qJ(6));
t136 = t203 * t263 + t204 * t267;
t196 = t260 * t238;
t185 = t196 * qJD(1);
t178 = Ifges(7,3) * t185;
t201 = t269 * t210;
t141 = t202 * t265 + t201;
t128 = pkin(9) * t260 + t141;
t152 = -t228 * pkin(4) - t229 * pkin(9) + t230;
t78 = t128 * t268 + t152 * t264;
t58 = pkin(10) * t203 + t78;
t335 = t263 * t58;
t77 = -t128 * t264 + t268 * t152;
t57 = -pkin(10) * t204 + t77;
t52 = pkin(5) * t225 + t57;
t19 = t267 * t52 - t335;
t333 = t267 * t58;
t20 = t263 * t52 + t333;
t299 = t267 * t203 - t204 * t263;
t354 = Ifges(7,4) * t136;
t219 = qJD(6) + t225;
t370 = -t219 / 0.2e1;
t376 = -t136 / 0.2e1;
t92 = -pkin(5) * t203 + t127;
t417 = t178 + (Ifges(7,5) * t299 - Ifges(7,6) * t136) * t370 + (t136 * t20 + t19 * t299) * mrSges(7,3) - t92 * (mrSges(7,1) * t136 + mrSges(7,2) * t299) + (Ifges(7,1) * t299 - t354) * t376;
t233 = t360 * t264;
t259 = t268 * pkin(10);
t234 = t254 * t268 + t259;
t186 = t233 * t267 - t234 * t263;
t416 = qJD(6) * t186 + t424 * t263 + t425 * t267;
t187 = t233 * t263 + t234 * t267;
t415 = -qJD(6) * t187 - t425 * t263 + t424 * t267;
t236 = t265 * t266 - t269 * t270;
t195 = t260 * t236;
t184 = t195 * qJD(1);
t124 = qJD(5) * t203 - t184 * t268;
t125 = -qJD(5) * t204 + t184 * t264;
t62 = -mrSges(6,1) * t125 + mrSges(6,2) * t124;
t252 = qJD(3) * t258;
t198 = -qJD(3) * t285 + t252;
t274 = qJD(3) * t210;
t73 = qJD(4) * t141 + t198 * t265 + t269 * t274;
t413 = m(6) * t73 + t62;
t237 = t263 * t268 + t264 * t267;
t158 = t237 * t228;
t394 = qJD(5) + qJD(6);
t194 = t394 * t237;
t412 = t158 - t194;
t283 = t263 * t264 - t267 * t268;
t159 = t283 * t228;
t193 = t394 * t283;
t411 = t159 - t193;
t145 = t209 * t265 + t201;
t216 = pkin(5) * t327;
t317 = qJD(5) * t264;
t311 = pkin(5) * t317;
t410 = t265 * t353 - t145 - t216 + t311;
t318 = qJD(3) * t266;
t313 = pkin(3) * t318;
t100 = pkin(4) * t185 + pkin(9) * t184 + qJD(1) * t313;
t316 = qJD(5) * t268;
t72 = qJD(4) * t140 + t269 * t198 - t265 * t274;
t16 = t264 * t100 - t128 * t317 + t152 * t316 + t268 * t72;
t17 = -qJD(5) * t78 + t268 * t100 - t264 * t72;
t409 = t16 * t268 - t17 * t264;
t11 = pkin(5) * t185 - pkin(10) * t124 + t17;
t13 = pkin(10) * t125 + t16;
t3 = qJD(6) * t19 + t11 * t263 + t13 * t267;
t39 = qJD(6) * t299 + t124 * t267 + t125 * t263;
t4 = -qJD(6) * t20 + t11 * t267 - t13 * t263;
t40 = -qJD(6) * t136 - t124 * t263 + t125 * t267;
t408 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,5) * t39 + Ifges(7,6) * t40;
t132 = Ifges(7,4) * t299;
t407 = -Ifges(7,2) * t136 + t132;
t153 = -mrSges(6,2) * t225 + mrSges(6,3) * t203;
t154 = mrSges(6,1) * t225 - mrSges(6,3) * t204;
t80 = mrSges(6,1) * t185 - mrSges(6,3) * t124;
t81 = -mrSges(6,2) * t185 + mrSges(6,3) * t125;
t286 = -t264 * t80 + t268 * t81;
t406 = (-t153 * t264 - t154 * t268) * qJD(5) + t286;
t290 = Ifges(6,5) * t268 - Ifges(6,6) * t264;
t355 = Ifges(6,4) * t268;
t292 = -Ifges(6,2) * t264 + t355;
t356 = Ifges(6,4) * t264;
t294 = Ifges(6,1) * t268 - t356;
t199 = Ifges(6,4) * t203;
t117 = Ifges(6,1) * t204 + Ifges(6,5) * t225 + t199;
t322 = t268 * t117;
t371 = t204 / 0.2e1;
t405 = t322 / 0.2e1 + t225 * t290 / 0.2e1 + t294 * t371 + t203 * t292 / 0.2e1 + t420;
t404 = -t230 * mrSges(5,1) - t77 * mrSges(6,1) - t19 * mrSges(7,1) + t78 * mrSges(6,2) + t20 * mrSges(7,2) + t418;
t389 = t39 / 0.2e1;
t388 = t40 / 0.2e1;
t374 = t185 / 0.2e1;
t248 = t383 * t264;
t249 = pkin(9) * t268 + t259;
t205 = t248 * t267 - t249 * t263;
t398 = qJD(6) * t205 + t421 * t263 + t422 * t267;
t206 = t248 * t263 + t249 * t267;
t397 = -qJD(6) * t206 - t422 * t263 + t421 * t267;
t175 = t283 * t238;
t361 = pkin(8) + t251;
t231 = t361 * t266;
t232 = t361 * t270;
t183 = -t231 * t265 + t232 * t269;
t170 = t268 * t183;
t171 = t236 * pkin(4) - t238 * pkin(9) + t242;
t104 = t264 * t171 + t170;
t396 = -t269 * t231 - t232 * t265;
t395 = -t264 * t77 + t268 * t78;
t392 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + Ifges(6,5) * t124 + Ifges(6,6) * t125 + t408;
t391 = Ifges(7,4) * t389 + Ifges(7,2) * t388 + Ifges(7,6) * t374;
t390 = Ifges(7,1) * t389 + Ifges(7,4) * t388 + Ifges(7,5) * t374;
t60 = Ifges(7,2) * t299 + Ifges(7,6) * t219 + t354;
t387 = -t60 / 0.2e1;
t386 = t60 / 0.2e1;
t61 = Ifges(7,1) * t136 + Ifges(7,5) * t219 + t132;
t385 = -t61 / 0.2e1;
t384 = t61 / 0.2e1;
t380 = t124 / 0.2e1;
t379 = t125 / 0.2e1;
t378 = -t299 / 0.2e1;
t377 = t299 / 0.2e1;
t375 = t136 / 0.2e1;
t373 = -t203 / 0.2e1;
t372 = -t204 / 0.2e1;
t369 = t219 / 0.2e1;
t368 = -t225 / 0.2e1;
t364 = t268 / 0.2e1;
t363 = m(5) * t230;
t362 = pkin(3) * t269;
t359 = mrSges(5,3) * t228;
t358 = Ifges(4,4) * t266;
t223 = Ifges(5,4) * t228;
t352 = t299 * Ifges(7,6);
t351 = t136 * Ifges(7,5);
t348 = t396 * t73;
t347 = t203 * Ifges(6,6);
t346 = t204 * Ifges(6,5);
t345 = t219 * Ifges(7,3);
t344 = t225 * Ifges(6,3);
t342 = t229 * mrSges(5,3);
t341 = t229 * Ifges(5,1);
t339 = t236 * t73;
t246 = t310 * qJD(1);
t338 = t246 * mrSges(4,2);
t331 = Ifges(4,5) * qJD(3);
t330 = Ifges(4,6) * qJD(3);
t324 = t238 * t264;
t321 = mrSges(5,1) * t260 + mrSges(6,1) * t203 - mrSges(6,2) * t204 - t342;
t257 = Ifges(4,4) * t319;
t256 = -pkin(5) * t268 - pkin(4);
t308 = t331 / 0.2e1;
t307 = -t330 / 0.2e1;
t304 = m(4) * t251 + mrSges(4,3);
t303 = qJD(3) * t361;
t221 = t266 * t303;
t222 = t270 * t303;
t107 = qJD(4) * t396 - t221 * t269 - t222 * t265;
t123 = pkin(4) * t196 + pkin(9) * t195 + t313;
t300 = -t107 * t264 + t268 * t123;
t103 = t268 * t171 - t183 * t264;
t298 = t330 / 0.2e1 + (Ifges(4,2) * t270 + t358) * qJD(1) / 0.2e1 - t246 * mrSges(4,1);
t296 = mrSges(6,1) * t268 - mrSges(6,2) * t264;
t293 = Ifges(6,1) * t264 + t355;
t291 = Ifges(6,2) * t268 + t356;
t289 = Ifges(6,5) * t264 + Ifges(6,6) * t268;
t79 = pkin(5) * t236 - t238 * t259 + t103;
t86 = -pkin(10) * t324 + t104;
t35 = -t263 * t86 + t267 * t79;
t36 = t263 * t79 + t267 * t86;
t288 = -t264 * t78 - t268 * t77;
t218 = t244 * t270 + t315;
t211 = -mrSges(5,2) * t260 + t359;
t282 = t153 * t268 - t154 * t264 + t211;
t281 = t288 * mrSges(6,3);
t280 = -t195 * t264 + t238 * t316;
t30 = t268 * t107 + t264 * t123 + t171 * t316 - t183 * t317;
t108 = qJD(4) * t183 - t221 * t265 + t269 * t222;
t273 = qJD(5) * t288 + t409;
t272 = m(6) * t273;
t115 = t344 + t346 + t347;
t168 = t223 + t337 + t341;
t41 = -pkin(5) * t125 + t73;
t46 = t124 * Ifges(6,4) + t125 * Ifges(6,2) + t185 * Ifges(6,6);
t47 = t124 * Ifges(6,1) + t125 * Ifges(6,4) + t185 * Ifges(6,5);
t59 = t345 + t351 + t352;
t271 = (-Ifges(7,5) * t159 - Ifges(7,6) * t158) * t370 + (-t296 - mrSges(5,1)) * t73 - (Ifges(5,1) * t228 + t115 - t340 + t59) * t229 / 0.2e1 - (-Ifges(5,2) * t229 + t168 + t223 + t322) * t228 / 0.2e1 + (Ifges(7,5) * t237 - Ifges(7,6) * t283 + t289) * t374 + (Ifges(7,4) * t237 - Ifges(7,2) * t283) * t388 + (Ifges(7,1) * t237 - Ifges(7,4) * t283) * t389 + t41 * (mrSges(7,1) * t283 + mrSges(7,2) * t237) - t283 * t391 + (-Ifges(7,1) * t193 - Ifges(7,4) * t194) * t375 + (-Ifges(7,4) * t193 - Ifges(7,2) * t194) * t377 + (-Ifges(7,5) * t193 - Ifges(7,6) * t194) * t369 + (t326 * t77 + t327 * t78 + t409) * mrSges(6,3) + (-t19 * t411 + t20 * t412 - t237 * t4 - t3 * t283) * mrSges(7,3) + (-mrSges(7,1) * t412 + mrSges(7,2) * t411) * t92 + (-Ifges(7,1) * t159 - Ifges(7,4) * t158) * t376 + (-Ifges(7,4) * t159 - Ifges(7,2) * t158) * t378 + t405 * qJD(5) - t158 * t387 + t237 * t390 + t291 * t379 + t293 * t380 - t193 * t384 - t159 * t385 - t194 * t386 + t46 * t364 + t140 * t359 + t264 * t47 / 0.2e1 - Ifges(5,6) * t185 - Ifges(5,5) * t184 - t72 * mrSges(5,2) + (t290 * t368 + t292 * t373 + t294 * t372 - t419 - t420) * t228 + (Ifges(6,5) * t372 + Ifges(7,5) * t376 + Ifges(6,6) * t373 + Ifges(7,6) * t378 + Ifges(6,3) * t368 + Ifges(7,3) * t370 + t404 + t423) * t229;
t247 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t319;
t245 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t320;
t243 = t256 - t362;
t227 = Ifges(4,1) * t320 + t257 + t331;
t217 = -t244 * t266 + t258;
t214 = t218 * qJD(3);
t213 = -t244 * t318 + t252;
t189 = -mrSges(5,1) * t228 + mrSges(5,2) * t229;
t179 = Ifges(6,3) * t185;
t174 = t237 * t238;
t148 = pkin(5) * t324 - t396;
t111 = t141 + t216;
t102 = mrSges(7,1) * t219 - mrSges(7,3) * t136;
t101 = -mrSges(7,2) * t219 + mrSges(7,3) * t299;
t74 = -mrSges(7,1) * t299 + mrSges(7,2) * t136;
t66 = pkin(5) * t280 + t108;
t56 = t175 * t394 + t237 * t195;
t55 = -t194 * t238 + t195 * t283;
t31 = -qJD(5) * t104 + t300;
t29 = -mrSges(7,2) * t185 + mrSges(7,3) * t40;
t28 = mrSges(7,1) * t185 - mrSges(7,3) * t39;
t23 = -pkin(10) * t280 + t30;
t22 = t267 * t57 - t335;
t21 = -t263 * t57 - t333;
t18 = t195 * t259 + pkin(5) * t196 + (-t170 + (pkin(10) * t238 - t171) * t264) * qJD(5) + t300;
t12 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t6 = -qJD(6) * t36 + t18 * t267 - t23 * t263;
t5 = qJD(6) * t35 + t18 * t263 + t23 * t267;
t1 = [-t396 * t62 + (t140 * t195 - t141 * t196 - t183 * t185 + t184 * t396) * mrSges(5,3) + (t290 * t374 + t292 * t379 + t294 * t380 - Ifges(5,4) * t185 - Ifges(5,1) * t184 + t46 * t365 + t47 * t364 + (mrSges(5,3) + t295) * t73 + (-t16 * t264 - t17 * t268) * mrSges(6,3) + (-t268 * t116 / 0.2e1 + t117 * t365 + t289 * t368 + t291 * t373 + t293 * t372 + t127 * t296 - t395 * mrSges(6,3)) * qJD(5)) * t238 + (Ifges(5,4) * t184 + t178 / 0.2e1 + t179 / 0.2e1 - t72 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) + Ifges(7,3) / 0.2e1) * t185 + t392) * t236 + t242 * (mrSges(5,1) * t185 - mrSges(5,2) * t184) + m(5) * (t107 * t141 - t108 * t140 + t183 * t72 - t348) + m(6) * (t103 * t17 + t104 * t16 + t108 * t127 + t30 * t78 + t31 * t77 - t348) + (-Ifges(7,4) * t175 - Ifges(7,2) * t174) * t388 + (-Ifges(7,1) * t175 - Ifges(7,4) * t174) * t389 + (-Ifges(7,5) * t175 - Ifges(7,6) * t174) * t374 + (-t174 * t3 + t175 * t4 - t19 * t55 + t20 * t56) * mrSges(7,3) + t41 * (mrSges(7,1) * t174 - mrSges(7,2) * t175) + (t304 * t214 + (-t251 * t247 + t307 - t304 * t218 + (t310 * mrSges(4,1) - 0.3e1 / 0.2e1 * t358 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t270) * qJD(1) + (t189 + 0.2e1 * t363 + qJD(1) * (mrSges(5,1) * t236 + mrSges(5,2) * t238)) * pkin(3) - t298) * qJD(3)) * t266 - t175 * t390 - t174 * t391 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t375 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t377 + t55 * t384 + t56 * t386 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t369 + m(7) * (t148 * t41 + t19 * t6 + t20 * t5 + t3 * t36 + t35 * t4 + t66 * t92) + t107 * t211 + t148 * t12 + t30 * t153 + t31 * t154 + t5 * t101 + t6 * t102 + t103 * t80 + t104 * t81 + t92 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t66 * t74 + t35 * t28 + t36 * t29 + (t344 / 0.2e1 + t347 / 0.2e1 + t346 / 0.2e1 + t345 / 0.2e1 + t351 / 0.2e1 + t115 / 0.2e1 + t59 / 0.2e1 + t352 / 0.2e1 - t404 - t418) * t196 - (t341 / 0.2e1 + t223 / 0.2e1 + t168 / 0.2e1 + t281 + t405 + t419) * t195 - t321 * t108 + (t304 * t213 + (-t251 * t245 + t227 / 0.2e1 + 0.3e1 / 0.2e1 * t257 + t308 - t304 * t217 + 0.2e1 * t338) * qJD(3)) * t270; t55 * t101 + t56 * t102 - t174 * t28 - t175 * t29 + (-t184 * mrSges(5,3) + t12 + t62) * t236 + (t74 - t321) * t196 - t282 * t195 + (-t266 * t245 + t270 * t247 + (-t266 ^ 2 - t270 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t185 * mrSges(5,3) + t406) * t238 + m(7) * (-t174 * t4 - t175 * t3 + t19 * t56 + t196 * t92 + t20 * t55 + t236 * t41) + m(5) * (-t140 * t196 - t141 * t195 + t238 * t72 + t339) + m(6) * (t127 * t196 - t195 * t395 + t238 * t273 + t339) + m(4) * (t213 * t266 - t214 * t270 + (-t217 * t266 + t218 * t270) * qJD(3)); t271 + (m(5) * (t265 * t72 - t269 * t73) + (t184 * t269 - t185 * t265) * mrSges(5,3) + ((-m(5) * t140 + m(6) * t127 - t321) * t265 + (m(5) * t141 + m(6) * t395 + t282) * t269) * qJD(4)) * pkin(3) + t281 * qJD(5) + t141 * t342 + t416 * t101 - m(5) * (-t140 * t145 + t141 * t146) - m(6) * (t127 * t145 + t77 * t82 + t78 * t83) + t410 * t74 + t415 * t102 + t243 * t12 + t218 * t245 - t217 * t247 - t213 * mrSges(4,2) - t214 * mrSges(4,1) - t146 * t211 + t186 * t28 + t187 * t29 - t83 * t153 - t82 * t154 + t321 * t145 + ((t217 * mrSges(4,3) + t308 - t338 - t227 / 0.2e1 - t257 / 0.2e1) * t270 + (t218 * mrSges(4,3) + t307 + (t358 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t270) * qJD(1) + (-t189 - t363) * pkin(3) + t298) * t266) * qJD(1) + t413 * (-pkin(4) - t362) + (t272 + t406) * t254 + (t186 * t4 + t187 * t3 + t415 * t19 + t20 * t416 + t243 * t41 + t410 * t92) * m(7); t271 + pkin(9) * t272 + ((-t77 * mrSges(6,3) - pkin(9) * t154) * t268 + (-t78 * mrSges(6,3) + pkin(5) * t74 - pkin(9) * t153) * t264) * qJD(5) + t286 * pkin(9) + (t321 + t342) * t141 - m(6) * (t127 * t141 + t77 * t84 + t78 * t85) + t397 * t102 + t398 * t101 + t256 * t12 + t205 * t28 + t206 * t29 - t140 * t211 - t85 * t153 - t84 * t154 - t111 * t74 - t413 * pkin(4) + (t205 * t4 + t206 * t3 + t256 * t41 + (-t111 + t311) * t92 + t398 * t20 + t397 * t19) * m(7); (-Ifges(6,2) * t204 + t117 + t199) * t373 - t136 * t387 + t299 * t385 + t116 * t371 + (Ifges(6,1) * t203 - t357) * t372 + t407 * t378 + (Ifges(6,5) * t203 - Ifges(6,6) * t204) * t368 + t392 - m(7) * (t19 * t21 + t20 * t22) + (t203 * t77 + t204 * t78) * mrSges(6,3) - t127 * (mrSges(6,1) * t204 + mrSges(6,2) * t203) - t77 * t153 + t78 * t154 + t179 - t22 * t101 - t21 * t102 + (-t204 * t74 + t263 * t29 + t267 * t28 + (t101 * t267 - t102 * t263) * qJD(6) + (-t204 * t92 + t263 * t3 + t267 * t4 + (-t19 * t263 + t20 * t267) * qJD(6)) * m(7)) * pkin(5) + t417; t60 * t375 - t19 * t101 + t20 * t102 + (t407 + t61) * t378 + t408 + t417;];
tauc  = t1(:);
