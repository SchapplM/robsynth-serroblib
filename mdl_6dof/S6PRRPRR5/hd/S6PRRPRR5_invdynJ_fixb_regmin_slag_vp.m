% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:51
% EndTime: 2019-03-08 22:21:06
% DurationCPUTime: 7.03s
% Computational Cost: add. (5744->541), mult. (13463->782), div. (0->0), fcn. (11063->18), ass. (0->252)
t236 = cos(qJ(3));
t320 = qJD(2) * t236;
t209 = -qJD(5) + t320;
t204 = -qJD(6) + t209;
t226 = sin(pkin(12));
t228 = cos(pkin(12));
t310 = t228 * qJD(3);
t232 = sin(qJ(3));
t321 = qJD(2) * t232;
t177 = t226 * t321 - t310;
t319 = qJD(3) * t226;
t179 = t228 * t321 + t319;
t231 = sin(qJ(5));
t235 = cos(qJ(5));
t102 = t177 * t231 - t179 * t235;
t103 = t235 * t177 + t179 * t231;
t230 = sin(qJ(6));
t234 = cos(qJ(6));
t370 = t102 * t230 - t234 * t103;
t384 = t204 * t370;
t272 = pkin(3) * t232 - qJ(4) * t236;
t156 = qJD(3) * t272 - t232 * qJD(4);
t233 = sin(qJ(2));
t318 = qJD(3) * t232;
t304 = pkin(8) * t318;
t227 = sin(pkin(6));
t325 = qJD(1) * t227;
t237 = cos(qJ(2));
t333 = t236 * t237;
t331 = t228 * t156 + t226 * t304 - (-t226 * t333 + t228 * t233) * t325;
t383 = t226 * t156 - (t226 * t233 + t228 * t333) * t325;
t335 = t228 * t236;
t263 = pkin(4) * t232 - pkin(9) * t335;
t382 = -qJD(3) * t263 - t331;
t336 = t228 * t232;
t339 = t226 * t236;
t381 = (-pkin(8) * t336 - pkin(9) * t339) * qJD(3) + t383;
t300 = t226 * t320;
t187 = t272 * qJD(2);
t189 = qJD(2) * pkin(8) + t233 * t325;
t229 = cos(pkin(6));
t334 = t229 * t236;
t367 = qJD(1) * t334 - t232 * t189;
t86 = t226 * t187 + t228 * t367;
t74 = -pkin(9) * t300 + t86;
t380 = qJD(4) * t228 - t74;
t267 = t102 * t234 + t103 * t230;
t379 = t267 * t370;
t378 = t102 * t209;
t377 = t103 * t209;
t376 = t204 * t267;
t341 = t226 * t231;
t185 = -t235 * t228 + t341;
t257 = t236 * t185;
t329 = qJD(2) * t257 - t185 * qJD(5);
t186 = t226 * t235 + t228 * t231;
t256 = t186 * t236;
t328 = -qJD(2) * t256 + t186 * qJD(5);
t375 = t267 ^ 2 - t370 ^ 2;
t311 = qJD(6) * t234;
t312 = qJD(6) * t230;
t216 = t228 * qJDD(3);
t308 = qJD(2) * qJD(3);
t292 = t236 * t308;
t305 = t232 * qJDD(2);
t255 = t292 + t305;
t140 = t226 * t255 - t216;
t306 = t226 * qJDD(3);
t141 = t228 * t255 + t306;
t313 = qJD(5) * t235;
t315 = qJD(5) * t231;
t35 = -t231 * t140 + t235 * t141 - t177 * t313 - t179 * t315;
t36 = -qJD(5) * t102 + t235 * t140 + t231 * t141;
t6 = t102 * t312 - t103 * t311 - t230 * t36 + t234 * t35;
t374 = t6 + t384;
t324 = qJD(1) * t232;
t205 = t229 * t324;
t138 = t236 * t189 + t205;
t129 = qJD(3) * qJ(4) + t138;
t262 = pkin(3) * t236 + qJ(4) * t232 + pkin(2);
t323 = qJD(1) * t237;
t297 = t227 * t323;
t139 = -qJD(2) * t262 - t297;
t65 = -t129 * t226 + t228 * t139;
t41 = -pkin(4) * t320 - pkin(9) * t179 + t65;
t66 = t228 * t129 + t226 * t139;
t49 = -pkin(9) * t177 + t66;
t21 = t231 * t41 + t235 * t49;
t13 = -pkin(10) * t103 + t21;
t11 = t13 * t312;
t348 = sin(pkin(11));
t280 = t348 * t237;
t349 = cos(pkin(11));
t283 = t349 * t233;
t160 = t229 * t283 + t280;
t285 = t227 * t349;
t119 = t160 * t236 - t232 * t285;
t281 = t348 * t233;
t282 = t349 * t237;
t162 = -t229 * t281 + t282;
t284 = t227 * t348;
t121 = t162 * t236 + t232 * t284;
t159 = -t229 * t282 + t281;
t161 = t229 * t280 + t283;
t338 = t227 * t233;
t167 = t229 * t232 + t236 * t338;
t223 = pkin(12) + qJ(5);
t221 = qJ(6) + t223;
t212 = sin(t221);
t213 = cos(t221);
t337 = t227 * t237;
t125 = -qJD(3) * pkin(3) + qJD(4) - t367;
t87 = pkin(4) * t177 + t125;
t40 = pkin(5) * t103 + t87;
t373 = -t40 * t370 - g(1) * (-t121 * t213 - t161 * t212) - g(2) * (-t119 * t213 - t159 * t212) - g(3) * (-t167 * t213 + t212 * t337) + t11;
t220 = t236 * qJDD(2);
t365 = t232 * t308 - t220;
t184 = qJDD(5) + t365;
t309 = qJD(1) * qJD(2);
t149 = qJDD(2) * pkin(8) + (qJDD(1) * t233 + t237 * t309) * t227;
t307 = qJDD(1) * t229;
t290 = t232 * t307;
t53 = t290 + qJDD(3) * qJ(4) + t236 * t149 + (qJD(4) + t367) * qJD(3);
t294 = t233 * t309;
t271 = -qJDD(1) * t337 + t227 * t294;
t73 = qJD(2) * t156 - qJDD(2) * t262 + t271;
t28 = -t226 * t53 + t228 * t73;
t19 = pkin(4) * t365 - t141 * pkin(9) + t28;
t29 = t226 * t73 + t228 * t53;
t23 = -pkin(9) * t140 + t29;
t288 = t235 * t19 - t231 * t23;
t245 = -qJD(5) * t21 + t288;
t2 = t184 * pkin(5) - t35 * pkin(10) + t245;
t260 = -t231 * t19 - t235 * t23 - t41 * t313 + t315 * t49;
t3 = -pkin(10) * t36 - t260;
t303 = t234 * t2 - t230 * t3;
t20 = -t231 * t49 + t235 * t41;
t12 = pkin(10) * t102 + t20;
t10 = -pkin(5) * t209 + t12;
t351 = t13 * t234;
t5 = t10 * t230 + t351;
t372 = t40 * t267 - g(1) * (-t121 * t212 + t161 * t213) - g(2) * (-t119 * t212 + t159 * t213) - g(3) * (-t167 * t212 - t213 * t337) - qJD(6) * t5 + t303;
t243 = qJD(6) * t267 - t230 * t35 - t234 * t36;
t371 = t243 + t376;
t176 = t228 * t262;
t109 = -pkin(9) * t336 - t176 + (-pkin(8) * t226 - pkin(4)) * t236;
t143 = pkin(8) * t335 - t226 * t262;
t340 = t226 * t232;
t124 = -pkin(9) * t340 + t143;
t369 = t109 * t313 - t124 * t315 - t231 * t382 + t381 * t235;
t368 = t381 * t231 + t235 * t382;
t356 = pkin(9) + qJ(4);
t193 = t356 * t226;
t194 = t356 * t228;
t327 = -t231 * t193 + t235 * t194;
t265 = qJD(4) * t226 + qJD(5) * t194;
t85 = t228 * t187 - t226 * t367;
t64 = qJD(2) * t263 + t85;
t366 = -t193 * t313 + t380 * t235 + (-t265 - t64) * t231;
t364 = -qJDD(3) * pkin(3) + qJDD(4);
t238 = qJD(3) ^ 2;
t358 = g(2) * t159;
t359 = g(1) * t161;
t274 = t358 + t359;
t363 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t238 + t227 * (-g(3) * t237 + t294) - t271 + t274;
t289 = qJD(6) * t10 + t3;
t362 = t230 * t2 + t234 * t289;
t350 = t231 * t109 + t235 * t124;
t314 = qJD(5) * t232;
t90 = -qJD(3) * t257 - t186 * t314;
t361 = -pkin(5) * t318 + t90 * pkin(10) + qJD(5) * t350 + t368;
t91 = qJD(3) * t256 + t313 * t336 - t314 * t341;
t360 = -pkin(10) * t91 + t369;
t107 = t234 * t185 + t186 * t230;
t355 = -qJD(6) * t107 - t230 * t328 + t234 * t329;
t108 = -t185 * t230 + t186 * t234;
t354 = qJD(6) * t108 + t230 * t329 + t234 * t328;
t352 = qJD(2) * pkin(2);
t345 = t212 * t236;
t344 = t213 * t236;
t218 = sin(t223);
t343 = t218 * t236;
t219 = cos(t223);
t342 = t219 * t236;
t332 = qJDD(1) - g(3);
t330 = -t228 * t304 + t383;
t317 = qJD(3) * t236;
t171 = (pkin(4) * t226 + pkin(8)) * t317;
t188 = pkin(4) * t340 + t232 * pkin(8);
t224 = t232 ^ 2;
t326 = -t236 ^ 2 + t224;
t322 = qJD(2) * t227;
t112 = pkin(4) * t300 + t138;
t214 = -pkin(4) * t228 - pkin(3);
t301 = t232 * t323;
t299 = t233 * t322;
t298 = t237 * t322;
t296 = pkin(5) * t328 - t112;
t295 = qJ(4) * t220;
t291 = t237 * t308;
t279 = t235 * t109 - t124 * t231;
t277 = -t235 * t193 - t194 * t231;
t58 = t235 * t64;
t89 = -pkin(10) * t185 + t327;
t276 = pkin(5) * t321 + pkin(10) * t329 + t186 * qJD(4) + t327 * qJD(5) + qJD(6) * t89 - t231 * t74 + t58;
t88 = -pkin(10) * t186 + t277;
t275 = -pkin(10) * t328 + qJD(6) * t88 + t366;
t273 = g(1) * t162 + g(2) * t160;
t154 = t185 * t232;
t30 = -pkin(5) * t236 + pkin(10) * t154 + t279;
t153 = t186 * t232;
t34 = -pkin(10) * t153 + t350;
t270 = t230 * t30 + t234 * t34;
t116 = -t167 * t226 - t228 * t337;
t117 = t167 * t228 - t226 * t337;
t51 = t116 * t235 - t117 * t231;
t52 = t116 * t231 + t117 * t235;
t269 = -t230 * t52 + t234 * t51;
t268 = t230 * t51 + t234 * t52;
t80 = t234 * t153 - t154 * t230;
t81 = -t153 * t230 - t154 * t234;
t239 = qJD(2) ^ 2;
t264 = qJDD(2) * t237 - t233 * t239;
t166 = t232 * t338 - t334;
t254 = g(1) * (t162 * t232 - t236 * t284) + g(2) * (t160 * t232 + t236 * t285) + g(3) * t166;
t253 = g(1) * t121 + g(2) * t119 + g(3) * t167;
t252 = qJD(3) * t205 + t232 * t149 + t189 * t317 - t236 * t307;
t56 = t252 + t364;
t250 = t254 - t56;
t249 = g(3) * t337 - t274;
t248 = -g(3) * t338 - t273;
t247 = -qJ(4) * t318 + (qJD(4) - t125) * t236;
t190 = -t297 - t352;
t242 = -pkin(8) * qJDD(3) + (t190 + t297 - t352) * qJD(3);
t37 = pkin(4) * t140 + t56;
t241 = -t252 + t254;
t170 = qJDD(6) + t184;
t147 = pkin(5) * t185 + t214;
t142 = -pkin(8) * t339 - t176;
t123 = qJD(3) * t167 + t232 * t298;
t122 = -qJD(3) * t166 + t236 * t298;
t111 = pkin(5) * t153 + t188;
t84 = t122 * t228 + t226 * t299;
t83 = -t122 * t226 + t228 * t299;
t63 = pkin(5) * t91 + t171;
t27 = qJD(6) * t81 + t230 * t90 + t234 * t91;
t26 = -qJD(6) * t80 - t230 * t91 + t234 * t90;
t18 = -qJD(5) * t52 - t231 * t84 + t235 * t83;
t17 = qJD(5) * t51 + t231 * t83 + t235 * t84;
t14 = pkin(5) * t36 + t37;
t4 = t10 * t234 - t13 * t230;
t1 = [t332, 0, t264 * t227 (-qJDD(2) * t233 - t237 * t239) * t227, 0, 0, 0, 0, 0, -t123 * qJD(3) - t166 * qJDD(3) + (-t232 * t291 + t236 * t264) * t227, -t122 * qJD(3) - t167 * qJDD(3) + (-t232 * t264 - t236 * t291) * t227, -t116 * t220 + t123 * t177 + t166 * t140 + (t116 * t318 - t236 * t83) * qJD(2), t117 * t220 + t123 * t179 + t166 * t141 + (-t117 * t318 + t236 * t84) * qJD(2), -t116 * t141 - t117 * t140 - t177 * t84 - t179 * t83, t116 * t28 + t117 * t29 + t123 * t125 + t166 * t56 + t65 * t83 + t66 * t84 - g(3), 0, 0, 0, 0, 0, t103 * t123 + t166 * t36 - t18 * t209 + t184 * t51, -t102 * t123 + t166 * t35 + t17 * t209 - t184 * t52, 0, 0, 0, 0, 0 -(-qJD(6) * t268 - t230 * t17 + t234 * t18) * t204 + t269 * t170 - t123 * t370 - t166 * t243 (qJD(6) * t269 + t234 * t17 + t230 * t18) * t204 - t268 * t170 - t123 * t267 + t166 * t6; 0, qJDD(2), t332 * t337 + t274, -t332 * t338 + t273, qJDD(2) * t224 + 0.2e1 * t232 * t292, 0.2e1 * t220 * t232 - 0.2e1 * t308 * t326, qJDD(3) * t232 + t236 * t238, qJDD(3) * t236 - t232 * t238, 0, t242 * t232 + t236 * t363, -t232 * t363 + t242 * t236, t248 * t226 + (-t177 * t297 + pkin(8) * t140 + t56 * t226 + (qJD(2) * t142 + t65) * qJD(3)) * t232 + (-t142 * qJDD(2) - t28 + (pkin(8) * t177 + t125 * t226) * qJD(3) - t331 * qJD(2) - t249 * t228) * t236, t248 * t228 + (-t179 * t297 + pkin(8) * t141 + t56 * t228 + (-qJD(2) * t143 - t66) * qJD(3)) * t232 + (t143 * qJDD(2) + t29 + (pkin(8) * t179 + t125 * t228) * qJD(3) + t330 * qJD(2) + t249 * t226) * t236, -t140 * t143 - t142 * t141 - t331 * t179 - t330 * t177 + (-t226 * t66 - t228 * t65) * t317 + (-t226 * t29 - t228 * t28 - t249) * t232, t28 * t142 + t29 * t143 + t330 * t66 + t331 * t65 + t262 * t359 + t262 * t358 + (t125 * t317 + t232 * t56 - t273) * pkin(8) + (-g(3) * pkin(8) * t233 + (-g(3) * t262 - t125 * t324) * t237) * t227, -t102 * t90 - t154 * t35, t102 * t91 - t103 * t90 - t153 * t35 + t154 * t36, -t102 * t318 - t154 * t184 - t209 * t90 - t236 * t35, -t103 * t318 - t153 * t184 + t209 * t91 + t236 * t36, -t184 * t236 - t209 * t318, t279 * t184 - t288 * t236 + t20 * t318 + t171 * t103 + t188 * t36 + t37 * t153 + t87 * t91 - g(1) * (-t161 * t342 + t162 * t218) - g(2) * (-t159 * t342 + t160 * t218) + t368 * t209 + (t209 * t350 + t21 * t236) * qJD(5) + (-t103 * t301 - g(3) * (t218 * t233 + t219 * t333)) * t227, -t350 * t184 - t260 * t236 - t21 * t318 - t171 * t102 + t188 * t35 - t37 * t154 + t87 * t90 - g(1) * (t161 * t343 + t162 * t219) - g(2) * (t159 * t343 + t160 * t219) + t369 * t209 + (t102 * t301 - g(3) * (-t218 * t333 + t219 * t233)) * t227, -t26 * t267 + t6 * t81, t243 * t81 + t26 * t370 + t267 * t27 - t6 * t80, t170 * t81 - t204 * t26 - t236 * t6 - t267 * t318, -t170 * t80 + t204 * t27 - t236 * t243 + t318 * t370, -t170 * t236 - t204 * t318 (-t230 * t34 + t234 * t30) * t170 - t303 * t236 + t4 * t318 - t63 * t370 - t111 * t243 + t14 * t80 + t40 * t27 - g(1) * (-t161 * t344 + t162 * t212) - g(2) * (-t159 * t344 + t160 * t212) + (t360 * t230 + t234 * t361) * t204 + (t204 * t270 + t236 * t5) * qJD(6) + (t370 * t301 - g(3) * (t212 * t233 + t213 * t333)) * t227, -t270 * t170 + (-t11 + t362) * t236 - t5 * t318 - t63 * t267 + t111 * t6 + t14 * t81 + t40 * t26 - g(1) * (t161 * t345 + t162 * t213) - g(2) * (t159 * t345 + t160 * t213) + ((qJD(6) * t30 + t360) * t234 + (-qJD(6) * t34 - t361) * t230) * t204 + (t267 * t301 - g(3) * (-t212 * t333 + t213 * t233)) * t227; 0, 0, 0, 0, -t232 * t239 * t236, t326 * t239, t305, t220, qJDD(3), qJD(3) * t138 - t190 * t321 + t241, -t290 + (-qJD(2) * t190 - t149) * t236 + t253, t226 * t295 - pkin(3) * t140 - t138 * t177 + t250 * t228 + (t226 * t247 - t65 * t232 + t236 * t85) * qJD(2), t228 * t295 - pkin(3) * t141 - t138 * t179 - t250 * t226 + (t228 * t247 + t66 * t232 - t236 * t86) * qJD(2), t86 * t177 + t85 * t179 + (-qJ(4) * t140 - qJD(4) * t177 + t320 * t65 + t29) * t228 + (qJ(4) * t141 + qJD(4) * t179 + t320 * t66 - t28) * t226 - t253, -t125 * t138 - t65 * t85 - t66 * t86 + (-t226 * t65 + t228 * t66) * qJD(4) + t250 * pkin(3) + (-t28 * t226 + t29 * t228 - t253) * qJ(4), -t102 * t329 + t35 * t186, t102 * t328 - t103 * t329 - t35 * t185 - t186 * t36, t102 * t321 + t186 * t184 - t209 * t329, t103 * t321 - t185 * t184 + t209 * t328, t209 * t321, t277 * t184 + t214 * t36 + t37 * t185 - t20 * t321 - t112 * t103 + t328 * t87 + (t58 + t265 * t235 + (-qJD(5) * t193 + t380) * t231) * t209 + t254 * t219, t102 * t112 - t327 * t184 + t37 * t186 + t209 * t366 + t21 * t321 + t214 * t35 - t254 * t218 + t329 * t87, t108 * t6 - t267 * t355, -t107 * t6 + t108 * t243 + t267 * t354 + t355 * t370, t108 * t170 - t204 * t355 + t267 * t321, -t107 * t170 + t204 * t354 - t321 * t370, t204 * t321 (-t230 * t89 + t234 * t88) * t170 - t147 * t243 + t14 * t107 - t4 * t321 - t296 * t370 + t354 * t40 + (t230 * t275 + t234 * t276) * t204 + t254 * t213 -(t230 * t88 + t234 * t89) * t170 + t147 * t6 + t14 * t108 + t5 * t321 - t296 * t267 + t355 * t40 + (-t230 * t276 + t234 * t275) * t204 - t254 * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226 * t305 - t216 + (-t179 + t319) * t320, t228 * t305 + t306 + (t177 + t310) * t320, -t177 ^ 2 - t179 ^ 2, t177 * t66 + t179 * t65 - t241 + t364, 0, 0, 0, 0, 0, t36 + t378, t35 + t377, 0, 0, 0, 0, 0, -t243 + t376, t6 - t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * t103, t102 ^ 2 - t103 ^ 2, t35 - t377, -t36 + t378, t184, -t21 * t209 + t87 * t102 - g(1) * (-t121 * t218 + t161 * t219) - g(2) * (-t119 * t218 + t159 * t219) - g(3) * (-t167 * t218 - t219 * t337) + t245, -t20 * t209 + t87 * t103 - g(1) * (-t121 * t219 - t161 * t218) - g(2) * (-t119 * t219 - t159 * t218) - g(3) * (-t167 * t219 + t218 * t337) + t260, t379, t375, t374, t371, t170 (-t12 * t230 - t351) * t204 + (-t102 * t370 + t170 * t234 + t204 * t312) * pkin(5) + t372 (t13 * t204 - t2) * t230 + (-t12 * t204 - t289) * t234 + (-t102 * t267 - t170 * t230 + t204 * t311) * pkin(5) + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t379, t375, t374, t371, t170, -t5 * t204 + t372, -t4 * t204 - t362 + t373;];
tau_reg  = t1;
