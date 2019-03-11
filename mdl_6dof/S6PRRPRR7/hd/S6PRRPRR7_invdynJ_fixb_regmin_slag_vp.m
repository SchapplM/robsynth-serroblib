% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:15
% EndTime: 2019-03-08 22:34:28
% DurationCPUTime: 5.69s
% Computational Cost: add. (3228->472), mult. (7153->663), div. (0->0), fcn. (5448->14), ass. (0->257)
t177 = sin(pkin(6));
t186 = cos(qJ(2));
t308 = t177 * t186;
t181 = sin(qJ(3));
t277 = t181 * qJD(4);
t185 = cos(qJ(3));
t284 = qJD(3) * t185;
t212 = -qJ(4) * t284 - t277;
t182 = sin(qJ(2));
t276 = qJD(1) * qJD(2);
t257 = t182 * t276;
t232 = -qJDD(1) * t308 + t177 * t257;
t275 = qJD(2) * qJD(3);
t254 = t181 * t275;
t217 = pkin(3) * t254 + t232;
t251 = -t181 * qJ(4) - pkin(2);
t219 = pkin(3) * t185 - t251;
t364 = qJDD(2) * t219;
t35 = qJD(2) * t212 + t217 - t364;
t371 = -g(3) * t308 - t35;
t256 = t186 * t276;
t178 = cos(pkin(6));
t288 = qJD(3) * t178;
t319 = qJDD(2) * pkin(8);
t369 = qJD(1) * t288 + t319 + (qJDD(1) * t182 + t256) * t177;
t272 = t185 * qJDD(2);
t367 = t254 - t272;
t324 = cos(pkin(11));
t246 = t324 * t182;
t176 = sin(pkin(11));
t311 = t176 * t186;
t101 = t178 * t311 + t246;
t245 = t324 * t186;
t312 = t176 * t182;
t99 = -t178 * t245 + t312;
t238 = g(1) * t101 + g(2) * t99;
t294 = qJD(1) * t186;
t310 = t177 * t182;
t266 = qJD(1) * t310;
t130 = qJD(2) * pkin(8) + t266;
t295 = qJD(1) * t178;
t298 = -t181 * t130 + t185 * t295;
t365 = qJD(4) - t298;
t342 = pkin(4) + pkin(8);
t179 = sin(qJ(6));
t180 = sin(qJ(5));
t183 = cos(qJ(6));
t184 = cos(qJ(5));
t123 = t179 * t184 + t180 * t183;
t209 = t123 * t181;
t352 = qJD(5) + qJD(6);
t332 = -qJD(2) * t209 - t352 * t123;
t287 = qJD(3) * t180;
t291 = qJD(2) * t185;
t120 = t184 * t291 + t287;
t258 = t180 * t291;
t285 = qJD(3) * t184;
t122 = -t258 + t285;
t227 = t120 * t179 - t183 * t122;
t50 = t183 * t120 + t122 * t179;
t368 = t227 * t50;
t75 = -qJD(3) * pkin(3) + t365;
t366 = qJD(2) * t219;
t363 = t227 ^ 2 - t50 ^ 2;
t293 = qJD(2) * t181;
t158 = qJD(5) + t293;
t149 = qJD(6) + t158;
t278 = qJD(6) * t183;
t279 = qJD(6) * t179;
t46 = -qJD(5) * t120 + t184 * qJDD(3) + t367 * t180;
t47 = -qJD(5) * t258 + qJDD(3) * t180 + (qJD(3) * qJD(5) - t367) * t184;
t7 = -t120 * t278 - t122 * t279 - t179 * t47 + t183 * t46;
t362 = t149 * t50 + t7;
t270 = t181 * t310;
t108 = -t178 * t185 + t270;
t187 = -pkin(3) - pkin(9);
t300 = pkin(4) * t293 + t365;
t48 = t187 * qJD(3) + t300;
t115 = t187 * t185 + t251;
t265 = t177 * t294;
t71 = qJD(2) * t115 - t265;
t22 = t180 * t48 + t184 * t71;
t14 = -pkin(10) * t120 + t22;
t12 = t14 * t279;
t175 = qJ(5) + qJ(6);
t166 = sin(t175);
t167 = cos(t175);
t172 = qJD(3) * qJ(4);
t83 = t185 * t130 + t181 * t295;
t70 = pkin(4) * t291 + t83;
t56 = t172 + t70;
t36 = pkin(5) * t120 + t56;
t100 = t178 * t246 + t311;
t247 = t177 * t324;
t61 = t100 * t181 + t185 * t247;
t102 = -t178 * t312 + t245;
t309 = t177 * t185;
t63 = t102 * t181 - t176 * t309;
t361 = t36 * t50 - g(1) * (-t101 * t167 - t166 * t63) - g(2) * (-t166 * t61 - t167 * t99) - g(3) * (-t108 * t166 + t167 * t308) + t12;
t253 = t185 * t275;
t273 = t181 * qJDD(2);
t208 = t253 + t273;
t119 = qJDD(5) + t208;
t274 = qJDD(1) * t178;
t351 = t130 * t284 + t181 * t369;
t206 = -t185 * t274 + t351;
t201 = qJDD(4) + t206;
t17 = t208 * pkin(4) + t187 * qJDD(3) + t201;
t16 = t184 * t17;
t322 = qJ(4) * t185;
t233 = pkin(9) * t181 - t322;
t200 = qJD(3) * t233 - t277;
t25 = qJD(2) * t200 + qJDD(2) * t115 + t217;
t198 = -qJD(5) * t22 - t180 * t25 + t16;
t2 = pkin(5) * t119 - pkin(10) * t46 + t198;
t282 = qJD(5) * t184;
t271 = -t180 * t17 - t184 * t25 - t48 * t282;
t283 = qJD(5) * t180;
t215 = t71 * t283 + t271;
t3 = -pkin(10) * t47 - t215;
t267 = -t179 * t3 + t183 * t2;
t21 = -t180 * t71 + t184 * t48;
t13 = -pkin(10) * t122 + t21;
t10 = pkin(5) * t158 + t13;
t328 = t14 * t183;
t5 = t10 * t179 + t328;
t360 = t36 * t227 - g(1) * (-t101 * t166 + t167 * t63) - g(2) * (-t166 * t99 + t167 * t61) - g(3) * (t108 * t167 + t166 * t308) - qJD(6) * t5 + t267;
t197 = qJD(6) * t227 - t179 * t46 - t183 * t47;
t359 = -t149 * t227 + t197;
t109 = t178 * t181 + t182 * t309;
t189 = qJD(2) ^ 2;
t222 = qJDD(2) * t186 - t182 * t189;
t252 = t186 * t275;
t262 = qJD(2) * t308;
t65 = -qJD(3) * t270 + (t262 + t288) * t185;
t358 = qJD(3) * t65 + qJDD(3) * t109 + t177 * (t181 * t222 + t185 * t252);
t66 = qJD(3) * t109 + t181 * t262;
t357 = -qJD(3) * t66 - qJDD(3) * t108 + t177 * (-t181 * t252 + t185 * t222);
t76 = -t172 - t83;
t129 = t342 * t284;
t141 = t342 * t181;
t305 = t181 * t186;
t202 = t177 * (t180 * t305 + t182 * t184);
t286 = qJD(3) * t181;
t163 = pkin(3) * t286;
t85 = t163 + t200;
t356 = qJD(1) * t202 + t115 * t283 - t180 * t129 - t141 * t282 - t184 * t85;
t302 = t184 * t186;
t203 = t177 * (-t180 * t182 + t181 * t302);
t355 = -qJD(1) * t203 + t184 * t129 - t180 * t85;
t98 = t184 * t119;
t354 = -t158 * t283 + t98;
t84 = -t265 - t366;
t353 = t84 * t293 + qJDD(4);
t350 = g(1) * t63 + g(2) * t61 + g(3) * t108;
t349 = t187 * t119 + t158 * t56;
t188 = qJD(3) ^ 2;
t207 = pkin(8) * t188 - t238;
t348 = 0.2e1 * qJDD(2) * pkin(2) + t177 * (-g(3) * t186 + t257) - t207 - t232;
t249 = qJD(6) * t10 + t3;
t346 = t179 * t2 + t183 * t249;
t97 = t163 + t212;
t344 = qJD(2) * (-t97 + t266) + t364 - t207 + t371;
t343 = t185 * t352;
t125 = t180 * t141;
t306 = t180 * t181;
t220 = pkin(5) * t185 - pkin(10) * t306;
t250 = pkin(10) * t185 - t115;
t338 = t220 * qJD(3) + (t250 * t184 - t125) * qJD(5) + t355;
t334 = pkin(10) - t187;
t281 = qJD(5) * t185;
t259 = t180 * t281;
t333 = -(t181 * t285 + t259) * pkin(10) + t356;
t261 = t184 * t293;
t304 = t183 * t184;
t307 = t179 * t180;
t331 = -t179 * t283 - t180 * t279 + t183 * t261 - t293 * t307 + t304 * t352;
t164 = pkin(3) * t293;
t92 = qJD(2) * t233 + t164;
t330 = t180 * t70 + t184 * t92;
t329 = qJD(2) * pkin(2);
t327 = t46 * t184;
t268 = -pkin(5) * t184 - pkin(4);
t326 = pkin(5) * t282 - t268 * t293 + t365;
t325 = t184 * t115 + t125;
t323 = pkin(8) * qJDD(3);
t321 = qJD(3) * t83;
t318 = qJDD(3) * pkin(3);
t317 = t119 * t180;
t316 = t120 * t158;
t315 = t122 * t158;
t314 = t166 * t181;
t313 = t167 * t181;
t303 = t184 * t185;
t299 = qJDD(1) - g(3);
t142 = t342 * t185;
t173 = t181 ^ 2;
t174 = t185 ^ 2;
t297 = t173 - t174;
t296 = t173 + t174;
t292 = qJD(2) * t182;
t290 = qJD(3) * t120;
t289 = qJD(3) * t122;
t280 = qJD(5) * t187;
t269 = t181 * t189 * t185;
t264 = t185 * t294;
t263 = t177 * t292;
t133 = t334 * t184;
t242 = -t130 * t286 + t181 * t274 + t185 * t369;
t112 = qJDD(6) + t119;
t226 = -t304 + t307;
t241 = -t226 * t112 + t332 * t149;
t237 = g(1) * t102 + g(2) * t100;
t132 = t334 * t180;
t58 = t184 * t70;
t235 = qJD(2) * t220 - qJD(6) * t132 - t180 * t92 - t334 * t283 + t58;
t234 = pkin(10) * t261 + t352 * t133 + t330;
t126 = t184 * t141;
t37 = t181 * pkin(5) + t250 * t180 + t126;
t44 = -pkin(10) * t303 + t325;
t231 = t179 * t37 + t183 * t44;
t67 = t108 * t184 + t180 * t308;
t68 = -t108 * t180 + t177 * t302;
t230 = t179 * t68 + t183 * t67;
t229 = t179 * t67 - t183 * t68;
t225 = t158 * t180;
t170 = qJDD(3) * qJ(4);
t171 = qJD(3) * qJD(4);
t26 = -t170 - t171 - t242;
t214 = -t158 * t282 - t317;
t213 = -t123 * t112 - t331 * t149;
t62 = t100 * t185 - t181 * t247;
t64 = t176 * t177 * t181 + t102 * t185;
t211 = -g(1) * t64 - g(2) * t62 - g(3) * t109;
t18 = -t367 * pkin(4) - t26;
t204 = t18 + t211;
t131 = -t265 - t329;
t196 = -t323 + (t131 + t265 - t329) * qJD(3);
t195 = t323 + (-t265 - t84 + t366) * qJD(3);
t194 = -qJD(3) * t298 + t211 + t242;
t193 = -t206 + t350;
t27 = t201 - t318;
t191 = t27 * t181 - t26 * t185 + (t181 * t76 + t185 * t75) * qJD(3) - t237;
t160 = pkin(5) * t180 + qJ(4);
t128 = t342 * t286;
t127 = -qJ(4) * t291 + t164;
t107 = pkin(5) * t303 + t142;
t94 = t123 * t185;
t93 = t226 * t185;
t72 = -pkin(5) * t259 + (-pkin(8) + t268) * t286;
t30 = t123 * t343 - t226 * t286;
t29 = qJD(3) * t209 + t226 * t343;
t20 = qJD(5) * t67 + t180 * t66 + t184 * t263;
t19 = qJD(5) * t68 - t180 * t263 + t184 * t66;
t6 = pkin(5) * t47 + t18;
t4 = t10 * t183 - t14 * t179;
t1 = [t299, 0, t222 * t177 (-qJDD(2) * t182 - t186 * t189) * t177, 0, 0, 0, 0, 0, t357, -t358 (t108 * t181 + t109 * t185) * qJDD(2) + (t181 * t66 + t185 * t65 + (t108 * t185 - t109 * t181) * qJD(3)) * qJD(2), -t357, t358, t108 * t27 - t109 * t26 - t65 * t76 + t66 * t75 - g(3) + (-t186 * t35 + t84 * t292) * t177, 0, 0, 0, 0, 0, t109 * t47 + t119 * t67 + t120 * t65 + t158 * t19, t109 * t46 + t119 * t68 + t122 * t65 - t158 * t20, 0, 0, 0, 0, 0 (-qJD(6) * t229 - t179 * t20 + t183 * t19) * t149 + t230 * t112 + t65 * t50 - t109 * t197 -(qJD(6) * t230 + t179 * t19 + t183 * t20) * t149 - t229 * t112 - t65 * t227 + t109 * t7; 0, qJDD(2), t299 * t308 + t238, -t299 * t310 + t237, qJDD(2) * t173 + 0.2e1 * t181 * t253, 0.2e1 * t181 * t272 - 0.2e1 * t297 * t275, qJDD(3) * t181 + t185 * t188, qJDD(3) * t185 - t181 * t188, 0, t196 * t181 + t185 * t348, -t181 * t348 + t196 * t185, t296 * t319 + (-g(3) * t182 - t296 * t256) * t177 + t191, t195 * t181 - t185 * t344, t181 * t344 + t195 * t185, t84 * t97 + t191 * pkin(8) + ((-pkin(8) * g(3) - qJD(1) * t84) * t182 + (-t181 * t75 + t185 * t76) * t294) * t177 + (t238 + t371) * t219, -t46 * t180 * t185 + (t180 * t286 - t184 * t281) * t122 (-t120 * t180 + t122 * t184) * t286 + (t180 * t47 - t327 + (t120 * t184 + t122 * t180) * qJD(5)) * t185 (t158 * t287 + t46) * t181 + (t214 + t289) * t185 (t158 * t285 - t47) * t181 + (-t290 - t354) * t185, t119 * t181 + t158 * t284 (-t115 * t180 + t126) * t119 - t128 * t120 + t142 * t47 - t237 * t184 + (-t56 * t285 + t16 + (t238 - t25) * t180) * t181 - g(3) * t202 + t355 * t158 + (-t325 * t158 - t22 * t181) * qJD(5) + (qJD(3) * t21 - t120 * t265 + t18 * t184 - t56 * t283) * t185, -t325 * t119 - t128 * t122 + t142 * t46 + t237 * t180 + (t238 * t184 + (qJD(3) * t56 + qJD(5) * t71) * t180 + t271) * t181 - g(3) * t203 + t356 * t158 + (-qJD(3) * t22 - t122 * t265 - t18 * t180 - t282 * t56) * t185, -t227 * t29 - t7 * t94, -t197 * t94 - t227 * t30 - t29 * t50 + t7 * t93, -t112 * t94 + t149 * t29 + t181 * t7 - t227 * t284, t112 * t93 + t149 * t30 + t181 * t197 - t284 * t50, t112 * t181 + t149 * t284 (-t179 * t44 + t183 * t37) * t112 + t267 * t181 + t4 * t284 + t72 * t50 - t107 * t197 - t6 * t93 - t36 * t30 - g(1) * (-t101 * t314 + t102 * t167) - g(2) * (t100 * t167 - t314 * t99) + (t333 * t179 + t338 * t183) * t149 + (-t149 * t231 - t181 * t5) * qJD(6) + (-t50 * t264 - g(3) * (t166 * t305 + t167 * t182)) * t177, -t231 * t112 - (-t12 + t346) * t181 - t5 * t284 - t72 * t227 + t107 * t7 - t6 * t94 + t36 * t29 - g(1) * (-t101 * t313 - t102 * t166) - g(2) * (-t100 * t166 - t313 * t99) + ((-qJD(6) * t37 + t333) * t183 + (qJD(6) * t44 - t338) * t179) * t149 + (t227 * t264 - g(3) * (-t166 * t182 + t167 * t305)) * t177; 0, 0, 0, 0, -t269, t297 * t189, t273, t272, qJDD(3), -t131 * t293 + t193 + t321, -t131 * t291 - t194 (-pkin(3) * t181 + t322) * qJDD(2), -0.2e1 * t318 - t321 + (-qJD(2) * t127 - t274) * t185 - t350 + t351 + t353, 0.2e1 * t170 + 0.2e1 * t171 + (t127 * t181 + t185 * t84) * qJD(2) + t194, -t26 * qJ(4) - t27 * pkin(3) - t84 * t127 - t75 * t83 - g(1) * (-pkin(3) * t63 + qJ(4) * t64) - g(2) * (-pkin(3) * t61 + qJ(4) * t62) - g(3) * (-pkin(3) * t108 + qJ(4) * t109) - t365 * t76, -t122 * t225 + t327 (-t47 - t315) * t184 + (-t46 + t316) * t180 (-t122 * t185 - t158 * t306) * qJD(2) + t354 (-t158 * t181 * t184 + t120 * t185) * qJD(2) + t214, -t158 * t291, -t21 * t291 + qJ(4) * t47 - t58 * t158 + t300 * t120 + t349 * t184 + ((t92 - t280) * t158 + t204) * t180, qJ(4) * t46 + t330 * t158 + t22 * t291 + t300 * t122 - t349 * t180 + (-t158 * t280 + t204) * t184, -t226 * t7 - t227 * t332, -t123 * t7 - t197 * t226 + t227 * t331 - t332 * t50, t227 * t291 + t241, t291 * t50 + t213, -t149 * t291 (t132 * t179 - t133 * t183) * t112 - t160 * t197 + t6 * t123 - t4 * t291 + t326 * t50 + t331 * t36 + (t179 * t234 - t183 * t235) * t149 + t211 * t166 -(-t132 * t183 - t133 * t179) * t112 + t160 * t7 - t6 * t226 + t5 * t291 - t326 * t227 + t332 * t36 + (t179 * t235 + t183 * t234) * t149 + t211 * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, qJDD(3) + t269, -t173 * t189 - t188, qJD(3) * t76 - t193 - t318 + t353, 0, 0, 0, 0, 0, -t158 * t225 - t290 + t98, -t158 ^ 2 * t184 - t289 - t317, 0, 0, 0, 0, 0, -qJD(3) * t50 + t241, qJD(3) * t227 + t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 * t120, -t120 ^ 2 + t122 ^ 2, t46 + t316, t315 - t47, t119, t22 * t158 - t56 * t122 - g(1) * (-t101 * t180 + t184 * t63) - g(2) * (-t180 * t99 + t184 * t61) - g(3) * t67 + t198, t21 * t158 + t56 * t120 - g(1) * (-t101 * t184 - t180 * t63) - g(2) * (-t180 * t61 - t184 * t99) - g(3) * t68 + t215, -t368, t363, t362, t359, t112 -(-t13 * t179 - t328) * t149 + (t112 * t183 - t122 * t50 - t149 * t279) * pkin(5) + t360 (-t14 * t149 - t2) * t179 + (t13 * t149 - t249) * t183 + (-t112 * t179 + t122 * t227 - t149 * t278) * pkin(5) + t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t368, t363, t362, t359, t112, t149 * t5 + t360, t149 * t4 - t346 + t361;];
tau_reg  = t1;
