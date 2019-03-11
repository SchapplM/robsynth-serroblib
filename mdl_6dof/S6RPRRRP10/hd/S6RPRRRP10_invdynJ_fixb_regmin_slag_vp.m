% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:36
% EndTime: 2019-03-09 06:32:50
% DurationCPUTime: 5.56s
% Computational Cost: add. (6608->539), mult. (12750->671), div. (0->0), fcn. (8393->10), ass. (0->262)
t182 = sin(qJ(4));
t183 = sin(qJ(3));
t293 = qJD(1) * t183;
t258 = t182 * t293;
t346 = pkin(9) + pkin(8);
t259 = qJD(4) * t346;
t187 = cos(qJ(3));
t233 = pkin(3) * t187 + pkin(8) * t183;
t130 = t233 * qJD(1);
t190 = -pkin(1) - pkin(7);
t155 = t190 * qJD(1) + qJD(2);
t186 = cos(qJ(4));
t304 = t186 * t187;
t298 = t182 * t130 + t155 * t304;
t371 = pkin(9) * t258 + t182 * t259 + t298;
t114 = t186 * t130;
t310 = t183 * t186;
t266 = pkin(9) * t310;
t313 = t182 * t187;
t370 = t186 * t259 - t155 * t313 + t114 + (pkin(4) * t187 + t266) * qJD(1);
t276 = qJD(1) * qJD(3);
t250 = t187 * t276;
t270 = t183 * qJDD(1);
t122 = qJDD(4) + t250 + t270;
t118 = qJDD(5) + t122;
t180 = qJ(4) + qJ(5);
t172 = sin(t180);
t184 = sin(qJ(1));
t188 = cos(qJ(1));
t366 = -g(1) * t184 + g(2) * t188;
t218 = t366 * t187;
t340 = g(3) * t183;
t144 = t346 * t182;
t145 = t346 * t186;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t82 = -t144 * t181 + t145 * t185;
t369 = t82 * t118 - t172 * (-t218 - t340);
t292 = qJD(1) * t187;
t252 = t182 * t292;
t289 = qJD(3) * t186;
t125 = -t252 + t289;
t291 = qJD(3) * t182;
t126 = t186 * t292 + t291;
t224 = t125 * t181 + t185 * t126;
t68 = -t185 * t125 + t126 * t181;
t319 = t155 * t187;
t331 = qJD(3) * pkin(3);
t116 = -t319 - t331;
t78 = -pkin(4) * t125 + t116;
t22 = pkin(5) * t68 - qJ(6) * t224 + t78;
t368 = t22 * t68;
t367 = t68 * t78;
t268 = qJD(4) + qJD(5);
t280 = qJD(5) * t185;
t283 = qJD(4) * t186;
t305 = t185 * t186;
t315 = t181 * t182;
t333 = t181 * t258 - t185 * t283 - t186 * t280 + t268 * t315 - t293 * t305;
t129 = t181 * t186 + t182 * t185;
t107 = t129 * qJD(1);
t76 = t268 * t129;
t332 = t183 * t107 + t76;
t338 = t224 * t68;
t290 = qJD(3) * t183;
t257 = t182 * t290;
t282 = qJD(4) * t187;
t364 = t186 * t282 - t257;
t269 = t187 * qJDD(1);
t221 = qJD(3) * qJD(4) + t269;
t347 = t224 ^ 2;
t363 = -t68 ^ 2 + t347;
t164 = qJD(4) + t293;
t153 = qJD(5) + t164;
t240 = qJD(4) - t293;
t207 = qJD(3) * t240 + t269;
t275 = qJD(1) * qJD(4);
t249 = t187 * t275;
t225 = -qJDD(3) + t249;
t214 = t225 * t182;
t195 = t207 * t186 - t214;
t261 = t182 * t221 + t186 * t249;
t219 = t186 * qJDD(3) - t261;
t251 = t183 * t276;
t200 = t182 * t251 + t219;
t281 = qJD(5) * t181;
t18 = -t125 * t280 + t126 * t281 - t181 * t200 - t185 * t195;
t10 = t153 * t68 - t18;
t32 = pkin(5) * t224 + qJ(6) * t68;
t361 = -qJD(5) * t82 + t371 * t181 - t370 * t185;
t223 = -t144 * t185 - t145 * t181;
t360 = -qJD(5) * t223 + t370 * t181 + t371 * t185;
t135 = pkin(3) * t183 - pkin(8) * t187 + qJ(2);
t108 = t135 * qJD(1);
t134 = t183 * t155;
t115 = qJD(3) * pkin(8) + t134;
t59 = t182 * t108 + t186 * t115;
t49 = pkin(9) * t125 + t59;
t285 = qJD(4) * t182;
t239 = -t134 + (t258 + t285) * pkin(4);
t100 = t118 * qJ(6);
t136 = t153 * qJD(6);
t359 = t100 + t136;
t192 = qJD(1) ^ 2;
t357 = qJ(2) * t192 - t366;
t303 = t186 * t188;
t311 = t183 * t184;
t109 = -t182 * t311 + t303;
t306 = t184 * t186;
t309 = t183 * t188;
t111 = t182 * t309 + t306;
t356 = -g(1) * t109 - g(2) * t111;
t179 = t187 ^ 2;
t355 = qJDD(1) * t179;
t103 = t118 * pkin(5);
t354 = t103 - qJDD(6);
t123 = qJD(3) * t233 + qJD(2);
t66 = qJD(1) * t123 + qJDD(1) * t135;
t151 = t190 * qJDD(1) + qJDD(2);
t288 = qJD(3) * t187;
t80 = qJDD(3) * pkin(8) + t151 * t183 + t155 * t288;
t267 = -t108 * t283 - t182 * t66 - t186 * t80;
t215 = -t115 * t285 - t267;
t13 = pkin(9) * t200 + t215;
t58 = t186 * t108 - t115 * t182;
t48 = -pkin(9) * t126 + t58;
t42 = pkin(4) * t164 + t48;
t271 = t182 * qJDD(3);
t63 = t186 * t66;
t9 = -t182 * t80 + t63 - (t271 + (-t251 + t269) * t186) * pkin(9) + t122 * pkin(4) - t49 * qJD(4);
t247 = t181 * t13 - t185 * t9 + t49 * t280 + t42 * t281;
t318 = t172 * t187;
t173 = cos(t180);
t300 = t188 * t173;
t307 = t184 * t172;
t90 = t183 * t307 - t300;
t317 = t173 * t184;
t92 = t172 * t309 + t317;
t209 = g(1) * t90 - g(2) * t92 + g(3) * t318 - t247;
t196 = t22 * t224 - t209 - t354;
t353 = -t224 * t78 + t209;
t19 = qJD(5) * t224 + t181 * t195 - t185 * t200;
t352 = t153 * t224 - t19;
t348 = t183 * t268;
t128 = -t305 + t315;
t98 = t128 * t187;
t325 = qJD(3) * t98 + t129 * t348 + t107;
t97 = t128 * t183;
t351 = t97 * t118 + t153 * t325 + t187 * t18 + t224 * t290;
t121 = t186 * t135;
t312 = t182 * t190;
t246 = pkin(4) - t312;
t65 = -pkin(9) * t304 + t183 * t246 + t121;
t302 = t186 * t190;
t152 = t183 * t302;
t297 = t182 * t135 + t152;
t77 = -pkin(9) * t313 + t297;
t226 = t181 * t65 + t185 * t77;
t102 = t186 * t123;
t30 = t102 + (-t152 + (pkin(9) * t187 - t135) * t182) * qJD(4) + (t187 * t246 + t266) * qJD(3);
t287 = qJD(3) * t190;
t262 = t182 * t123 + t135 * t283 + t287 * t304;
t284 = qJD(4) * t183;
t33 = -pkin(9) * t364 - t284 * t312 + t262;
t349 = -qJD(5) * t226 - t181 * t33 + t185 * t30;
t345 = pkin(8) * t122;
t339 = g(3) * t187;
t337 = pkin(5) * t332 + qJ(6) * t333 - qJD(6) * t129 + t239;
t336 = -qJ(6) * t292 - t360;
t335 = pkin(5) * t292 - t361;
t328 = t185 * t49;
t17 = t181 * t42 + t328;
t330 = t153 * t17;
t329 = t181 * t49;
t96 = t129 * t187;
t326 = -qJD(3) * t96 + (qJD(1) + t348) * t128;
t21 = t185 * t48 - t329;
t324 = pkin(4) * t280 + qJD(6) - t21;
t323 = pkin(1) * qJDD(1);
t322 = qJDD(3) * pkin(3);
t321 = t122 * t183;
t320 = t126 * t164;
t316 = t173 * t187;
t314 = t182 * t122;
t301 = t187 * t346;
t16 = t185 * t42 - t329;
t299 = qJD(6) - t16;
t296 = t188 * pkin(1) + t184 * qJ(2);
t295 = t183 ^ 2 - t179;
t191 = qJD(3) ^ 2;
t294 = -t191 - t192;
t286 = qJD(4) * t125;
t279 = t116 * qJD(4);
t278 = t125 * qJD(3);
t277 = t126 * qJD(4);
t273 = qJDD(1) * qJ(2);
t272 = qJDD(3) * t183;
t265 = g(1) * t317;
t263 = 0.2e1 * qJD(1) * qJD(2);
t170 = pkin(4) * t186 + pkin(3);
t260 = pkin(4) * t182 + pkin(7);
t256 = t126 * t290;
t255 = t183 * t289;
t254 = t182 * t282;
t248 = t185 * t13 + t181 * t9 + t42 * t280 - t49 * t281;
t244 = -qJD(4) * t108 - t80;
t243 = t164 * t190 + t115;
t165 = pkin(4) * t313;
t124 = -t187 * t190 + t165;
t242 = qJDD(2) - t323;
t241 = qJD(1) + t284;
t20 = t181 * t48 + t328;
t238 = pkin(4) * t281 - t20;
t237 = -pkin(5) * t318 + qJ(6) * t316;
t236 = g(1) * t92 + g(2) * t90;
t91 = t172 * t188 + t173 * t311;
t93 = t183 * t300 - t307;
t235 = -g(1) * t93 - g(2) * t91;
t234 = g(2) * t187 * t300 + t118 * t223 + t173 * t340;
t232 = g(1) * t188 + g(2) * t184;
t227 = -t181 * t77 + t185 * t65;
t222 = t170 * t183 - t301;
t220 = pkin(5) * t173 + qJ(6) * t172 + t170;
t79 = -t187 * t151 + t155 * t290 - t322;
t217 = t181 * t30 + t185 * t33 + t65 * t280 - t77 * t281;
t216 = t164 * t283 + t314;
t83 = pkin(4) * t364 + t183 * t287;
t213 = t225 * t187;
t212 = 0.2e1 * qJ(2) * t276 + qJDD(3) * t190;
t211 = -t151 + t357;
t208 = g(1) * t91 - g(2) * t93 + g(3) * t316 - t248;
t206 = -qJD(4) * pkin(8) * t164 + t322 + t340 - t79;
t205 = t221 - t251;
t95 = t129 * t183;
t204 = -t95 * t118 + t326 * t153 - t187 * t19 + t68 * t290;
t203 = -g(1) * (-t90 * pkin(5) + t91 * qJ(6)) - g(2) * (t92 * pkin(5) - t93 * qJ(6));
t202 = -t232 + t263 + 0.2e1 * t273;
t198 = t153 * t16 + t208;
t197 = -t190 * t191 + t202;
t34 = -pkin(4) * t200 + t79;
t175 = t188 * qJ(2);
t171 = qJDD(3) * t187;
t169 = -pkin(4) * t185 - pkin(5);
t166 = pkin(4) * t181 + qJ(6);
t112 = -t182 * t184 + t183 * t303;
t110 = t182 * t188 + t183 * t306;
t64 = pkin(5) * t128 - qJ(6) * t129 - t170;
t47 = pkin(5) * t96 + qJ(6) * t98 + t124;
t41 = -t281 * t313 + (t268 * t304 - t257) * t185 + (-t254 - t255) * t181;
t39 = -t181 * t257 + t185 * t255 + t187 * t76;
t29 = -pkin(5) * t183 - t227;
t28 = qJ(6) * t183 + t226;
t26 = pkin(4) * t126 + t32;
t15 = qJ(6) * t153 + t17;
t14 = -pkin(5) * t153 + t299;
t6 = pkin(5) * t41 + qJ(6) * t39 + qJD(6) * t98 + t83;
t5 = -pkin(5) * t288 - t349;
t4 = qJ(6) * t288 + qJD(6) * t183 + t217;
t3 = t19 * pkin(5) + t18 * qJ(6) - qJD(6) * t224 + t34;
t2 = t247 - t354;
t1 = t248 + t359;
t7 = [qJDD(1), -t366, t232, qJDD(2) + t366 - 0.2e1 * t323, t202, -t242 * pkin(1) - g(1) * (-pkin(1) * t184 + t175) - g(2) * t296 + (t263 + t273) * qJ(2), -0.2e1 * t183 * t250 + t355, -0.2e1 * t183 * t269 + 0.2e1 * t295 * t276, -t183 * t191 + t171, -t187 * t191 - t272, 0, t183 * t197 + t187 * t212, -t183 * t212 + t187 * t197, -t126 * t254 + (-t182 * t213 + t205 * t304 - t256) * t186 (-t125 * t186 + t126 * t182) * t290 + ((t219 - t277) * t186 + (-t286 + t214 + (-t269 + (-qJD(4) + 0.2e1 * t293) * qJD(3)) * t186) * t182) * t187, t126 * t288 + (-t164 * t282 - t183 * t225) * t182 + ((t122 + t270) * t187 + (-t164 + t240) * t290) * t186 ((t164 + t293) * t291 + t219) * t183 + (-t216 + t278) * t187, t164 * t288 + t321 (-t135 * t285 + t102) * t164 + t121 * t122 - g(1) * t112 - g(2) * t110 + (-t190 * t278 + t63 - t243 * t283 + (-qJD(3) * t116 - t122 * t190 + t244) * t182) * t183 + (t190 * t219 + t79 * t182 + t186 * t279 + (t58 + (-t164 + t293) * t312) * qJD(3)) * t187, -t262 * t164 - t297 * t122 + g(1) * t111 - g(2) * t109 + (t243 * t285 + (-t116 * t186 + t126 * t190) * qJD(3) + t267) * t183 + (-t182 * t279 + t79 * t186 + (-t271 - (qJDD(1) * t186 - t182 * t275) * t187) * t190 + (-t240 * t302 - t59) * qJD(3)) * t187, t18 * t98 - t224 * t39, t18 * t96 + t19 * t98 - t224 * t41 + t39 * t68, -t118 * t98 - t153 * t39 - t18 * t183 + t224 * t288, -t118 * t96 - t153 * t41 - t183 * t19 - t288 * t68, t118 * t183 + t153 * t288, t227 * t118 + t124 * t19 + t153 * t349 + t16 * t288 - t247 * t183 + t34 * t96 + t78 * t41 + t83 * t68 + t235, -t118 * t226 - t124 * t18 - t153 * t217 - t17 * t288 - t183 * t248 + t224 * t83 - t34 * t98 - t78 * t39 + t236, -t118 * t29 - t14 * t288 - t153 * t5 - t183 * t2 + t19 * t47 + t22 * t41 + t3 * t96 + t6 * t68 + t235, -t1 * t96 - t14 * t39 - t15 * t41 - t18 * t29 + t187 * t232 - t19 * t28 - t2 * t98 + t224 * t5 - t4 * t68, t1 * t183 + t118 * t28 + t15 * t288 + t153 * t4 + t18 * t47 + t22 * t39 - t224 * t6 + t3 * t98 - t236, t1 * t28 + t15 * t4 + t3 * t47 + t22 * t6 + t2 * t29 + t14 * t5 - g(1) * (pkin(5) * t93 + qJ(6) * t92 + t175) - g(2) * (pkin(5) * t91 + qJ(6) * t90 + t296) + (-g(1) * t222 - g(2) * t260) * t188 + (-g(1) * (-pkin(1) - t260) - g(2) * t222) * t184; 0, 0, 0, qJDD(1), -t192, t242 - t357, 0, 0, 0, 0, 0, t294 * t183 + t171, t294 * t187 - t272, 0, 0, 0, 0, 0, t187 * t219 + (-t314 + (-t125 + t252) * qJD(3)) * t183 + (-t182 * t288 - t186 * t241) * t164, t256 + (t164 * t241 + t213) * t182 + (-t355 - t321 + (-t164 - t240) * t288) * t186, 0, 0, 0, 0, 0, t204, t351, t204, -t18 * t95 + t19 * t97 - t224 * t326 + t325 * t68, -t351, -t1 * t97 - t14 * t326 - t15 * t325 - t187 * t3 + t2 * t95 + t22 * t290 + t366; 0, 0, 0, 0, 0, 0, t187 * t192 * t183, -t295 * t192, t269, -t270, qJDD(3), -t187 * t211 + t340, t183 * t211 + t339, -t225 * t182 ^ 2 + (t182 * t207 + t320) * t186 (-t277 + (-t126 + t291) * t293 - t261) * t182 + (t286 + 0.2e1 * t271 + t221 * t186 + (-t254 + (t125 - t289) * t183) * qJD(1)) * t186 (-t126 * t187 + t164 * t310) * qJD(1) + t216, -t164 * t285 + t186 * t122 + (-t182 * t183 * t164 - t125 * t187) * qJD(1), -t164 * t292, -pkin(3) * t261 - t114 * t164 - t58 * t292 + t125 * t134 + (t164 * t319 - t345 + t279 + (t116 + t331) * t293) * t182 + (t218 + t206) * t186, t298 * t164 + t59 * t292 - t126 * t134 + (-pkin(3) * t205 + t116 * t164 - t345) * t186 + ((pkin(3) * t275 - t366) * t187 - t206) * t182, -t18 * t129 - t224 * t333, t128 * t18 - t129 * t19 - t224 * t332 + t333 * t68, t129 * t118 - t153 * t333 - t224 * t292, -t128 * t118 - t153 * t332 + t292 * t68, -t153 * t292, t34 * t128 - t170 * t19 + t332 * t78 + t239 * t68 + (-qJD(1) * t16 - t265) * t187 + t361 * t153 + t234, t34 * t129 + t360 * t153 + t17 * t292 + t170 * t18 + t239 * t224 - t333 * t78 - t369, t128 * t3 + t19 * t64 + t337 * t68 + t332 * t22 + (qJD(1) * t14 - t265) * t187 - t335 * t153 + t234, -t1 * t128 + t129 * t2 - t14 * t333 - t15 * t332 + t18 * t223 + t183 * t366 - t19 * t82 + t224 * t335 - t336 * t68 - t339, -t129 * t3 - t15 * t292 + t336 * t153 + t18 * t64 + t333 * t22 - t337 * t224 + t369, -g(3) * t301 + t1 * t82 + t335 * t14 + t336 * t15 - t2 * t223 + t337 * t22 + t220 * t340 + t3 * t64 + t366 * (t183 * t346 + t187 * t220); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t125, -t125 ^ 2 + t126 ^ 2, -t125 * t164 + t195, t200 + t320, t122, -t115 * t283 - t116 * t126 + t164 * t59 + t63 + (t244 + t339) * t182 + t356, g(1) * t110 - g(2) * t112 + g(3) * t304 - t116 * t125 + t164 * t58 - t215, t338, t363, t10, t352, t118, t153 * t20 + (t118 * t185 - t126 * t68 - t153 * t281) * pkin(4) + t353, t153 * t21 + t367 + (-t118 * t181 - t126 * t224 - t153 * t280) * pkin(4) + t208, -t118 * t169 - t153 * t238 - t26 * t68 - t196, -t166 * t19 - t169 * t18 + (t15 + t238) * t224 + (t14 - t324) * t68, t118 * t166 + t153 * t324 + t224 * t26 - t208 + t359 - t368, t1 * t166 + t2 * t169 - t22 * t26 - t14 * t20 - g(3) * (-t165 + t237) + t324 * t15 + (t14 * t281 + t356) * pkin(4) + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t363, t10, t352, t118, t330 + t353, t198 + t367, -t32 * t68 + t103 - t196 + t330, pkin(5) * t18 - qJ(6) * t19 + (t15 - t17) * t224 + (t14 - t299) * t68, t224 * t32 + 0.2e1 * t100 + 0.2e1 * t136 - t198 - t368, -t2 * pkin(5) - g(3) * t237 + t1 * qJ(6) - t14 * t17 + t15 * t299 - t22 * t32 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118 + t338, t10, -t153 ^ 2 - t347, -t15 * t153 + t196;];
tau_reg  = t7;
