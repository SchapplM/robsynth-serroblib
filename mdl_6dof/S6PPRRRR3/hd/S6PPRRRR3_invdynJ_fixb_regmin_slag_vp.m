% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:44
% EndTime: 2019-03-08 19:12:00
% DurationCPUTime: 7.03s
% Computational Cost: add. (6810->521), mult. (18210->800), div. (0->0), fcn. (17674->18), ass. (0->262)
t199 = sin(pkin(6));
t195 = sin(pkin(14));
t208 = sin(qJ(3));
t212 = cos(qJ(3));
t200 = cos(pkin(14));
t203 = cos(pkin(7));
t324 = t200 * t203;
t249 = -t195 * t208 + t212 * t324;
t364 = t199 * t249;
t204 = cos(pkin(6));
t185 = qJD(1) * t204 + qJD(2);
t198 = sin(pkin(7));
t326 = t198 * t212;
t167 = t185 * t326;
t111 = qJD(1) * t364 + t167;
t327 = t198 * t208;
t248 = t195 * t212 + t208 * t324;
t360 = t248 * t199;
t112 = -qJD(1) * t360 - t185 * t327;
t211 = cos(qJ(4));
t202 = cos(pkin(8));
t207 = sin(qJ(4));
t322 = t202 * t207;
t197 = sin(pkin(8));
t331 = t197 * t207;
t186 = pkin(10) * t331;
t321 = t202 * t211;
t359 = pkin(3) * t321 - t186;
t366 = t359 * qJD(4) - t111 * t211 - t112 * t322;
t184 = t204 * qJDD(1) + qJDD(2);
t107 = qJD(3) * pkin(3) + t111;
t344 = t107 * t202;
t357 = qJD(1) * t249;
t365 = -qJD(4) * t344 - qJDD(3) * t197 * pkin(10) - (qJD(3) * t185 * t212 + t184 * t208) * t198 - (qJD(3) * t357 + qJDD(1) * t248) * t199;
t119 = t204 * t327 + t360;
t311 = qJD(3) * t197;
t287 = t211 * t311;
t363 = qJD(5) - t287;
t210 = cos(qJ(5));
t206 = sin(qJ(5));
t288 = t207 * t311;
t267 = t206 * t288;
t310 = qJD(3) * t202;
t273 = qJD(4) + t310;
t144 = -t210 * t273 + t267;
t138 = qJD(6) + t144;
t299 = qJDD(3) * t202;
t183 = qJDD(4) + t299;
t254 = qJD(5) * t273;
t308 = qJD(4) * t211;
t283 = t206 * t308;
t297 = t207 * qJDD(3);
t304 = qJD(5) * t210;
t88 = -t210 * t183 + t197 * (qJD(3) * (t207 * t304 + t283) + t206 * t297) + t206 * t254;
t329 = t197 * t211;
t313 = pkin(3) * t322 + pkin(10) * t329;
t347 = t313 * qJD(4) - t111 * t207 + t112 * t321;
t151 = pkin(11) * t202 + t313;
t263 = -pkin(4) * t211 - pkin(11) * t207;
t152 = (-pkin(3) + t263) * t197;
t252 = (pkin(4) * t207 - pkin(11) * t211) * t197;
t155 = qJD(4) * t252;
t306 = qJD(5) * t206;
t332 = t197 * t206;
t362 = t112 * t332 - t151 * t306 + t152 * t304 + t206 * t155 + t210 * t366;
t314 = t210 * t151 + t206 * t152;
t330 = t197 * t210;
t361 = qJD(5) * t314 - t112 * t330 - t210 * t155 + t206 * t366;
t201 = cos(pkin(13));
t196 = sin(pkin(13));
t323 = t201 * t204;
t251 = -t195 * t196 + t200 * t323;
t328 = t198 * t199;
t358 = -t201 * t328 + t251 * t203;
t106 = pkin(10) * t311 - t112;
t293 = t200 * t328;
t137 = -qJDD(1) * t293 + t184 * t203;
t141 = -qJD(1) * t293 + t185 * t203;
t284 = t197 * t308;
t309 = qJD(4) * t207;
t235 = qJDD(1) * t364 + t184 * t326;
t301 = t112 * qJD(3);
t73 = qJDD(3) * pkin(3) + t235 + t301;
t226 = t106 * t309 - t137 * t331 - t141 * t284 + t211 * t365 - t73 * t322;
t12 = pkin(11) * t183 - t226;
t49 = t211 * t106 + t107 * t322 + t141 * t331;
t44 = pkin(11) * t273 + t49;
t132 = t202 * t141;
t69 = t132 + (qJD(3) * t263 - t107) * t197;
t23 = t206 * t69 + t210 * t44;
t131 = t202 * t137;
t300 = qJD(3) * qJD(4);
t281 = t207 * t300;
t298 = qJDD(3) * t211;
t237 = t281 - t298;
t280 = t211 * t300;
t238 = t280 + t297;
t39 = t131 + (pkin(4) * t237 - pkin(11) * t238 - t73) * t197;
t222 = -qJD(5) * t23 - t206 * t12 + t210 * t39;
t48 = -t207 * t106 + t211 * (t141 * t197 + t344);
t146 = t206 * t273 + t210 * t288;
t325 = t199 * t203;
t220 = -t198 * t251 - t201 * t325;
t215 = t220 * t197;
t158 = t195 * t323 + t196 * t200;
t89 = -t158 * t208 + t212 * t358;
t90 = t158 * t212 + t208 * t358;
t30 = t90 * t211 + (t89 * t202 + t215) * t207;
t333 = t196 * t204;
t250 = -t195 * t201 - t200 * t333;
t219 = t196 * t325 - t198 * t250;
t214 = t219 * t197;
t159 = -t195 * t333 + t200 * t201;
t218 = t196 * t328 + t203 * t250;
t91 = -t159 * t208 + t212 * t218;
t92 = t159 * t212 + t208 * t218;
t32 = t92 * t211 + (t91 * t202 + t214) * t207;
t118 = t204 * t326 + t364;
t247 = t203 * t204 - t293;
t228 = t247 * t197;
t221 = t118 * t202 + t228;
t61 = t119 * t211 + t207 * t221;
t93 = -t118 * t197 + t202 * t247;
t33 = t206 * t61 - t93 * t210;
t62 = -t89 * t197 + t202 * t220;
t63 = -t91 * t197 + t202 * t219;
t243 = g(1) * (-t206 * t32 + t210 * t63) + g(2) * (-t206 * t30 + t210 * t62) - g(3) * t33;
t182 = t197 * t298;
t153 = t197 * t281 + qJDD(5) - t182;
t4 = -t153 * pkin(5) - t222;
t356 = t138 * (pkin(5) * t146 + pkin(12) * t138) + t243 + t4;
t180 = -pkin(5) * t210 - pkin(12) * t206 - pkin(4);
t84 = qJDD(6) + t88;
t355 = t138 * (-t49 + t363 * (pkin(5) * t206 - pkin(12) * t210)) + t180 * t84;
t205 = sin(qJ(6));
t209 = cos(qJ(6));
t110 = t146 * t209 + t205 * t363;
t266 = t210 * t284;
t279 = t197 * t297;
t87 = qJD(3) * t266 - qJD(5) * t267 + t206 * t183 + (t254 + t279) * t210;
t51 = qJD(6) * t110 - t209 * t153 + t205 * t87;
t21 = pkin(12) * t363 + t23;
t43 = -pkin(4) * t273 - t48;
t28 = t144 * pkin(5) - t146 * pkin(12) + t43;
t259 = t205 * t21 - t209 * t28;
t244 = t210 * t12 + t206 * t39 + t69 * t304 - t306 * t44;
t3 = pkin(12) * t153 + t244;
t285 = t197 * t309;
t253 = -t106 * t308 + t137 * t329 - t141 * t285 + t207 * t365 + t73 * t321;
t13 = -pkin(4) * t183 - t253;
t6 = pkin(5) * t88 - pkin(12) * t87 + t13;
t1 = -t259 * qJD(6) + t205 * t6 + t209 * t3;
t213 = qJD(3) ^ 2;
t353 = -pkin(5) * t285 + t361;
t352 = pkin(11) * qJD(5);
t350 = t205 * t84;
t349 = t209 * t84;
t302 = qJD(6) * t209;
t303 = qJD(6) * t205;
t50 = -t146 * t303 + t205 * t153 + t209 * t87 + t302 * t363;
t348 = t50 * t205;
t154 = qJD(3) * t252;
t345 = t206 * t154 + t210 * t48;
t108 = t146 * t205 - t209 * t363;
t343 = t108 * t138;
t342 = t110 * t138;
t117 = t119 * qJD(3);
t191 = t197 ^ 2;
t341 = t117 * t191;
t340 = t119 * t207;
t339 = t144 * t363;
t338 = t146 * t363;
t160 = -t197 * t326 + t202 * t203;
t337 = t160 * t197;
t336 = t363 * t206;
t335 = t191 * t213;
t320 = t207 * t208;
t319 = t207 * t212;
t318 = t208 * t211;
t317 = t208 * t213;
t316 = t210 * t211;
t315 = t211 * t212;
t193 = t207 ^ 2;
t312 = -t211 ^ 2 + t193;
t307 = qJD(5) * t205;
t305 = qJD(5) * t209;
t294 = t205 * t329;
t282 = t197 * t202 * t213;
t274 = t138 * t209;
t272 = qJD(4) + 0.2e1 * t310;
t271 = t183 + t299;
t270 = t191 * t198 * t317;
t268 = t311 * t327;
t150 = t186 + (-pkin(3) * t211 - pkin(4)) * t202;
t161 = -t210 * t202 + t206 * t331;
t162 = t202 * t206 + t207 * t330;
t96 = t161 * pkin(5) - t162 * pkin(12) + t150;
t265 = -pkin(12) * t285 - qJD(6) * t96 - t362;
t126 = -qJD(5) * t161 + t266;
t127 = qJD(5) * t162 + t197 * t283;
t98 = -pkin(12) * t329 + t314;
t264 = -t127 * pkin(5) + t126 * pkin(12) + qJD(6) * t98 - t347;
t135 = (t205 * t207 + t209 * t316) * t311;
t260 = t209 * t304 - t135;
t8 = t205 * t28 + t209 * t21;
t34 = t206 * t93 + t210 * t61;
t60 = -t118 * t321 - t211 * t228 + t340;
t15 = t205 * t60 + t209 * t34;
t14 = -t205 * t34 + t209 * t60;
t22 = -t206 * t44 + t210 * t69;
t246 = t202 * t315 - t320;
t120 = -t198 * t246 - t203 * t329;
t245 = t202 * t319 + t318;
t121 = t198 * t245 + t203 * t331;
t95 = t121 * t210 + t160 * t206;
t258 = t120 * t209 - t205 * t95;
t257 = t120 * t205 + t209 * t95;
t94 = t121 * t206 - t160 * t210;
t255 = -t151 * t206 + t152 * t210;
t128 = t162 * t205 + t209 * t329;
t29 = t207 * t90 - t211 * t215 - t321 * t89;
t31 = t207 * t92 - t211 * t214 - t321 * t91;
t242 = g(1) * t31 + g(2) * t29 + g(3) * t60;
t241 = -g(1) * t32 - g(2) * t30 - g(3) * t61;
t53 = t211 * t89 - t322 * t90;
t55 = t211 * t91 - t322 * t92;
t78 = t118 * t211 - t119 * t322;
t240 = g(1) * t55 + g(2) * t53 + g(3) * t78;
t239 = g(1) * t92 + g(2) * t90 + g(3) * t119;
t233 = -t13 + t242;
t20 = -pkin(5) * t363 - t22;
t225 = -pkin(12) * t84 + (t20 + t22) * t138;
t224 = -pkin(11) * t153 + t363 * t43;
t223 = pkin(11) * qJD(6) * t138 - t242;
t2 = -qJD(6) * t8 - t205 * t3 + t209 * t6;
t217 = (pkin(12) * t288 - qJD(6) * t180 + t345) * t138 + t241;
t134 = t205 * t210 * t287 - t209 * t288;
t129 = t162 * t209 - t294;
t116 = t118 * qJD(3);
t97 = pkin(5) * t329 - t255;
t86 = t203 * t285 + (t245 * qJD(4) + (t202 * t318 + t319) * qJD(3)) * t198;
t85 = t203 * t284 + (t246 * qJD(4) + (-t202 * t320 + t315) * qJD(3)) * t198;
t79 = -t197 * t107 + t132;
t77 = t118 * t207 + t119 * t321;
t76 = -qJD(6) * t294 + t126 * t205 + t162 * t302 - t209 * t285;
t75 = -qJD(6) * t128 + t209 * t126 + t205 * t285;
t57 = t119 * t332 + t210 * t78;
t56 = -t197 * t73 + t131;
t54 = t207 * t91 + t321 * t92;
t52 = t207 * t89 + t321 * t90;
t46 = qJD(5) * t95 + t206 * t85 - t210 * t268;
t45 = -qJD(5) * t94 + t206 * t268 + t210 * t85;
t35 = -pkin(5) * t288 - t154 * t210 + t206 * t48;
t27 = -t117 * t322 + t116 * t211 + (t211 * t221 - t340) * qJD(4);
t26 = qJD(4) * t61 + t116 * t207 + t117 * t321;
t25 = t210 * t55 + t332 * t92;
t24 = t210 * t53 + t332 * t90;
t19 = t206 * t63 + t210 * t32;
t17 = t206 * t62 + t210 * t30;
t10 = -qJD(5) * t33 + t117 * t332 + t27 * t210;
t9 = qJD(5) * t34 - t117 * t330 + t27 * t206;
t5 = [qJDD(1) - g(3), t184 * t204 - g(3) + (t195 ^ 2 + t200 ^ 2) * t199 ^ 2 * qJDD(1), 0, -qJD(3) * t117 + qJDD(3) * t118, -qJD(3) * t116 - qJDD(3) * t119, 0, 0, 0, 0, 0, -t93 * t182 - t26 * qJD(4) - t60 * t183 + (-t202 * t26 - t211 * t341 + t285 * t93) * qJD(3), t93 * t279 - t27 * qJD(4) - t61 * t183 + (-t202 * t27 + t207 * t341 + t284 * t93) * qJD(3), 0, 0, 0, 0, 0, t144 * t26 - t153 * t33 - t363 * t9 + t60 * t88, -t10 * t363 + t146 * t26 - t153 * t34 + t60 * t87, 0, 0, 0, 0, 0 (-qJD(6) * t15 - t10 * t205 + t26 * t209) * t138 + t14 * t84 + t9 * t108 + t33 * t51 -(qJD(6) * t14 + t10 * t209 + t26 * t205) * t138 - t15 * t84 + t9 * t110 + t33 * t50; 0, -g(3) * t204 + (-g(1) * t196 + g(2) * t201) * t199 + t184, 0 (qJDD(3) * t212 - t317) * t198 (-qJDD(3) * t208 - t212 * t213) * t198, 0, 0, 0, 0, 0, -t120 * t183 - t211 * t270 + t237 * t337 - t273 * t86, -t121 * t183 + t207 * t270 + t238 * t337 - t273 * t85, 0, 0, 0, 0, 0, t120 * t88 + t144 * t86 - t153 * t94 - t363 * t46, t120 * t87 + t146 * t86 - t153 * t95 - t363 * t45, 0, 0, 0, 0, 0 (-qJD(6) * t257 - t205 * t45 + t209 * t86) * t138 + t258 * t84 + t46 * t108 + t94 * t51 -(qJD(6) * t258 + t205 * t86 + t209 * t45) * t138 - t257 * t84 + t46 * t110 + t94 * t50; 0, 0, qJDD(3), -g(1) * t91 - g(2) * t89 - g(3) * t118 + t235, -t184 * t327 - qJDD(1) * t360 + (-t199 * t357 + t111 - t167) * qJD(3) + t239 (qJDD(3) * t193 + 0.2e1 * t207 * t280) * t191, 0.2e1 * (t211 * t297 - t300 * t312) * t191 (t207 * t271 + t272 * t308) * t197 (t211 * t271 - t272 * t309) * t197, t183 * t202, t359 * t183 + t253 * t202 + (-t211 * t56 + t309 * t79) * t197 + (-pkin(3) * t237 - t211 * t301) * t191 - t240 - t347 * t273, -t313 * t183 + t226 * t202 + g(1) * t54 + g(2) * t52 + g(3) * t77 + (t207 * t56 + t308 * t79) * t197 + (-pkin(3) * t238 + t207 * t301) * t191 - t366 * t273, t126 * t146 + t162 * t87, -t126 * t144 - t127 * t146 - t161 * t87 - t162 * t88, t126 * t363 + t162 * t153 + (t146 * t309 - t211 * t87) * t197, -t127 * t363 - t161 * t153 + (-t144 * t309 + t211 * t88) * t197 (-t153 * t211 + t309 * t363) * t197, t255 * t153 + t150 * t88 + t13 * t161 + t43 * t127 - g(1) * t25 - g(2) * t24 - g(3) * t57 + (-t211 * t222 + t22 * t309) * t197 - t361 * t363 + t347 * t144, -t314 * t153 + t150 * t87 + t13 * t162 + t43 * t126 + t240 * t206 + (-t210 * t239 + t211 * t244 - t23 * t309) * t197 - t362 * t363 + t347 * t146, t110 * t75 + t129 * t50, -t108 * t75 - t110 * t76 - t128 * t50 - t129 * t51, t110 * t127 + t129 * t84 + t138 * t75 + t161 * t50, -t108 * t127 - t128 * t84 - t138 * t76 - t161 * t51, t127 * t138 + t161 * t84 (-t205 * t98 + t209 * t96) * t84 + t2 * t161 - t259 * t127 + t97 * t51 + t4 * t128 + t20 * t76 - g(1) * (t205 * t54 + t209 * t25) - g(2) * (t205 * t52 + t209 * t24) - g(3) * (t205 * t77 + t209 * t57) + (t205 * t265 - t209 * t264) * t138 + t353 * t108 -(t205 * t96 + t209 * t98) * t84 - t1 * t161 - t8 * t127 + t97 * t50 + t4 * t129 + t20 * t75 - g(1) * (-t205 * t25 + t209 * t54) - g(2) * (-t205 * t24 + t209 * t52) - g(3) * (-t205 * t57 + t209 * t77) + (t205 * t264 + t209 * t265) * t138 + t353 * t110; 0, 0, 0, 0, 0, -t211 * t207 * t335, t312 * t335, -t211 * t282 + t279, t207 * t282 + t182, t183, t273 * t49 - t288 * t79 + t242 + t253, t273 * t48 - t287 * t79 + t226 - t241, t87 * t206 + t210 * t338 (t87 - t339) * t210 + (-t338 - t88) * t206, t363 * t304 + t206 * t153 + (-t146 * t207 - t316 * t363) * t311, -t363 * t306 + t210 * t153 + (t144 * t207 + t211 * t336) * t311, -t363 * t288, -t22 * t288 - pkin(4) * t88 - t49 * t144 + (t363 * t48 + t224) * t206 + (-(t154 + t352) * t363 + t233) * t210, -pkin(4) * t87 + t345 * t363 + t23 * t288 - t49 * t146 + t224 * t210 + (t352 * t363 - t233) * t206, t50 * t209 * t206 + (-t206 * t303 + t260) * t110, t135 * t108 + t110 * t134 + (-t108 * t209 - t110 * t205) * t304 + (-t348 - t209 * t51 + (t108 * t205 - t110 * t209) * qJD(6)) * t206, -t50 * t210 + t260 * t138 + (t110 * t363 - t138 * t303 + t349) * t206, t51 * t210 + (-t205 * t304 + t134) * t138 + (-t108 * t363 - t138 * t302 - t350) * t206, t138 * t336 - t84 * t210, -t35 * t108 - t20 * t134 + t355 * t209 + t217 * t205 + (t20 * t307 - t2 + (qJD(5) * t108 - t350) * pkin(11) - t223 * t209) * t210 + (t20 * t302 + t4 * t205 - t363 * t259 + (t138 * t307 + t51) * pkin(11)) * t206, -t35 * t110 - t20 * t135 - t355 * t205 + t217 * t209 + (t20 * t305 + t1 + (qJD(5) * t110 - t349) * pkin(11) + t223 * t205) * t210 + (-t20 * t303 + t4 * t209 - t363 * t8 + (t138 * t305 + t50) * pkin(11)) * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146 * t144, -t144 ^ 2 + t146 ^ 2, t87 + t339, -t88 + t338, t153, -t43 * t146 + t23 * t363 + t222 - t243, g(1) * t19 + g(2) * t17 + g(3) * t34 + t144 * t43 + t22 * t363 - t244, t110 * t274 + t348 (t50 - t343) * t209 + (-t51 - t342) * t205, -t110 * t146 + t138 * t274 + t350, -t138 ^ 2 * t205 + t108 * t146 + t349, -t138 * t146, -pkin(5) * t51 - t23 * t108 + t146 * t259 + t225 * t205 - t209 * t356, -pkin(5) * t50 - t23 * t110 + t8 * t146 + t205 * t356 + t225 * t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t108, -t108 ^ 2 + t110 ^ 2, t50 + t343, t342 - t51, t84, t8 * t138 - t20 * t110 - g(1) * (-t19 * t205 + t209 * t31) - g(2) * (-t17 * t205 + t209 * t29) - g(3) * t14 + t2, -t259 * t138 + t20 * t108 - g(1) * (-t19 * t209 - t205 * t31) - g(2) * (-t17 * t209 - t205 * t29) + g(3) * t15 - t1;];
tau_reg  = t5;
