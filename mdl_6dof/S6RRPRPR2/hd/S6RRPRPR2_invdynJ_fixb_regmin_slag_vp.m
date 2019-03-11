% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:59
% EndTime: 2019-03-09 10:14:09
% DurationCPUTime: 4.45s
% Computational Cost: add. (6785->431), mult. (16144->519), div. (0->0), fcn. (12510->14), ass. (0->243)
t189 = qJD(2) + qJD(4);
t195 = sin(qJ(6));
t199 = cos(qJ(6));
t197 = sin(qJ(2));
t200 = cos(qJ(2));
t294 = sin(pkin(10));
t295 = cos(pkin(10));
t228 = t294 * t197 - t295 * t200;
t141 = t228 * qJD(1);
t227 = -t295 * t197 - t294 * t200;
t142 = t227 * qJD(1);
t196 = sin(qJ(4));
t323 = cos(qJ(4));
t96 = -t323 * t141 + t142 * t196;
t78 = t189 * t195 + t199 * t96;
t312 = t96 * t78;
t207 = qJD(2) * t142 - qJDD(1) * t228;
t274 = qJD(1) * qJD(2);
t210 = t227 * qJDD(1) + t228 * t274;
t268 = qJD(4) * t323;
t277 = qJD(4) * t196;
t236 = -t141 * t268 + t142 * t277 + t196 * t207 - t210 * t323;
t42 = -qJDD(6) - t236;
t304 = t195 * t42;
t239 = -t196 * t141 - t323 * t142;
t332 = qJD(6) + t239;
t349 = t199 * t332;
t355 = -t332 * t349 + t304;
t360 = t355 + t312;
t305 = t189 * t96;
t223 = t236 - t305;
t350 = t195 * t332;
t36 = t199 * t42;
t242 = -t332 * t350 - t36;
t80 = t189 * t199 - t195 * t96;
t315 = t80 * t96;
t359 = t242 - t315;
t313 = t96 ^ 2;
t341 = t239 ^ 2;
t358 = -t313 + t341;
t187 = qJDD(2) + qJDD(4);
t275 = qJD(6) * t199;
t276 = qJD(6) * t195;
t208 = -t196 * t210 - t323 * t207;
t44 = t239 * qJD(4) + t208;
t20 = t199 * t187 - t189 * t276 + t195 * t44 - t275 * t96;
t19 = t20 * t199;
t357 = -t350 * t80 + t19;
t256 = t195 * t187 - t199 * t44;
t21 = t80 * qJD(6) + t256;
t353 = t332 * t78;
t356 = -t349 * t80 + (-t20 + t353) * t195 - t199 * t21;
t325 = t96 * pkin(5);
t311 = t239 * t96;
t310 = qJ(3) + pkin(7);
t166 = t310 * t200;
t155 = qJD(1) * t166;
t144 = t294 * t155;
t165 = t310 * t197;
t154 = qJD(1) * t165;
t309 = qJD(2) * pkin(2);
t148 = -t154 + t309;
t104 = t295 * t148 - t144;
t319 = t142 * pkin(8);
t75 = qJD(2) * pkin(3) + t104 + t319;
t261 = t295 * t155;
t105 = t294 * t148 + t261;
t320 = t141 * pkin(8);
t77 = t105 - t320;
t49 = t196 * t75 + t323 * t77;
t37 = -qJ(5) * t189 - t49;
t24 = -t37 + t325;
t354 = t24 * t332;
t352 = t332 * t96;
t300 = t96 * qJ(5);
t270 = pkin(2) * t294;
t172 = t196 * t270;
t269 = t295 * pkin(2);
t178 = t269 + pkin(3);
t109 = t294 * t154 - t261;
t81 = t109 + t320;
t110 = -t295 * t154 - t144;
t82 = t110 + t319;
t351 = qJD(4) * t172 - t178 * t268 + t196 * t81 + t323 * t82;
t48 = t196 * t77 - t323 * t75;
t283 = qJD(5) + t48;
t206 = -t141 * t277 - t142 * t268 + t208;
t301 = t239 * t189;
t348 = -t206 + t301;
t185 = t200 * pkin(2);
t191 = qJ(2) + pkin(10);
t184 = qJ(4) + t191;
t176 = sin(t184);
t177 = cos(t184);
t280 = t177 * pkin(4) + t176 * qJ(5);
t347 = pkin(3) * cos(t191) + t185 + t280;
t339 = pkin(5) * t239;
t284 = t339 + t283;
t326 = pkin(4) + pkin(9);
t22 = -t326 * t189 + t284;
t179 = t185 + pkin(1);
t158 = -qJD(1) * t179 + qJD(3);
t111 = t141 * pkin(3) + t158;
t214 = -qJ(5) * t239 + t111;
t27 = -t326 * t96 + t214;
t11 = t195 * t22 + t199 * t27;
t262 = qJD(2) * t310;
t139 = -t197 * qJD(3) - t200 * t262;
t103 = qJDD(2) * pkin(2) + t139 * qJD(1) - qJDD(1) * t165;
t138 = t200 * qJD(3) - t197 * t262;
t108 = t138 * qJD(1) + qJDD(1) * t166;
t59 = t295 * t103 - t294 * t108;
t47 = qJDD(2) * pkin(3) + pkin(8) * t210 + t59;
t60 = t294 * t103 + t295 * t108;
t50 = pkin(8) * t207 + t60;
t257 = -t196 * t47 - t75 * t268 + t77 * t277 - t323 * t50;
t180 = t187 * qJ(5);
t334 = -t189 * qJD(5) - t180;
t7 = t257 + t334;
t5 = -pkin(5) * t44 - t7;
t346 = t11 * t96 + t5 * t199;
t10 = -t195 * t27 + t199 * t22;
t345 = -t10 * t96 + t5 * t195 + t24 * t275;
t170 = g(3) * t176;
t201 = cos(qJ(1));
t289 = t177 * t201;
t198 = sin(qJ(1));
t290 = t177 * t198;
t237 = -g(1) * t289 - g(2) * t290 - t170 - t257;
t51 = -pkin(4) * t96 + t214;
t344 = t51 * t96 + t237 - t334;
t343 = -t111 * t96 - t237;
t340 = pkin(4) * t239;
t338 = t239 * t24;
t298 = -qJD(5) + t351;
t226 = t196 * t178 + t323 * t270;
t296 = t226 * qJD(4) - t196 * t82 + t323 * t81;
t337 = t239 * t326;
t336 = t296 * t189;
t335 = -qJD(6) + t332;
t333 = g(1) * t198 - g(2) * t201;
t171 = g(3) * t177;
t258 = t196 * t50 + t77 * t268 + t75 * t277 - t323 * t47;
t291 = t176 * t201;
t292 = t176 * t198;
t238 = -g(1) * t291 - g(2) * t292 + t171 + t258;
t216 = t51 * t239 + qJDD(5) + t238;
t248 = t323 * t178 - t172;
t136 = -pkin(4) - t248;
t128 = -pkin(9) + t136;
t331 = (t296 - t325) * t332 - t128 * t42;
t330 = -t111 * t239 - t238;
t218 = qJD(2) * t228;
t86 = -t294 * t138 + t295 * t139;
t67 = pkin(8) * t218 + t86;
t217 = qJD(2) * t227;
t87 = t295 * t138 + t294 * t139;
t68 = pkin(8) * t217 + t87;
t112 = -t295 * t165 - t294 * t166;
t88 = pkin(8) * t227 + t112;
t113 = -t294 * t165 + t295 * t166;
t89 = -pkin(8) * t228 + t113;
t15 = -t196 * t67 - t88 * t268 + t89 * t277 - t323 * t68;
t57 = t196 * t88 + t323 * t89;
t329 = t15 * t189 - t176 * t333 - t57 * t187;
t16 = t57 * qJD(4) + t196 * t68 - t323 * t67;
t56 = t196 * t89 - t323 * t88;
t328 = -t16 * t189 + t177 * t333 - t56 * t187;
t321 = g(3) * t200;
t318 = t187 * pkin(4);
t317 = t197 * pkin(2);
t215 = t323 * t228;
t106 = -t196 * t227 + t215;
t219 = t196 * t228;
t107 = -t227 * t323 - t219;
t118 = pkin(3) * t228 - t179;
t213 = -t107 * qJ(5) + t118;
t31 = t326 * t106 + t213;
t316 = t31 * t42;
t308 = t106 * t24;
t306 = t189 * t49;
t303 = t195 * t80;
t299 = -t298 + t339;
t293 = qJDD(1) * pkin(1);
t288 = t195 * t198;
t287 = t195 * t201;
t286 = t198 * t199;
t285 = t199 * t201;
t192 = t197 ^ 2;
t278 = -t200 ^ 2 + t192;
t273 = t200 * qJDD(1);
t267 = t197 * t274;
t272 = pkin(2) * t267 + qJDD(3);
t182 = t197 * t309;
t114 = -t142 * pkin(3) + qJD(1) * t317;
t266 = -pkin(4) * t176 - pkin(3) * sin(t191) - t317;
t76 = -pkin(2) * t273 - t207 * pkin(3) + t272 - t293;
t9 = t44 * pkin(4) - qJ(5) * t236 - qJD(5) * t239 + t76;
t6 = t44 * pkin(9) + t9;
t264 = qJD(6) * t22 + t6;
t247 = qJDD(5) + t258;
t4 = pkin(5) * t236 - t326 * t187 + t247;
t263 = -qJD(6) * t27 + t4;
t250 = g(1) * t201 + g(2) * t198;
t246 = -t326 * t42 + (t49 + t325) * t332;
t243 = t114 - t300;
t241 = pkin(1) + t347;
t240 = -0.2e1 * pkin(1) * t274 - pkin(7) * qJDD(2);
t32 = t107 * pkin(5) + t56;
t62 = -qJD(4) * t219 - t196 * t218 - t323 * t217 - t227 * t268;
t234 = t5 * t106 + t24 * t62 + t32 * t42;
t230 = -qJDD(1) * t179 + t272;
t224 = -t250 * t177 - t170;
t222 = t236 + t305;
t203 = qJD(2) ^ 2;
t221 = -pkin(7) * t203 + 0.2e1 * t293 + t333;
t204 = qJD(1) ^ 2;
t220 = pkin(1) * t204 - pkin(7) * qJDD(1) + t250;
t115 = -pkin(3) * t217 + t182;
t212 = (-qJD(6) * t128 + t243 + t337) * t332 + t224;
t211 = (qJD(6) * t326 - t300 + t337) * t332 + t224;
t61 = t189 * t215 - t196 * t217 - t227 * t277;
t17 = t62 * pkin(4) + t61 * qJ(5) - t107 * qJD(5) + t115;
t205 = -t206 - t301;
t188 = -pkin(8) - t310;
t160 = qJ(5) * t289;
t159 = qJ(5) * t290;
t135 = qJ(5) + t226;
t134 = -t176 * t288 + t285;
t133 = t176 * t286 + t287;
t132 = t176 * t287 + t286;
t131 = t176 * t285 - t288;
t58 = -t300 + t340;
t55 = t106 * pkin(4) + t213;
t54 = t243 + t340;
t35 = -pkin(4) * t189 + t283;
t33 = -t106 * pkin(5) + t57;
t14 = t62 * pkin(9) + t17;
t13 = -t61 * pkin(5) + t16;
t12 = -pkin(5) * t62 - t15;
t8 = t247 - t318;
t1 = t199 * t4;
t2 = [qJDD(1), t333, t250, qJDD(1) * t192 + 0.2e1 * t200 * t267, 0.2e1 * t197 * t273 - 0.2e1 * t278 * t274, qJDD(2) * t197 + t200 * t203, qJDD(2) * t200 - t197 * t203, 0, t240 * t197 + t221 * t200, -t221 * t197 + t240 * t200, t105 * t217 + t112 * t210 + t113 * t207 - t87 * t141 + t86 * t142 + t59 * t227 - t250 + (qJD(2) * t104 - t60) * t228, t60 * t113 + t105 * t87 + t59 * t112 + t104 * t86 - t230 * t179 + t158 * t182 - g(1) * (-t179 * t198 + t201 * t310) - g(2) * (t179 * t201 + t198 * t310) t107 * t236 - t239 * t61, -t106 * t236 - t107 * t44 - t239 * t62 - t61 * t96, t107 * t187 - t189 * t61, -t106 * t187 - t189 * t62, 0, t76 * t106 + t111 * t62 - t115 * t96 + t118 * t44 + t328, t76 * t107 - t111 * t61 + t115 * t239 + t118 * t236 + t329, t106 * t7 + t107 * t8 - t15 * t96 + t16 * t239 + t236 * t56 - t35 * t61 + t37 * t62 - t44 * t57 - t250, -t9 * t106 + t17 * t96 - t55 * t44 - t51 * t62 - t328, -t9 * t107 - t17 * t239 - t236 * t55 + t51 * t61 - t329, t37 * t15 + t35 * t16 + t51 * t17 + t9 * t55 + t8 * t56 - t7 * t57 + (g(1) * t188 - g(2) * t241) * t201 + (g(1) * t241 + g(2) * t188) * t198, t62 * t303 + (t20 * t195 + t275 * t80) * t106 (-t195 * t78 + t199 * t80) * t62 + (-t195 * t21 + t19 + (-t199 * t78 - t303) * qJD(6)) * t106, t62 * t350 + t20 * t107 - t80 * t61 + (t275 * t332 - t304) * t106, t62 * t349 - t21 * t107 + t78 * t61 + (-t276 * t332 - t36) * t106, -t107 * t42 - t332 * t61, -g(1) * t134 - g(2) * t132 + t1 * t107 - t10 * t61 + t12 * t78 + t33 * t21 + (-t6 * t107 - t14 * t332 + t316) * t195 + (t13 * t332 - t234) * t199 + ((-t195 * t32 - t199 * t31) * t332 - t11 * t107 + t195 * t308) * qJD(6), g(1) * t133 - g(2) * t131 + t11 * t61 + t12 * t80 + t33 * t20 + (-(qJD(6) * t32 + t14) * t332 + t316 - t264 * t107 + qJD(6) * t308) * t199 + (-(-qJD(6) * t31 + t13) * t332 - t263 * t107 + t234) * t195; 0, 0, 0, -t197 * t204 * t200, t278 * t204, t197 * qJDD(1), t273, qJDD(2), t197 * t220 - t321, g(3) * t197 + t200 * t220, t207 * t270 + t210 * t269 + (-t105 - t109) * t142 - (-t110 + t104) * t141, -t104 * t109 - t105 * t110 + (t294 * t60 + t295 * t59 - t321 + (-qJD(1) * t158 + t250) * t197) * pkin(2), -t311, t358, t223, t348, t187, t114 * t96 + t248 * t187 + t330 - t336, -t114 * t239 - t226 * t187 + t351 * t189 + t343, -t135 * t44 + t136 * t236 + (t296 - t37) * t239 + (-t298 - t35) * t96, -t54 * t96 + t336 + (-pkin(4) + t136) * t187 + t216, t135 * t187 - t298 * t189 + t239 * t54 + t344, -t7 * t135 + t8 * t136 - t51 * t54 - g(1) * (t201 * t266 + t160) - g(2) * (t198 * t266 + t159) - g(3) * t347 + t298 * t37 + t296 * t35, t357, t356, t359, t360, -t352, t135 * t21 + t299 * t78 + (t331 + t338) * t199 + t212 * t195 + t345, t135 * t20 + t299 * t80 + t212 * t199 + (-t331 - t354) * t195 + t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 ^ 2 - t142 ^ 2, -t104 * t142 + t105 * t141 + t230 - t333, 0, 0, 0, 0, 0, -t205, t222, -t313 - t341, t205, -t222, -t239 * t35 + t37 * t96 - t333 + t9, 0, 0, 0, 0, 0, t355 - t312, -t242 - t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t311, t358, t223, t348, t187, t306 + t330, -t189 * t48 + t343, -pkin(4) * t236 - qJ(5) * t44 + (-t37 - t49) * t239 - (t35 - t283) * t96, -t58 * t96 + t216 - t306 - 0.2e1 * t318, t189 * t283 + t239 * t58 + t180 + t344, -t7 * qJ(5) - t8 * pkin(4) - t51 * t58 - t35 * t49 - g(1) * (-pkin(4) * t291 + t160) - g(2) * (-pkin(4) * t292 + t159) - g(3) * t280 - t283 * t37, t357, t356, t359, t360, -t352, qJ(5) * t21 + t284 * t78 + (-t246 + t338) * t199 + t211 * t195 + t345, qJ(5) * t20 + t284 * t80 + (t246 - t354) * t195 + t211 * t199 + t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t187 + t311, -t189 ^ 2 - t341, t189 * t37 + t216 - t318, 0, 0, 0, 0, 0, -t189 * t78 + t242, -t189 * t80 + t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t78, -t78 ^ 2 + t80 ^ 2, t20 + t353, t335 * t80 - t256, -t42, -g(1) * t131 - g(2) * t133 + t335 * t11 + t199 * t171 - t195 * t6 - t24 * t80 + t1, g(1) * t132 - g(2) * t134 + t10 * t332 + t24 * t78 - t264 * t199 + (-t263 - t171) * t195;];
tau_reg  = t2;
