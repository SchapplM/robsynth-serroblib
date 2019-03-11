% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:55
% EndTime: 2019-03-09 02:24:04
% DurationCPUTime: 5.74s
% Computational Cost: add. (6951->523), mult. (13012->667), div. (0->0), fcn. (8158->14), ass. (0->267)
t193 = sin(qJ(4));
t302 = qJD(1) * t193;
t158 = qJD(5) + t302;
t197 = cos(qJ(4));
t348 = g(3) * t193;
t183 = qJ(1) + pkin(10);
t174 = cos(t183);
t164 = g(2) * t174;
t173 = sin(t183);
t165 = g(1) * t173;
t361 = t165 - t164;
t215 = t197 * t361 - t348;
t190 = cos(pkin(10));
t166 = -pkin(1) * t190 - pkin(2);
t157 = -pkin(7) + t166;
t119 = qJDD(1) * t157 + qJDD(3);
t122 = qJD(1) * t157 + qJD(3);
t286 = qJD(2) * qJD(4);
t297 = qJD(4) * t193;
t245 = t193 * qJDD(2) + t122 * t297 + (-t119 + t286) * t197;
t327 = qJDD(4) * pkin(4);
t44 = t245 - t327;
t211 = -t44 - t215;
t372 = qJD(5) * pkin(8) * t158 - t211;
t192 = sin(qJ(5));
t301 = qJD(1) * t197;
t265 = t192 * t301;
t196 = cos(qJ(5));
t288 = t196 * qJD(4);
t127 = t265 - t288;
t298 = qJD(4) * t192;
t129 = t196 * t301 + t298;
t191 = sin(qJ(6));
t195 = cos(qJ(6));
t227 = t127 * t191 - t195 * t129;
t65 = t195 * t127 + t129 * t191;
t346 = t65 * t227;
t266 = t192 * t302;
t199 = -pkin(9) - pkin(8);
t267 = qJD(5) * t199;
t240 = pkin(4) * t197 + pkin(8) * t193;
t132 = t240 * qJD(1);
t86 = -t193 * qJD(2) + t122 * t197;
t53 = t192 * t132 + t196 * t86;
t371 = -pkin(9) * t266 + t192 * t267 - t53;
t312 = t193 * t196;
t276 = pkin(9) * t312;
t52 = t196 * t132 - t192 * t86;
t370 = t196 * t267 - (pkin(5) * t197 + t276) * qJD(1) - t52;
t278 = t197 * qJDD(1);
t369 = qJD(4) * qJD(5) + t278;
t368 = t227 ^ 2 - t65 ^ 2;
t152 = qJD(6) + t158;
t293 = qJD(5) * t197;
t258 = t196 * t293;
t270 = qJD(1) * t258 + t192 * t369;
t222 = t196 * qJDD(4) - t270;
t287 = qJD(1) * qJD(4);
t254 = t193 * t287;
t209 = t192 * t254 + t222;
t291 = qJD(6) * t195;
t292 = qJD(6) * t191;
t261 = t193 * t288;
t217 = t192 * t293 + t261;
t62 = qJD(1) * t217 - t192 * qJDD(4) - t196 * t369;
t16 = t127 * t291 + t129 * t292 - t191 * t209 + t195 * t62;
t367 = t152 * t65 - t16;
t300 = qJD(2) * t197;
t87 = t122 * t193 + t300;
t78 = qJD(4) * pkin(8) + t87;
t189 = sin(pkin(10));
t161 = pkin(1) * t189 + qJ(3);
t239 = pkin(4) * t193 - pkin(8) * t197;
t118 = t161 + t239;
t88 = t118 * qJD(1);
t41 = -t192 * t78 + t196 * t88;
t29 = -pkin(9) * t129 + t41;
t26 = pkin(5) * t158 + t29;
t42 = t192 * t88 + t196 * t78;
t30 = -pkin(9) * t127 + t42;
t253 = t197 * t287;
t279 = t193 * qJDD(1);
t124 = qJDD(5) + t253 + t279;
t296 = qJD(4) * t197;
t274 = -t197 * qJDD(2) - t193 * t119 - t122 * t296;
t46 = -t193 * t286 - t274;
t43 = qJDD(4) * pkin(8) + t46;
t125 = qJD(4) * t240 + qJD(3);
t328 = pkin(1) * qJDD(1);
t167 = t189 * t328;
t306 = qJDD(1) * qJ(3) + t167;
t59 = qJD(1) * t125 + qJDD(1) * t239 + t306;
t51 = t196 * t59;
t9 = -qJD(5) * t42 - t192 * t43 + t51;
t6 = pkin(5) * t124 + pkin(9) * t62 + t9;
t294 = qJD(5) * t196;
t295 = qJD(5) * t192;
t8 = t192 * t59 + t196 * t43 + t88 * t294 - t295 * t78;
t7 = pkin(9) * t209 + t8;
t1 = t195 * (qJD(6) * t26 + t7) + t191 * t6 - t30 * t292;
t188 = qJ(5) + qJ(6);
t181 = cos(t188);
t347 = g(3) * t197;
t338 = qJD(4) * pkin(4);
t77 = -t86 - t338;
t54 = pkin(5) * t127 + t77;
t180 = sin(t188);
t317 = t181 * t193;
t83 = t173 * t317 + t174 * t180;
t85 = -t173 * t180 + t174 * t317;
t366 = g(1) * t83 - g(2) * t85 + t181 * t347 + t54 * t65 - t1;
t334 = t195 * t30;
t11 = t191 * t26 + t334;
t2 = -qJD(6) * t11 - t191 * t7 + t195 * t6;
t318 = t180 * t193;
t82 = -t173 * t318 + t174 * t181;
t84 = t173 * t181 + t174 * t318;
t365 = -g(1) * t82 - g(2) * t84 + t180 * t347 + t54 * t227 + t2;
t210 = qJD(6) * t227 + t191 * t62 + t195 * t209;
t364 = -t152 * t227 + t210;
t363 = -t158 * t41 + t8;
t311 = t195 * t196;
t316 = t191 * t192;
t130 = -t311 + t316;
t104 = t130 * t193;
t131 = t191 * t196 + t192 * t195;
t102 = t131 * t193;
t362 = t209 * t193;
t314 = t192 * t193;
t91 = -t173 * t314 + t174 * t196;
t93 = t173 * t196 + t174 * t314;
t360 = -g(1) * t91 - g(2) * t93;
t359 = qJDD(1) * t166;
t277 = qJD(5) + qJD(6);
t232 = t192 * t41 - t196 * t42;
t358 = qJD(4) * t232 + t44;
t310 = t196 * t197;
t357 = t124 * t310 - t158 * t217;
t244 = qJD(5) * t193 + qJD(1);
t257 = t192 * t296;
t356 = t196 * t244 + t257;
t207 = -(t193 * t86 - t197 * t87) * qJD(4) + t46 * t193 - t245 * t197;
t137 = qJD(1) * t161;
t355 = 0.2e1 * qJD(4) * t137 + qJDD(4) * t157;
t354 = -pkin(2) - pkin(7);
t351 = pkin(5) * t192;
t350 = pkin(8) * t124;
t349 = g(2) * t173;
t105 = t130 * t197;
t263 = t192 * t297;
t73 = t277 * t131;
t34 = -t191 * t263 + t195 * t261 + t197 * t73;
t345 = -t105 * t210 + t34 * t65;
t103 = t131 * t197;
t117 = qJDD(6) + t124;
t313 = t192 * t197;
t36 = -t292 * t313 + (t277 * t310 - t263) * t195 - t217 * t191;
t344 = -t103 * t117 - t36 * t152;
t144 = t199 * t192;
t145 = t199 * t196;
t79 = t144 * t195 + t145 * t191;
t343 = qJD(6) * t79 + t191 * t370 + t195 * t371;
t80 = t144 * t191 - t145 * t195;
t342 = -qJD(6) * t80 - t191 * t371 + t195 * t370;
t341 = -t16 * t193 - t227 * t296;
t340 = t191 * t266 - t195 * t294 - t196 * t291 + t277 * t316 - t302 * t311;
t114 = t131 * qJD(1);
t339 = t193 * t114 + t73;
t336 = t158 * t42;
t335 = t191 * t30;
t333 = t196 * t41;
t332 = t62 * t192;
t331 = t129 * t296 - t62 * t193;
t330 = -t130 * qJD(1) - t104 * t277 + t131 * t296;
t329 = qJD(4) * t105 + t102 * t277 + t114;
t126 = t157 * t312;
t70 = t192 * t118 + t126;
t326 = t124 * t196;
t325 = t127 * t158;
t324 = t127 * t192;
t323 = t129 * t127;
t322 = t129 * t196;
t321 = t174 * t192;
t320 = t174 * t193;
t319 = t174 * t197;
t315 = t192 * t124;
t309 = t197 * t199;
t308 = g(1) * t319 + t197 * t349;
t185 = t193 ^ 2;
t186 = t197 ^ 2;
t305 = t185 - t186;
t200 = qJD(4) ^ 2;
t201 = qJD(1) ^ 2;
t304 = -t200 - t201;
t303 = qJD(1) * t137;
t299 = qJD(4) * t127;
t290 = t127 * qJD(5);
t289 = t129 * qJD(5);
t285 = qJD(3) * qJD(1);
t283 = qJDD(1) * t161;
t281 = qJDD(4) * t193;
t280 = qJDD(4) * t197;
t272 = t157 * t314;
t271 = t197 * t201 * t193;
t198 = cos(qJ(1));
t268 = t198 * pkin(1) + t174 * pkin(2) + t173 * qJ(3);
t264 = t129 * t297;
t262 = t127 * t297;
t260 = t197 * t289;
t256 = t197 * t288;
t194 = sin(qJ(1));
t252 = -pkin(1) * t194 + t174 * qJ(3);
t251 = -t157 * t192 + pkin(5);
t250 = -g(2) * t320 + t347;
t246 = (-t185 - t186) * qJDD(1);
t243 = t174 * pkin(7) + t268;
t242 = t193 * t253;
t241 = pkin(5) * t295 - t300 - (-qJD(1) * t351 + t122) * t193;
t238 = g(1) * t174 + t349;
t236 = g(1) * t194 - g(2) * t198;
t235 = -t103 * t16 - t227 * t36;
t234 = -t303 - t361;
t99 = t196 * t118;
t48 = -pkin(9) * t310 + t193 * t251 + t99;
t60 = -pkin(9) * t313 + t70;
t22 = -t191 * t60 + t195 * t48;
t23 = t191 * t48 + t195 * t60;
t233 = -t192 * t42 - t333;
t229 = -t303 - t165;
t228 = t105 * t117 + t152 * t34;
t172 = pkin(5) * t196 + pkin(4);
t225 = t172 * t193 + t309;
t224 = g(2) * t243;
t123 = t285 + t306;
t223 = t137 * qJD(3) + t123 * t161;
t221 = t193 * t210 - t296 * t65;
t220 = t158 * t294 + t315;
t219 = -t165 * t193 - t250;
t218 = qJDD(3) + t359;
t216 = -t258 + t263;
t31 = -qJD(5) * t272 + t118 * t294 + t192 * t125 + t157 * t256;
t213 = -t157 * t200 + t123 + t283 + t285;
t212 = t192 * t244 - t256;
t208 = qJD(5) * t233 - t9 * t192 + t8 * t196;
t206 = t209 * t196;
t204 = qJD(4) * t77 + t208;
t187 = qJDD(2) - g(3);
t136 = -t193 * t200 + t280;
t135 = -t197 * t200 - t281;
t121 = t158 * t263;
t110 = t196 * t125;
t106 = (-t157 + t351) * t197;
t94 = -t173 * t192 + t174 * t312;
t92 = t173 * t312 + t321;
t74 = -pkin(5) * t216 + t157 * t297;
t69 = t99 - t272;
t49 = t197 * t206;
t32 = -qJD(5) * t70 - t157 * t257 + t110;
t25 = pkin(9) * t216 + t31;
t24 = t110 + (-t126 + (pkin(9) * t197 - t118) * t192) * qJD(5) + (t197 * t251 + t276) * qJD(4);
t21 = -pkin(5) * t209 + t44;
t13 = t195 * t29 - t335;
t12 = -t191 * t29 - t334;
t10 = t195 * t26 - t335;
t4 = -qJD(6) * t23 - t191 * t25 + t195 * t24;
t3 = qJD(6) * t22 + t191 * t24 + t195 * t25;
t5 = [0, 0, 0, 0, 0, qJDD(1), t236, g(1) * t198 + g(2) * t194, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t190 * t328 + t361, -0.2e1 * t167 + t238, 0 (t236 + (t189 ^ 2 + t190 ^ 2) * t328) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t361 + 0.2e1 * t359, -t238 + t283 + 0.2e1 * t285 + t306, t218 * t166 - g(1) * (-pkin(2) * t173 + t252) - g(2) * t268 + t223, qJDD(1) * t186 - 0.2e1 * t242, -0.2e1 * t193 * t278 + 0.2e1 * t287 * t305, t136, qJDD(1) * t185 + 0.2e1 * t242, t135, 0, t355 * t197 + (t213 - t238) * t193, -t193 * t355 + t213 * t197 - t308, t157 * t246 - t207 + t361, -g(1) * (t173 * t354 + t252) - t224 + t207 * t157 + t223, -t129 * t217 - t310 * t62, t49 + (-t260 + t262) * t196 + (t264 + (t62 + t290) * t197) * t192, t331 + t357, -t127 * t216 - t209 * t313, t121 + t362 + (-t220 - t299) * t197, t124 * t193 + t158 * t296, -g(1) * t94 - g(2) * t92 + t69 * t124 + t32 * t158 + (t9 + (t127 * t157 - t192 * t77) * qJD(4)) * t193 + (t41 * qJD(4) + t157 * t209 + t44 * t192 + t294 * t77) * t197, g(1) * t93 - g(2) * t91 - t124 * t70 - t158 * t31 + (-t8 + (t129 * t157 - t196 * t77) * qJD(4)) * t193 + (-qJD(4) * t42 + t157 * t62 + t44 * t196 - t295 * t77) * t197, -t31 * t127 + t70 * t222 - t32 * t129 + t69 * t62 + (t333 + (qJD(1) * t70 + t42) * t192) * t297 + (qJD(5) * t232 - t8 * t192 - t9 * t196) * t197 + t308, t8 * t70 + t42 * t31 + t9 * t69 + t41 * t32 - g(1) * (pkin(4) * t320 - pkin(8) * t319 + t252) - t224 + (-g(1) * t354 - g(2) * t239) * t173 + (-t197 * t44 + t297 * t77) * t157, t105 * t16 + t227 * t34, -t235 + t345, -t228 + t341, -t103 * t210 + t36 * t65, t221 + t344, t117 * t193 + t152 * t296, -g(1) * t85 - g(2) * t83 + t10 * t296 + t103 * t21 - t106 * t210 + t117 * t22 + t152 * t4 + t193 * t2 + t36 * t54 + t65 * t74, g(1) * t84 - g(2) * t82 - t1 * t193 - t105 * t21 - t106 * t16 - t11 * t296 - t117 * t23 - t152 * t3 - t227 * t74 - t34 * t54, -t1 * t103 + t10 * t34 + t105 * t2 - t11 * t36 + t16 * t22 + t210 * t23 + t227 * t4 - t3 * t65 + t308, t1 * t23 + t11 * t3 + t2 * t22 + t10 * t4 + t21 * t106 + t54 * t74 - g(1) * (t172 * t320 + t174 * t309 + t252) - g(2) * (pkin(5) * t321 + t243) + (-g(1) * (-t351 + t354) - g(2) * t225) * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, 0, 0, 0, 0, 0, 0, t135, -t136, 0, t193 * t245 + t197 * t46 - g(3) + (-t193 * t87 - t197 * t86) * qJD(4), 0, 0, 0, 0, 0, 0, t121 - t362 + (-t220 + t299) * t197, t331 - t357, t49 + (t260 + t262) * t196 + (-t264 + (-t62 + t290) * t197) * t192, t193 * t358 + t204 * t197 - g(3), 0, 0, 0, 0, 0, 0, -t221 + t344, t228 + t341, t235 + t345, -t1 * t105 - t10 * t36 - t103 * t2 - t11 * t34 + t193 * t21 + t296 * t54 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t201, t218 + t234, 0, 0, 0, 0, 0, 0, t193 * t304 + t280, t197 * t304 - t281, t246, t207 + t234, 0, 0, 0, 0, 0, 0, t197 * t222 + (-t315 + (t127 + t265) * qJD(4)) * t193 - t356 * t158, t197 * t62 + (qJD(4) * t129 - t326) * t193 + t212 * t158 (t222 * t196 + (t196 * t254 - t62) * t192) * t193 + t356 * t129 + t212 * t127, t233 * qJD(1) + t204 * t193 - t197 * t358 - t361, 0, 0, 0, 0, 0, 0, -t102 * t117 - t152 * t330 + t197 * t210 + t297 * t65, t104 * t117 + t152 * t329 + t16 * t197 - t227 * t297, -t102 * t16 - t104 * t210 - t227 * t330 + t329 * t65, -t1 * t104 - t10 * t330 - t102 * t2 - t11 * t329 - t197 * t21 + t297 * t54 - t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, -t305 * t201, t278, -t271, -t279, qJDD(4), t348 + qJD(4) * t87 + (t229 + t164) * t197 - t245, qJD(4) * t86 + (-t229 + t286) * t193 + t250 + t274, 0, 0, t158 * t322 - t332 (-t62 - t325) * t196 + (-t289 + (-t129 + t298) * t302 + t222) * t192 (-t129 * t197 + t158 * t312) * qJD(1) + t220, t158 * t324 + t206, -t158 * t295 + t326 + (t127 * t197 - t158 * t314) * qJD(1), -t158 * t301, -pkin(4) * t270 - t52 * t158 - t41 * t301 - t87 * t127 + (-t350 + t77 * qJD(5) + (t77 + t338) * t302) * t192 + (t327 - t372) * t196, t42 * t301 + pkin(4) * t62 - t129 * t87 + t158 * t53 + (t158 * t77 - t350) * t196 + t372 * t192, t53 * t127 + t52 * t129 + t363 * t196 + (-t9 - t336) * t192 + (-t332 + t206 + (t322 + t324) * qJD(5)) * pkin(8) + t219, -t41 * t52 - t42 * t53 - t77 * t87 + t211 * pkin(4) + (-t193 * t361 + t208 - t347) * pkin(8), -t131 * t16 + t227 * t340, t130 * t16 + t131 * t210 + t227 * t339 + t340 * t65, t117 * t131 - t152 * t340 + t227 * t301, -t130 * t210 + t339 * t65, -t117 * t130 - t152 * t339 + t301 * t65, -t152 * t301, -t10 * t301 + t117 * t79 + t130 * t21 + t152 * t342 + t172 * t210 - t181 * t215 + t241 * t65 + t339 * t54, t11 * t301 - t117 * t80 + t131 * t21 - t152 * t343 + t16 * t172 + t180 * t215 - t227 * t241 - t340 * t54, -t1 * t130 + t10 * t340 - t11 * t339 - t131 * t2 + t16 * t79 + t210 * t80 + t227 * t342 - t343 * t65 + t219, g(3) * t225 + t1 * t80 + t342 * t10 + t343 * t11 - t21 * t172 + t2 * t79 + t241 * t54 - t361 * (t172 * t197 - t193 * t199); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, -t127 ^ 2 + t129 ^ 2, t325 - t62, -t323, t129 * t158 + t209, t124, -t78 * t294 - t129 * t77 + t336 + t51 + (-qJD(5) * t88 + t347 - t43) * t192 + t360, g(1) * t92 - g(2) * t94 + g(3) * t310 + t127 * t77 - t363, 0, 0, -t346, t368, t367, t346, t364, t117, -t12 * t152 + (t117 * t195 - t129 * t65 - t152 * t292) * pkin(5) + t365, t13 * t152 + (-t117 * t191 + t129 * t227 - t152 * t291) * pkin(5) + t366, -t10 * t65 - t11 * t227 - t12 * t227 + t13 * t65 + (t16 * t195 + t210 * t191 + (-t191 * t227 - t195 * t65) * qJD(6)) * pkin(5), -t10 * t12 - t11 * t13 + (t1 * t191 + t2 * t195 - t54 * t129 + g(3) * t313 + (-t10 * t191 + t11 * t195) * qJD(6) + t360) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t346, t368, t367, t346, t364, t117, t11 * t152 + t365, t10 * t152 + t366, 0, 0;];
tau_reg  = t5;
