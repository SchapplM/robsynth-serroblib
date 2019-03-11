% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:51
% EndTime: 2019-03-09 08:28:06
% DurationCPUTime: 6.76s
% Computational Cost: add. (8777->505), mult. (20471->632), div. (0->0), fcn. (15452->14), ass. (0->241)
t206 = sin(pkin(9));
t211 = sin(qJ(2));
t213 = cos(qJ(2));
t311 = cos(pkin(9));
t169 = t206 * t213 + t311 * t211;
t154 = t169 * qJD(1);
t205 = sin(pkin(10));
t207 = cos(pkin(10));
t133 = qJD(2) * t205 + t154 * t207;
t210 = sin(qJ(5));
t335 = cos(qJ(5));
t295 = t205 * t154;
t350 = qJD(2) * t207 - t295;
t236 = t335 * t350;
t78 = -t133 * t210 + t236;
t358 = t78 ^ 2;
t266 = t311 * t213;
t185 = qJD(1) * t266;
t288 = qJD(1) * t211;
t151 = t206 * t288 - t185;
t147 = qJD(5) + t151;
t357 = t147 * t78;
t248 = t210 * t350;
t355 = t335 * t133 + t248;
t336 = t355 ^ 2;
t356 = t147 * t355;
t202 = qJ(2) + pkin(9);
t198 = cos(t202);
t189 = g(3) * t198;
t196 = sin(t202);
t212 = sin(qJ(1));
t214 = cos(qJ(1));
t255 = g(1) * t214 + g(2) * t212;
t228 = t255 * t196 - t189;
t323 = qJ(3) + pkin(7);
t268 = qJD(2) * t323;
t149 = -t211 * qJD(3) - t213 * t268;
t177 = t323 * t211;
t115 = qJDD(2) * pkin(2) + t149 * qJD(1) - qJDD(1) * t177;
t148 = t213 * qJD(3) - t211 * t268;
t178 = t323 * t213;
t124 = t148 * qJD(1) + qJDD(1) * t178;
t64 = t311 * t115 - t206 * t124;
t63 = -qJDD(2) * pkin(3) + qJDD(4) - t64;
t354 = t228 - t63;
t153 = t169 * qJD(2);
t265 = qJDD(1) * t311;
t283 = t211 * qJDD(1);
t252 = t206 * t283 - t213 * t265;
t121 = qJD(1) * t153 + t252;
t116 = qJDD(5) + t121;
t170 = t335 * t205 + t210 * t207;
t237 = -t210 * t205 + t335 * t207;
t273 = qJD(5) * t335;
t285 = qJD(5) * t210;
t339 = -t205 * t285 + t207 * t273;
t349 = t237 * t151 + t339;
t259 = t170 * t116 + t349 * t147;
t314 = t355 * t154;
t353 = t259 + t314;
t334 = pkin(2) * t206;
t188 = qJ(4) + t334;
t324 = pkin(8) + t188;
t162 = t324 * t205;
t163 = t324 * t207;
t112 = -t210 * t162 + t335 * t163;
t201 = pkin(10) + qJ(5);
t195 = sin(t201);
t352 = t112 * t116 + t195 * t228;
t158 = t170 * qJD(5);
t312 = t170 * t151 + t158;
t258 = t237 * t116 - t147 * t312;
t316 = t154 * t78;
t351 = t258 + t316;
t348 = t133 * t151;
t298 = t196 * t214;
t299 = t196 * t212;
t345 = g(1) * t298 + g(2) * t299 - t189;
t199 = t213 * pkin(2);
t194 = t199 + pkin(1);
t173 = qJD(1) * t178;
t159 = t206 * t173;
t172 = qJD(1) * t177;
t126 = -t311 * t172 - t159;
t99 = pkin(2) * t288 + pkin(3) * t154 + qJ(4) * t151;
t62 = t207 * t126 + t205 * t99;
t344 = qJD(4) * t207 - t62;
t238 = -t335 * t162 - t210 * t163;
t306 = t151 * t207;
t61 = -t126 * t205 + t207 * t99;
t41 = pkin(4) * t154 + pkin(8) * t306 + t61;
t307 = t151 * t205;
t52 = pkin(8) * t307 + t62;
t343 = -t237 * qJD(4) - t238 * qJD(5) + t210 * t41 + t335 * t52;
t342 = -t170 * qJD(4) - t112 * qJD(5) + t210 * t52 - t335 * t41;
t341 = g(1) * t212 - g(2) * t214;
t284 = qJD(1) * qJD(2);
t272 = t211 * t284;
t233 = qJD(2) * t185 - t206 * t272;
t221 = t169 * qJDD(1) + t233;
t220 = t205 * qJDD(2) + t207 * t221;
t110 = t205 * t221;
t263 = qJDD(2) * t207 - t110;
t26 = -qJD(5) * t236 + t133 * t285 - t210 * t263 - t335 * t220;
t338 = -t237 * t26 - t312 * t355;
t235 = -t206 * t211 + t266;
t302 = t169 * t207;
t118 = -pkin(3) * t235 - qJ(4) * t169 - t194;
t131 = -t206 * t177 + t311 * t178;
t67 = t207 * t118 - t131 * t205;
t51 = -pkin(4) * t235 - pkin(8) * t302 + t67;
t303 = t169 * t205;
t68 = t205 * t118 + t207 * t131;
t56 = -pkin(8) * t303 + t68;
t243 = t210 * t51 + t335 * t56;
t156 = t235 * qJD(2);
t304 = t156 * t207;
t101 = t311 * t148 + t206 * t149;
t318 = qJD(2) * pkin(2);
t278 = t211 * t318;
t79 = pkin(3) * t153 - qJ(4) * t156 - qJD(4) * t169 + t278;
t47 = -t101 * t205 + t207 * t79;
t30 = pkin(4) * t153 - pkin(8) * t304 + t47;
t305 = t156 * t205;
t48 = t207 * t101 + t205 * t79;
t36 = -pkin(8) * t305 + t48;
t337 = -t243 * qJD(5) - t210 * t36 + t335 * t30;
t150 = t151 ^ 2;
t333 = pkin(2) * t211;
t332 = pkin(5) * t116;
t328 = g(3) * t196;
t327 = g(3) * t213;
t326 = t207 * pkin(4);
t325 = t355 * t78;
t281 = pkin(2) * t272 + qJDD(3);
t282 = t213 * qJDD(1);
t309 = qJDD(1) * pkin(1);
t45 = -pkin(2) * t282 + t121 * pkin(3) - qJ(4) * t221 - t154 * qJD(4) + t281 - t309;
t65 = t206 * t115 + t311 * t124;
t60 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t65;
t23 = t205 * t45 + t207 * t60;
t321 = -qJ(6) * t154 - t343;
t320 = t154 * pkin(5) - t342;
t267 = t311 * t173;
t125 = -t172 * t206 + t267;
t82 = -pkin(4) * t307 + t125;
t319 = t312 * pkin(5) - t349 * qJ(6) - qJD(6) * t170 - t82;
t164 = -t172 + t318;
t120 = t206 * t164 + t267;
t109 = qJD(2) * qJ(4) + t120;
t176 = -qJD(1) * t194 + qJD(3);
t91 = t151 * pkin(3) - t154 * qJ(4) + t176;
t54 = t207 * t109 + t205 * t91;
t53 = -t109 * t205 + t207 * t91;
t33 = pkin(4) * t151 - pkin(8) * t133 + t53;
t42 = pkin(8) * t350 + t54;
t13 = t210 * t33 + t335 * t42;
t317 = t13 * t147;
t310 = qJ(6) * t116;
t301 = t188 * t205;
t209 = -pkin(8) - qJ(4);
t300 = t196 * t209;
t297 = t198 * t212;
t296 = t198 * t214;
t294 = t207 * t121;
t197 = cos(t201);
t292 = t212 * t197;
t291 = t214 * t195;
t12 = -t210 * t42 + t335 * t33;
t290 = qJD(6) - t12;
t203 = t211 ^ 2;
t289 = -t213 ^ 2 + t203;
t275 = t311 * pkin(2);
t22 = -t205 * t60 + t207 * t45;
t10 = t121 * pkin(4) - pkin(8) * t220 + t22;
t18 = t263 * pkin(8) + t23;
t271 = -t335 * t10 + t210 * t18 + t42 * t273 + t33 * t285;
t270 = pkin(4) * t205 + t323;
t100 = t148 * t206 - t311 * t149;
t130 = t311 * t177 + t178 * t206;
t219 = t210 * t220 - t335 * t263;
t27 = qJD(5) * t355 + t219;
t262 = -t170 * t27 + t349 * t78;
t260 = g(1) * t299 - g(2) * t298;
t193 = -t275 - pkin(3);
t136 = t195 * t297 + t197 * t214;
t138 = t198 * t291 - t292;
t257 = -g(1) * t136 + g(2) * t138;
t137 = t198 * t292 - t291;
t139 = t195 * t212 + t197 * t296;
t256 = g(1) * t137 - g(2) * t139;
t70 = pkin(4) * t305 + t100;
t97 = pkin(4) * t303 + t130;
t253 = pkin(3) * t198 + qJ(4) * t196;
t251 = -t22 * t205 + t23 * t207;
t250 = -t205 * t53 + t207 * t54;
t119 = t311 * t164 - t159;
t192 = pkin(3) + t326;
t249 = t192 * t198 - t300;
t247 = pkin(5) * t197 + qJ(6) * t195 + t192;
t245 = -t210 * t56 + t335 * t51;
t242 = t210 * t10 + t335 * t18 + t33 * t273 - t42 * t285;
t240 = -0.2e1 * pkin(1) * t284 - pkin(7) * qJDD(2);
t239 = t210 * t30 + t51 * t273 - t56 * t285 + t335 * t36;
t175 = t193 - t326;
t232 = -qJDD(1) * t194 + t281;
t231 = t238 * t116 + t345 * t197;
t103 = -qJD(2) * pkin(3) + qJD(4) - t119;
t230 = t26 - t357;
t229 = g(1) * t138 + g(2) * t136 + t195 * t328 - t271;
t215 = qJD(2) ^ 2;
t227 = -pkin(7) * t215 + 0.2e1 * t309 + t341;
t216 = qJD(1) ^ 2;
t226 = pkin(1) * t216 - pkin(7) * qJDD(1) + t255;
t69 = -pkin(4) * t350 + t103;
t24 = -pkin(5) * t78 - qJ(6) * t355 + t69;
t225 = t24 * t355 + qJDD(6) - t229;
t223 = -g(1) * t139 - g(2) * t137 - t197 * t328 + t242;
t37 = -t263 * pkin(4) + t63;
t4 = t27 * pkin(5) + t26 * qJ(6) - qJD(6) * t355 + t37;
t218 = -qJD(5) * t248 - t133 * t273 - t219;
t217 = -t218 + t356;
t181 = t214 * t194;
t105 = t237 * t169;
t104 = t170 * t169;
t102 = -pkin(5) * t237 - t170 * qJ(6) + t175;
t59 = t170 * t156 + t339 * t169;
t58 = -t237 * t156 + t169 * t158;
t40 = pkin(5) * t355 - qJ(6) * t78;
t35 = pkin(5) * t104 - qJ(6) * t105 + t97;
t21 = pkin(5) * t235 - t245;
t20 = -qJ(6) * t235 + t243;
t19 = -t26 - t357;
t11 = pkin(5) * t59 + qJ(6) * t58 - qJD(6) * t105 + t70;
t9 = t147 * qJ(6) + t13;
t8 = -t147 * pkin(5) + t290;
t5 = -t153 * pkin(5) - t337;
t3 = qJ(6) * t153 - qJD(6) * t235 + t239;
t2 = qJDD(6) + t271 - t332;
t1 = qJD(6) * t147 + t242 + t310;
t6 = [qJDD(1), t341, t255, qJDD(1) * t203 + 0.2e1 * t213 * t272, 0.2e1 * t211 * t282 - 0.2e1 * t289 * t284, qJDD(2) * t211 + t213 * t215, qJDD(2) * t213 - t211 * t215, 0, t240 * t211 + t227 * t213, -t227 * t211 + t240 * t213, t100 * t154 - t101 * t151 - t119 * t156 - t120 * t153 - t131 * t121 + t130 * t221 - t64 * t169 + t235 * t65 - t255, t65 * t131 + t120 * t101 - t64 * t130 - t119 * t100 - t232 * t194 + t176 * t278 - g(1) * (-t194 * t212 + t214 * t323) - g(2) * (t212 * t323 + t181) t47 * t151 + t67 * t121 - t22 * t235 + t53 * t153 - t100 * t350 - t130 * t263 + t341 * t207 * t198 + (t103 * t156 + t63 * t169 - t255) * t205, -t48 * t151 - t68 * t121 + t23 * t235 - t54 * t153 + t100 * t133 + t130 * t220 + t63 * t302 + t103 * t304 - g(1) * (t205 * t297 + t207 * t214) - g(2) * (-t205 * t296 + t207 * t212) -t47 * t133 - t22 * t302 - t67 * t220 - t23 * t303 + t68 * t263 - t53 * t304 - t54 * t305 + t350 * t48 + t260, -g(2) * t181 + t103 * t100 + t63 * t130 + t22 * t67 + t23 * t68 + t53 * t47 + t54 * t48 + (-g(1) * t323 - g(2) * t253) * t214 + (-g(1) * (-t194 - t253) - g(2) * t323) * t212, -t105 * t26 - t355 * t58, t104 * t26 - t105 * t27 - t355 * t59 - t58 * t78, t105 * t116 - t147 * t58 + t153 * t355 + t235 * t26, -t104 * t116 - t147 * t59 + t153 * t78 + t235 * t27, -t116 * t235 + t147 * t153, t37 * t104 + t245 * t116 + t12 * t153 + t337 * t147 + t235 * t271 + t97 * t27 + t69 * t59 - t70 * t78 + t256, t37 * t105 - t116 * t243 - t13 * t153 - t147 * t239 + t235 * t242 - t97 * t26 + t355 * t70 - t69 * t58 + t257, t104 * t4 - t11 * t78 - t116 * t21 - t147 * t5 - t153 * t8 + t2 * t235 + t24 * t59 + t27 * t35 + t256, -t1 * t104 + t105 * t2 - t20 * t27 - t21 * t26 + t3 * t78 + t355 * t5 - t58 * t8 - t59 * t9 + t260, -t1 * t235 - t105 * t4 - t11 * t355 + t116 * t20 + t147 * t3 + t153 * t9 + t24 * t58 + t26 * t35 - t257, t1 * t20 + t9 * t3 + t4 * t35 + t24 * t11 + t2 * t21 + t8 * t5 - g(1) * (-pkin(5) * t137 - qJ(6) * t136) - g(2) * (pkin(5) * t139 + qJ(6) * t138 + t181) + (-g(1) * t270 - g(2) * t249) * t214 + (-g(1) * (-t194 - t249) - g(2) * t270) * t212; 0, 0, 0, -t211 * t216 * t213, t289 * t216, t283, t282, qJDD(2), t211 * t226 - t327, g(3) * t211 + t213 * t226, -t121 * t334 - t221 * t275 - (-t120 + t125) * t154 + (t126 - t119) * t151, t119 * t125 - t120 * t126 + (t311 * t64 - t327 + t206 * t65 + (-qJD(1) * t176 + t255) * t211) * pkin(2), -t121 * t301 + t193 * t110 - t125 * t295 - t53 * t154 + (-t61 + (-qJD(4) + t103) * t205) * t151 + (t125 * qJD(2) - t193 * qJDD(2) + t354) * t207, t103 * t306 - t125 * t133 + t54 * t154 - t188 * t294 + t193 * t220 - t344 * t151 + (t63 - t345) * t205, t207 * t188 * t263 - g(1) * t296 - g(2) * t297 + t220 * t301 - t53 * t306 - t54 * t307 + t251 - t328 + t344 * t350 + (qJD(4) * t205 + t61) * t133, t63 * t193 - t54 * t62 - t53 * t61 - t103 * t125 - g(3) * (t199 + t253) + t251 * t188 + t250 * qJD(4) + t255 * (pkin(3) * t196 - qJ(4) * t198 + t333) -t26 * t170 + t349 * t355, t262 + t338, t259 - t314, t258 - t316, -t147 * t154, -t12 * t154 + t342 * t147 + t175 * t27 - t237 * t37 + t312 * t69 + t78 * t82 + t231, t13 * t154 + t343 * t147 + t37 * t170 - t175 * t26 + t349 * t69 - t355 * t82 - t352, t102 * t27 - t147 * t320 + t154 * t8 - t237 * t4 + t24 * t312 - t319 * t78 + t231, t1 * t237 - t112 * t27 + t170 * t2 - t198 * t255 + t238 * t26 - t312 * t9 + t320 * t355 + t321 * t78 + t349 * t8 - t328, t102 * t26 + t147 * t321 - t154 * t9 - t170 * t4 - t24 * t349 - t319 * t355 + t352, t1 * t112 + t4 * t102 - t2 * t238 - g(3) * (t199 - t300) + t321 * t9 + t320 * t8 + t319 * t24 - t247 * t189 + t255 * (t247 * t196 + t198 * t209 + t333); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154 ^ 2 - t150, t119 * t154 + t120 * t151 + t232 - t341, -t205 * t150 + t154 * t350 + t294, -t121 * t205 - t133 * t154 - t150 * t207 (-t151 * t295 + (t151 * qJD(2) - t206 * t282 - t211 * t265 - t233) * t207) * t207 + (-t110 + t348) * t205, -t103 * t154 + t250 * t151 + t205 * t23 + t207 * t22 - t341, 0, 0, 0, 0, 0, t351, -t353, t351, t262 - t338, t353, t1 * t170 - t154 * t24 - t2 * t237 + t312 * t8 + t349 * t9 - t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263 + t348, t151 * t350 + t220, -t133 ^ 2 - t350 ^ 2, t133 * t53 - t350 * t54 - t354, 0, 0, 0, 0, 0, t217, -t230, t217, -t336 - t358, t230, -t355 * t8 - t9 * t78 - t228 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t325, t336 - t358, t19, t218 + t356, t116, -t355 * t69 + t229 + t317, t12 * t147 - t69 * t78 - t223, t40 * t78 - t225 + t317 + 0.2e1 * t332, pkin(5) * t26 - qJ(6) * t27 + (-t13 + t9) * t355 - (t8 - t290) * t78, 0.2e1 * t310 + t24 * t78 + t40 * t355 + (0.2e1 * qJD(6) - t12) * t147 + t223, t1 * qJ(6) - t2 * pkin(5) - t24 * t40 - t8 * t13 - g(1) * (-pkin(5) * t138 + qJ(6) * t139) - g(2) * (-pkin(5) * t136 + qJ(6) * t137) + t290 * t9 - (-pkin(5) * t195 + qJ(6) * t197) * t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t154 - qJDD(5) - t252 - t325, t19, -t147 ^ 2 - t336, -t147 * t9 + t225 - t332;];
tau_reg  = t6;
