% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:10
% EndTime: 2019-03-09 18:04:26
% DurationCPUTime: 7.46s
% Computational Cost: add. (13285->448), mult. (32819->585), div. (0->0), fcn. (26066->16), ass. (0->266)
t259 = cos(qJ(5));
t260 = cos(qJ(3));
t261 = cos(qJ(2));
t338 = qJD(1) * t261;
t321 = t260 * t338;
t255 = sin(qJ(3));
t256 = sin(qJ(2));
t339 = qJD(1) * t256;
t322 = t255 * t339;
t176 = -t321 + t322;
t178 = -t255 * t338 - t260 * t339;
t251 = sin(pkin(11));
t252 = cos(pkin(11));
t291 = t176 * t252 - t178 * t251;
t141 = t259 * t291;
t145 = t176 * t251 + t178 * t252;
t254 = sin(qJ(5));
t102 = -t145 * t254 + t141;
t258 = cos(qJ(6));
t333 = qJD(6) * t258;
t410 = t102 * t258 + t333;
t245 = qJDD(2) + qJDD(3);
t240 = qJDD(5) + t245;
t247 = qJD(2) + qJD(3);
t241 = qJD(5) + t247;
t253 = sin(qJ(6));
t334 = qJD(6) * t253;
t335 = qJD(5) * t254;
t329 = t261 * qJDD(1);
t331 = qJD(1) * qJD(2);
t319 = t261 * t331;
t330 = t256 * qJDD(1);
t396 = t319 + t330;
t129 = qJD(3) * t321 - t247 * t322 + t255 * t329 + t260 * t396;
t190 = t255 * t261 + t256 * t260;
t153 = t247 * t190;
t298 = t255 * t330 - t260 * t329;
t130 = qJD(1) * t153 + t298;
t80 = -t129 * t251 - t130 * t252;
t81 = t129 * t252 - t130 * t251;
t38 = -qJD(5) * t141 + t145 * t335 + t254 * t80 + t259 * t81;
t394 = -t259 * t145 - t254 * t291;
t21 = t253 * t240 + t241 * t333 + t258 * t38 - t334 * t394;
t19 = t21 * t253;
t94 = t241 * t253 + t258 * t394;
t11 = t410 * t94 + t19;
t39 = qJD(5) * t394 + t254 * t81 - t259 * t80;
t37 = qJDD(6) + t39;
t34 = t253 * t37;
t366 = t94 * t394;
t395 = qJD(6) + t102;
t10 = t395 * t410 + t34 - t366;
t142 = t145 * pkin(9);
t173 = t178 * qJ(4);
t386 = pkin(7) + pkin(8);
t210 = t386 * t261;
t197 = qJD(1) * t210;
t179 = t255 * t197;
t209 = t386 * t256;
t195 = qJD(1) * t209;
t371 = qJD(2) * pkin(2);
t185 = -t195 + t371;
t305 = t260 * t185 - t179;
t127 = t173 + t305;
t118 = pkin(3) * t247 + t127;
t183 = t260 * t197;
t290 = -t185 * t255 - t183;
t362 = qJ(4) * t176;
t128 = -t290 - t362;
t119 = t251 * t128;
t73 = t252 * t118 - t119;
t64 = pkin(4) * t247 + t142 + t73;
t383 = pkin(9) * t291;
t351 = t252 * t128;
t74 = t251 * t118 + t351;
t66 = t74 - t383;
t32 = -t254 * t66 + t259 * t64;
t30 = -pkin(5) * t241 - t32;
t369 = t102 * t30;
t356 = t102 * t241;
t23 = t38 + t356;
t354 = t394 * t102;
t27 = -t102 ^ 2 + t394 ^ 2;
t244 = t261 * pkin(2);
t377 = pkin(1) + t244;
t208 = t377 * qJD(1);
t155 = t176 * pkin(3) + qJD(4) - t208;
t113 = pkin(4) * t291 + t155;
t250 = qJ(2) + qJ(3);
t234 = pkin(11) + qJ(5) + t250;
t224 = sin(t234);
t211 = g(3) * t224;
t225 = cos(t234);
t257 = sin(qJ(1));
t262 = cos(qJ(1));
t300 = g(1) * t262 + g(2) * t257;
t154 = qJDD(2) * pkin(2) - t386 * t396;
t320 = t256 * t331;
t156 = t386 * (-t320 + t329);
t274 = qJD(3) * t290 + t260 * t154 - t255 * t156;
t59 = pkin(3) * t245 - qJ(4) * t129 + qJD(4) * t178 + t274;
t337 = qJD(3) * t255;
t388 = (qJD(3) * t185 + t156) * t260 + t255 * t154 - t197 * t337;
t63 = -qJ(4) * t130 - qJD(4) * t176 + t388;
t25 = -t251 * t63 + t252 * t59;
t16 = pkin(4) * t245 - pkin(9) * t81 + t25;
t26 = t251 * t59 + t252 * t63;
t17 = pkin(9) * t80 + t26;
t389 = (qJD(5) * t64 + t17) * t259 + t254 * t16 - t66 * t335;
t270 = t102 * t113 + t225 * t300 + t211 - t389;
t22 = qJD(6) * t94 - t258 * t240 + t253 * t38;
t92 = -t258 * t241 + t253 * t394;
t294 = t253 * t94 + t258 * t92;
t392 = t21 * t258 - t94 * t334;
t1 = -t102 * t294 - t253 * t22 - t92 * t333 + t392;
t408 = t300 * t224;
t405 = pkin(5) * t394;
t370 = t394 * t92;
t304 = t195 * t255 - t183;
t131 = t304 + t362;
t344 = -t260 * t195 - t179;
t132 = t173 + t344;
t350 = t252 * t255;
t372 = pkin(2) * qJD(3);
t363 = -t252 * t131 + t132 * t251 + (-t251 * t260 - t350) * t372;
t352 = t251 * t255;
t393 = -t251 * t131 - t252 * t132 + (t252 * t260 - t352) * t372;
t358 = t394 * t241;
t24 = -t39 + t358;
t360 = t395 * t394;
t242 = sin(t250);
t243 = cos(t250);
t404 = -g(3) * t243 + t242 * t300;
t33 = t254 * t64 + t259 * t66;
t31 = pkin(10) * t241 + t33;
t48 = t102 * pkin(5) - pkin(10) * t394 + t113;
t12 = -t253 * t31 + t258 * t48;
t403 = -t12 * t394 + t258 * t408 + t30 * t334;
t13 = t253 * t48 + t258 * t31;
t380 = g(3) * t225;
t390 = qJD(5) * t33 - t259 * t16 + t254 * t17;
t4 = -t240 * pkin(5) + t390;
t324 = -t4 - t380;
t402 = t13 * t394 - t253 * t324 + t30 * t333;
t266 = -t113 * t394 - t380 - t390 + t408;
t401 = t383 - t363;
t400 = -t142 + t393;
t391 = t395 * (pkin(10) * t395 + t405);
t236 = pkin(2) * t260 + pkin(3);
t171 = -pkin(2) * t352 + t252 * t236;
t165 = pkin(4) + t171;
t172 = pkin(2) * t350 + t236 * t251;
t345 = t254 * t165 + t259 * t172;
t343 = -t255 * t209 + t260 * t210;
t231 = pkin(3) * t252 + pkin(4);
t384 = pkin(3) * t251;
t342 = t254 * t231 + t259 * t384;
t189 = t255 * t256 - t260 * t261;
t149 = -t189 * t252 - t190 * t251;
t150 = -t189 * t251 + t190 * t252;
t107 = t149 * t254 + t150 * t259;
t293 = t259 * t149 - t150 * t254;
t316 = t240 * pkin(10) + qJD(6) * t48 + t389;
t303 = -t260 * t209 - t210 * t255;
t143 = -qJ(4) * t190 + t303;
t144 = -qJ(4) * t189 + t343;
t96 = t252 * t143 - t144 * t251;
t71 = -pkin(9) * t150 + t96;
t97 = t251 * t143 + t252 * t144;
t72 = pkin(9) * t149 + t97;
t47 = t254 * t71 + t259 * t72;
t152 = t247 * t189;
t112 = -t152 * t252 - t153 * t251;
t323 = qJD(2) * t386;
t196 = t256 * t323;
t198 = t261 * t323;
t336 = qJD(3) * t260;
t281 = -t260 * t196 - t255 * t198 - t209 * t336 - t210 * t337;
t86 = -qJ(4) * t153 - qJD(4) * t189 + t281;
t273 = -qJD(3) * t343 + t255 * t196 - t260 * t198;
t87 = qJ(4) * t152 - qJD(4) * t190 + t273;
t53 = -t251 * t86 + t252 * t87;
t44 = -pkin(9) * t112 + t53;
t111 = t152 * t251 - t153 * t252;
t54 = t251 * t87 + t252 * t86;
t45 = pkin(9) * t111 + t54;
t46 = t254 * t72 - t259 * t71;
t5 = -qJD(5) * t46 + t254 * t44 + t259 * t45;
t51 = qJD(5) * t293 + t254 * t111 + t259 * t112;
t301 = pkin(3) * t189 - t377;
t123 = -pkin(4) * t149 + t301;
t55 = -pkin(5) * t293 - pkin(10) * t107 + t123;
t387 = -(qJD(6) * t55 + t5) * t395 + t293 * t316 + t4 * t107 + t30 * t51 - t47 * t37;
t385 = pkin(3) * t178;
t378 = t55 * t37;
t292 = t165 * t259 - t172 * t254;
t374 = -qJD(5) * t292 + t254 * t401 - t259 * t400;
t373 = qJD(5) * t345 + t254 * t400 + t259 * t401;
t368 = t107 * t30;
t35 = t258 * t37;
t367 = t258 * t94;
t287 = t231 * t259 - t254 * t384;
t78 = -t127 * t251 - t351;
t67 = t78 + t383;
t79 = t252 * t127 - t119;
t68 = t142 + t79;
t365 = -t287 * qJD(5) + t254 * t67 + t259 * t68;
t364 = qJD(5) * t342 - t254 * t68 + t259 * t67;
t359 = t395 * t253;
t353 = t178 * t176;
t349 = t253 * t257;
t348 = t253 * t262;
t347 = t257 * t258;
t346 = t258 * t262;
t341 = pkin(3) * t243 + t244;
t248 = t256 ^ 2;
t340 = -t261 ^ 2 + t248;
t239 = t256 * t371;
t318 = pkin(3) * t153 + t239;
t174 = pkin(2) * t320 - qJDD(1) * t377;
t276 = t130 * pkin(3) + qJDD(4) + t174;
t56 = -pkin(4) * t80 + t276;
t8 = pkin(5) * t39 - pkin(10) * t38 + t56;
t315 = qJD(6) * t31 - t8;
t135 = pkin(10) + t345;
t238 = pkin(2) * t339;
t117 = -pkin(4) * t145 - t385;
t50 = pkin(10) * t102 + t117 + t405;
t308 = qJD(6) * t135 + t238 + t50;
t167 = pkin(10) + t342;
t307 = qJD(6) * t167 + t50;
t299 = g(1) * t257 - g(2) * t262;
t297 = -t135 * t37 + t369;
t296 = -t167 * t37 + t369;
t295 = -t74 * t145 - t291 * t73;
t289 = -t102 * t359 - t334 * t395 + t35;
t85 = -pkin(4) * t111 + t318;
t288 = -t316 + t211;
t285 = -0.2e1 * pkin(1) * t331 - pkin(7) * qJDD(2);
t279 = -pkin(10) * t37 + t32 * t395 + t369;
t263 = qJD(2) ^ 2;
t278 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t263 + t299;
t264 = qJD(1) ^ 2;
t277 = pkin(1) * t264 - pkin(7) * qJDD(1) + t300;
t275 = t258 * t324 + t403;
t271 = -t253 * t408 + t402;
t268 = g(3) * t242 - t208 * t176 + t243 * t300 - t388;
t265 = -t208 * t178 + t274 + t404;
t246 = -qJ(4) - t386;
t194 = pkin(1) + t341;
t166 = -pkin(5) - t287;
t164 = t225 * t346 + t349;
t163 = -t225 * t348 + t347;
t162 = -t225 * t347 + t348;
t161 = t225 * t349 + t346;
t134 = -pkin(5) - t292;
t133 = -t176 ^ 2 + t178 ^ 2;
t114 = t117 + t238;
t110 = -t298 + (-qJD(1) * t190 - t178) * t247;
t109 = t176 * t247 + t129;
t52 = qJD(5) * t107 - t259 * t111 + t254 * t112;
t14 = pkin(5) * t52 - pkin(10) * t51 + t85;
t9 = t289 + t370;
t7 = t258 * t8;
t6 = qJD(5) * t47 + t254 * t45 - t259 * t44;
t2 = [qJDD(1), t299, t300, qJDD(1) * t248 + 0.2e1 * t256 * t319, 0.2e1 * t256 * t329 - 0.2e1 * t331 * t340, qJDD(2) * t256 + t261 * t263, qJDD(2) * t261 - t256 * t263, 0, t256 * t285 + t261 * t278, -t256 * t278 + t261 * t285, t129 * t190 + t152 * t178, -t129 * t189 - t130 * t190 + t152 * t176 + t153 * t178, -t152 * t247 + t190 * t245, -t153 * t247 - t189 * t245, 0, -t130 * t377 - t208 * t153 + t174 * t189 + t176 * t239 + t243 * t299 + t245 * t303 + t247 * t273, -t129 * t377 + t208 * t152 + t174 * t190 - t178 * t239 - t242 * t299 - t245 * t343 - t247 * t281, t74 * t111 - t73 * t112 + t145 * t53 + t26 * t149 - t25 * t150 - t291 * t54 + t97 * t80 - t96 * t81 - t300, t26 * t97 + t74 * t54 + t25 * t96 + t73 * t53 + t276 * t301 + t155 * t318 - g(1) * (-t194 * t257 - t246 * t262) - g(2) * (t194 * t262 - t246 * t257) t107 * t38 + t394 * t51, -t102 * t51 - t107 * t39 + t293 * t38 - t394 * t52, t107 * t240 + t241 * t51, t240 * t293 - t241 * t52, 0, t102 * t85 + t113 * t52 + t123 * t39 + t225 * t299 - t240 * t46 - t241 * t6 - t293 * t56, t107 * t56 + t113 * t51 + t123 * t38 - t224 * t299 - t240 * t47 - t241 * t5 + t394 * t85, t107 * t392 + t51 * t367, -t294 * t51 + (-t19 - t22 * t258 + (t253 * t92 - t367) * qJD(6)) * t107, t107 * t35 - t21 * t293 + t94 * t52 + (-t107 * t334 + t258 * t51) * t395, -t107 * t34 + t22 * t293 - t92 * t52 + (-t107 * t333 - t253 * t51) * t395, -t293 * t37 + t395 * t52, -g(1) * t162 - g(2) * t164 - t7 * t293 + t12 * t52 + t46 * t22 + t6 * t92 + (t14 * t395 + t378 + (t293 * t31 - t395 * t47 + t368) * qJD(6)) * t258 + t387 * t253, -g(1) * t161 - g(2) * t163 - t13 * t52 + t46 * t21 + t6 * t94 + (-(-qJD(6) * t47 + t14) * t395 - t378 - t315 * t293 - qJD(6) * t368) * t253 + t387 * t258; 0, 0, 0, -t256 * t264 * t261, t340 * t264, t330, t329, qJDD(2), -g(3) * t261 + t256 * t277, g(3) * t256 + t261 * t277, -t353, t133, t109, t110, t245, -t304 * t247 + (-t176 * t339 + t245 * t260 - t247 * t337) * pkin(2) + t265, t344 * t247 + (t178 * t339 - t245 * t255 - t247 * t336) * pkin(2) + t268, t145 * t363 - t171 * t81 + t172 * t80 - t291 * t393 + t295, t26 * t172 + t25 * t171 - t155 * (t238 - t385) - g(3) * t341 + t393 * t74 + t363 * t73 - t300 * (-pkin(2) * t256 - pkin(3) * t242) t354, t27, t23, t24, t240, -t114 * t102 + t240 * t292 - t241 * t373 + t266, -t114 * t394 - t240 * t345 + t241 * t374 + t270, t11, t1, t10, t9, -t360, t134 * t22 + t373 * t92 + t297 * t253 + (t253 * t374 - t258 * t308) * t395 + t275, t134 * t21 + t373 * t94 + t297 * t258 + (t253 * t308 + t258 * t374) * t395 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, t133, t109, t110, t245, -t247 * t290 + t265, t247 * t305 + t268, t79 * t291 - t78 * t145 + (t251 * t80 - t252 * t81) * pkin(3) + t295, -t73 * t78 - t74 * t79 + (t155 * t178 + t25 * t252 + t251 * t26 + t404) * pkin(3), t354, t27, t23, t24, t240, -t117 * t102 + t240 * t287 - t241 * t364 + t266, -t117 * t394 - t240 * t342 + t241 * t365 + t270, t11, t1, t10, t9, -t360, t166 * t22 + t364 * t92 + t296 * t253 + (t253 * t365 - t258 * t307) * t395 + t275, t166 * t21 + t364 * t94 + t296 * t258 + (t253 * t307 + t258 * t365) * t395 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145 ^ 2 - t291 ^ 2, -t145 * t73 + t291 * t74 + t276 - t299, 0, 0, 0, 0, 0, t39 + t358, t38 - t356, 0, 0, 0, 0, 0, t289 - t370, -t395 ^ 2 * t258 - t34 - t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, t27, t23, t24, t240, t241 * t33 + t266, t241 * t32 + t270, t11, t1, t10, -t359 * t395 + t35 + t370, -t360, -pkin(5) * t22 - t33 * t92 + t279 * t253 + (t324 - t391) * t258 + t403, -pkin(5) * t21 - t33 * t94 + t279 * t258 + (-t408 + t391) * t253 + t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t92, -t92 ^ 2 + t94 ^ 2, t395 * t92 + t21, t395 * t94 - t22, t37, -g(1) * t163 + g(2) * t161 + t13 * t395 + t253 * t288 - t30 * t94 - t31 * t333 + t7, g(1) * t164 - g(2) * t162 + t12 * t395 + t253 * t315 + t258 * t288 + t30 * t92;];
tau_reg  = t2;
