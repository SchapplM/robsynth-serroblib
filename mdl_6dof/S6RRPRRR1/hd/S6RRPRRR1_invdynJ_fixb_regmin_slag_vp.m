% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:33
% EndTime: 2019-03-09 13:14:46
% DurationCPUTime: 5.78s
% Computational Cost: add. (12716->426), mult. (31475->554), div. (0->0), fcn. (25747->16), ass. (0->247)
t235 = cos(qJ(6));
t303 = qJD(6) * t235;
t231 = sin(qJ(5));
t236 = cos(qJ(5));
t227 = sin(pkin(11));
t228 = cos(pkin(11));
t233 = sin(qJ(2));
t238 = cos(qJ(2));
t181 = -t227 * t233 + t228 * t238;
t169 = t181 * qJD(1);
t182 = t227 * t238 + t228 * t233;
t171 = t182 * qJD(1);
t232 = sin(qJ(4));
t237 = cos(qJ(4));
t263 = t169 * t232 + t237 * t171;
t264 = t169 * t237 - t232 * t171;
t76 = t231 * t263 - t236 * t264;
t391 = t235 * t76;
t396 = t303 + t391;
t223 = qJDD(2) + qJDD(4);
t220 = qJDD(5) + t223;
t224 = qJD(2) + qJD(4);
t221 = qJD(5) + t224;
t230 = sin(qJ(6));
t305 = qJD(5) * t236;
t306 = qJD(5) * t231;
t170 = t182 * qJD(2);
t134 = -qJD(1) * t170 + qJDD(1) * t181;
t302 = qJD(1) * qJD(2);
t293 = t238 * t302;
t294 = t233 * t302;
t135 = qJDD(1) * t182 - t227 * t294 + t228 * t293;
t307 = qJD(4) * t237;
t308 = qJD(4) * t232;
t60 = t232 * t134 + t237 * t135 + t169 * t307 - t171 * t308;
t61 = qJD(4) * t263 - t237 * t134 + t232 * t135;
t27 = -t231 * t61 + t236 * t60 - t263 * t306 + t264 * t305;
t304 = qJD(6) * t230;
t349 = t231 * t264 + t236 * t263;
t17 = t230 * t220 + t221 * t303 + t235 * t27 - t304 * t349;
t69 = t221 * t230 + t235 * t349;
t18 = qJD(6) * t69 - t235 * t220 + t230 * t27;
t67 = -t235 * t221 + t230 * t349;
t387 = t17 * t235 - t230 * t18 - t396 * t67;
t378 = qJD(6) + t76;
t393 = t230 * t378;
t395 = -t393 * t69 + t387;
t28 = qJD(5) * t349 + t231 * t60 + t236 * t61;
t26 = qJDD(6) + t28;
t22 = t235 * t26;
t338 = t67 * t349;
t394 = -t378 * t393 + t22 + t338;
t330 = t221 * t76;
t374 = t27 + t330;
t222 = qJ(2) + pkin(11) + qJ(4);
t214 = qJ(5) + t222;
t208 = sin(t214);
t234 = sin(qJ(1));
t239 = cos(qJ(1));
t270 = g(1) * t239 + g(2) * t234;
t392 = t270 * t208;
t15 = t17 * t230;
t386 = t396 * t69 + t15;
t21 = t230 * t26;
t337 = t69 * t349;
t72 = t378 * t303;
t385 = t378 * t391 + t21 - t337 + t72;
t336 = qJ(3) + pkin(7);
t199 = t336 * t233;
t186 = qJD(1) * t199;
t331 = qJD(2) * pkin(2);
t178 = -t186 + t331;
t200 = t336 * t238;
t187 = qJD(1) * t200;
t319 = t228 * t187;
t133 = t227 * t178 + t319;
t346 = pkin(8) * t169;
t103 = t133 + t346;
t174 = t227 * t187;
t132 = t228 * t178 - t174;
t345 = pkin(8) * t171;
t97 = qJD(2) * pkin(3) + t132 - t345;
t267 = -t103 * t237 - t232 * t97;
t368 = pkin(9) * t264;
t51 = -t267 + t368;
t328 = t231 * t51;
t281 = -t103 * t232 + t237 * t97;
t369 = pkin(9) * t263;
t50 = t281 - t369;
t48 = pkin(4) * t224 + t50;
t29 = t236 * t48 - t328;
t24 = -pkin(5) * t221 - t29;
t384 = t24 * t76;
t381 = t349 * t76;
t323 = t349 * t221;
t352 = -t28 + t323;
t139 = t186 * t227 - t319;
t104 = t139 - t346;
t140 = -t228 * t186 - t174;
t105 = t140 - t345;
t213 = pkin(2) * t228 + pkin(3);
t347 = pkin(2) * t227;
t272 = t237 * t213 - t232 * t347;
t390 = -t272 * qJD(4) + t232 * t104 + t237 * t105;
t165 = t213 * t232 + t237 * t347;
t389 = -t165 * qJD(4) - t237 * t104 + t105 * t232;
t375 = t349 ^ 2 - t76 ^ 2;
t370 = pkin(5) * t349;
t388 = pkin(10) * t76 + t370;
t203 = g(3) * t208;
t209 = cos(t214);
t284 = qJD(2) * t336;
t167 = -t233 * qJD(3) - t238 * t284;
t131 = qJDD(2) * pkin(2) + qJD(1) * t167 - qJDD(1) * t199;
t166 = t238 * qJD(3) - t233 * t284;
t138 = qJD(1) * t166 + qJDD(1) * t200;
t82 = t228 * t131 - t138 * t227;
t62 = qJDD(2) * pkin(3) - pkin(8) * t135 + t82;
t83 = t227 * t131 + t228 * t138;
t63 = pkin(8) * t134 + t83;
t250 = qJD(4) * t267 - t232 * t63 + t237 * t62;
t12 = t223 * pkin(4) - t60 * pkin(9) + t250;
t350 = (qJD(4) * t97 + t63) * t237 - t103 * t308 + t232 * t62;
t13 = -t61 * pkin(9) + t350;
t351 = (qJD(5) * t48 + t13) * t236 + t231 * t12 - t51 * t306;
t217 = pkin(2) * t238 + pkin(1);
t193 = -qJD(1) * t217 + qJD(3);
t141 = -t169 * pkin(3) + t193;
t88 = -pkin(4) * t264 + t141;
t373 = t270 * t209 + t88 * t76 + t203 - t351;
t324 = t236 * t51;
t30 = t231 * t48 + t324;
t355 = qJD(5) * t30 - t236 * t12 + t231 * t13;
t3 = -t220 * pkin(5) + t355;
t342 = g(3) * t209;
t383 = t3 + t342;
t380 = t369 - t390;
t379 = -t368 - t389;
t364 = t378 * t349;
t321 = t263 * t224;
t376 = -t61 + t321;
t25 = pkin(10) * t221 + t30;
t39 = t76 * pkin(5) - pkin(10) * t349 + t88;
t8 = -t230 * t25 + t235 * t39;
t359 = t235 * t392 + t24 * t304 - t8 * t349;
t9 = t230 * t39 + t235 * t25;
t358 = t383 * t230 + t24 * t303 + t9 * t349;
t357 = -t349 * t88 - t342 - t355 + t392;
t371 = pkin(4) * t263;
t215 = pkin(4) * t231 + pkin(10);
t367 = (qJD(6) * t215 + t371 + t388) * t378;
t164 = pkin(4) + t272;
t310 = t231 * t164 + t236 * t165;
t112 = pkin(10) + t310;
t144 = t233 * qJD(1) * pkin(2) + pkin(3) * t171;
t89 = t144 + t371;
t366 = (qJD(6) * t112 + t388 + t89) * t378;
t365 = (t378 * pkin(10) + t370) * t378;
t320 = t264 * t224;
t363 = t60 - t320;
t362 = t263 * t264;
t142 = -t228 * t199 - t200 * t227;
t118 = -pkin(8) * t182 + t142;
t143 = -t227 * t199 + t228 * t200;
t119 = pkin(8) * t181 + t143;
t311 = t232 * t118 + t237 * t119;
t361 = t263 ^ 2 - t264 ^ 2;
t211 = sin(t222);
t212 = cos(t222);
t354 = -g(3) * t212 - t141 * t263 + t270 * t211 + t250;
t353 = g(3) * t211 - t141 * t264 + t270 * t212 - t350;
t291 = t220 * pkin(10) + qJD(6) * t39 + t351;
t137 = t181 * t232 + t182 * t237;
t274 = t237 * t118 - t119 * t232;
t54 = -pkin(9) * t137 + t274;
t262 = t237 * t181 - t182 * t232;
t55 = pkin(9) * t262 + t311;
t38 = t231 * t54 + t236 * t55;
t116 = -t166 * t227 + t228 * t167;
t173 = t181 * qJD(2);
t93 = -pkin(8) * t173 + t116;
t117 = t228 * t166 + t227 * t167;
t94 = -pkin(8) * t170 + t117;
t257 = t118 * t307 - t119 * t308 + t232 * t93 + t237 * t94;
t85 = qJD(4) * t137 + t237 * t170 + t232 * t173;
t35 = -pkin(9) * t85 + t257;
t249 = -qJD(4) * t311 - t232 * t94 + t237 * t93;
t84 = qJD(4) * t262 - t232 * t170 + t237 * t173;
t36 = -t84 * pkin(9) + t249;
t37 = t231 * t55 - t236 * t54;
t4 = -qJD(5) * t37 + t231 * t36 + t236 * t35;
t86 = t137 * t231 - t236 * t262;
t41 = -qJD(5) * t86 - t231 * t85 + t236 * t84;
t87 = t137 * t236 + t231 * t262;
t148 = -pkin(3) * t181 - t217;
t98 = -pkin(4) * t262 + t148;
t44 = pkin(5) * t86 - pkin(10) * t87 + t98;
t348 = t24 * t41 - t38 * t26 - (qJD(6) * t44 + t4) * t378 - t291 * t86 + t3 * t87;
t341 = g(3) * t238;
t340 = t24 * t87;
t339 = t44 * t26;
t265 = t164 * t236 - t165 * t231;
t333 = -qJD(5) * t265 + t379 * t231 - t380 * t236;
t332 = qJD(5) * t310 + t380 * t231 + t379 * t236;
t329 = t230 * t69;
t318 = t230 * t234;
t317 = t230 * t239;
t316 = t234 * t235;
t315 = t235 * t239;
t225 = t233 ^ 2;
t309 = -t238 ^ 2 + t225;
t301 = t238 * qJDD(1);
t219 = t233 * t331;
t296 = t87 * t304;
t145 = pkin(3) * t170 + t219;
t254 = pkin(2) * t294 - qJDD(1) * t217 + qJDD(3);
t99 = -t134 * pkin(3) + t254;
t47 = t61 * pkin(4) + t99;
t7 = t28 * pkin(5) - t27 * pkin(10) + t47;
t290 = qJD(6) * t25 - t7;
t31 = t231 * t50 + t324;
t271 = pkin(4) * t306 - t31;
t269 = g(1) * t234 - g(2) * t239;
t268 = t26 * t87 + t378 * t41;
t261 = t22 - (t230 * t76 + t304) * t378;
t66 = pkin(4) * t85 + t145;
t260 = -t291 + t203;
t258 = -0.2e1 * pkin(1) * t302 - pkin(7) * qJDD(2);
t256 = -pkin(10) * t26 + t29 * t378 + t384;
t253 = -t112 * t26 + t333 * t378 + t384;
t240 = qJD(2) ^ 2;
t252 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t240 + t269;
t241 = qJD(1) ^ 2;
t251 = pkin(1) * t241 - pkin(7) * qJDD(1) + t270;
t32 = t236 * t50 - t328;
t246 = -t215 * t26 + t384 + (-pkin(4) * t305 + t32) * t378;
t216 = -pkin(4) * t236 - pkin(5);
t158 = t209 * t315 + t318;
t157 = -t209 * t317 + t316;
t156 = -t209 * t316 + t317;
t155 = t209 * t318 + t315;
t111 = -pkin(5) - t265;
t42 = qJD(5) * t87 + t231 * t84 + t236 * t85;
t10 = pkin(5) * t42 - pkin(10) * t41 + t66;
t6 = t235 * t7;
t5 = qJD(5) * t38 + t231 * t35 - t236 * t36;
t1 = [qJDD(1), t269, t270, qJDD(1) * t225 + 0.2e1 * t233 * t293, 0.2e1 * t233 * t301 - 0.2e1 * t309 * t302, qJDD(2) * t233 + t238 * t240, qJDD(2) * t238 - t233 * t240, 0, t233 * t258 + t238 * t252, -t233 * t252 + t238 * t258, -t116 * t171 + t117 * t169 - t132 * t173 - t133 * t170 + t134 * t143 - t135 * t142 + t181 * t83 - t182 * t82 - t270, t83 * t143 + t133 * t117 + t82 * t142 + t132 * t116 - t254 * t217 + t193 * t219 - g(1) * (-t217 * t234 + t239 * t336) - g(2) * (t217 * t239 + t234 * t336) t137 * t60 + t263 * t84, -t137 * t61 + t262 * t60 - t263 * t85 + t264 * t84, t137 * t223 + t224 * t84, t223 * t262 - t224 * t85, 0, t141 * t85 - t145 * t264 + t148 * t61 + t212 * t269 + t223 * t274 + t224 * t249 - t262 * t99, t99 * t137 + t141 * t84 + t145 * t263 + t148 * t60 - t269 * t211 - t311 * t223 - t257 * t224, t27 * t87 + t349 * t41, -t27 * t86 - t28 * t87 - t349 * t42 - t41 * t76, t220 * t87 + t221 * t41, -t220 * t86 - t221 * t42, 0, t209 * t269 - t37 * t220 - t5 * t221 + t98 * t28 + t88 * t42 + t47 * t86 + t66 * t76, -t208 * t269 - t38 * t220 - t4 * t221 + t98 * t27 + t349 * t66 + t88 * t41 + t47 * t87, -t69 * t296 + (t17 * t87 + t41 * t69) * t235 (-t235 * t67 - t329) * t41 + (-t15 - t18 * t235 + (t230 * t67 - t235 * t69) * qJD(6)) * t87, t17 * t86 + t235 * t268 - t296 * t378 + t69 * t42, -t18 * t86 - t230 * t268 - t67 * t42 - t72 * t87, t26 * t86 + t378 * t42, -g(1) * t156 - g(2) * t158 + t37 * t18 + t8 * t42 + t5 * t67 + t6 * t86 + (t10 * t378 + t339 + (-t25 * t86 - t378 * t38 + t340) * qJD(6)) * t235 + t348 * t230, -g(1) * t155 - g(2) * t157 + t37 * t17 - t9 * t42 + t5 * t69 + (-(-qJD(6) * t38 + t10) * t378 - t339 + t290 * t86 - qJD(6) * t340) * t230 + t348 * t235; 0, 0, 0, -t233 * t241 * t238, t309 * t241, t233 * qJDD(1), t301, qJDD(2), t233 * t251 - t341, g(3) * t233 + t238 * t251 (t133 + t139) * t171 + (t132 - t140) * t169 + (t134 * t227 - t135 * t228) * pkin(2), -t132 * t139 - t133 * t140 + (-t341 + t227 * t83 + t228 * t82 + (-qJD(1) * t193 + t270) * t233) * pkin(2), -t362, t361, t363, t376, t223, t144 * t264 + t272 * t223 + t389 * t224 + t354, -t144 * t263 - t165 * t223 + t390 * t224 + t353, t381, t375, t374, t352, t220, t265 * t220 - t332 * t221 - t89 * t76 + t357, -t310 * t220 + t333 * t221 - t349 * t89 + t373, t386, -t329 * t378 + t387, t385, t261 + t338, -t364, t111 * t18 + t332 * t67 + (-t383 - t366) * t235 + t253 * t230 + t359, t111 * t17 + t332 * t69 + t253 * t235 + (-t392 + t366) * t230 + t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169 ^ 2 - t171 ^ 2, t132 * t171 - t133 * t169 + t254 - t269, 0, 0, 0, 0, 0, t61 + t321, t60 + t320, 0, 0, 0, 0, 0, t28 + t323, t27 - t330, 0, 0, 0, 0, 0, t261 - t338, -t235 * t378 ^ 2 - t21 - t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t362, t361, t363, t376, t223, -t224 * t267 + t354, t224 * t281 + t353, t381, t375, t374, t352, t220, t31 * t221 + (t220 * t236 - t221 * t306 - t263 * t76) * pkin(4) + t357, t32 * t221 + (-t220 * t231 - t221 * t305 - t263 * t349) * pkin(4) + t373, t386, t395, t385, t394, -t364, t216 * t18 + t271 * t67 + (-t383 - t367) * t235 + t246 * t230 + t359, t216 * t17 + t271 * t69 + t246 * t235 + (-t392 + t367) * t230 + t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t381, t375, t374, t352, t220, t30 * t221 + t357, t29 * t221 + t373, t386, t395, t385, t394, -t364, -pkin(5) * t18 - t30 * t67 + t256 * t230 + (-t383 - t365) * t235 + t359, -pkin(5) * t17 - t30 * t69 + t256 * t235 + (-t392 + t365) * t230 + t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t67, -t67 ^ 2 + t69 ^ 2, t378 * t67 + t17, t378 * t69 - t18, t26, -g(1) * t157 + g(2) * t155 + t230 * t260 - t24 * t69 - t25 * t303 + t378 * t9 + t6, g(1) * t158 - g(2) * t156 + t230 * t290 + t235 * t260 + t24 * t67 + t378 * t8;];
tau_reg  = t1;
