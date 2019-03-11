% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:43
% EndTime: 2019-03-09 11:56:00
% DurationCPUTime: 7.30s
% Computational Cost: add. (11279->541), mult. (25846->669), div. (0->0), fcn. (19439->14), ass. (0->274)
t218 = qJ(2) + pkin(10);
t214 = cos(t218);
t363 = g(3) * t214;
t213 = sin(t218);
t227 = sin(qJ(1));
t230 = cos(qJ(1));
t273 = g(1) * t230 + g(2) * t227;
t402 = t213 * t273;
t247 = -t363 + t402;
t225 = sin(qJ(4));
t222 = sin(pkin(10));
t370 = pkin(2) * t222;
t205 = pkin(8) + t370;
t357 = pkin(9) + t205;
t284 = qJD(4) * t357;
t229 = cos(qJ(2));
t342 = cos(pkin(10));
t282 = t342 * t229;
t200 = qJD(1) * t282;
t226 = sin(qJ(2));
t314 = qJD(1) * t226;
t167 = -t222 * t314 + t200;
t336 = t167 * t225;
t223 = -qJ(3) - pkin(7);
t193 = t223 * t229;
t186 = qJD(1) * t193;
t172 = t222 * t186;
t192 = t223 * t226;
t185 = qJD(1) * t192;
t128 = t185 * t342 + t172;
t228 = cos(qJ(4));
t181 = t222 * t229 + t226 * t342;
t169 = t181 * qJD(1);
t98 = pkin(2) * t314 + pkin(3) * t169 - pkin(8) * t167;
t346 = t128 * t228 + t225 * t98;
t401 = -pkin(9) * t336 + t225 * t284 + t346;
t335 = t167 * t228;
t91 = t228 * t98;
t400 = pkin(4) * t169 - pkin(9) * t335 - t128 * t225 + t228 * t284 + t91;
t168 = t181 * qJD(2);
t308 = t226 * qJDD(1);
t270 = qJDD(1) * t282 - t222 * t308;
t118 = qJD(1) * t168 + qJDD(4) - t270;
t113 = qJDD(5) + t118;
t157 = qJD(4) - t167;
t154 = qJD(5) + t157;
t224 = sin(qJ(5));
t371 = cos(qJ(5));
t299 = t371 * t225;
t184 = t224 * t228 + t299;
t325 = t224 * t225;
t256 = t228 * t371 - t325;
t379 = qJD(4) + qJD(5);
t291 = t371 * qJD(5);
t382 = qJD(4) * t371 + t291;
t345 = t167 * t256 - t228 * t382 + t325 * t379;
t278 = t184 * t113 - t154 * t345;
t137 = qJD(2) * t228 - t169 * t225;
t138 = qJD(2) * t225 + t169 * t228;
t258 = t137 * t224 + t138 * t371;
t349 = t169 * t258;
t399 = t278 + t349;
t131 = t379 * t184;
t344 = -t167 * t184 + t131;
t279 = t113 * t256 - t154 * t344;
t75 = -t137 * t371 + t138 * t224;
t350 = t169 * t75;
t398 = t279 - t350;
t313 = qJD(4) * t225;
t397 = t313 - t336;
t352 = qJD(2) * pkin(2);
t175 = t185 + t352;
t123 = t175 * t342 + t172;
t110 = -qJD(2) * pkin(3) - t123;
t73 = -pkin(4) * t137 + t110;
t31 = pkin(5) * t75 - qJ(6) * t258 + t73;
t396 = t31 * t75;
t395 = t73 * t75;
t358 = t258 * t75;
t252 = -t222 * t226 + t282;
t171 = t252 * qJD(2);
t312 = qJD(4) * t228;
t394 = t171 * t225 + t181 * t312;
t176 = t357 * t225;
t177 = t357 * t228;
t120 = -t176 * t224 + t177 * t371;
t221 = qJ(4) + qJ(5);
t215 = sin(t221);
t391 = t120 * t113 + t215 * t247;
t372 = t258 ^ 2;
t390 = -t75 ^ 2 + t372;
t310 = qJD(1) * qJD(2);
t290 = t226 * t310;
t244 = qJDD(1) * t181 - t222 * t290;
t238 = qJD(2) * t200 + t244;
t236 = t225 * qJDD(2) + t228 * t238;
t235 = qJD(4) * t137 + t236;
t309 = qJD(2) * qJD(4);
t300 = t169 * t312 + (t238 + t309) * t225;
t265 = t228 * qJDD(2) - t300;
t311 = qJD(5) * t224;
t27 = -t137 * t291 + t138 * t311 - t224 * t265 - t235 * t371;
t14 = t154 * t75 - t27;
t47 = pkin(5) * t258 + qJ(6) * t75;
t285 = qJD(2) * t223;
t166 = -qJD(3) * t226 + t229 * t285;
t117 = qJDD(2) * pkin(2) + qJD(1) * t166 + qJDD(1) * t192;
t165 = qJD(3) * t229 + t226 * t285;
t126 = qJD(1) * t165 - qJDD(1) * t193;
t68 = t117 * t342 - t126 * t222;
t64 = -qJDD(2) * pkin(3) - t68;
t388 = qJD(4) * t205 * t157 + t64;
t257 = -t176 * t371 - t177 * t224;
t387 = -t257 * qJD(5) + t224 * t400 + t371 * t401;
t386 = -t120 * qJD(5) + t224 * t401 - t371 * t400;
t283 = t342 * t186;
t127 = t185 * t222 - t283;
t277 = pkin(4) * t397 - t127;
t324 = t225 * t118;
t102 = t113 * qJ(6);
t146 = t154 * qJD(6);
t384 = t102 + t146;
t383 = g(1) * t227 - g(2) * t230;
t296 = t181 * t313;
t333 = t171 * t228;
t381 = -t296 + t333;
t319 = t228 * t230;
t323 = t225 * t227;
t158 = t214 * t323 + t319;
t320 = t227 * t228;
t322 = t225 * t230;
t160 = -t214 * t322 + t320;
t380 = -g(1) * t160 + g(2) * t158;
t104 = t113 * pkin(5);
t378 = t104 - qJDD(6);
t28 = t137 * t311 + t138 * t291 + t224 * t235 - t265 * t371;
t377 = t154 * t258 - t28;
t216 = cos(t221);
t326 = t216 * t230;
t327 = t215 * t227;
t142 = t214 * t327 + t326;
t318 = t230 * t215;
t321 = t227 * t216;
t144 = t214 * t318 - t321;
t124 = t175 * t222 - t283;
t111 = qJD(2) * pkin(8) + t124;
t125 = -qJD(2) * t169 + t270;
t306 = pkin(2) * t290 + qJDD(3);
t307 = t229 * qJDD(1);
t341 = qJDD(1) * pkin(1);
t57 = -pkin(2) * t307 - t125 * pkin(3) - pkin(8) * t238 + t306 - t341;
t69 = t117 * t222 + t126 * t342;
t65 = qJDD(2) * pkin(8) + t69;
t359 = t229 * pkin(2);
t212 = pkin(1) + t359;
t190 = -qJD(1) * t212 + qJD(3);
t87 = -pkin(3) * t167 - pkin(8) * t169 + t190;
t254 = -t111 * t313 + t225 * t57 + t228 * t65 + t312 * t87;
t12 = pkin(9) * t265 + t254;
t55 = -t111 * t225 + t228 * t87;
t45 = -pkin(9) * t138 + t55;
t37 = pkin(4) * t157 + t45;
t56 = t228 * t111 + t225 * t87;
t46 = pkin(9) * t137 + t56;
t271 = -t111 * t312 + t228 * t57;
t8 = t118 * pkin(4) - pkin(9) * t235 - t225 * t65 - t313 * t87 + t271;
t288 = t12 * t224 + t291 * t46 + t311 * t37 - t371 * t8;
t330 = t213 * t215;
t249 = g(1) * t144 + g(2) * t142 + g(3) * t330 - t288;
t240 = t258 * t31 - t249 - t378;
t376 = -t258 * t73 + t249;
t375 = -t157 ^ 2 * t228 - t324;
t374 = -t256 * t27 - t258 * t344;
t122 = -pkin(3) * t252 - pkin(8) * t181 - t212;
t136 = t192 * t222 - t193 * t342;
t129 = t228 * t136;
t303 = t226 * t352;
t99 = pkin(3) * t168 - pkin(8) * t171 + t303;
t92 = t228 * t99;
t97 = t165 * t342 + t166 * t222;
t24 = -pkin(9) * t333 + pkin(4) * t168 - t225 * t97 + t92 + (-t129 + (pkin(9) * t181 - t122) * t225) * qJD(4);
t109 = t228 * t122;
t331 = t181 * t228;
t51 = -pkin(4) * t252 - pkin(9) * t331 - t136 * t225 + t109;
t316 = t122 * t225 + t129;
t332 = t181 * t225;
t60 = -pkin(9) * t332 + t316;
t261 = t224 * t51 + t371 * t60;
t253 = t122 * t312 - t136 * t313 + t225 * t99 + t228 * t97;
t30 = -pkin(9) * t394 + t253;
t373 = -qJD(5) * t261 - t224 * t30 + t24 * t371;
t364 = g(3) * t213;
t362 = g(3) * t225;
t361 = g(3) * t229;
t360 = t228 * pkin(4);
t355 = pkin(5) * t344 + qJ(6) * t345 - qJD(6) * t184 + t277;
t354 = -qJ(6) * t169 - t387;
t353 = pkin(5) * t169 - t386;
t301 = t371 * t46;
t17 = t224 * t37 + t301;
t351 = t154 * t17;
t347 = t224 * t46;
t21 = t371 * t45 - t347;
t343 = pkin(4) * t291 + qJD(6) - t21;
t339 = t137 * t169;
t338 = t138 * t157;
t337 = t138 * t169;
t329 = t213 * t216;
t231 = -pkin(9) - pkin(8);
t328 = t213 * t231;
t103 = t228 * t118;
t16 = t37 * t371 - t347;
t317 = qJD(6) - t16;
t219 = t226 ^ 2;
t315 = -t229 ^ 2 + t219;
t297 = t342 * pkin(2);
t293 = t110 * t312;
t289 = t12 * t371 + t224 * t8 + t291 * t37 - t311 * t46;
t287 = pkin(4) * t225 - t223;
t281 = -qJD(4) * t87 - t65;
t96 = t165 * t222 - t166 * t342;
t135 = -t192 * t342 - t193 * t222;
t280 = -t184 * t28 + t345 * t75;
t20 = t224 * t45 + t301;
t276 = pkin(4) * t311 - t20;
t207 = -t297 - pkin(3);
t275 = -g(1) * t142 + g(2) * t144;
t143 = t214 * t321 - t318;
t145 = t214 * t326 + t327;
t274 = g(1) * t143 - g(2) * t145;
t95 = pkin(4) * t332 + t135;
t211 = pkin(3) + t360;
t269 = t211 * t214 - t328;
t63 = pkin(4) * t394 + t96;
t267 = -t157 * t397 + t103;
t266 = pkin(5) * t216 + qJ(6) * t215 + t211;
t263 = -t224 * t60 + t371 * t51;
t260 = -0.2e1 * pkin(1) * t310 - pkin(7) * qJDD(2);
t259 = t224 * t24 + t291 * t51 + t30 * t371 - t311 * t60;
t191 = t207 - t360;
t251 = -qJDD(1) * t212 + t306;
t250 = t113 * t257 + t216 * t247;
t248 = g(1) * t145 + g(2) * t143 + g(3) * t329 - t289;
t232 = qJD(2) ^ 2;
t246 = -pkin(7) * t232 + 0.2e1 * t341 + t383;
t233 = qJD(1) ^ 2;
t245 = pkin(1) * t233 - pkin(7) * qJDD(1) + t273;
t242 = t154 * t16 + t248;
t32 = -pkin(4) * t265 + t64;
t237 = -g(1) * (-pkin(5) * t144 + qJ(6) * t145) - g(2) * (-pkin(5) * t142 + qJ(6) * t143) - g(3) * (-pkin(5) * t330 + qJ(6) * t329);
t234 = t235 * t228;
t210 = -pkin(4) * t371 - pkin(5);
t206 = pkin(4) * t224 + qJ(6);
t195 = t230 * t212;
t161 = t214 * t319 + t323;
t159 = -t214 * t320 + t322;
t107 = -pkin(5) * t256 - qJ(6) * t184 + t191;
t106 = t256 * t181;
t105 = t184 * t181;
t40 = pkin(5) * t105 - qJ(6) * t106 + t95;
t39 = t171 * t299 - t224 * t296 - t311 * t332 + (t171 * t224 + t181 * t382) * t228;
t38 = t131 * t181 - t171 * t256;
t34 = pkin(4) * t138 + t47;
t26 = pkin(5) * t252 - t263;
t25 = -qJ(6) * t252 + t261;
t15 = qJ(6) * t154 + t17;
t13 = -pkin(5) * t154 + t317;
t9 = pkin(5) * t39 + qJ(6) * t38 - qJD(6) * t106 + t63;
t5 = t28 * pkin(5) + t27 * qJ(6) - qJD(6) * t258 + t32;
t4 = -pkin(5) * t168 - t373;
t3 = qJ(6) * t168 - qJD(6) * t252 + t259;
t2 = t288 - t378;
t1 = t289 + t384;
t6 = [qJDD(1), t383, t273, qJDD(1) * t219 + 0.2e1 * t229 * t290, 0.2e1 * t226 * t307 - 0.2e1 * t310 * t315, qJDD(2) * t226 + t229 * t232, qJDD(2) * t229 - t226 * t232, 0, t226 * t260 + t229 * t246, -t226 * t246 + t229 * t260, -t123 * t171 - t124 * t168 + t136 * t125 + t135 * t238 + t97 * t167 + t96 * t169 - t68 * t181 + t252 * t69 - t273, t69 * t136 + t124 * t97 - t68 * t135 - t123 * t96 - t251 * t212 + t190 * t303 - g(1) * (-t212 * t227 - t223 * t230) - g(2) * (-t223 * t227 + t195) t138 * t381 + t181 * t234, t137 * t381 - t138 * t394 - t235 * t332 + t265 * t331, t103 * t181 + t138 * t168 + t157 * t381 - t235 * t252, t137 * t168 - t157 * t394 - t181 * t324 - t252 * t265, -t118 * t252 + t157 * t168 (-t136 * t312 + t92) * t157 + t109 * t118 - t271 * t252 + t55 * t168 - t96 * t137 - t135 * t265 + t181 * t293 - g(1) * t159 - g(2) * t161 + ((-qJD(4) * t122 - t97) * t157 - t136 * t118 - t281 * t252 + t64 * t181 + t110 * t171) * t225, -g(1) * t158 - g(2) * t160 + t110 * t381 - t118 * t316 + t135 * t235 + t96 * t138 - t157 * t253 - t56 * t168 + t252 * t254 + t331 * t64, -t106 * t27 - t258 * t38, t105 * t27 - t106 * t28 - t258 * t39 + t38 * t75, t106 * t113 - t154 * t38 + t168 * t258 + t252 * t27, -t105 * t113 - t154 * t39 - t168 * t75 + t252 * t28, -t113 * t252 + t154 * t168, t32 * t105 + t113 * t263 + t154 * t373 + t16 * t168 + t252 * t288 + t95 * t28 + t73 * t39 + t63 * t75 + t274, t106 * t32 - t113 * t261 - t154 * t259 - t168 * t17 + t252 * t289 + t258 * t63 - t27 * t95 - t38 * t73 + t275, t105 * t5 - t113 * t26 - t13 * t168 - t154 * t4 + t2 * t252 + t28 * t40 + t31 * t39 + t75 * t9 + t274, -t1 * t105 + t106 * t2 - t13 * t38 - t15 * t39 + t213 * t383 - t25 * t28 + t258 * t4 - t26 * t27 - t3 * t75, -t1 * t252 - t106 * t5 + t113 * t25 + t15 * t168 + t154 * t3 - t258 * t9 + t27 * t40 + t31 * t38 - t275, t1 * t25 + t15 * t3 + t5 * t40 + t31 * t9 + t2 * t26 + t13 * t4 - g(1) * (-pkin(5) * t143 - qJ(6) * t142) - g(2) * (pkin(5) * t145 + qJ(6) * t144 + t195) + (-g(1) * t287 - g(2) * t269) * t230 + (-g(1) * (-t212 - t269) - g(2) * t287) * t227; 0, 0, 0, -t226 * t233 * t229, t315 * t233, t308, t307, qJDD(2), t226 * t245 - t361, g(3) * t226 + t229 * t245, t125 * t370 - t238 * t297 - (-t124 + t127) * t169 + (-t128 + t123) * t167, t123 * t127 - t124 * t128 + (t342 * t68 - t361 + t222 * t69 + (-qJD(1) * t190 + t273) * t226) * pkin(2) (-qJD(4) * t169 + qJDD(2)) * t225 ^ 2 + (((t200 + qJD(4)) * qJD(2) + t244) * t225 + t338) * t228, t225 * t265 + t234 - t397 * t138 + (t312 - t335) * t137, -t337 - t375, t267 - t339, -t157 * t169, t207 * t300 - t55 * t169 + t127 * t137 - t205 * t324 + (-t207 * qJDD(2) + t247 - t388) * t228 + (-t91 + (t110 + t128) * t225) * t157, -t103 * t205 - t110 * t335 - t127 * t138 + t157 * t346 + t56 * t169 + t207 * t235 + t214 * t362 + t293 + (t388 - t402) * t225, -t184 * t27 - t258 * t345, t280 + t374, t278 - t349, t279 + t350, -t154 * t169, t154 * t386 - t16 * t169 + t191 * t28 - t256 * t32 + t277 * t75 + t344 * t73 + t250, t154 * t387 + t17 * t169 + t32 * t184 - t191 * t27 + t258 * t277 - t345 * t73 - t391, t107 * t28 + t13 * t169 - t154 * t353 - t256 * t5 + t31 * t344 + t355 * t75 + t250, t1 * t256 - t120 * t28 - t13 * t345 - t15 * t344 + t184 * t2 - t214 * t273 + t257 * t27 + t258 * t353 - t354 * t75 - t364, t107 * t27 - t15 * t169 + t154 * t354 - t184 * t5 - t258 * t355 + t31 * t345 + t391, t1 * t120 + t5 * t107 - t2 * t257 - g(3) * (-t328 + t359) + t355 * t31 - t266 * t363 + t354 * t15 + t353 * t13 + t273 * (pkin(2) * t226 + t213 * t266 + t214 * t231); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167 ^ 2 - t169 ^ 2, t123 * t169 - t124 * t167 + t251 - t383, 0, 0, 0, 0, 0, t267 + t339, -t337 + t375, 0, 0, 0, 0, 0, t398, -t399, t398, t280 - t374, t399, t1 * t184 + t13 * t344 - t15 * t345 - t169 * t31 - t2 * t256 - t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138 * t137, -t137 ^ 2 + t138 ^ 2, -t137 * t157 - t169 * t313 + t228 * t309 + t236, t265 + t338, t118, -t110 * t138 + t157 * t56 + (t281 + t364) * t225 + t271 + t380, g(1) * t161 - g(2) * t159 - t110 * t137 + t157 * t55 + t228 * t364 - t254, t358, t390, t14, t377, t113, t20 * t154 + (t113 * t371 - t138 * t75 - t154 * t311) * pkin(4) + t376, t21 * t154 + t395 + (-t113 * t224 - t138 * t258 - t154 * t291) * pkin(4) + t248, -t113 * t210 - t154 * t276 - t34 * t75 - t240, -t206 * t28 - t210 * t27 + (t15 + t276) * t258 + (t13 - t343) * t75, t113 * t206 + t154 * t343 + t258 * t34 - t248 + t384 - t396, t1 * t206 + t2 * t210 - t31 * t34 - t13 * t20 + t343 * t15 + (t13 * t311 + t213 * t362 + t380) * pkin(4) + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t358, t390, t14, t377, t113, t351 + t376, t242 + t395, -t47 * t75 + t104 - t240 + t351, pkin(5) * t27 - qJ(6) * t28 + (t15 - t17) * t258 + (t13 - t317) * t75, t258 * t47 + 0.2e1 * t102 + 0.2e1 * t146 - t242 - t396, -t2 * pkin(5) + t1 * qJ(6) - t13 * t17 + t15 * t317 - t31 * t47 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - qJDD(5) + t125 + t358, t14, -t154 ^ 2 - t372, -t15 * t154 + t240;];
tau_reg  = t6;
