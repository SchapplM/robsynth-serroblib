% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:04
% EndTime: 2019-03-09 22:04:16
% DurationCPUTime: 5.81s
% Computational Cost: add. (9559->493), mult. (22201->611), div. (0->0), fcn. (16729->14), ass. (0->269)
t219 = sin(qJ(2));
t222 = cos(qJ(2));
t376 = sin(qJ(3));
t311 = qJD(1) * t376;
t378 = cos(qJ(3));
t312 = qJD(1) * t378;
t145 = -t219 * t311 + t222 * t312;
t146 = -t219 * t312 - t222 * t311;
t218 = sin(qJ(4));
t377 = cos(qJ(4));
t112 = t377 * t145 + t146 * t218;
t217 = sin(qJ(6));
t273 = -t218 * t145 + t146 * t377;
t394 = qJD(6) - t273;
t221 = cos(qJ(6));
t405 = t394 * t221;
t155 = t219 * t378 + t222 * t376;
t212 = qJD(2) + qJD(3);
t396 = t212 * t155;
t233 = t396 * qJD(1);
t265 = t219 * t376 - t222 * t378;
t229 = -t265 * qJDD(1) - t233;
t234 = t212 * t265;
t230 = -t234 * qJD(1) + t155 * qJDD(1);
t309 = qJD(4) * t377;
t325 = qJD(4) * t218;
t270 = t145 * t309 + t146 * t325 + t218 * t229 + t377 * t230;
t44 = -qJDD(6) - t270;
t264 = t217 * t44 - t394 * t405;
t205 = qJD(4) + t212;
t96 = t112 * t221 + t205 * t217;
t15 = t112 * t96 + t264;
t40 = t221 * t44;
t406 = t394 * t217;
t275 = -t394 * t406 - t40;
t98 = -t112 * t217 + t205 * t221;
t14 = -t112 * t98 + t275;
t211 = qJDD(2) + qJDD(3);
t204 = qJDD(4) + t211;
t227 = qJD(4) * t273 - t230 * t218 + t377 * t229;
t323 = qJD(6) * t221;
t324 = qJD(6) * t217;
t24 = -t112 * t323 + t221 * t204 - t205 * t324 - t217 * t227;
t23 = t24 * t221;
t16 = -t406 * t98 + t23;
t300 = t204 * t217 + t221 * t227;
t25 = qJD(6) * t98 + t300;
t407 = t394 * t96;
t1 = -t405 * t98 + (-t24 + t407) * t217 - t221 * t25;
t401 = t112 * pkin(5);
t141 = t146 * pkin(9);
t380 = pkin(7) + pkin(8);
t174 = t380 * t222;
t162 = qJD(1) * t174;
t147 = t376 * t162;
t173 = t380 * t219;
t160 = qJD(1) * t173;
t364 = qJD(2) * pkin(2);
t153 = -t160 + t364;
t284 = t378 * t153 - t147;
t94 = t141 + t284;
t86 = t212 * pkin(3) + t94;
t151 = t378 * t162;
t266 = -t153 * t376 - t151;
t367 = t145 * pkin(9);
t95 = -t266 + t367;
t90 = t377 * t95;
t61 = t218 * t86 + t90;
t56 = -qJ(5) * t205 - t61;
t31 = -t56 + t401;
t299 = t394 * t31;
t352 = qJ(5) * t112;
t78 = -pkin(4) * t273 - t352;
t382 = t273 ^ 2;
t37 = -t112 ^ 2 + t382;
t369 = t273 * pkin(5);
t210 = t222 * pkin(2);
t365 = pkin(1) + t210;
t398 = t273 * t31;
t285 = t160 * t376 - t151;
t100 = t285 - t367;
t331 = -t378 * t160 - t147;
t101 = t141 + t331;
t200 = pkin(2) * t378 + pkin(3);
t308 = t376 * qJD(3);
t310 = qJD(3) * t378;
t375 = pkin(2) * t218;
t397 = t218 * t100 + (qJD(4) * t376 + t308) * t375 - t200 * t309 + (-pkin(2) * t310 + t101) * t377;
t350 = t394 * t112;
t395 = t273 * t112;
t89 = t218 * t95;
t60 = -t377 * t86 + t89;
t332 = qJD(5) + t60;
t216 = qJ(2) + qJ(3);
t208 = cos(t216);
t209 = qJ(4) + t216;
t197 = sin(t209);
t198 = cos(t209);
t328 = t198 * pkin(4) + t197 * qJ(5);
t316 = pkin(3) * t208 + t328;
t333 = t332 - t369;
t381 = pkin(4) + pkin(10);
t29 = -t205 * t381 + t333;
t172 = t365 * qJD(1);
t126 = -pkin(3) * t145 - t172;
t247 = qJ(5) * t273 + t126;
t49 = -t112 * t381 + t247;
t18 = t217 * t29 + t221 * t49;
t320 = qJD(1) * qJD(2);
t306 = t222 * t320;
t319 = t219 * qJDD(1);
t124 = qJDD(2) * pkin(2) + t380 * (-t306 - t319);
t307 = t219 * t320;
t318 = t222 * qJDD(1);
t125 = t380 * (-t307 + t318);
t241 = qJD(3) * t266 + t378 * t124 - t125 * t376;
t34 = t211 * pkin(3) - t230 * pkin(9) + t241;
t243 = t124 * t376 + t125 * t378 + t153 * t310 - t162 * t308;
t43 = t229 * pkin(9) + t243;
t301 = -t218 * t34 - t86 * t309 + t95 * t325 - t377 * t43;
t187 = t204 * qJ(5);
t389 = -t205 * qJD(5) - t187;
t8 = t301 + t389;
t6 = pkin(5) * t227 - t8;
t313 = t18 * t112 + t6 * t221;
t17 = -t217 * t49 + t221 * t29;
t393 = -t17 * t112 + t6 * t217 + t31 * t323;
t188 = g(3) * t197;
t223 = cos(qJ(1));
t342 = t198 * t223;
t220 = sin(qJ(1));
t343 = t198 * t220;
t271 = -g(1) * t342 - g(2) * t343 - t188 - t301;
t64 = -pkin(4) * t112 + t247;
t244 = t112 * t64 + t271 - t389;
t189 = g(3) * t198;
t302 = t218 * t43 + t95 * t309 + t86 * t325 - t377 * t34;
t344 = t197 * t223;
t345 = t197 * t220;
t272 = -g(1) * t344 - g(2) * t345 + t189 + t302;
t392 = -t273 * t64 + qJDD(5) + t272;
t254 = -t126 * t112 - t271;
t391 = -t112 * t205 + t270;
t255 = t126 * t273 - t272;
t28 = -t205 * t273 + t227;
t357 = -qJD(5) + t397;
t293 = t376 * t377;
t355 = t100 * t377 - t218 * t101 + t200 * t325 + (qJD(4) * t293 + (t218 * t378 + t293) * qJD(3)) * pkin(2);
t63 = t377 * t94 - t89;
t353 = pkin(3) * t309 + qJD(5) - t63;
t62 = t218 * t94 + t90;
t291 = pkin(3) * t325 - t62;
t390 = t355 * t205;
t129 = pkin(3) * t265 - t365;
t207 = sin(t216);
t290 = -pkin(3) * t207 - pkin(4) * t197;
t388 = -qJD(6) + t394;
t199 = -pkin(3) * t377 - pkin(4);
t193 = -pkin(10) + t199;
t387 = t394 * (-t401 + t291) - t193 * t44;
t286 = t200 * t377 - t376 * t375;
t139 = -pkin(4) - t286;
t132 = -pkin(10) + t139;
t386 = t394 * (-t401 + t355) - t132 * t44;
t283 = -t378 * t173 - t174 * t376;
t105 = -t155 * pkin(9) + t283;
t330 = -t376 * t173 + t378 * t174;
t106 = -pkin(9) * t265 + t330;
t314 = qJD(2) * t380;
t161 = t219 * t314;
t163 = t222 * t314;
t261 = -t378 * t161 - t376 * t163 - t173 * t310 - t174 * t308;
t73 = -pkin(9) * t396 + t261;
t287 = t161 * t376 - t378 * t163;
t74 = pkin(9) * t234 + t173 * t308 - t174 * t310 + t287;
t20 = -t105 * t309 + t106 * t325 - t218 * t74 - t377 * t73;
t288 = g(1) * t220 - g(2) * t223;
t77 = t218 * t105 + t106 * t377;
t385 = t197 * t288 - t20 * t205 + t204 * t77;
t21 = qJD(4) * t77 + t218 * t73 - t377 * t74;
t76 = -t105 * t377 + t218 * t106;
t384 = -t198 * t288 + t204 * t76 + t205 * t21;
t374 = pkin(3) * t146;
t371 = pkin(4) * t204;
t248 = t377 * t265;
t122 = t155 * t218 + t248;
t258 = t218 * t265;
t123 = t155 * t377 - t258;
t242 = -t123 * qJ(5) + t129;
t54 = t122 * t381 + t242;
t366 = t54 * t44;
t361 = t205 * t61;
t359 = t217 * t98;
t358 = -t369 - t357;
t354 = -t369 + t353;
t351 = qJDD(1) * pkin(1);
t347 = t122 * t217;
t346 = t146 * t145;
t341 = t207 * t220;
t340 = t207 * t223;
t339 = t208 * t220;
t338 = t208 * t223;
t337 = t217 * t220;
t336 = t217 * t223;
t335 = t220 * t221;
t334 = t221 * t223;
t214 = t219 ^ 2;
t327 = -t222 ^ 2 + t214;
t326 = qJD(1) * t219;
t203 = t219 * t364;
t305 = -pkin(2) * t219 + t290;
t190 = pkin(2) * t307;
t10 = -pkin(2) * t318 - pkin(3) * t229 - pkin(4) * t227 - qJ(5) * t270 + qJD(5) * t273 + t190 - t351;
t7 = -pkin(10) * t227 + t10;
t304 = qJD(6) * t29 + t7;
t282 = qJDD(5) + t302;
t5 = pkin(5) * t270 - t204 * t381 + t282;
t303 = -qJD(6) * t49 + t5;
t289 = g(1) * t223 + g(2) * t220;
t281 = (t61 + t401) * t394 - t381 * t44;
t55 = -pkin(4) * t205 + t332;
t280 = -t112 * t55 + t273 * t56;
t277 = -t221 * t398 + t393;
t276 = t365 + t316;
t274 = -0.2e1 * pkin(1) * t320 - pkin(7) * qJDD(2);
t57 = t123 * pkin(5) + t76;
t69 = -qJD(4) * t258 + t155 * t309 - t218 * t234 + t377 * t396;
t267 = t6 * t122 + t31 * t69 + t57 * t44;
t70 = -t374 + t78;
t260 = pkin(2) * t293 + t218 * t200;
t202 = pkin(2) * t326;
t67 = t202 + t70;
t257 = -t198 * t289 - t188;
t225 = qJD(2) ^ 2;
t251 = -pkin(7) * t225 + t288 + 0.2e1 * t351;
t226 = qJD(1) ^ 2;
t250 = pkin(1) * t226 - pkin(7) * qJDD(1) + t289;
t107 = t273 * pkin(10);
t240 = (-qJD(6) * t132 - t107 + t67) * t394 + t257;
t239 = (-qJD(6) * t193 - t107 + t70) * t394 + t257;
t238 = (t394 * t381 - t352) * t394 + t257;
t237 = g(1) * t338 + g(2) * t339 + g(3) * t207 + t172 * t145 - t243;
t231 = g(1) * t340 + g(2) * t341 - g(3) * t208 - t172 * t146 + t241;
t115 = pkin(3) * t396 + t203;
t68 = qJD(4) * t248 + t155 * t325 + t218 * t396 + t234 * t377;
t19 = t69 * pkin(4) + t68 * qJ(5) - t123 * qJD(5) + t115;
t213 = -pkin(9) - t380;
t194 = pkin(3) * t218 + qJ(5);
t167 = qJ(5) * t342;
t166 = qJ(5) * t343;
t142 = -qJDD(1) * t365 + t190;
t138 = qJ(5) + t260;
t136 = -t197 * t337 + t334;
t135 = t197 * t335 + t336;
t134 = t197 * t336 + t335;
t133 = t197 * t334 - t337;
t127 = t202 - t374;
t99 = -t145 ^ 2 + t146 ^ 2;
t81 = pkin(3) * t233 + t129 * qJDD(1) + t190;
t80 = -t146 * t212 + t229;
t79 = -t145 * t212 + t230;
t75 = t122 * pkin(4) + t242;
t58 = -t122 * pkin(5) + t77;
t13 = -t68 * pkin(5) + t21;
t12 = -pkin(5) * t69 - t20;
t11 = t69 * pkin(10) + t19;
t9 = t282 - t371;
t2 = t221 * t5;
t3 = [qJDD(1), t288, t289, qJDD(1) * t214 + 0.2e1 * t219 * t306, 0.2e1 * t219 * t318 - 0.2e1 * t320 * t327, qJDD(2) * t219 + t222 * t225, qJDD(2) * t222 - t219 * t225, 0, t219 * t274 + t222 * t251, -t219 * t251 + t222 * t274, t146 * t234 + t230 * t155, -t145 * t234 + t146 * t396 + t155 * t229 - t230 * t265, t155 * t211 - t212 * t234, -t211 * t265 - t212 * t396, 0, -t145 * t203 + t365 * t229 + t142 * t265 - t172 * t396 + (-qJD(3) * t330 + t287) * t212 + t283 * t211 + g(1) * t339 - g(2) * t338, -g(1) * t341 + g(2) * t340 + t142 * t155 - t146 * t203 + t172 * t234 - t211 * t330 - t212 * t261 - t230 * t365, t123 * t270 + t273 * t68, -t112 * t68 - t122 * t270 + t123 * t227 + t273 * t69, t123 * t204 - t205 * t68, -t122 * t204 - t205 * t69, 0, -t112 * t115 + t122 * t81 + t126 * t69 - t129 * t227 - t384, -t115 * t273 + t123 * t81 - t126 * t68 + t129 * t270 - t385, -t112 * t20 + t122 * t8 + t123 * t9 - t21 * t273 + t227 * t77 + t270 * t76 - t55 * t68 + t56 * t69 - t289, -t10 * t122 + t112 * t19 + t227 * t75 - t64 * t69 + t384, -t10 * t123 + t19 * t273 - t270 * t75 + t64 * t68 + t385, t10 * t75 + t64 * t19 + t56 * t20 + t55 * t21 + t9 * t76 - t8 * t77 + (g(1) * t213 - g(2) * t276) * t223 + (g(1) * t276 + g(2) * t213) * t220, t69 * t359 + (t217 * t24 + t323 * t98) * t122 (-t217 * t96 + t221 * t98) * t69 + (-t217 * t25 + t23 + (-t221 * t96 - t359) * qJD(6)) * t122, -t44 * t347 + t123 * t24 - t68 * t98 + (t122 * t323 + t217 * t69) * t394, -t122 * t40 - t123 * t25 + t68 * t96 + (-t122 * t324 + t221 * t69) * t394, -t123 * t44 - t394 * t68, -g(1) * t136 - g(2) * t134 + t12 * t96 + t2 * t123 - t17 * t68 + t58 * t25 + (-t11 * t394 - t7 * t123 + t366) * t217 + (t13 * t394 - t267) * t221 + ((-t217 * t57 - t221 * t54) * t394 - t18 * t123 + t31 * t347) * qJD(6), g(1) * t135 - g(2) * t133 + t12 * t98 + t18 * t68 + t58 * t24 + (-(qJD(6) * t57 + t11) * t394 + t366 - t304 * t123 + t31 * qJD(6) * t122) * t221 + (-(-qJD(6) * t54 + t13) * t394 - t303 * t123 + t267) * t217; 0, 0, 0, -t219 * t226 * t222, t327 * t226, t319, t318, qJDD(2), -g(3) * t222 + t219 * t250, g(3) * t219 + t222 * t250, t346, t99, t79, t80, t211, -t285 * t212 + (t145 * t326 + t211 * t378 - t212 * t308) * pkin(2) + t231, t331 * t212 + (t146 * t326 - t211 * t376 - t212 * t310) * pkin(2) + t237, t395, t37, t391, t28, t204, t112 * t127 + t204 * t286 + t255 - t390, t127 * t273 - t260 * t204 + t397 * t205 + t254, -t112 * t357 + t138 * t227 + t139 * t270 - t273 * t355 + t280, -t112 * t67 + t390 + (-pkin(4) + t139) * t204 + t392, t138 * t204 - t205 * t357 - t273 * t67 + t244, -t8 * t138 + t9 * t139 - t64 * t67 - g(1) * (t223 * t305 + t167) - g(2) * (t220 * t305 + t166) - g(3) * (t210 + t316) + t357 * t56 + t355 * t55, t16, t1, t14, t15, -t350, t138 * t25 + t240 * t217 + t221 * t386 + t358 * t96 + t277, t138 * t24 + t358 * t98 + t240 * t221 + (-t299 - t386) * t217 + t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t346, t99, t79, t80, t211, -t212 * t266 + t231, t212 * t284 + t237, t395, t37, t391, t28, t204, t62 * t205 + (-t112 * t146 + t204 * t377 - t205 * t325) * pkin(3) + t255, t63 * t205 + (-t146 * t273 - t204 * t218 - t205 * t309) * pkin(3) + t254, t112 * t353 + t194 * t227 + t199 * t270 - t273 * t291 + t280, -t112 * t70 + t291 * t205 + (-pkin(4) + t199) * t204 + t392, t194 * t204 + t205 * t353 - t273 * t70 + t244, -t8 * t194 + t9 * t199 - t64 * t70 - g(1) * (t223 * t290 + t167) - g(2) * (t220 * t290 + t166) - g(3) * t316 - t353 * t56 + t291 * t55, t16, t1, t14, t15, -t350, t194 * t25 + t239 * t217 + t221 * t387 + t354 * t96 + t277, t194 * t24 + t354 * t98 + t239 * t221 + (-t299 - t387) * t217 + t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t395, t37, t391, t28, t204, t255 + t361, -t205 * t60 + t254, -pkin(4) * t270 + qJ(5) * t227 - (-t56 - t61) * t273 - (t55 - t332) * t112, -t112 * t78 - t361 - 0.2e1 * t371 + t392, t205 * t332 - t273 * t78 + t187 + t244, -t8 * qJ(5) - t9 * pkin(4) - t64 * t78 - t55 * t61 - g(1) * (-pkin(4) * t344 + t167) - g(2) * (-pkin(4) * t345 + t166) - g(3) * t328 - t332 * t56, t16, t1, t14, t15, -t350, qJ(5) * t25 + t333 * t96 + (-t281 - t398) * t221 + t238 * t217 + t393, qJ(5) * t24 + t333 * t98 + (t281 - t299) * t217 + t238 * t221 + t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391, t204 - t395, -t205 ^ 2 - t382, t205 * t56 - t371 + t392, 0, 0, 0, 0, 0, -t205 * t96 + t275, -t205 * t98 + t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t96, -t96 ^ 2 + t98 ^ 2, t24 + t407, t388 * t98 - t300, -t44, -g(1) * t133 - g(2) * t135 + t18 * t388 + t189 * t221 - t217 * t7 - t31 * t98 + t2, g(1) * t134 - g(2) * t136 + t394 * t17 + t31 * t96 - t304 * t221 + (-t303 - t189) * t217;];
tau_reg  = t3;
