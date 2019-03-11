% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:17
% EndTime: 2019-03-09 11:45:31
% DurationCPUTime: 6.63s
% Computational Cost: add. (12246->527), mult. (29125->646), div. (0->0), fcn. (22460->14), ass. (0->273)
t227 = cos(qJ(5));
t307 = qJD(5) * t227;
t220 = sin(pkin(10));
t221 = cos(pkin(10));
t225 = sin(qJ(2));
t228 = cos(qJ(2));
t164 = -t220 * t225 + t221 * t228;
t154 = t164 * qJD(1);
t165 = t220 * t228 + t221 * t225;
t156 = t165 * qJD(1);
t224 = sin(qJ(4));
t245 = t224 * t156;
t367 = cos(qJ(4));
t374 = t154 * t367 - t245;
t393 = t374 * t227;
t402 = t307 - t393;
t214 = t228 * pkin(2);
t209 = t214 + pkin(1);
t133 = -pkin(3) * t164 - t209;
t107 = qJD(5) - t374;
t223 = sin(qJ(5));
t304 = qJD(1) * qJD(2);
t284 = t228 * t304;
t285 = t225 * t304;
t121 = qJDD(1) * t165 - t220 * t285 + t221 * t284;
t155 = t165 * qJD(2);
t244 = qJD(1) * t155;
t248 = t164 * qJDD(1);
t234 = t248 - t244;
t279 = t224 * t121 - t234 * t367;
t250 = -t154 * t224 - t156 * t367;
t389 = qJD(4) * t250;
t55 = t279 - t389;
t54 = qJDD(5) + t55;
t47 = t223 * t54;
t259 = t107 * t402 + t47;
t301 = qJD(2) + qJD(4);
t98 = t223 * t301 - t227 * t250;
t341 = t250 * t98;
t401 = t259 + t341;
t299 = qJDD(2) + qJDD(4);
t308 = qJD(5) * t223;
t291 = t367 * t121;
t388 = t374 * qJD(4);
t232 = t224 * t234 + t291 + t388;
t391 = qJD(5) * t301 + t232;
t39 = -t223 * t299 - t227 * t391 - t250 * t308;
t37 = t39 * t223;
t400 = t402 * t98 - t37;
t399 = pkin(5) * t250;
t96 = -t223 * t250 - t227 * t301;
t342 = t250 * t96;
t222 = -qJ(3) - pkin(7);
t190 = t222 * t228;
t175 = qJD(1) * t190;
t159 = t220 * t175;
t189 = t222 * t225;
t174 = qJD(1) * t189;
t345 = qJD(2) * pkin(2);
t163 = t174 + t345;
t119 = t163 * t221 + t159;
t363 = pkin(8) * t156;
t91 = qJD(2) * pkin(3) + t119 - t363;
t320 = t221 * t175;
t120 = t163 * t220 - t320;
t364 = pkin(8) * t154;
t95 = t120 + t364;
t58 = -t224 * t95 + t367 * t91;
t49 = -pkin(4) * t301 - t58;
t28 = pkin(5) * t96 - qJ(6) * t98 + t49;
t398 = t28 * t374;
t397 = t374 * t49;
t396 = qJ(6) * t250;
t395 = t107 * t250;
t394 = t250 * t374;
t217 = qJ(2) + pkin(10);
t213 = qJ(4) + t217;
t204 = sin(t213);
t229 = cos(qJ(1));
t323 = t204 * t229;
t226 = sin(qJ(1));
t324 = t204 * t226;
t392 = g(1) * t323 + g(2) * t324;
t306 = t250 * qJD(2);
t390 = -t306 - t279;
t268 = pkin(5) * t227 + qJ(6) * t223;
t387 = t250 ^ 2 - t374 ^ 2;
t59 = t224 * t91 + t367 * t95;
t50 = pkin(9) * t301 + t59;
t177 = -qJD(1) * t209 + qJD(3);
t127 = -pkin(3) * t154 + t177;
t64 = -pkin(4) * t374 + pkin(9) * t250 + t127;
t25 = -t223 * t50 + t227 * t64;
t314 = qJD(6) - t25;
t18 = -pkin(5) * t107 + t314;
t27 = t28 * t308;
t313 = t392 * t227;
t386 = -t18 * t250 + t27 + t313;
t26 = t223 * t64 + t227 * t50;
t19 = qJ(6) * t107 + t26;
t205 = cos(t213);
t358 = g(3) * t223;
t295 = -t205 * t358 + t223 * t392;
t286 = qJD(4) * t367;
t309 = qJD(4) * t224;
t281 = qJD(2) * t222;
t151 = -qJD(3) * t225 + t228 * t281;
t118 = qJDD(2) * pkin(2) + qJD(1) * t151 + qJDD(1) * t189;
t150 = qJD(3) * t228 + t225 * t281;
t124 = qJD(1) * t150 - qJDD(1) * t190;
t79 = t118 * t221 - t124 * t220;
t57 = qJDD(2) * pkin(3) - pkin(8) * t121 + t79;
t80 = t118 * t220 + t124 * t221;
t62 = pkin(8) * t234 + t80;
t280 = t224 * t62 + t286 * t95 + t309 * t91 - t367 * t57;
t13 = -pkin(4) * t299 + t280;
t40 = t223 * t391 - t227 * t299 - t250 * t307;
t5 = pkin(5) * t40 + qJ(6) * t39 - qJD(6) * t98 + t13;
t356 = t223 * t5;
t385 = t19 * t250 + t295 - t356;
t44 = t49 * t308;
t384 = t25 * t250 + t313 + t44;
t383 = t13 * t223 - t250 * t26 + t307 * t49 - t295;
t240 = t224 * t57 + t286 * t91 - t309 * t95 + t367 * t62;
t196 = g(3) * t204;
t321 = t205 * t229;
t322 = t205 * t226;
t294 = g(1) * t321 + g(2) * t322 + t196;
t382 = -t127 * t374 - t240 + t294;
t359 = g(3) * t205;
t381 = t127 * t250 - t280 - t359 + t392;
t267 = pkin(5) * t223 - qJ(6) * t227;
t380 = pkin(5) * t308 - qJ(6) * t307 - qJD(6) * t223 - t267 * t374;
t334 = t227 * t96;
t337 = t223 * t98;
t263 = t334 + t337;
t275 = -t227 * t39 - t308 * t98;
t352 = -t223 * t40 - t307 * t96;
t379 = t263 * t374 + t275 + t352;
t78 = -pkin(4) * t250 - pkin(9) * t374;
t238 = t224 * t248 + t291;
t378 = (-t374 - t245) * qJD(2) + t238;
t12 = pkin(9) * t299 + t240;
t300 = pkin(2) * t285 + qJDD(3);
t302 = t228 * qJDD(1);
t328 = qJDD(1) * pkin(1);
t16 = -pkin(2) * t302 - pkin(3) * t234 + t55 * pkin(4) - pkin(9) * t232 + t300 - t328;
t255 = t12 * t227 + t16 * t223 + t307 * t64 - t308 * t50;
t346 = qJ(6) * t54;
t2 = qJD(6) * t107 + t255 + t346;
t283 = t12 * t223 - t16 * t227 + t307 * t50 + t308 * t64;
t369 = pkin(5) * t54;
t4 = qJDD(6) + t283 - t369;
t377 = t2 * t227 + t4 * t223;
t206 = pkin(2) * t221 + pkin(3);
t366 = pkin(2) * t220;
t270 = t206 * t367 - t224 * t366;
t137 = t270 * qJD(4);
t125 = -t174 * t220 + t320;
t100 = t125 - t364;
t126 = t174 * t221 + t159;
t101 = t126 - t363;
t67 = t100 * t224 + t101 * t367;
t376 = -t137 + t67;
t312 = t206 * t224 + t366 * t367;
t330 = qJD(4) * t312 + t100 * t367 - t101 * t224;
t48 = t227 * t54;
t375 = t107 * t308 - t48;
t373 = pkin(4) * t205 + pkin(9) * t204;
t372 = g(1) * t226 - g(2) * t229;
t272 = g(1) * t229 + g(2) * t226;
t371 = t98 ^ 2;
t370 = t107 ^ 2;
t368 = pkin(9) * t54;
t357 = g(3) * t228;
t355 = t225 * pkin(2);
t353 = t98 * t96;
t351 = t223 * t78 + t227 * t58;
t131 = pkin(3) * t156 + qJD(1) * t355;
t68 = t131 + t78;
t350 = t223 * t68 + t227 * t67;
t123 = t164 * t224 + t165 * t367;
t249 = t164 * t367 - t165 * t224;
t73 = -pkin(4) * t249 - pkin(9) * t123 + t133;
t129 = t189 * t221 + t190 * t220;
t105 = -pkin(8) * t165 + t129;
t130 = t189 * t220 - t190 * t221;
t106 = pkin(8) * t164 + t130;
t75 = t105 * t224 + t106 * t367;
t349 = t223 * t73 + t227 * t75;
t347 = pkin(9) * qJD(5);
t344 = t107 * t19;
t343 = t107 * t26;
t149 = pkin(9) + t312;
t340 = t149 * t54;
t158 = t164 * qJD(2);
t81 = qJD(4) * t249 - t155 * t224 + t158 * t367;
t339 = t223 * t81;
t338 = t223 * t96;
t336 = t227 * t40;
t335 = t227 * t81;
t333 = t227 * t98;
t331 = t330 + t380;
t329 = t380 - t59;
t327 = t107 * t223;
t318 = t223 * t226;
t317 = t226 * t227;
t316 = t227 * t229;
t315 = t229 * t223;
t104 = t150 * t221 + t151 * t220;
t311 = pkin(3) * cos(t217) + t214;
t218 = t225 ^ 2;
t310 = -t228 ^ 2 + t218;
t303 = t225 * qJDD(1);
t212 = t225 * t345;
t292 = -t5 - t359;
t288 = t149 * t308;
t287 = -t13 - t359;
t132 = pkin(3) * t155 + t212;
t103 = -t150 * t220 + t151 * t221;
t276 = t205 * t268 + t373;
t144 = t205 * t318 + t316;
t146 = t205 * t315 - t317;
t274 = -g(1) * t144 + g(2) * t146;
t145 = t205 * t317 - t315;
t147 = t205 * t316 + t318;
t273 = g(1) * t145 - g(2) * t147;
t266 = -t340 - t397;
t265 = t18 * t227 - t19 * t223;
t264 = t18 * t223 + t19 * t227;
t262 = t18 * t307 - t294 + t377;
t260 = t327 * t374 - t375;
t258 = pkin(4) + t268;
t257 = pkin(1) + t311 + t373;
t256 = -0.2e1 * pkin(1) * t304 - pkin(7) * qJDD(2);
t252 = t105 * t367 - t106 * t224;
t84 = -pkin(8) * t158 + t103;
t85 = -pkin(8) * t155 + t104;
t32 = qJD(4) * t252 + t224 * t84 + t367 * t85;
t82 = qJD(4) * t123 + t155 * t367 + t224 * t158;
t42 = pkin(4) * t82 - pkin(9) * t81 + t132;
t254 = t223 * t42 + t227 * t32 + t307 * t73 - t308 * t75;
t247 = -qJDD(1) * t209 + t300;
t243 = g(1) * t146 + g(2) * t144 + t204 * t358 - t283;
t230 = qJD(2) ^ 2;
t242 = -pkin(7) * t230 + 0.2e1 * t328 + t372;
t231 = qJD(1) ^ 2;
t241 = pkin(1) * t231 - pkin(7) * qJDD(1) + t272;
t239 = qJD(5) * t265 + t377;
t237 = t28 * t98 + qJDD(6) - t243;
t236 = -g(1) * t147 - g(2) * t145 - t196 * t227 + t255;
t33 = qJD(4) * t75 + t224 * t85 - t367 * t84;
t233 = t272 * t258 * t204;
t216 = -pkin(8) + t222;
t182 = pkin(9) * t321;
t180 = pkin(9) * t322;
t176 = -pkin(3) * sin(t217) - t355;
t148 = -pkin(4) - t270;
t128 = -t258 - t270;
t94 = pkin(3) * t244 + qJDD(1) * t133 + t300;
t69 = pkin(5) * t98 + qJ(6) * t96;
t43 = t123 * t267 - t252;
t30 = pkin(5) * t249 + t223 * t75 - t227 * t73;
t29 = -qJ(6) * t249 + t349;
t24 = t223 * t58 - t227 * t78 + t399;
t23 = t351 - t396;
t22 = t223 * t67 - t227 * t68 + t399;
t21 = t350 - t396;
t20 = t107 * t96 - t39;
t8 = t267 * t81 + (qJD(5) * t268 - qJD(6) * t227) * t123 + t33;
t7 = -pkin(5) * t82 + qJD(5) * t349 + t223 * t32 - t227 * t42;
t6 = qJ(6) * t82 - qJD(6) * t249 + t254;
t1 = [qJDD(1), t372, t272, qJDD(1) * t218 + 0.2e1 * t225 * t284, 0.2e1 * t225 * t302 - 0.2e1 * t304 * t310, qJDD(2) * t225 + t228 * t230, qJDD(2) * t228 - t225 * t230, 0, t225 * t256 + t228 * t242, -t225 * t242 + t228 * t256, -t103 * t156 + t104 * t154 - t119 * t158 - t120 * t155 - t129 * t121 + t130 * t234 + t80 * t164 - t79 * t165 - t272, t80 * t130 + t120 * t104 + t79 * t129 + t119 * t103 - t247 * t209 + t177 * t212 - g(1) * (-t209 * t226 - t222 * t229) - g(2) * (t209 * t229 - t222 * t226) t123 * t232 - t250 * t81, -t123 * t55 + t232 * t249 + t250 * t82 + t374 * t81, t123 * t299 + t301 * t81, t249 * t299 - t301 * t82, 0, t127 * t82 - t132 * t374 + t133 * t55 + t205 * t372 - t249 * t94 + t252 * t299 - t301 * t33, -g(1) * t324 + g(2) * t323 + t94 * t123 + t127 * t81 - t132 * t250 + t133 * t232 - t299 * t75 - t301 * t32, t123 * t275 + t333 * t81, -t263 * t81 + (t37 - t336 + (-t333 + t338) * qJD(5)) * t123, t123 * t48 + t249 * t39 + t82 * t98 + (-t123 * t308 + t335) * t107, -t123 * t47 + t249 * t40 - t82 * t96 + (-t123 * t307 - t339) * t107, t107 * t82 - t249 * t54, t283 * t249 + t25 * t82 + t33 * t96 - t252 * t40 + ((-qJD(5) * t75 + t42) * t107 + t73 * t54 + t49 * qJD(5) * t123) * t227 + ((-qJD(5) * t73 - t32) * t107 - t75 * t54 + t13 * t123 + t49 * t81) * t223 + t273, -t254 * t107 - t349 * t54 + t255 * t249 - t26 * t82 + t33 * t98 + t252 * t39 + t49 * t335 + (t13 * t227 - t44) * t123 + t274, t28 * t339 - t107 * t7 + t249 * t4 - t18 * t82 - t30 * t54 + t40 * t43 + t8 * t96 + (t28 * t307 + t356) * t123 + t273, -t29 * t40 - t30 * t39 - t6 * t96 + t7 * t98 + t265 * t81 + t372 * t204 + (-qJD(5) * t264 - t2 * t223 + t227 * t4) * t123, -t28 * t335 + t107 * t6 - t249 * t2 + t19 * t82 + t29 * t54 + t39 * t43 - t8 * t98 + (-t227 * t5 + t27) * t123 - t274, t2 * t29 + t19 * t6 + t5 * t43 + t28 * t8 + t4 * t30 + t18 * t7 - g(1) * (-pkin(5) * t145 - qJ(6) * t144) - g(2) * (pkin(5) * t147 + qJ(6) * t146) + (g(1) * t216 - g(2) * t257) * t229 + (g(1) * t257 + g(2) * t216) * t226; 0, 0, 0, -t225 * t231 * t228, t310 * t231, t303, t302, qJDD(2), t225 * t241 - t357, g(3) * t225 + t228 * t241 (t120 + t125) * t156 + (-t126 + t119) * t154 + (-t221 * t121 + ((-t284 - t303) * t220 + (-t285 + t302) * t221) * t220) * pkin(2), -t119 * t125 - t120 * t126 + (-t357 + t220 * t80 + t221 * t79 + (-qJD(1) * t177 + t272) * t225) * pkin(2), t394, t387, t378, t390, t299, t131 * t374 + t270 * t299 - t301 * t330 + t381, t131 * t250 - t299 * t312 + t301 * t376 + t382, t400, t379, t401, t260 - t342, t395, t148 * t40 + t330 * t96 + t287 * t227 + t266 * t223 + ((-qJD(5) * t149 - t68) * t227 + t376 * t223) * t107 + t384, -t148 * t39 + t330 * t98 + t266 * t227 + (-t137 * t227 + t288 + t350) * t107 + t383, t128 * t40 + t331 * t96 + t292 * t227 + (-t340 - t398) * t223 + (-t223 * t137 - t149 * t307 + t22) * t107 + t386, t21 * t96 - t22 * t98 + (-t374 * t18 - t137 * t96 + (qJD(5) * t98 - t40) * t149) * t227 + (t374 * t19 + t137 * t98 - t149 * t39 + (t149 * t96 - t19) * qJD(5)) * t223 + t262, t128 * t39 + t340 * t227 - t331 * t98 + (-t21 - t288 + (t137 - t28) * t227) * t107 + t385, t5 * t128 - t19 * t21 - t18 * t22 - g(1) * (t176 * t229 + t182) - g(2) * (t176 * t226 + t180) - g(3) * (t276 + t311) + t331 * t28 + t264 * t137 + t233 + t239 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154 ^ 2 - t156 ^ 2, t119 * t156 - t120 * t154 + t247 - t372, 0, 0, 0, 0, 0, t279 - t306 - 0.2e1 * t389 (t374 - t245) * qJD(2) + 0.2e1 * t388 + t238, 0, 0, 0, 0, 0, t260 + t342, -t227 * t370 + t341 - t47, -t223 * t370 + t342 + t48 (t334 - t337) * t374 - t275 + t352, t259 - t341, t250 * t28 + (-t4 + t344) * t227 + (t107 * t18 + t2) * t223 - t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t394, t387, t378, t390, t299, t301 * t59 + t381, t301 * t58 + t382, t400, t379, t401, -t107 * t327 - t342 + t48, t395, -pkin(4) * t40 - t59 * t96 + (t107 * t58 - t368 - t397) * t223 + ((-t78 - t347) * t107 + t287) * t227 + t384, pkin(4) * t39 + pkin(9) * t375 + t107 * t351 - t393 * t49 - t59 * t98 + t383, t107 * t24 - t258 * t40 + t329 * t96 + (-t368 - t398) * t223 + (-t107 * t347 + t292) * t227 + t386, -t19 * t308 + t23 * t96 - t24 * t98 - t265 * t374 + (-t37 - t336 + (t333 + t338) * qJD(5)) * pkin(9) + t262, -t258 * t39 - t329 * t98 + (-pkin(9) * t308 - t23) * t107 + (-t107 * t28 + t368) * t227 + t385, pkin(9) * t239 - g(1) * t182 - g(2) * t180 - g(3) * t276 - t18 * t24 - t19 * t23 - t258 * t5 + t28 * t329 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, -t96 ^ 2 + t371, t20, t98 * t107 - t40, t54, -t49 * t98 + t243 + t343, t107 * t25 + t49 * t96 - t236, -t69 * t96 - t237 + t343 + 0.2e1 * t369, pkin(5) * t39 - qJ(6) * t40 + (t19 - t26) * t98 + (t18 - t314) * t96, 0.2e1 * t346 - t28 * t96 + t69 * t98 + (0.2e1 * qJD(6) - t25) * t107 + t236, t2 * qJ(6) - t4 * pkin(5) - t28 * t69 - t18 * t26 - g(1) * (-pkin(5) * t146 + qJ(6) * t147) - g(2) * (-pkin(5) * t144 + qJ(6) * t145) + t267 * t196 + t314 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353 - t54, t20, -t370 - t371, t237 - t344 - t369;];
tau_reg  = t1;
