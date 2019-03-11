% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:11
% EndTime: 2019-03-09 12:20:24
% DurationCPUTime: 5.42s
% Computational Cost: add. (7029->549), mult. (14666->659), div. (0->0), fcn. (9917->8), ass. (0->264)
t185 = sin(qJ(4));
t186 = sin(qJ(2));
t189 = cos(qJ(4));
t190 = cos(qJ(2));
t296 = t189 * t190;
t114 = t185 * t186 + t296;
t107 = t114 * qJD(1);
t370 = qJD(5) + t107;
t381 = t370 ^ 2;
t184 = sin(qJ(5));
t188 = cos(qJ(5));
t285 = qJD(1) * t190;
t286 = qJD(1) * t186;
t110 = -t185 * t285 + t189 * t286;
t112 = -qJD(1) * pkin(1) - pkin(2) * t285 - qJ(3) * t286;
t87 = pkin(3) * t285 - t112;
t40 = pkin(4) * t107 - pkin(9) * t110 + t87;
t272 = qJD(2) - qJD(4);
t163 = pkin(7) * t285;
t121 = -pkin(8) * t285 + t163;
t179 = qJD(2) * qJ(3);
t111 = t121 + t179;
t350 = pkin(2) + pkin(3);
t267 = t350 * qJD(2);
t162 = pkin(7) * t286;
t119 = pkin(8) * t286 - t162;
t375 = qJD(3) - t119;
t91 = -t267 + t375;
t54 = t189 * t111 + t185 * t91;
t49 = -t272 * pkin(9) + t54;
t18 = t184 * t40 + t188 * t49;
t12 = qJ(6) * t370 + t18;
t380 = t370 * t12;
t251 = t188 * t272;
t77 = t110 * t184 + t251;
t379 = t370 * t77;
t79 = t188 * t110 - t184 * t272;
t378 = t370 * t79;
t274 = qJD(1) * qJD(2);
t264 = t186 * t274;
t280 = qJD(4) * t189;
t281 = qJD(4) * t185;
t282 = qJD(2) * t190;
t376 = t185 * t282 + t186 * t280 - t190 * t281;
t45 = qJD(1) * t376 + t114 * qJDD(1) - t189 * t264;
t377 = t110 * t272 + t45;
t42 = qJDD(5) + t45;
t241 = pkin(5) * t188 + qJ(6) * t184;
t227 = pkin(4) + t241;
t191 = cos(qJ(1));
t293 = t190 * t191;
t299 = t186 * t191;
t100 = t185 * t293 - t189 * t299;
t187 = sin(qJ(1));
t230 = t190 * t185 - t186 * t189;
t98 = t230 * t187;
t360 = -g(1) * t100 - g(2) * t98 - g(3) * t114;
t374 = t360 * t227;
t166 = t186 * qJDD(1);
t263 = t190 * t274;
t373 = -t263 - t166;
t257 = -t185 * qJ(3) - t189 * t350;
t92 = qJD(3) * t189 + qJD(4) * t257;
t312 = t92 * t370;
t289 = t189 * qJ(3) - t185 * t350;
t118 = -pkin(9) + t289;
t326 = t118 * t42;
t53 = -t185 * t111 + t189 * t91;
t48 = t272 * pkin(4) - t53;
t21 = t77 * pkin(5) - t79 * qJ(6) + t48;
t363 = t370 * t21;
t372 = t363 + t312 + t326;
t347 = pkin(9) * t42;
t371 = t363 - t347;
t215 = t230 * qJDD(1);
t216 = t114 * qJD(4);
t69 = t114 * qJD(2) - t216;
t196 = t69 * qJD(1) - t215;
t271 = qJDD(2) - qJDD(4);
t278 = qJD(5) * t184;
t24 = qJD(5) * t251 + t110 * t278 + t184 * t271 - t188 * t196;
t323 = t184 * t24;
t361 = t188 * t370;
t369 = t361 * t79 - t323;
t322 = t184 * t42;
t368 = -t110 * t79 + t361 * t370 + t322;
t283 = qJD(2) * t189;
t109 = t184 * t286 + t188 * t283;
t317 = t188 * t42;
t367 = t185 * (t272 * t79 - t278 * t370 + t317) + (t188 * t280 - t109) * t370 - t189 * t24;
t308 = qJD(5) * t79;
t25 = t184 * t196 + t188 * t271 + t308;
t366 = (t24 + t379) * t188 + (t25 + t378) * t184;
t288 = t190 * pkin(2) + t186 * qJ(3);
t365 = -pkin(1) - t288;
t66 = t119 * t189 + t121 * t185;
t364 = t92 - t66;
t332 = qJD(4) * t289 + t121 * t189 + t185 * t375;
t362 = t370 * t48;
t243 = g(1) * t191 + g(2) * t187;
t355 = -t110 * t77 + t184 * t381 - t317;
t310 = pkin(7) * qJDD(2);
t354 = (qJD(1) * t365 + t112) * qJD(2) - t310;
t353 = t272 ^ 2;
t240 = pkin(5) * t184 - qJ(6) * t188;
t352 = -qJD(6) * t184 + t370 * t240;
t351 = t79 ^ 2;
t349 = pkin(7) - pkin(8);
t348 = pkin(5) * t42;
t345 = pkin(5) * t110;
t343 = g(1) * t187;
t340 = g(2) * t191;
t338 = g(3) * t230;
t172 = t190 * pkin(3);
t337 = t79 * t77;
t61 = pkin(4) * t110 + pkin(9) * t107;
t336 = t184 * t61 + t188 * t53;
t335 = -t352 + t332;
t152 = qJ(3) * t285;
t97 = -t350 * t286 + t152;
t43 = -t61 + t97;
t334 = t184 * t43 + t188 * t66;
t269 = t172 + t288;
t229 = pkin(9) * t230 + t269;
t52 = pkin(4) * t114 + pkin(1) + t229;
t127 = t349 * t186;
t128 = t349 * t190;
t82 = t127 * t185 + t128 * t189;
t333 = t184 * t52 + t188 * t82;
t331 = pkin(9) * qJD(5);
t330 = qJ(6) * t42;
t329 = t370 * t18;
t17 = -t184 * t49 + t188 * t40;
t292 = qJD(6) - t17;
t11 = -pkin(5) * t370 + t292;
t328 = t11 * t110;
t327 = t110 * t12;
t325 = t17 * t110;
t324 = t18 * t110;
t321 = t184 * t69;
t320 = t184 * t77;
t319 = t184 * t79;
t318 = t188 * t25;
t316 = t188 * t69;
t315 = t188 * t77;
t314 = t188 * t79;
t311 = t352 - t54;
t309 = qJ(6) * t110;
t182 = qJDD(1) * pkin(1);
t307 = qJDD(2) * pkin(2);
t306 = t370 * t110;
t304 = t110 * t107;
t303 = t230 * t184;
t302 = t230 * t188;
t301 = t186 * t187;
t194 = qJD(1) ^ 2;
t298 = t186 * t194;
t297 = t187 * t190;
t167 = t186 * qJD(3);
t290 = qJ(3) * t282 + t167;
t180 = t186 ^ 2;
t181 = t190 ^ 2;
t287 = t180 - t181;
t284 = qJD(2) * t186;
t279 = qJD(5) * t118;
t277 = qJD(5) * t188;
t273 = t190 * qJDD(1);
t270 = t190 * t298;
t268 = -g(1) * t299 - g(2) * t301 + g(3) * t190;
t265 = t107 ^ 2 - t110 ^ 2;
t157 = pkin(7) * t166;
t262 = pkin(7) * t263 + qJDD(3) + t157;
t63 = t373 * pkin(8) - t350 * qJDD(2) + t262;
t158 = pkin(7) * t273;
t177 = qJDD(2) * qJ(3);
t178 = qJD(2) * qJD(3);
t89 = -pkin(7) * t264 + t158 + t177 + t178;
t64 = (t264 - t273) * pkin(8) + t89;
t220 = -t111 * t281 + t185 * t63 + t189 * t64 + t91 * t280;
t15 = -t271 * pkin(9) + t220;
t234 = pkin(2) * t273 - t373 * qJ(3) + qJD(1) * t167 + t182;
t226 = pkin(3) * t273 + t234;
t9 = t45 * pkin(4) + pkin(9) * t215 + (pkin(9) * t216 + (-pkin(9) * t296 + (-pkin(9) * t185 - t350) * t186) * qJD(2)) * qJD(1) + t226;
t261 = t184 * t15 - t188 * t9 + t49 * t277 + t40 * t278;
t99 = t114 * t187;
t71 = t184 * t99 - t191 * t188;
t258 = -qJD(2) * pkin(2) + qJD(3);
t250 = t191 * pkin(1) + pkin(2) * t293 + t187 * pkin(7) + qJ(3) * t299;
t249 = -t157 - t268;
t248 = t186 * t267;
t101 = t114 * t191;
t75 = t101 * t184 + t187 * t188;
t247 = -g(1) * t71 + g(2) * t75;
t72 = t184 * t191 + t188 * t99;
t76 = t101 * t188 - t187 * t184;
t246 = g(1) * t72 - g(2) * t76;
t245 = g(1) * t98 - g(2) * t100;
t193 = qJD(2) ^ 2;
t244 = pkin(7) * t193 + t340;
t238 = t11 * t188 - t12 * t184;
t237 = t11 * t184 + t12 * t188;
t233 = t111 * t280 + t185 * t64 - t189 * t63 + t91 * t281;
t123 = t162 + t258;
t126 = t163 + t179;
t232 = t123 * t190 - t126 * t186;
t231 = t127 * t189 - t128 * t185;
t96 = t262 - t307;
t225 = t188 * t15 + t184 * t9 + t40 * t277 - t49 * t278;
t224 = -0.2e1 * pkin(1) * t274 - t310;
t68 = -t186 * t283 + t376;
t85 = -t248 + t290;
t29 = pkin(4) * t68 - pkin(9) * t69 + t85;
t120 = t349 * t284;
t122 = qJD(2) * t128;
t36 = t231 * qJD(4) - t120 * t189 + t122 * t185;
t222 = t184 * t29 + t188 * t36 + t52 * t277 - t82 * t278;
t221 = -t347 + t362;
t219 = -t326 - t362;
t217 = g(1) * t101 + g(2) * t99 - t338;
t214 = -t244 + 0.2e1 * t182;
t16 = t271 * pkin(4) + t233;
t212 = -t16 - t360;
t211 = g(1) * t75 + g(2) * t71 - g(3) * t303 - t261;
t210 = t331 * t370 + t360;
t209 = t279 * t370 - t360;
t3 = t25 * pkin(5) + t24 * qJ(6) - t79 * qJD(6) + t16;
t208 = -t210 - t3;
t103 = pkin(2) * t284 - t290;
t59 = pkin(2) * t264 - t234;
t207 = -qJD(1) * t103 - qJDD(1) * t365 - t244 - t59;
t206 = -t209 + t3;
t1 = qJD(6) * t370 + t225 + t330;
t2 = qJDD(6) + t261 - t348;
t205 = t238 * qJD(5) + t1 * t188 + t2 * t184;
t204 = t21 * t79 + qJDD(6) - t211;
t203 = t232 * qJD(2) + t96 * t186 + t89 * t190;
t37 = t82 * qJD(4) - t120 * t185 - t122 * t189;
t202 = -g(1) * t76 - g(2) * t72 + g(3) * t302 + t225;
t201 = -t87 * t110 - t233 - t360;
t200 = t87 * t107 + t217 - t220;
t106 = t184 * t283 - t188 * t286;
t197 = -t189 * t25 + t77 * t281 + (-qJD(2) * t77 - t322) * t185 + (-t184 * t280 - t185 * t277 + t106) * t370;
t174 = t191 * pkin(7);
t146 = g(1) * t297;
t140 = qJ(3) * t293;
t138 = qJ(3) * t297;
t117 = pkin(4) - t257;
t116 = pkin(2) * t286 - t152;
t113 = pkin(1) + t269;
t88 = t227 - t257;
t46 = -qJD(1) * t248 + t226;
t38 = pkin(5) * t79 + qJ(6) * t77;
t33 = -t230 * t240 - t231;
t27 = -pkin(5) * t114 + t184 * t82 - t188 * t52;
t26 = qJ(6) * t114 + t333;
t23 = t184 * t53 - t188 * t61 - t345;
t22 = t309 + t336;
t20 = t184 * t66 - t188 * t43 + t345;
t19 = -t309 + t334;
t10 = -t24 + t379;
t6 = t240 * t69 - (t241 * qJD(5) - qJD(6) * t188) * t230 + t37;
t5 = -pkin(5) * t68 + t333 * qJD(5) + t184 * t36 - t188 * t29;
t4 = qJ(6) * t68 + qJD(6) * t114 + t222;
t7 = [qJDD(1), -t340 + t343, t243, qJDD(1) * t180 + 0.2e1 * t186 * t263, 0.2e1 * t186 * t273 - 0.2e1 * t274 * t287, qJDD(2) * t186 + t190 * t193, qJDD(2) * t190 - t186 * t193, 0, t186 * t224 + t190 * t214 + t146, t224 * t190 + (-t214 - t343) * t186, t186 * t354 + t207 * t190 + t146 (t180 + t181) * qJDD(1) * pkin(7) + t203 - t243, -t354 * t190 + (t207 + t343) * t186, pkin(7) * t203 - g(1) * t174 - g(2) * t250 + t112 * t103 + (-t343 + t59) * t365, t110 * t69 - t196 * t230, -t69 * t107 - t110 * t68 - t114 * t196 + t230 * t45, t230 * t271 - t272 * t69, t114 * t271 + t272 * t68, 0, g(1) * t99 - g(2) * t101 + t85 * t107 + t113 * t45 + t46 * t114 - t231 * t271 + t272 * t37 + t87 * t68, t85 * t110 + t113 * t196 - t230 * t46 + t271 * t82 + t272 * t36 + t87 * t69 - t245, t69 * t314 - (-t188 * t24 - t278 * t79) * t230 (-t315 - t319) * t69 - (t323 - t318 + (-t314 + t320) * qJD(5)) * t230, -t42 * t302 - t114 * t24 + t68 * t79 + (t230 * t278 + t316) * t370, t42 * t303 - t114 * t25 - t68 * t77 + (t230 * t277 - t321) * t370, t114 * t42 + t370 * t68, -t261 * t114 + t17 * t68 + t37 * t77 - t231 * t25 + ((-qJD(5) * t82 + t29) * t370 + t52 * t42 - t48 * qJD(5) * t230) * t188 + ((-qJD(5) * t52 - t36) * t370 - t82 * t42 - t16 * t230 + t48 * t69) * t184 + t246, -t222 * t370 - t333 * t42 - t225 * t114 - t18 * t68 + t37 * t79 + t231 * t24 + t48 * t316 - (t16 * t188 - t278 * t48) * t230 + t247, t21 * t321 - t370 * t5 - t11 * t68 - t114 * t2 + t25 * t33 - t27 * t42 + t6 * t77 - (t184 * t3 + t21 * t277) * t230 + t246, -t24 * t27 - t25 * t26 - t4 * t77 + t5 * t79 + t238 * t69 - (-qJD(5) * t237 - t1 * t184 + t188 * t2) * t230 + t245, -t21 * t316 + t1 * t114 + t370 * t4 + t12 * t68 + t24 * t33 + t26 * t42 - t6 * t79 - (-t188 * t3 + t21 * t278) * t230 - t247, t1 * t26 + t12 * t4 + t3 * t33 + t21 * t6 + t2 * t27 + t11 * t5 - g(1) * (-pkin(4) * t99 - pkin(5) * t72 - pkin(8) * t191 - pkin(9) * t98 - qJ(6) * t71 + t174) - g(2) * (pkin(3) * t293 + pkin(4) * t101 + pkin(5) * t76 + pkin(9) * t100 + qJ(6) * t75 + t250) + (-g(1) * (t365 - t172) + g(2) * pkin(8)) * t187; 0, 0, 0, -t270, t287 * t194, t166, t273, qJDD(2), pkin(1) * t298 + t249, g(3) * t186 - t158 + (pkin(1) * t194 + t243) * t190, 0.2e1 * t307 - qJDD(3) + (-t112 * t186 + t116 * t190) * qJD(1) + t249 (-t186 * pkin(2) + qJ(3) * t190) * qJDD(1) + ((t126 - t179) * t186 + (-t123 + t258) * t190) * qJD(1), t158 + 0.2e1 * t177 + 0.2e1 * t178 + (qJD(1) * t116 - g(3)) * t186 + (qJD(1) * t112 - t243) * t190, t89 * qJ(3) + t126 * qJD(3) - t96 * pkin(2) - t112 * t116 - g(1) * (-pkin(2) * t299 + t140) - g(2) * (-pkin(2) * t301 + t138) - g(3) * t288 - t232 * qJD(1) * pkin(7), -t304, t265, t215, t377, t271, -t97 * t107 - t257 * t271 + t272 * t332 - t201, -t97 * t110 + t271 * t289 + t272 * t364 - t200, -t369, t366, -t368, t355, t306, t325 + t117 * t25 + t332 * t77 + (-t364 * t370 + t219) * t184 + ((-t43 - t279) * t370 - t212) * t188, -t117 * t24 + t334 * t370 - t324 + t332 * t79 + (t219 - t312) * t188 + (-t16 + t209) * t184, -t372 * t184 + t206 * t188 + t20 * t370 + t25 * t88 + t335 * t77 - t328, t19 * t77 - t20 * t79 + (-t107 * t11 - t118 * t25 - t77 * t92 - t1 + (t118 * t79 - t11) * qJD(5)) * t188 + (t107 * t12 - t118 * t24 + t79 * t92 - t2 + (t118 * t77 + t12) * qJD(5)) * t184 + t217, t206 * t184 + t372 * t188 - t19 * t370 + t24 * t88 - t335 * t79 + t327, t3 * t88 - g(1) * (-pkin(9) * t101 + t140) - g(2) * (-pkin(9) * t99 + t138) - g(3) * t229 + t335 * t21 + (t188 * t92 - t19) * t12 + (t184 * t92 - t20) * t11 + t205 * t118 + t243 * t186 * t350 + t374; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t270, t166, -t180 * t194 - t193, -qJD(2) * t126 + t112 * t286 + t268 + t96, 0, 0, 0, 0, 0, -t107 * t286 - t185 * t353 - t189 * t271, -t110 * t286 + t185 * t271 - t189 * t353, 0, 0, 0, 0, 0, t197, -t367, t197, -t106 * t79 + t109 * t77 + (-t315 + t319) * t280 + (-t323 - t318 + (t314 + t320) * qJD(5)) * t185, t367, -t106 * t11 - t109 * t12 + (qJD(4) * t237 - t3) * t189 + (-t21 * t272 + t205) * t185 + t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, -t265, -t107 * t272 + t196, -t377, -t271, -t272 * t54 + t201, -t272 * t53 + t200, t369, -t366, t368, -t355, -t306, -pkin(4) * t25 - t325 - t54 * t77 + (t370 * t53 + t221) * t184 + ((-t61 - t331) * t370 + t212) * t188, pkin(4) * t24 + t336 * t370 + t324 - t54 * t79 + t221 * t188 + (t16 + t210) * t184, t371 * t184 + t208 * t188 - t227 * t25 + t23 * t370 + t311 * t77 + t328, t22 * t77 - t23 * t79 + (t1 + t370 * t11 + (-t25 + t308) * pkin(9)) * t188 + (t2 - t380 + (qJD(5) * t77 - t24) * pkin(9)) * t184 - t217, t208 * t184 - t371 * t188 - t22 * t370 - t227 * t24 - t311 * t79 - t327, -t11 * t23 - t12 * t22 - t3 * t227 + t311 * t21 + (t205 - t217) * pkin(9) - t374; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, -t77 ^ 2 + t351, t10, -t25 + t378, t42, -t48 * t79 + t211 + t329, t17 * t370 + t48 * t77 - t202, -t38 * t77 - t204 + t329 + 0.2e1 * t348, pkin(5) * t24 - qJ(6) * t25 + (t12 - t18) * t79 + (t11 - t292) * t77, 0.2e1 * t330 - t21 * t77 + t38 * t79 + (0.2e1 * qJD(6) - t17) * t370 + t202, t1 * qJ(6) - t2 * pkin(5) - t21 * t38 - t11 * t18 - g(1) * (-pkin(5) * t75 + qJ(6) * t76) - g(2) * (-pkin(5) * t71 + qJ(6) * t72) + t292 * t12 - t240 * t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337 - t42, t10, -t351 - t381, t204 - t348 - t380;];
tau_reg  = t7;
