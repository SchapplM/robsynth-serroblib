% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:46
% EndTime: 2019-03-09 08:15:55
% DurationCPUTime: 5.82s
% Computational Cost: add. (5629->582), mult. (11256->670), div. (0->0), fcn. (6697->10), ass. (0->297)
t211 = cos(qJ(2));
t209 = sin(qJ(2));
t362 = g(3) * t209;
t210 = sin(qJ(1));
t190 = g(2) * t210;
t212 = cos(qJ(1));
t377 = g(1) * t212 + t190;
t223 = t211 * t377 + t362;
t376 = qJDD(2) * pkin(4) + qJDD(5);
t174 = t211 * qJDD(1);
t314 = t211 * qJD(4);
t318 = qJD(2) * t209;
t244 = pkin(7) * t318 + t314;
t312 = qJD(1) * qJD(2);
t291 = t209 * t312;
t139 = qJ(4) * t291;
t163 = pkin(7) * t174;
t195 = qJDD(2) * qJ(3);
t196 = qJD(2) * qJD(3);
t303 = t163 + t195 + t196;
t273 = -t139 - t303;
t49 = qJ(4) * t174 + qJD(1) * t244 + t273;
t41 = -t49 + t376;
t383 = t41 - t223;
t173 = t209 * qJDD(1);
t290 = t211 * t312;
t236 = t290 + t173;
t114 = qJDD(6) + t236;
t202 = sin(pkin(9));
t203 = cos(pkin(9));
t208 = sin(qJ(6));
t368 = cos(qJ(6));
t116 = t202 * t368 + t208 * t203;
t322 = qJD(1) * t209;
t146 = qJD(6) + t322;
t293 = qJD(6) * t368;
t295 = t203 * t322;
t297 = t202 * t322;
t316 = qJD(6) * t208;
t348 = t202 * t316 - t203 * t293 + t208 * t297 - t295 * t368;
t269 = -t116 * t114 + t146 * t348;
t321 = qJD(1) * t211;
t296 = t202 * t321;
t315 = t203 * qJD(2);
t112 = t296 - t315;
t294 = t203 * t321;
t319 = qJD(2) * t202;
t113 = t294 + t319;
t58 = t368 * t112 + t113 * t208;
t382 = t58 ^ 2;
t59 = t208 * t112 - t113 * t368;
t381 = t59 ^ 2;
t380 = t146 * t58;
t370 = pkin(2) + pkin(3);
t379 = t209 * t370;
t241 = -t208 * t202 + t368 * t203;
t95 = t241 * t211;
t299 = t370 * qJD(2);
t292 = t370 * qJDD(2);
t166 = pkin(7) * t321;
t122 = -qJ(4) * t321 + t166;
t197 = qJD(2) * qJ(3);
t107 = -t122 - t197;
t183 = t209 * pkin(4);
t260 = t211 * qJ(5) + t183;
t191 = g(1) * t210;
t363 = g(2) * t212;
t378 = t191 - t363;
t198 = t209 ^ 2;
t199 = t211 ^ 2;
t277 = qJD(1) * (t198 - t199);
t310 = t203 * qJDD(2);
t78 = t310 + (t291 - t174) * t202;
t134 = t203 * t174;
t329 = -t202 * qJDD(2) - t134;
t79 = t203 * t291 + t329;
t282 = t208 * t79 + t368 * t78;
t17 = qJD(6) * t59 + t282;
t105 = t116 * qJD(6);
t347 = t116 * t322 + t105;
t375 = t17 * t241 + t347 * t58;
t374 = -t116 * t17 - t348 * t58;
t239 = t112 * t293 + t113 * t316 - t208 * t78 + t368 * t79;
t373 = -t116 * t239 + t348 * t59;
t237 = -t114 * t241 + t146 * t347;
t108 = -qJD(1) * pkin(1) - pkin(2) * t321 - qJ(3) * t322;
t177 = t209 * qJ(3);
t186 = t211 * pkin(2);
t327 = t186 + t177;
t124 = -pkin(1) - t327;
t346 = pkin(7) * qJDD(2);
t372 = (qJD(1) * t124 + t108) * qJD(2) - t346;
t81 = pkin(3) * t321 + qJD(4) - t108;
t60 = qJD(1) * t260 + t81;
t157 = qJ(4) * t322;
t307 = pkin(7) * t322;
t264 = qJD(3) + t307;
t248 = -t157 + t264;
t309 = qJ(5) + t370;
t77 = -qJD(2) * t309 + t248;
t26 = -t202 * t77 + t203 * t60;
t20 = pkin(5) * t322 + pkin(8) * t113 + t26;
t27 = t202 * t60 + t203 * t77;
t22 = pkin(8) * t112 + t27;
t246 = -t20 * t368 + t208 * t22;
t235 = pkin(4) * t211 - t209 * t309;
t220 = qJD(2) * t235 + t211 * qJD(5);
t176 = t209 * qJD(3);
t201 = qJDD(1) * pkin(1);
t256 = pkin(2) * t174 + qJ(3) * t236 + qJD(1) * t176 + t201;
t231 = pkin(3) * t174 + qJDD(4) + t256;
t21 = qJD(1) * t220 + qJDD(1) * t260 + t231;
t145 = pkin(7) * t290;
t162 = pkin(7) * t173;
t286 = qJDD(3) + t145 + t162;
t311 = qJD(1) * qJD(4);
t218 = -qJ(4) * t236 - t209 * t311 + t286;
t39 = -qJD(2) * qJD(5) - qJDD(2) * t309 + t218;
t10 = -t202 * t39 + t203 * t21;
t6 = pkin(5) * t236 - t79 * pkin(8) + t10;
t11 = t202 * t21 + t203 * t39;
t9 = -pkin(8) * t78 + t11;
t1 = -qJD(6) * t246 + t208 * t6 + t368 * t9;
t8 = t208 * t20 + t22 * t368;
t2 = -qJD(6) * t8 - t208 * t9 + t368 * t6;
t189 = g(3) * t211;
t336 = t209 * t212;
t337 = t209 * t210;
t301 = -g(1) * t336 - g(2) * t337 + t189;
t371 = t1 * t241 - t116 * t2 - t246 * t348 - t347 * t8 + t301;
t369 = pkin(5) * t78;
t367 = pkin(5) * t202;
t366 = pkin(8) * t211;
t187 = t212 * pkin(7);
t365 = g(1) * t187;
t364 = g(2) * qJ(4);
t185 = t211 * pkin(3);
t361 = t59 * t58;
t360 = pkin(7) - qJ(4);
t359 = pkin(8) + t309;
t154 = t203 * pkin(5) + pkin(4);
t250 = -pkin(8) * t203 * t209 + pkin(5) * t211;
t159 = qJ(3) * t321;
t65 = qJD(1) * t235 + t159;
t37 = -t202 * t122 + t203 * t65;
t28 = qJD(1) * t250 + t37;
t38 = t203 * t122 + t202 * t65;
t31 = -pkin(8) * t297 + t38;
t118 = t359 * t202;
t119 = t359 * t203;
t63 = t118 * t368 + t208 * t119;
t358 = -qJD(5) * t241 + qJD(6) * t63 - t208 * t28 - t31 * t368;
t64 = t208 * t118 - t119 * t368;
t357 = qJD(5) * t116 - qJD(6) * t64 + t208 * t31 - t28 * t368;
t317 = qJD(2) * t211;
t100 = -t209 * qJD(4) + t317 * t360;
t328 = qJ(3) * t317 + t176;
t50 = t220 + t328;
t30 = t203 * t100 + t202 * t50;
t356 = qJD(2) * pkin(2);
t352 = t202 * t78;
t351 = t78 * t203;
t350 = t79 * t202;
t349 = t79 * t203;
t131 = t360 * t209;
t302 = t185 + t327;
t254 = t302 + t260;
t75 = pkin(1) + t254;
t48 = t203 * t131 + t202 * t75;
t345 = qJ(3) * t211;
t344 = qJDD(2) * pkin(2);
t343 = t112 * t203;
t342 = t113 * t202;
t215 = qJD(1) ^ 2;
t175 = t198 * t215;
t340 = t199 * t215;
t338 = t209 * t154;
t335 = t209 * t215;
t334 = t210 * t211;
t207 = -pkin(8) - qJ(5);
t333 = t211 * t207;
t332 = t211 * t212;
t300 = -pkin(7) + t367;
t267 = t209 * t300;
t331 = -qJD(1) * t267 + qJD(3) - t157;
t330 = -qJD(4) - t81;
t326 = t212 * pkin(1) + t210 * pkin(7);
t323 = t198 + t199;
t320 = qJD(2) * t107;
t120 = -t157 + t307;
t313 = qJD(3) + t120;
t308 = -t207 + t370;
t306 = t58 * t321;
t305 = t59 * t321;
t304 = t163 + 0.2e1 * t195 + 0.2e1 * t196;
t289 = t202 * t173;
t288 = t202 * t174;
t287 = t203 * t173;
t285 = -pkin(1) - t177;
t284 = qJ(4) + t367;
t29 = -t202 * t100 + t203 * t50;
t47 = -t131 * t202 + t203 * t75;
t281 = t330 * t209;
t109 = pkin(1) + t302;
t280 = qJD(1) * t109 + t81;
t279 = g(1) * t309;
t278 = g(1) * t308;
t275 = t112 + t315;
t274 = -t113 + t319;
t272 = pkin(2) * t332 + qJ(3) * t336 + t326;
t271 = t162 + t301;
t270 = t209 * t299;
t268 = t209 * t290;
t266 = t323 * qJDD(1) * pkin(7);
t214 = qJD(2) ^ 2;
t265 = pkin(7) * t214 + t363;
t261 = g(1) * (-t212 * qJ(4) + t187);
t259 = t202 * t175 - t287;
t258 = -t10 * t202 + t11 * t203;
t257 = t202 * t26 - t203 * t27;
t255 = pkin(3) * t332 + t272;
t253 = t342 + t343;
t123 = t264 - t356;
t129 = t166 + t197;
t252 = t123 * t211 - t129 * t209;
t251 = -qJDD(3) - t271;
t249 = t285 - t186;
t32 = pkin(5) * t209 + t203 * t366 + t47;
t40 = t202 * t366 + t48;
t14 = -t208 * t40 + t32 * t368;
t15 = t208 * t32 + t368 * t40;
t245 = -0.2e1 * pkin(1) * t312 - t346;
t243 = g(2) * t255;
t242 = -t239 * t241 + t347 * t59;
t240 = t145 - t251;
t238 = (-t202 * t27 - t203 * t26) * t209;
t234 = -t211 * t41 + t377;
t233 = -t265 + 0.2e1 * t201;
t232 = t10 * t203 + t11 * t202 - t363;
t230 = -t175 * t203 - t289;
t229 = -qJ(4) * qJDD(1) - t377;
t88 = qJD(2) * pkin(4) + qJD(5) - t107;
t227 = t309 * t317 + (qJD(5) - t88) * t209;
t36 = -qJD(1) * t270 + t231;
t80 = -t270 + t328;
t225 = -qJD(1) * t80 - qJDD(1) * t109 - t36 + t363;
t222 = -qJ(4) * t173 + t240;
t101 = pkin(2) * t318 - t328;
t53 = pkin(2) * t291 - t256;
t221 = -qJD(1) * t101 - qJDD(1) * t124 - t265 - t53;
t86 = -pkin(7) * t291 + t303;
t96 = t286 - t344;
t219 = qJD(2) * t252 + t96 * t209 + t86 * t211;
t217 = (-pkin(7) * t312 - g(3)) * t209 + (t229 - t311) * t211 - t273 + t376;
t206 = qJ(3) + pkin(4);
t194 = pkin(9) + qJ(6);
t180 = t211 * qJ(4);
t172 = cos(t194);
t171 = sin(t194);
t160 = qJ(4) * t318;
t151 = g(1) * t334;
t150 = g(1) * t337;
t144 = qJ(3) * t332;
t142 = qJ(3) * t334;
t141 = t211 * t335;
t138 = qJ(3) + t154;
t137 = -t175 - t214;
t132 = pkin(7) * t211 - t180;
t130 = qJDD(2) + t141;
t128 = -t175 + t340;
t126 = qJDD(2) * t211 - t209 * t214;
t125 = qJDD(2) * t209 + t211 * t214;
t121 = pkin(2) * t322 - t159;
t111 = qJDD(1) * t199 - 0.2e1 * t268;
t110 = qJDD(1) * t198 + 0.2e1 * t268;
t103 = -t211 * t300 - t180;
t99 = -t160 + t244;
t98 = -t322 * t370 + t159;
t94 = t116 * t211;
t93 = -t299 + t248;
t92 = -t171 * t210 + t172 * t336;
t91 = -t171 * t336 - t172 * t210;
t90 = -t171 * t212 - t172 * t337;
t89 = t171 * t337 - t172 * t212;
t87 = -qJD(2) * t277 + t174 * t209;
t76 = 0.2e1 * t87;
t73 = qJD(2) * t267 + t160 - t314;
t54 = -pkin(5) * t112 + t88;
t46 = -t292 + t218;
t45 = qJD(6) * t95 - t116 * t318;
t44 = t105 * t211 + t241 * t318;
t25 = t41 + t369;
t24 = -pkin(8) * t202 * t318 + t30;
t23 = qJD(2) * t250 + t29;
t4 = -qJD(6) * t15 - t208 * t24 + t23 * t368;
t3 = qJD(6) * t14 + t208 * t23 + t24 * t368;
t5 = [0, 0, 0, 0, 0, qJDD(1), t378, t377, 0, 0, t110, t76, t125, t111, t126, 0, t209 * t245 + t211 * t233 + t151, -t209 * t233 + t211 * t245 - t150, 0.2e1 * t266 - t377, -g(1) * (-t210 * pkin(1) + t187) - g(2) * t326 + (pkin(7) ^ 2 * t323 + pkin(1) ^ 2) * qJDD(1), t110, t125, -0.2e1 * t87, 0, -t126, t111, t209 * t372 + t221 * t211 + t151, t266 + t219 - t377, t221 * t209 - t211 * t372 + t150, pkin(7) * t219 - g(2) * t272 + t108 * t101 + t53 * t124 - t191 * t249 - t365, t111, t76, t126, t110, t125, 0, qJDD(2) * t132 + t150 + (t211 * t280 - t99) * qJD(2) - t225 * t209, qJDD(2) * t131 - t151 + (t209 * t280 + t100) * qJD(2) + t225 * t211 (-qJD(2) * t93 - qJDD(1) * t132 + t49 + (-qJD(2) * t131 + t99) * qJD(1)) * t211 + (-t320 - qJDD(1) * t131 - t46 + (qJD(2) * t132 - t100) * qJD(1)) * t209 + t377, t46 * t131 + t93 * t100 - t49 * t132 + t107 * t99 + t36 * t109 + t81 * t80 - t261 - t243 + (-g(1) * (t249 - t185) + t364) * t210 (-t113 * t318 - t211 * t79) * t203 (t350 + t351) * t211 + t253 * t318 (t79 - t134) * t209 + (-t113 * t211 + t203 * t277) * qJD(2) (-t112 * t318 - t211 * t78) * t202 (-t78 + t288) * t209 + (t112 * t211 - t202 * t277) * qJD(2), t110, t99 * t112 + t132 * t78 + (qJD(1) * t47 + t26) * t317 + t234 * t202 + (t29 * qJD(1) + t47 * qJDD(1) + t203 * t378 + t319 * t88 + t10) * t209, t99 * t113 + t132 * t79 + (-qJD(1) * t48 - t27) * t317 + t234 * t203 + (-t30 * qJD(1) - t48 * qJDD(1) - t202 * t378 + t315 * t88 - t11) * t209, qJD(2) * t238 + t112 * t30 + t113 * t29 + t211 * t232 - t47 * t79 - t48 * t78 + t151, t11 * t48 + t27 * t30 + t10 * t47 + t26 * t29 + t41 * t132 - t88 * t99 - t261 - g(2) * (pkin(4) * t336 + qJ(5) * t332 + t255) + (-g(1) * (t285 - t183) + t364 + t211 * t279) * t210, -t239 * t95 + t44 * t59, t17 * t95 + t239 * t94 + t44 * t58 + t45 * t59, -t114 * t95 + t146 * t44 + t209 * t239 + t317 * t59, -t17 * t94 + t45 * t58, t114 * t94 + t146 * t45 - t17 * t209 + t317 * t58, t114 * t209 + t146 * t317, -g(1) * t90 - g(2) * t92 + t103 * t17 + t114 * t14 + t146 * t4 + t2 * t209 - t246 * t317 - t25 * t94 - t45 * t54 - t58 * t73, -g(1) * t89 - g(2) * t91 - t1 * t209 + t103 * t239 - t114 * t15 - t146 * t3 - t25 * t95 - t317 * t8 + t44 * t54 + t59 * t73, -g(2) * t332 + t1 * t94 - t14 * t239 - t15 * t17 + t2 * t95 + t246 * t44 + t3 * t58 - t4 * t59 + t45 * t8 + t151, t1 * t15 + t8 * t3 + t2 * t14 - t246 * t4 + t25 * t103 + t54 * t73 - t365 - t243 + (g(1) * t284 - g(2) * (-t333 + t338)) * t212 + (-g(1) * (t285 - t338) + g(2) * t284 + t211 * t278) * t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t128, t173, t141, t174, qJDD(2), pkin(1) * t335 - t271, t362 - t163 + (pkin(1) * t215 + t377) * t211, 0, 0, -t141, t173, t128, qJDD(2), -t174, t141, 0.2e1 * t344 + (-t108 * t209 + t121 * t211) * qJD(1) + t251 (-pkin(2) * t209 + t345) * qJDD(1) + ((t129 - t197) * t209 + (qJD(3) - t123 - t356) * t211) * qJD(1) (qJD(1) * t121 - g(3)) * t209 + (qJD(1) * t108 - t377) * t211 + t304, t86 * qJ(3) + t129 * qJD(3) - t96 * pkin(2) - t108 * t121 - g(1) * (-pkin(2) * t336 + t144) - g(2) * (-pkin(2) * t337 + t142) - g(3) * t327 - t252 * qJD(1) * pkin(7), t141, -t128, t174, -t141, t173, qJDD(2), qJD(2) * t120 + t139 + (-g(3) + (-pkin(7) * qJD(2) - t98) * qJD(1)) * t209 + (qJD(1) * t330 + t229) * t211 + t304, -qJD(2) * t122 - 0.2e1 * t292 + ((-qJ(4) * qJD(2) + t98) * t211 + t281) * qJD(1) + t222 (-t345 + t379) * qJDD(1) + (-t313 + t93 + t299) * t321, -g(1) * t144 - g(2) * t142 - g(3) * t302 - t49 * qJ(3) - t313 * t107 - t93 * t122 - t370 * t46 + t377 * t379 - t81 * t98, t113 * t295 - t350, -t253 * t322 - t349 + t352, -t274 * t321 + t230, t112 * t297 + t351, -t275 * t321 + t259, -t141, t309 * t289 + t206 * t78 - t313 * t112 + t383 * t203 + (t202 * t227 - t209 * t37 - t211 * t26) * qJD(1), t309 * t287 + t206 * t79 - t313 * t113 - t383 * t202 + (t203 * t227 + t209 * t38 + t211 * t27) * qJD(1), -t112 * t38 - t113 * t37 + (-qJD(5) * t112 + t26 * t322 + t309 * t78 - t11) * t203 + (qJD(5) * t113 + t27 * t322 - t309 * t79 + t10) * t202 - t301, t41 * t206 - t27 * t38 - t26 * t37 - g(1) * (pkin(4) * t332 + t144) - g(2) * (pkin(4) * t334 + t142) - g(3) * t254 + t313 * t88 - t258 * t309 + t257 * qJD(5) + (t190 * t309 + t212 * t279) * t209, t373, t242 - t374, t269 - t305, t375, t237 - t306, -t146 * t321, t114 * t63 + t138 * t17 + t146 * t357 - t172 * t223 + t241 * t25 + t246 * t321 - t331 * t58 - t347 * t54, -t114 * t64 - t116 * t25 + t138 * t239 - t146 * t358 + t171 * t223 + t321 * t8 + t331 * t59 + t348 * t54, -t17 * t64 - t239 * t63 - t357 * t59 + t358 * t58 - t371, t1 * t64 + t2 * t63 + t25 * t138 - g(1) * (t154 * t332 + t144) - g(2) * (t154 * t334 + t142) - g(3) * (t302 - t333) + t358 * t8 - t357 * t246 + t331 * t54 + (-g(3) * t154 + t190 * t308 + t212 * t278) * t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t173, t137, -qJD(2) * t129 + t108 * t322 + t240 - t344, 0, 0, 0, 0, 0, 0, t137, t130, -t173, t320 - t292 + (-qJ(4) * t317 + t281) * qJD(1) + t222, 0, 0, 0, 0, 0, 0 (t112 - t296) * qJD(2) + t230 (t113 - t294) * qJD(2) + t259, t350 - t351 + (-t112 * t202 - t113 * t203) * t322, qJD(1) * t238 - t88 * qJD(2) + t258 + t301, 0, 0, 0, 0, 0, 0, qJD(2) * t58 + t269, -qJD(2) * t59 + t237, -t373 - t375, -qJD(2) * t54 + t371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173 + 0.2e1 * t290, -t174 + 0.2e1 * t291, -t175 - t340 (-t107 * t211 + (t93 - t299) * t209) * qJD(1) + t231 + t378, 0, 0, 0, 0, 0, 0 (-t112 + t315) * t321 - t259 (-t113 - t319) * t321 + t230, -t352 - t349 + (-t342 + t343) * t322, t191 + (-t209 * t257 + t211 * t88) * qJD(1) + t232, 0, 0, 0, 0, 0, 0, -t237 - t306, t269 + t305, t242 + t374, t1 * t116 + t2 * t241 + t246 * t347 + t321 * t54 - t348 * t8 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274 * t322 - t288 + t310, t275 * t322 + t329, -t112 ^ 2 - t113 ^ 2, -t112 * t27 - t113 * t26 + t217, 0, 0, 0, 0, 0, 0, t146 * t59 + t17, t239 + t380, -t381 - t382, -t246 * t59 - t8 * t58 + t217 + t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361, t381 - t382, t239 - t380, t361, -t282 + (-qJD(6) + t146) * t59, t114, -g(1) * t91 + g(2) * t89 + t8 * t146 - t171 * t189 - t54 * t59 + t2, g(1) * t92 - g(2) * t90 - t146 * t246 - t172 * t189 - t54 * t58 - t1, 0, 0;];
tau_reg  = t5;
