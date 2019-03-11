% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:46:00
% EndTime: 2019-03-09 16:46:14
% DurationCPUTime: 6.26s
% Computational Cost: add. (7825->560), mult. (16952->661), div. (0->0), fcn. (11921->10), ass. (0->284)
t206 = sin(qJ(3));
t207 = sin(qJ(2));
t210 = cos(qJ(2));
t367 = cos(qJ(3));
t141 = t206 * t210 + t207 * t367;
t131 = t141 * qJD(1);
t384 = qJD(5) + t131;
t279 = t384 ^ 2;
t291 = t367 * t210;
t270 = qJD(1) * t291;
t310 = qJD(1) * t207;
t290 = t206 * t310;
t129 = -t270 + t290;
t200 = qJD(2) + qJD(3);
t205 = sin(qJ(5));
t209 = cos(qJ(5));
t105 = -t209 * t129 + t200 * t205;
t379 = t105 * t384;
t199 = qJDD(2) + qJDD(3);
t100 = t200 * t141;
t283 = qJDD(1) * t367;
t173 = t210 * t283;
t302 = t207 * qJDD(1);
t258 = t206 * t302 - t173;
t216 = qJD(1) * t100 + t258;
t306 = qJD(5) * t209;
t307 = qJD(5) * t205;
t47 = -t129 * t306 - t209 * t199 + t200 * t307 - t205 * t216;
t387 = t47 - t379;
t321 = t206 * t207;
t261 = t200 * t321;
t301 = t210 * qJDD(1);
t272 = t200 * t270 + t206 * t301 + t207 * t283;
t78 = qJD(1) * t261 - t272;
t72 = -qJDD(5) + t78;
t67 = t209 * t72;
t386 = -t205 * t279 - t67;
t286 = qJD(3) * t367;
t212 = -pkin(8) - pkin(7);
t154 = t212 * t210;
t144 = qJD(1) * t154;
t133 = t206 * t144;
t153 = t212 * t207;
t142 = qJD(1) * t153;
t98 = t142 * t367 + t133;
t385 = pkin(2) * t286 - t98;
t107 = t129 * t205 + t200 * t209;
t195 = t200 * qJ(4);
t359 = t129 * pkin(4);
t136 = t367 * t144;
t351 = qJD(2) * pkin(2);
t138 = t142 + t351;
t94 = t206 * t138 - t136;
t77 = t94 - t359;
t66 = t195 + t77;
t38 = pkin(5) * t105 - qJ(6) * t107 + t66;
t303 = qJD(1) * qJD(2);
t284 = t210 * t303;
t103 = qJDD(2) * pkin(2) - t212 * (-t284 - t302);
t285 = t207 * t303;
t104 = t212 * (-t285 + t301);
t309 = qJD(3) * t206;
t275 = -t206 * t103 + t367 * t104 - t138 * t286 - t144 * t309;
t190 = t199 * qJ(4);
t376 = -t200 * qJD(4) - t190;
t39 = t275 + t376;
t19 = -pkin(4) * t216 - t39;
t281 = t199 * t205 - t209 * t216;
t48 = qJD(5) * t107 + t281;
t5 = t48 * pkin(5) + t47 * qJ(6) - t107 * qJD(6) + t19;
t383 = -t5 * t205 - t38 * t306;
t382 = t19 * t205 + t66 * t306;
t340 = qJD(4) + t385;
t381 = -t209 * t5 + t38 * t307;
t380 = -t105 * t200 + t386;
t298 = pkin(2) * t309;
t97 = t206 * t142 - t136;
t269 = -t97 + t298;
t378 = t131 * t200;
t277 = t384 * t107;
t93 = -t367 * t138 - t133;
t316 = qJD(4) + t93;
t260 = -pkin(5) * t209 - qJ(6) * t205;
t249 = -pkin(4) + t260;
t377 = pkin(5) * t306 + qJ(6) * t307 - t209 * qJD(6) - t249 * t131 + qJD(4);
t204 = qJ(2) + qJ(3);
t196 = sin(t204);
t197 = cos(t204);
t312 = t197 * pkin(3) + t196 * qJ(4);
t217 = -t378 - t258;
t198 = t205 * pkin(5);
t339 = qJ(6) * t209;
t375 = -t339 + t198;
t208 = sin(qJ(1));
t211 = cos(qJ(1));
t264 = g(1) * t211 + g(2) * t208;
t356 = t210 * pkin(2);
t189 = pkin(1) + t356;
t374 = qJDD(1) * t189;
t110 = t206 * t153 - t154 * t367;
t263 = g(1) * t208 - g(2) * t211;
t292 = qJD(2) * t212;
t143 = t207 * t292;
t145 = t210 * t292;
t62 = -t367 * t143 - t206 * t145 - t153 * t286 - t154 * t309;
t373 = -t110 * t199 - t196 * t263 + t200 * t62;
t320 = t208 * t209;
t322 = t205 * t211;
t121 = t196 * t322 + t320;
t319 = t209 * t211;
t323 = t205 * t208;
t123 = -t196 * t323 + t319;
t185 = g(3) * t197;
t178 = pkin(2) * t285;
t246 = t78 * qJ(4) - t131 * qJD(4) + t178;
t338 = qJDD(1) * pkin(1);
t12 = -pkin(2) * t301 - pkin(3) * t217 + pkin(9) * t216 + t246 - t338;
t274 = t367 * t103 + t206 * t104 - t138 * t309 + t144 * t286;
t251 = qJDD(4) - t274;
t369 = pkin(3) + pkin(9);
t18 = -pkin(4) * t78 - t199 * t369 + t251;
t358 = t131 * pkin(4);
t317 = t358 + t316;
t58 = -t200 * t369 + t317;
t300 = t209 * t12 + t205 * t18 + t58 * t306;
t152 = t189 * qJD(1);
t233 = -qJ(4) * t131 - t152;
t61 = t129 * t369 + t233;
t372 = -g(1) * t121 + g(2) * t123 + (-qJD(5) * t61 + t185) * t205 + t300;
t371 = t107 ^ 2;
t370 = t131 ^ 2;
t368 = pkin(5) * t72;
t366 = pkin(2) * t206;
t365 = pkin(2) * t207;
t364 = pkin(3) * t199;
t363 = pkin(5) * t129;
t184 = g(3) * t196;
t360 = g(3) * t210;
t182 = t197 * pkin(9);
t282 = t205 * t12 - t209 * t18 + t61 * t306 + t58 * t307;
t3 = qJDD(6) + t282 + t368;
t2 = t3 * t209;
t124 = t131 * pkin(9);
t91 = pkin(3) * t131 + qJ(4) * t129;
t87 = pkin(2) * t310 + t91;
t64 = t124 + t87;
t81 = t97 - t359;
t355 = t205 * t81 + t209 * t64;
t73 = t124 + t91;
t354 = t205 * t77 + t209 * t73;
t140 = -t291 + t321;
t250 = -qJ(4) * t141 - t189;
t79 = t140 * t369 + t250;
t109 = -t153 * t367 - t206 * t154;
t88 = t141 * pkin(4) + t109;
t353 = t205 * t88 + t209 * t79;
t352 = qJ(6) * t72;
t29 = t205 * t58 + t209 * t61;
t350 = t384 * t29;
t23 = qJ(6) * t384 + t29;
t349 = t131 * t23;
t188 = -pkin(2) * t367 - pkin(3);
t179 = -pkin(9) + t188;
t348 = t179 * t72;
t347 = t200 * t94;
t346 = t205 * t72;
t345 = t209 * t47;
t344 = t369 * t72;
t343 = t377 + t385;
t342 = t93 + t377;
t341 = t358 + t340;
t337 = t100 * t205;
t336 = t100 * t209;
t335 = t107 * t105;
t334 = t107 * t200;
t333 = t384 * t129;
t332 = t131 * t129;
t331 = t131 * t209;
t330 = t140 * t205;
t329 = t196 * t208;
t328 = t196 * t209;
t327 = t196 * t211;
t326 = t197 * t208;
t325 = t197 * t209;
t324 = t197 * t211;
t318 = t211 * t212;
t28 = -t205 * t61 + t209 * t58;
t315 = qJD(6) - t28;
t155 = qJ(4) * t326;
t299 = t197 * t198;
t314 = t208 * t299 + t155;
t157 = qJ(4) * t324;
t313 = t211 * t299 + t157;
t202 = t207 ^ 2;
t311 = -t210 ^ 2 + t202;
t308 = qJD(5) * t179;
t305 = qJD(5) * t369;
t193 = t207 * t351;
t297 = t197 * t320;
t296 = t197 * t319;
t295 = -g(1) * t296 - g(2) * t297 - g(3) * t328;
t294 = pkin(3) * t324 + qJ(4) * t327 + t211 * t189;
t293 = -g(1) * t327 - g(2) * t329 + t185;
t289 = t179 * t306;
t288 = t209 * t305;
t244 = -t307 * t61 + t300;
t1 = qJD(6) * t384 + t244 - t352;
t22 = -pkin(5) * t384 + t315;
t287 = -t131 * t22 - t1;
t280 = t384 * t66;
t278 = t209 * t384;
t276 = -t2 + t293;
t271 = t196 * t198 + t182 + t312;
t268 = -g(1) * t326 + g(2) * t324;
t267 = -pkin(3) * t196 - t365;
t120 = -t196 * t319 + t323;
t122 = t196 * t320 + t322;
t266 = g(1) * t122 + g(2) * t120;
t265 = -g(1) * t123 - g(2) * t121;
t150 = qJ(4) + t375;
t86 = -pkin(3) * t200 + t316;
t90 = -t195 - t94;
t257 = t86 * t129 - t90 * t131;
t256 = t205 * t23 - t209 * t22;
t255 = t205 * t22 + t209 * t23;
t254 = -t129 * t22 + t38 * t331 - t383;
t252 = t28 * t129 + t66 * t331 + t382;
t248 = -t189 - t312;
t247 = -0.2e1 * pkin(1) * t303 - pkin(7) * qJDD(2);
t99 = -qJD(2) * t291 - t210 * t286 + t261;
t245 = qJ(4) * t99 - qJD(4) * t141 + t193;
t37 = t100 * t369 + t245;
t63 = qJD(3) * t110 + t206 * t143 - t145 * t367;
t46 = -t99 * pkin(4) + t63;
t243 = t205 * t46 + t209 * t37 + t88 * t306 - t307 * t79;
t242 = t140 * t306 + t337;
t241 = -t29 * t129 + t19 * t209 + t295;
t238 = g(1) * t324 + g(2) * t326 + t184 + t275;
t237 = t274 - t293;
t234 = -t209 * t279 + t346;
t231 = -t197 * t264 - t184;
t230 = g(1) * t120 - g(2) * t122 + g(3) * t325 - t282;
t229 = t109 * t199 + t63 * t200 + t268;
t214 = qJD(2) ^ 2;
t227 = -pkin(7) * t214 + t263 + 0.2e1 * t338;
t215 = qJD(1) ^ 2;
t226 = pkin(1) * t215 - pkin(7) * qJDD(1) + t264;
t225 = t38 * t205 * t131 + t129 * t23 - t295 + t381;
t224 = t152 * t131 + t237;
t223 = -t152 * t129 + t238;
t84 = pkin(3) * t129 + t233;
t222 = t131 * t84 + qJDD(4) - t237;
t221 = -t129 * t84 - t238 - t376;
t220 = qJD(5) * t255 + t1 * t205 - t2;
t219 = t107 * t38 + qJDD(6) - t230;
t53 = t272 + (t129 - t290) * t200;
t181 = qJ(4) + t366;
t137 = t150 + t366;
t126 = t178 - t374;
t119 = t129 * qJ(6);
t92 = pkin(3) * t140 + t250;
t89 = -t140 * pkin(4) + t110;
t80 = -t129 ^ 2 + t370;
t65 = pkin(5) * t107 + qJ(6) * t105;
t54 = t217 + t378;
t52 = t140 * t249 + t110;
t50 = pkin(3) * t100 + t245;
t45 = -pkin(4) * t100 - t62;
t42 = t251 - t364;
t41 = -pkin(5) * t141 + t205 * t79 - t209 * t88;
t40 = qJ(6) * t141 + t353;
t33 = t205 * t73 - t209 * t77 + t363;
t32 = -t119 + t354;
t31 = t205 * t64 - t209 * t81 + t363;
t30 = -t119 + t355;
t24 = -t374 + t246 + (qJDD(1) * t321 - t173 + t378) * pkin(3);
t21 = -t105 * t129 + t234;
t20 = t107 * t129 + t386;
t13 = -t205 * t277 - t345;
t9 = (qJD(5) * t375 - qJD(6) * t205) * t140 + t249 * t100 - t62;
t8 = (-t48 - t277) * t209 + (t47 + t379) * t205;
t7 = pkin(5) * t99 + qJD(5) * t353 + t205 * t37 - t209 * t46;
t6 = -qJ(6) * t99 + qJD(6) * t141 + t243;
t4 = [qJDD(1), t263, t264, qJDD(1) * t202 + 0.2e1 * t207 * t284, 0.2e1 * t207 * t301 - 0.2e1 * t303 * t311, qJDD(2) * t207 + t210 * t214, qJDD(2) * t210 - t207 * t214, 0, t207 * t247 + t210 * t227, -t207 * t227 + t210 * t247, -t131 * t99 - t141 * t78, -t131 * t100 + t99 * t129 + t78 * t140 + t141 * t217, t141 * t199 - t200 * t99, -t100 * t200 - t140 * t199, 0, -t152 * t100 + t126 * t140 + t129 * t193 - t189 * t216 - t229, t126 * t141 + t131 * t193 + t152 * t99 + t189 * t78 + t373, t90 * t100 - t109 * t78 - t110 * t216 + t62 * t129 + t63 * t131 + t39 * t140 + t42 * t141 - t86 * t99 - t264, -t84 * t100 - t50 * t129 - t24 * t140 + t217 * t92 + t229, -t131 * t50 - t141 * t24 + t78 * t92 + t84 * t99 - t373, t24 * t92 + t84 * t50 - t39 * t110 + t90 * t62 + t42 * t109 + t86 * t63 + g(1) * t318 - g(2) * t294 + (-g(1) * t248 + g(2) * t212) * t208, t107 * t242 - t330 * t47 (-t105 * t205 + t107 * t209) * t100 + (-t205 * t48 - t345 + (-t105 * t209 - t107 * t205) * qJD(5)) * t140, -t107 * t99 - t141 * t47 + t242 * t384 - t330 * t72, -t140 * t67 + t105 * t99 - t141 * t48 + (-t140 * t307 + t336) * t384, -t141 * t72 - t384 * t99, -t282 * t141 - t28 * t99 + t45 * t105 + t89 * t48 + ((-qJD(5) * t88 - t37) * t384 + t79 * t72 + t66 * qJD(5) * t140) * t205 + ((-qJD(5) * t79 + t46) * t384 - t88 * t72 - t19 * t140 - t66 * t100) * t209 + t265, t45 * t107 + t140 * t382 - t244 * t141 - t243 * t384 + t29 * t99 + t66 * t337 + t353 * t72 - t89 * t47 + t266, t105 * t9 + t140 * t381 - t141 * t3 + t22 * t99 - t38 * t336 - t384 * t7 + t41 * t72 + t48 * t52 + t265, -t105 * t6 + t107 * t7 - t40 * t48 - t41 * t47 + t255 * t100 + (-qJD(5) * t256 + t1 * t209 + t205 * t3) * t140 - t268, t1 * t141 - t107 * t9 + t140 * t383 - t23 * t99 - t38 * t337 + t384 * t6 - t40 * t72 + t47 * t52 - t266, t1 * t40 + t23 * t6 + t5 * t52 + t38 * t9 + t3 * t41 + t22 * t7 - g(1) * (t211 * pkin(4) + t123 * pkin(5) + t122 * qJ(6) - t318) - g(2) * (t121 * pkin(5) + pkin(9) * t324 + t120 * qJ(6) + t294) + (-g(1) * (t248 - t182) - g(2) * (pkin(4) - t212)) * t208; 0, 0, 0, -t207 * t215 * t210, t311 * t215, t302, t301, qJDD(2), t207 * t226 - t360, g(3) * t207 + t210 * t226, t332, t80, t53, t54, t199, t97 * t200 + (-t129 * t310 + t199 * t367 - t200 * t309) * pkin(2) + t224, t98 * t200 + (-t131 * t310 - t199 * t206 - t200 * t286) * pkin(2) + t223, -t129 * t340 + t131 * t269 - t181 * t216 - t188 * t78 + t257, t129 * t87 + t269 * t200 + (-pkin(3) + t188) * t199 + t222, t131 * t87 + t181 * t199 + t200 * t340 + t221, -t39 * t181 + t42 * t188 - t84 * t87 - g(1) * (t211 * t267 + t157) - g(2) * (t208 * t267 + t155) - g(3) * (t312 + t356) - t340 * t90 + t269 * t86, t13, t8, t20, t21, t333, t181 * t48 + t341 * t105 + (-t348 + (-t81 + t298) * t384) * t209 + ((t64 - t308) * t384 + t231) * t205 + t252, -t181 * t47 + (-t289 + t355) * t384 + t341 * t107 + (-t298 * t384 - t280 + t348) * t205 + t241, -t179 * t67 + t137 * t48 + (t209 * t298 + t31) * t384 + t343 * t105 + (-t308 * t384 + t231) * t205 + t254, t105 * t30 - t107 * t31 + (-t107 * t298 - t349 + t179 * t47 + (-t105 * t179 - t23) * qJD(5)) * t209 + (-t105 * t298 - t179 * t48 + (t107 * t179 - t22) * qJD(5) + t287) * t205 - t276, -t179 * t346 + t137 * t47 - t343 * t107 + (t205 * t298 + t289 - t30) * t384 + t225, t5 * t137 - t23 * t30 - t22 * t31 - g(1) * t313 - g(2) * t314 - g(3) * (-qJ(6) * t328 + t271) + t343 * t38 + (t256 * t309 - t360) * pkin(2) + t220 * t179 + t264 * (qJ(6) * t325 + t196 * t369 + t365); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t80, t53, t54, t199, t224 + t347, -t200 * t93 + t223, pkin(3) * t78 - qJ(4) * t216 - t129 * t316 - t94 * t131 + t257, t129 * t91 + t222 - t347 - 0.2e1 * t364, t131 * t91 + t200 * t316 + t190 + t221, -t39 * qJ(4) - t42 * pkin(3) - t84 * t91 - t86 * t94 - g(1) * (-pkin(3) * t327 + t157) - g(2) * (-pkin(3) * t329 + t155) - g(3) * t312 - t316 * t90, t13, t8, t20, t21, t333, qJ(4) * t48 + (-t384 * t77 + t344) * t209 + t317 * t105 + ((t73 + t305) * t384 + t231) * t205 + t252, -qJ(4) * t47 + (t288 + t354) * t384 + t317 * t107 + (-t280 - t344) * t205 + t241, t209 * t344 + t384 * t33 + t150 * t48 + t342 * t105 + (t305 * t384 + t231) * t205 + t254, t105 * t32 - t107 * t33 + (-t349 - t369 * t47 + (t105 * t369 - t23) * qJD(5)) * t209 + (t369 * t48 + (-t107 * t369 - t22) * qJD(5) + t287) * t205 - t276, t205 * t344 + t150 * t47 + (-t32 - t288) * t384 - t342 * t107 + t225, t5 * t150 - t23 * t32 - t22 * t33 - g(1) * (-qJ(6) * t296 + t313) - g(2) * (-qJ(6) * t297 + t314) - g(3) * t271 + t342 * t38 + (g(3) * t339 + t264 * t369) * t196 - t220 * t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t199 - t332, -t200 ^ 2 - t370, t200 * t90 + t222 - t364, 0, 0, 0, 0, 0, t380, t234 - t334, t380, t387 * t209 + (-t48 + t277) * t205, t278 * t384 + t334 - t346, -t200 * t38 + t23 * t278 + (t22 * t384 + t1) * t205 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, -t105 ^ 2 + t371, -t387, -t281 + (-qJD(5) + t384) * t107, -t72, -t107 * t66 + t230 + t350, t105 * t66 + t28 * t384 - t372, -t105 * t65 - t219 + t350 - 0.2e1 * t368, pkin(5) * t47 - qJ(6) * t48 + (t23 - t29) * t107 + (t22 - t315) * t105, -0.2e1 * t352 - t105 * t38 + t107 * t65 + (0.2e1 * qJD(6) - t28) * t384 + t372, t1 * qJ(6) - t3 * pkin(5) - t38 * t65 - t22 * t29 - g(1) * (-pkin(5) * t120 + qJ(6) * t121) - g(2) * (pkin(5) * t122 - qJ(6) * t123) + t315 * t23 - t260 * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 + t335, -t387, -t371 - t279, -t23 * t384 + t219 + t368;];
tau_reg  = t4;
