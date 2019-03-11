% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:40
% EndTime: 2019-03-09 08:35:48
% DurationCPUTime: 5.26s
% Computational Cost: add. (4460->560), mult. (8794->622), div. (0->0), fcn. (4868->6), ass. (0->303)
t195 = cos(qJ(5));
t192 = sin(qJ(5));
t300 = qJD(2) * t192;
t196 = cos(qJ(2));
t303 = qJD(1) * t196;
t102 = t195 * t303 + t300;
t193 = sin(qJ(2));
t304 = qJD(1) * t193;
t132 = qJD(5) + t304;
t332 = t102 * t132;
t160 = t196 * qJDD(1);
t290 = qJD(1) * qJD(2);
t271 = t193 * t290;
t55 = qJD(5) * t102 - t195 * qJDD(2) + (t160 - t271) * t192;
t375 = -t332 - t55;
t294 = t195 * qJD(2);
t101 = t192 * t303 - t294;
t234 = t101 * t132;
t295 = qJD(5) * t196;
t217 = t192 * t295 + t193 * t294;
t54 = qJD(1) * t217 - qJD(5) * t294 - t192 * qJDD(2) - t195 * t160;
t36 = t54 - t234;
t374 = t54 + t234;
t373 = t55 - t332;
t361 = pkin(2) + pkin(3);
t372 = t193 * t361;
t270 = t196 * t290;
t288 = t193 * qJDD(1);
t371 = t270 + t288;
t275 = t361 * qJD(2);
t272 = t361 * qJDD(2);
t342 = t54 * qJ(6);
t100 = qJDD(5) + t371;
t352 = t100 * pkin(5);
t291 = pkin(8) + t361;
t222 = pkin(4) * t196 - t193 * t291;
t213 = t222 * qJD(2);
t162 = t193 * qJD(3);
t187 = qJDD(1) * pkin(1);
t241 = pkin(2) * t160 + qJ(3) * t371 + qJD(1) * t162 + t187;
t218 = pkin(3) * t160 + qJDD(4) + t241;
t168 = t193 * pkin(4);
t249 = t196 * pkin(8) + t168;
t25 = qJD(1) * t213 + qJDD(1) * t249 + t218;
t20 = t195 * t25;
t91 = -qJD(1) * pkin(1) - pkin(2) * t303 - qJ(3) * t304;
t74 = pkin(3) * t303 + qJD(4) - t91;
t60 = qJD(1) * t249 + t74;
t143 = qJ(4) * t304;
t285 = pkin(7) * t304;
t247 = qJD(3) + t285;
t232 = -t143 + t247;
t72 = -qJD(2) * t291 + t232;
t27 = t192 * t60 + t195 * t72;
t131 = pkin(7) * t270;
t151 = pkin(7) * t288;
t269 = qJDD(3) + t131 + t151;
t289 = qJD(1) * qJD(4);
t203 = -qJ(4) * t371 - t193 * t289 + t269;
t44 = -qJDD(2) * t291 + t203;
t4 = -qJD(5) * t27 - t192 * t44 + t20;
t1 = t102 * qJD(6) - t342 + t352 + t4;
t26 = -t192 * t72 + t195 * t60;
t23 = qJ(6) * t102 + t26;
t15 = pkin(5) * t132 + t23;
t296 = qJD(5) * t195;
t297 = qJD(5) * t192;
t3 = t192 * t25 + t195 * t44 + t60 * t296 - t72 * t297;
t340 = t55 * qJ(6);
t2 = qJD(6) * t101 + t3 + t340;
t176 = g(3) * t196;
t197 = cos(qJ(1));
t321 = t193 * t197;
t194 = sin(qJ(1));
t323 = t193 * t194;
t277 = -g(1) * t321 - g(2) * t323 + t176;
t24 = qJ(6) * t101 + t27;
t345 = t24 * t132;
t370 = -t192 * (t1 + t345) - t195 * (t132 * t15 - t2) + t277;
t369 = t132 ^ 2;
t351 = pkin(7) - qJ(4);
t116 = t351 * t193;
t103 = t195 * t116;
t172 = t196 * pkin(3);
t164 = t193 * qJ(3);
t173 = t196 * pkin(2);
t308 = t173 + t164;
t278 = t172 + t308;
t239 = t278 + t249;
t71 = pkin(1) + t239;
t50 = t192 * t71 + t103;
t155 = pkin(7) * t303;
t106 = -qJ(4) * t303 + t155;
t183 = qJD(2) * qJ(3);
t89 = -t106 - t183;
t177 = g(2) * t194;
t367 = g(1) * t197 + t177;
t313 = t197 * t195;
t319 = t194 * t192;
t92 = t193 * t319 - t313;
t318 = t194 * t195;
t94 = -t192 * t321 - t318;
t366 = -g(1) * t94 + g(2) * t92;
t365 = -t55 * pkin(5) + qJDD(6);
t299 = qJD(2) * t193;
t227 = pkin(7) * t299 + t196 * qJD(4);
t124 = qJ(4) * t271;
t152 = pkin(7) * t160;
t180 = qJDD(2) * qJ(3);
t181 = qJD(2) * qJD(3);
t279 = t152 + t180 + t181;
t257 = -t124 - t279;
t48 = qJ(4) * t160 + qJD(1) * t227 + t257;
t110 = -pkin(1) - t308;
t337 = pkin(7) * qJDD(2);
t363 = (qJD(1) * t110 + t91) * qJD(2) - t337;
t362 = t102 ^ 2;
t357 = pkin(5) * t192;
t174 = t197 * pkin(7);
t356 = g(1) * t174;
t178 = g(1) * t194;
t355 = g(2) * qJ(4);
t354 = g(2) * t197;
t353 = g(3) * t193;
t150 = t195 * pkin(5) + pkin(4);
t350 = -t23 + t15;
t311 = qJ(6) + t291;
t258 = qJD(5) * t311;
t293 = t195 * qJD(6);
t145 = qJ(3) * t303;
t64 = qJD(1) * t222 + t145;
t41 = t195 * t106 + t192 * t64;
t349 = -t293 - t41 + (qJ(6) * t304 + t258) * t192;
t322 = t193 * t195;
t231 = pkin(5) * t196 - qJ(6) * t322;
t40 = -t192 * t106 + t195 * t64;
t348 = -qJD(1) * t231 + t192 * qJD(6) + t195 * t258 - t40;
t347 = qJD(2) * pkin(2);
t346 = t132 * t26;
t344 = t27 * t132;
t186 = qJDD(2) * pkin(4);
t45 = t186 - t48;
t343 = t45 * t192;
t341 = t54 * t192;
t339 = t55 * t195;
t276 = -pkin(7) + t357;
t338 = -pkin(5) * t297 - t276 * t304 + qJD(3) - t143;
t336 = qJ(3) * t196;
t335 = qJDD(2) * pkin(2);
t334 = t101 * t102;
t333 = t101 * t196;
t331 = t102 * t195;
t330 = t102 * t196;
t329 = t150 * t193;
t185 = t196 ^ 2;
t200 = qJD(1) ^ 2;
t328 = t185 * t200;
t190 = -qJ(6) - pkin(8);
t327 = t190 * t196;
t326 = t192 * t100;
t325 = t192 * t193;
t324 = t192 * t196;
t320 = t193 * t200;
t317 = t194 * t196;
t316 = t195 * t100;
t315 = t195 * t196;
t314 = t196 * t197;
t312 = t89 * qJD(2);
t310 = -qJD(4) - t74;
t298 = qJD(2) * t196;
t309 = qJ(3) * t298 + t162;
t307 = t197 * pkin(1) + t194 * pkin(7);
t184 = t193 ^ 2;
t305 = t184 + t185;
t302 = qJD(2) * t101;
t301 = qJD(2) * t102;
t104 = -t143 + t285;
t292 = qJD(3) + t104;
t287 = -t190 + t361;
t59 = t213 + t309;
t86 = -t193 * qJD(4) + t298 * t351;
t286 = t192 * t59 + t195 * t86 + t71 * t296;
t284 = t132 * t325;
t283 = t132 * t322;
t263 = -pkin(5) * t101 + qJD(6);
t82 = qJD(2) * pkin(4) - t89;
t58 = t263 + t82;
t282 = t58 * t297;
t281 = g(3) * t325 + (g(1) * t314 + g(2) * t317) * t192;
t280 = t152 + 0.2e1 * t180 + 0.2e1 * t181;
t274 = t195 * t295;
t273 = t132 * t303;
t268 = -pkin(1) - t164;
t267 = qJ(4) + t357;
t266 = t178 - t354;
t265 = -t192 * t86 + t195 * t59;
t264 = g(1) * t291;
t96 = pkin(1) + t278;
t262 = qJD(1) * t96 + t74;
t49 = -t116 * t192 + t195 * t71;
t261 = t310 * t193;
t259 = g(1) * t287;
t256 = pkin(2) * t314 + qJ(3) * t321 + t307;
t255 = t151 + t277;
t254 = t193 * t275;
t253 = t193 * t270;
t252 = -g(1) * t92 - g(2) * t94;
t93 = -t192 * t197 - t193 * t318;
t95 = t193 * t313 - t319;
t251 = -g(1) * t93 - g(2) * t95;
t250 = t305 * qJDD(1) * pkin(7);
t199 = qJD(2) ^ 2;
t248 = pkin(7) * t199 + t354;
t245 = g(1) * (-t197 * qJ(4) + t174);
t244 = -t15 * t192 + t195 * t24;
t243 = -t192 * t27 - t195 * t26;
t242 = -t192 * t26 + t195 * t27;
t240 = pkin(3) * t314 + t256;
t109 = t247 - t347;
t114 = t155 + t183;
t238 = t109 * t196 - t114 * t193;
t237 = -qJDD(3) - t255;
t233 = t268 - t173;
t18 = t45 + t365;
t229 = -t18 * t192 - t296 * t58;
t228 = -0.2e1 * pkin(1) * t290 - t337;
t226 = g(2) * t240;
t225 = t131 - t237;
t224 = -t132 * t296 - t326;
t223 = t132 * t297 - t316;
t219 = -t248 + 0.2e1 * t187;
t85 = -qJ(4) * t299 + t227;
t216 = t192 * t299 - t274;
t215 = -qJ(4) * qJDD(1) - t367;
t212 = t100 * t291 - t132 * t82;
t39 = -qJD(1) * t254 + t218;
t73 = -t254 + t309;
t211 = -qJD(1) * t73 - qJDD(1) * t96 + t354 - t39;
t210 = -t196 * t367 - t353;
t209 = -qJ(4) * t288 + t225;
t57 = pkin(2) * t271 - t241;
t87 = pkin(2) * t299 - t309;
t208 = -qJD(1) * t87 - qJDD(1) * t110 - t248 - t57;
t207 = g(1) * t95 - g(2) * t93 - g(3) * t315 - t3;
t75 = -pkin(7) * t271 + t279;
t83 = t269 - t335;
t206 = qJD(2) * t238 + t83 * t193 + t75 * t196;
t205 = qJD(5) * t244 + t1 * t195 + t2 * t192 - t354;
t204 = qJD(5) * t242 + t3 * t192 + t4 * t195 - t354;
t202 = t20 + (-qJD(5) * t60 - t176 - t44) * t192 - t72 * t296 + t366;
t191 = qJ(3) + pkin(4);
t165 = t196 * qJ(4);
t161 = t184 * t200;
t138 = g(1) * t317;
t137 = g(1) * t323;
t130 = qJ(3) + t150;
t129 = qJ(3) * t314;
t127 = qJ(3) * t317;
t126 = t196 * t320;
t122 = -t161 - t199;
t117 = pkin(7) * t196 - t165;
t115 = qJDD(2) + t126;
t113 = -t161 + t328;
t112 = qJDD(2) * t196 - t193 * t199;
t111 = qJDD(2) * t193 + t196 * t199;
t108 = t311 * t195;
t107 = t311 * t192;
t105 = pkin(2) * t304 - t145;
t99 = t101 ^ 2;
t98 = qJDD(1) * t185 - 0.2e1 * t253;
t97 = qJDD(1) * t184 + 0.2e1 * t253;
t88 = -t196 * t276 - t165;
t84 = -t304 * t361 + t145;
t79 = -t275 + t232;
t77 = t193 * t160 + (-t184 + t185) * t290;
t70 = 0.2e1 * t77;
t65 = t100 * t193 + t132 * t298;
t61 = pkin(5) * t216 - t85;
t47 = -t99 + t362;
t46 = -t272 + t203;
t43 = qJ(6) * t324 + t50;
t38 = pkin(5) * t193 + qJ(6) * t315 + t49;
t34 = -t195 * t369 + t302 - t326;
t33 = t192 * t369 + t301 - t316;
t32 = (-t283 - t330) * qJD(1) + t224;
t31 = (-t283 + t330) * qJD(1) + t224;
t30 = (-t284 - t333) * qJD(1) - t223;
t29 = (t284 - t333) * qJD(1) + t223;
t22 = t192 * t234 - t339;
t21 = t132 * t331 - t341;
t17 = -t101 * t216 + t324 * t55;
t16 = -t102 * t217 - t315 * t54;
t14 = -qJD(5) * t50 + t265;
t13 = -t116 * t297 + t286;
t12 = (-t132 * t300 + t55) * t193 + (-t224 + t302) * t196;
t11 = (t132 * t294 + t54) * t193 + (t223 - t301) * t196;
t10 = qJ(6) * t274 + (-qJ(6) * t299 - qJD(5) * t116 + qJD(6) * t196) * t192 + t286;
t9 = t196 * t293 + t231 * qJD(2) + (-t103 + (-qJ(6) * t196 - t71) * t192) * qJD(5) + t265;
t8 = t192 * t36 + t195 * t373;
t7 = t192 * t375 - t374 * t195;
t6 = t192 * t373 - t195 * t36;
t5 = (t101 * t195 + t102 * t192) * t299 + (t341 - t339 + (t101 * t192 - t331) * qJD(5)) * t196;
t19 = [0, 0, 0, 0, 0, qJDD(1), t266, t367, 0, 0, t97, t70, t111, t98, t112, 0, t193 * t228 + t196 * t219 + t138, -t193 * t219 + t196 * t228 - t137, 0.2e1 * t250 - t367, -g(1) * (-t194 * pkin(1) + t174) - g(2) * t307 + (pkin(7) ^ 2 * t305 + pkin(1) ^ 2) * qJDD(1), t97, t111, -0.2e1 * t77, 0, -t112, t98, t193 * t363 + t208 * t196 + t138, t250 + t206 - t367, t208 * t193 - t196 * t363 + t137, pkin(7) * t206 - g(2) * t256 + t57 * t110 - t178 * t233 + t91 * t87 - t356, t98, t70, t112, t97, t111, 0, t117 * qJDD(2) + t137 + (t196 * t262 - t85) * qJD(2) - t211 * t193, t116 * qJDD(2) - t138 + (t193 * t262 + t86) * qJD(2) + t211 * t196 (-qJD(2) * t79 - qJDD(1) * t117 + t48 + (-qJD(2) * t116 + t85) * qJD(1)) * t196 + (-t312 - qJDD(1) * t116 - t46 + (qJD(2) * t117 - t86) * qJD(1)) * t193 + t367, t46 * t116 + t79 * t86 - t48 * t117 + t89 * t85 + t39 * t96 + t74 * t73 - t245 - t226 + (-g(1) * (t233 - t172) + t355) * t194, t16, t5, t11, t17, t12, t65, t49 * t100 + t85 * t101 - t117 * t55 + t14 * t132 + (t300 * t82 + t4) * t193 + (qJD(2) * t26 - t296 * t82 - t343) * t196 + t251, -t50 * t100 + t85 * t102 + t117 * t54 - t13 * t132 + (t294 * t82 - t3) * t193 + (-qJD(2) * t27 - t45 * t195 + t297 * t82) * t196 + t252, t13 * t101 + t14 * t102 + t196 * t204 + t243 * t299 - t49 * t54 + t50 * t55 + t138, t3 * t50 + t27 * t13 + t4 * t49 + t26 * t14 + t45 * t117 - t82 * t85 - t245 - g(2) * (pkin(4) * t321 + pkin(8) * t314 + t240) + (-g(1) * (t268 - t168) + t355 + t196 * t264) * t194, t16, t5, t11, t17, t12, t65, t38 * t100 - t61 * t101 + t9 * t132 - t88 * t55 + (t300 * t58 + t1) * t193 + (qJD(2) * t15 + t229) * t196 + t251, -t10 * t132 - t43 * t100 - t61 * t102 + t88 * t54 + (t294 * t58 - t2) * t193 + (-qJD(2) * t24 - t18 * t195 + t282) * t196 + t252, t10 * t101 + t9 * t102 - t38 * t54 + t43 * t55 + t138 + (-t15 * t195 - t192 * t24) * t299 + t205 * t196, t2 * t43 + t24 * t10 + t1 * t38 + t15 * t9 + t18 * t88 + t58 * t61 - t356 - t226 + (g(1) * t267 - g(2) * (-t327 + t329)) * t197 + (-g(1) * (t268 - t329) + g(2) * t267 + t196 * t259) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t113, t288, t126, t160, qJDD(2), pkin(1) * t320 - t255, t353 - t152 + (pkin(1) * t200 + t367) * t196, 0, 0, -t126, t288, t113, qJDD(2), -t160, t126, 0.2e1 * t335 + (t105 * t196 - t193 * t91) * qJD(1) + t237 (-pkin(2) * t193 + t336) * qJDD(1) + ((t114 - t183) * t193 + (qJD(3) - t109 - t347) * t196) * qJD(1) (qJD(1) * t105 - g(3)) * t193 + (qJD(1) * t91 - t367) * t196 + t280, t75 * qJ(3) + t114 * qJD(3) - t83 * pkin(2) - t91 * t105 - g(1) * (-pkin(2) * t321 + t129) - g(2) * (-pkin(2) * t323 + t127) - g(3) * t308 - t238 * qJD(1) * pkin(7), t126, -t113, t160, -t126, t288, qJDD(2), t104 * qJD(2) + t124 + (-g(3) + (-pkin(7) * qJD(2) - t84) * qJD(1)) * t193 + (qJD(1) * t310 + t215) * t196 + t280, -t106 * qJD(2) - 0.2e1 * t272 + ((-qJ(4) * qJD(2) + t84) * t196 + t261) * qJD(1) + t209 (-t336 + t372) * qJDD(1) + (-t292 + t79 + t275) * t303, -g(1) * t129 - g(2) * t127 - g(3) * t278 - t48 * qJ(3) - t79 * t106 - t292 * t89 - t361 * t46 + t367 * t372 - t74 * t84, t21, t7, t31, t22, t29, -t273, -t26 * t303 - t40 * t132 - t191 * t55 - t292 * t101 + t212 * t192 + (qJD(5) * t132 * t291 + t210 + t45) * t195, t27 * t303 + t191 * t54 - t343 + (-t291 * t297 + t41) * t132 - t292 * t102 + t212 * t195 + t281, -t41 * t101 - t40 * t102 + (t26 * t304 - t291 * t55 - t3 + (t102 * t291 + t26) * qJD(5)) * t195 + (t27 * t304 - t291 * t54 + t4 + (t101 * t291 + t27) * qJD(5)) * t192 - t277, t45 * t191 - t27 * t41 - t26 * t40 - g(1) * (pkin(4) * t314 + t129) - g(2) * (pkin(4) * t317 + t127) - g(3) * t239 + t292 * t82 + (t177 * t291 + t197 * t264) * t193 - (qJD(5) * t243 - t4 * t192 + t3 * t195) * t291, t21, t7, t31, t22, t29, -t273, -t282 + t107 * t100 - t130 * t55 + t348 * t132 - t338 * t101 + (-t15 * t196 - t325 * t58) * qJD(1) + (t18 + t210) * t195, t108 * t100 + t130 * t54 - t349 * t132 - t338 * t102 + (t196 * t24 - t322 * t58) * qJD(1) + t229 + t281, t349 * t101 + t348 * t102 - t107 * t54 - t108 * t55 - t370, -t2 * t108 + t1 * t107 + t18 * t130 - g(1) * (t150 * t314 + t129) - g(2) * (t150 * t317 + t127) - g(3) * (t278 - t327) + t338 * t58 + t349 * t24 + t348 * t15 + (-g(3) * t150 + t177 * t287 + t197 * t259) * t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t288, t122, -qJD(2) * t114 + t304 * t91 + t225 - t335, 0, 0, 0, 0, 0, 0, t122, t115, -t288, t312 - t272 + (-qJ(4) * t298 + t261) * qJD(1) + t209, 0, 0, 0, 0, 0, 0, t34, t33, t8, -t82 * qJD(2) + (t3 - t346) * t195 + (-t4 - t344) * t192 + t277, 0, 0, 0, 0, 0, 0, t34, t33, t8, -t58 * qJD(2) + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t270 + t288, -t160 + 0.2e1 * t271, -t161 - t328 (-t196 * t89 + (t79 - t275) * t193) * qJD(1) + t218 + t266, 0, 0, 0, 0, 0, 0, t30, t32, t6, t178 + (t193 * t242 + t196 * t82) * qJD(1) + t204, 0, 0, 0, 0, 0, 0, t30, t32, t6, t178 + (t193 * t244 + t196 * t58) * qJD(1) + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t47, t36, -t334, t373, t100, t82 * t102 + t202 + t344, -t101 * t82 + t207 + t346, 0, 0, t334, t47, t36, -t334, t373, t100, 0.2e1 * t352 - t342 + t345 + (t263 + t58) * t102 + t202, -t362 * pkin(5) - t340 + t23 * t132 + (-qJD(6) - t58) * t101 + t207, -t54 * pkin(5) + t101 * t350, t350 * t24 + (-g(3) * t324 + t58 * t102 + t1 + t366) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, t374, -t99 - t362, -t24 * t101 - t15 * t102 + t186 + (-pkin(7) * t290 - g(3)) * t193 + (t215 - t289) * t196 - t257 + t365;];
tau_reg  = t19;
