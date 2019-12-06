% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:41
% EndTime: 2019-12-05 19:00:49
% DurationCPUTime: 3.96s
% Computational Cost: add. (8658->427), mult. (13431->543), div. (0->0), fcn. (9125->16), ass. (0->252)
t249 = cos(qJ(4));
t250 = cos(qJ(3));
t331 = t249 * t250;
t244 = sin(qJ(4));
t245 = sin(qJ(3));
t335 = t244 * t245;
t168 = -t331 + t335;
t253 = -pkin(8) - pkin(7);
t308 = qJD(3) * t253;
t174 = t245 * t308;
t175 = t250 * t308;
t197 = t253 * t245;
t232 = t250 * pkin(8);
t198 = t250 * pkin(7) + t232;
t251 = cos(qJ(2));
t324 = qJD(1) * t251;
t312 = pkin(1) * t324;
t318 = qJD(4) * t249;
t319 = qJD(4) * t244;
t345 = t168 * t312 + t249 * t174 + t244 * t175 + t197 * t318 - t198 * t319;
t127 = t244 * t197 + t249 * t198;
t332 = t249 * t245;
t169 = t244 * t250 + t332;
t344 = -t127 * qJD(4) + t169 * t312 - t244 * t174 + t249 * t175;
t241 = qJ(3) + qJ(4);
t227 = sin(t241);
t229 = cos(t241);
t242 = qJ(1) + qJ(2);
t228 = sin(t242);
t211 = g(2) * t228;
t230 = cos(t242);
t368 = -g(3) * t230 + t211;
t377 = -g(1) * t229 - t368 * t227;
t237 = qJD(1) + qJD(2);
t311 = t237 * t335;
t143 = -t237 * t331 + t311;
t145 = t169 * t237;
t243 = sin(qJ(5));
t248 = cos(qJ(5));
t274 = t243 * t143 - t248 * t145;
t92 = t248 * t143 + t243 * t145;
t376 = t92 * t274;
t236 = qJD(3) + qJD(4);
t276 = t236 * t335;
t320 = qJD(3) * t250;
t371 = -t249 * t320 - t250 * t318;
t113 = t276 + t371;
t365 = t113 * pkin(9);
t375 = t365 + t344;
t114 = t236 * t169;
t111 = t114 * pkin(9);
t374 = -t111 + t345;
t239 = t245 ^ 2;
t240 = t250 ^ 2;
t325 = t239 + t240;
t369 = t325 * t251;
t373 = t237 * t369;
t235 = qJDD(1) + qJDD(2);
t246 = sin(qJ(2));
t315 = qJDD(1) * t246;
t322 = qJD(2) * t251;
t148 = t235 * pkin(7) + (qJD(1) * t322 + t315) * pkin(1);
t295 = t325 * t148;
t281 = g(2) * t230 + g(3) * t228;
t22 = t274 ^ 2 - t92 ^ 2;
t226 = qJD(5) + t236;
t316 = qJD(5) * t248;
t317 = qJD(5) * t243;
t330 = t250 * t235;
t285 = -t235 * t332 + t371 * t237 - t244 * t330;
t67 = t237 * t276 + t285;
t333 = t245 * t235;
t278 = t244 * t333 - t249 * t330;
t68 = t114 * t237 + t278;
t23 = t143 * t316 + t145 * t317 + t243 * t68 + t248 * t67;
t18 = t92 * t226 - t23;
t231 = qJ(5) + t241;
t216 = sin(t231);
t217 = cos(t231);
t340 = t217 * t230;
t341 = t217 * t228;
t234 = qJDD(3) + qJDD(4);
t350 = pkin(1) * qJD(1);
t313 = t246 * t350;
t181 = t237 * pkin(7) + t313;
t305 = t237 * t320;
t73 = -t181 * t320 + qJDD(3) * pkin(3) - t245 * t148 + (-t305 - t333) * pkin(8);
t301 = pkin(8) * t237 + t181;
t133 = t301 * t250;
t124 = t249 * t133;
t132 = t301 * t245;
t125 = qJD(3) * pkin(3) - t132;
t77 = t244 * t125 + t124;
t321 = qJD(3) * t245;
t306 = t237 * t321;
t78 = -t181 * t321 + t250 * t148 + (-t306 + t330) * pkin(8);
t26 = -t77 * qJD(4) - t244 * t78 + t249 * t73;
t16 = t234 * pkin(4) + t67 * pkin(9) + t26;
t290 = -t125 * t318 + t133 * t319 - t244 * t73 - t249 * t78;
t17 = -t68 * pkin(9) - t290;
t141 = t145 * pkin(9);
t122 = t244 * t133;
t76 = t249 * t125 - t122;
t60 = -t141 + t76;
t58 = t236 * pkin(4) + t60;
t364 = t143 * pkin(9);
t61 = t77 - t364;
t4 = (qJD(5) * t58 + t17) * t248 + t243 * t16 - t61 * t317;
t359 = t250 * pkin(3);
t220 = pkin(2) + t359;
t146 = -t220 * t237 - t312;
t97 = t143 * pkin(4) + t146;
t259 = g(1) * t216 - g(2) * t341 + g(3) * t340 + t97 * t92 - t4;
t348 = t248 * t61;
t28 = t243 * t58 + t348;
t5 = -t28 * qJD(5) + t248 * t16 - t243 * t17;
t257 = -g(1) * t217 - t368 * t216 + t274 * t97 + t5;
t361 = t237 * pkin(2);
t182 = -t312 - t361;
t370 = t181 * t369 + t182 * t246;
t263 = t274 * qJD(5) + t243 * t67 - t248 * t68;
t19 = -t226 * t274 + t263;
t218 = t246 * pkin(1) + pkin(7);
t353 = -pkin(8) - t218;
t165 = t353 * t245;
t166 = t250 * t218 + t232;
t105 = t244 * t165 + t249 * t166;
t366 = g(1) * t250;
t363 = t169 * pkin(9);
t362 = t235 * pkin(2);
t247 = sin(qJ(1));
t360 = t247 * pkin(1);
t358 = t251 * pkin(1);
t252 = cos(qJ(1));
t357 = t252 * pkin(1);
t126 = t249 * t197 - t244 * t198;
t98 = t126 - t363;
t164 = t168 * pkin(9);
t99 = -t164 + t127;
t56 = -t243 * t99 + t248 * t98;
t352 = t56 * qJD(5) + t375 * t243 + t374 * t248;
t57 = t243 * t98 + t248 * t99;
t351 = -t57 * qJD(5) - t374 * t243 + t375 * t248;
t349 = t243 * t61;
t219 = t249 * pkin(3) + pkin(4);
t336 = t243 * t244;
t79 = t244 * t132 - t124;
t62 = t79 + t364;
t80 = -t249 * t132 - t122;
t63 = -t141 + t80;
t347 = -t243 * t62 - t248 * t63 + t219 * t316 + (-t244 * t317 + (t248 * t249 - t336) * qJD(4)) * pkin(3);
t334 = t244 * t248;
t346 = t243 * t63 - t248 * t62 - t219 * t317 + (-t244 * t316 + (-t243 * t249 - t334) * qJD(4)) * pkin(3);
t343 = t145 * t143;
t339 = t228 * t229;
t338 = t229 * t230;
t337 = t237 * t245;
t327 = -qJD(2) * t313 + qJDD(1) * t358;
t147 = -t327 - t362;
t329 = t147 * t245 + t182 * t320;
t326 = t239 - t240;
t323 = qJD(2) * t246;
t314 = pkin(1) * t322;
t223 = pkin(3) * t321;
t233 = t237 ^ 2;
t310 = t245 * t233 * t250;
t309 = t182 * t321 + t281 * t250;
t307 = t237 * t323;
t302 = pkin(4) * t229 + t359;
t100 = t114 * pkin(4) + t223;
t298 = qJD(3) * t353;
t294 = t325 * t235;
t104 = t249 * t165 - t244 * t166;
t173 = pkin(2) + t302;
t238 = -pkin(9) + t253;
t292 = -t230 * t173 + t228 * t238;
t291 = -t230 * t220 + t228 * t253;
t109 = t248 * t168 + t243 * t169;
t110 = -t243 * t168 + t248 * t169;
t46 = t110 * qJD(5) - t243 * t113 + t248 * t114;
t101 = pkin(3) * t306 - t220 * t235 - t327;
t47 = t68 * pkin(4) + t101;
t289 = g(2) * t340 + g(3) * t341 + t47 * t109 + t97 * t46;
t288 = g(2) * t338 + g(3) * t339 + t101 * t168 + t146 * t114;
t287 = t237 * t313;
t286 = t368 + t295;
t284 = t327 + t281;
t283 = t245 * t305;
t282 = t100 - t313;
t279 = g(2) * t252 + g(3) * t247;
t27 = t248 * t58 - t349;
t277 = -t27 * t92 - t274 * t28;
t87 = t104 - t363;
t88 = -t164 + t105;
t43 = -t243 * t88 + t248 * t87;
t44 = t243 * t87 + t248 * t88;
t45 = t248 * t113 + t243 * t114 + t168 * t316 + t169 * t317;
t275 = -t4 * t109 - t5 * t110 + t27 * t45 - t28 * t46 + t368;
t273 = -t228 * t173 - t230 * t238;
t272 = -t228 * t220 - t230 * t253;
t142 = t168 * pkin(4) - t220;
t271 = t76 * t113 - t77 * t114 + t168 * t290 - t26 * t169 + t368;
t128 = t245 * t298 + t250 * t314;
t129 = -t245 * t314 + t250 * t298;
t48 = t249 * t128 + t244 * t129 + t165 * t318 - t166 * t319;
t270 = -t313 + t223;
t269 = -t182 * t237 - t148 - t368;
t268 = t47 * t110 - t281 * t216 - t97 * t45;
t267 = t101 * t169 - t146 * t113 - t281 * t227;
t254 = qJD(3) ^ 2;
t266 = -pkin(7) * t254 + t287 + t362;
t221 = -pkin(2) - t358;
t265 = -pkin(1) * t307 - t218 * t254 - t221 * t235;
t264 = -pkin(7) * qJDD(3) + (t312 - t361) * qJD(3);
t49 = -t105 * qJD(4) - t244 * t128 + t249 * t129;
t262 = -qJDD(3) * t218 + (t221 * t237 - t314) * qJD(3);
t261 = g(1) * t227 - g(2) * t339 + g(3) * t338 + t146 * t143 + t290;
t256 = -t146 * t145 + t26 + t377;
t225 = qJDD(5) + t234;
t224 = pkin(1) * t323;
t208 = t230 * pkin(7);
t192 = -t220 - t358;
t191 = qJDD(3) * t250 - t254 * t245;
t190 = qJDD(3) * t245 + t254 * t250;
t176 = t224 + t223;
t157 = pkin(3) * t334 + t243 * t219;
t156 = -pkin(3) * t336 + t248 * t219;
t150 = t240 * t235 - 0.2e1 * t283;
t149 = t239 * t235 + 0.2e1 * t283;
t131 = t142 - t358;
t116 = -0.2e1 * t326 * t237 * qJD(3) + 0.2e1 * t245 * t330;
t115 = pkin(3) * t337 + t145 * pkin(4);
t96 = t100 + t224;
t84 = -t114 * t236 - t168 * t234;
t83 = -t113 * t236 + t169 * t234;
t71 = -t143 ^ 2 + t145 ^ 2;
t52 = -t285 + (t143 - t311) * t236;
t42 = t49 + t365;
t41 = -t111 + t48;
t40 = t143 * t114 + t68 * t168;
t39 = -t145 * t113 - t67 * t169;
t34 = -t109 * t225 - t46 * t226;
t33 = t110 * t225 - t45 * t226;
t30 = t248 * t60 - t349;
t29 = -t243 * t60 - t348;
t14 = t113 * t143 - t145 * t114 + t67 * t168 - t169 * t68;
t9 = -t109 * t263 + t92 * t46;
t8 = -t23 * t110 + t274 * t45;
t7 = -t44 * qJD(5) - t243 * t41 + t248 * t42;
t6 = t43 * qJD(5) + t243 * t42 + t248 * t41;
t1 = t23 * t109 + t110 * t263 + t274 * t46 + t45 * t92;
t2 = [0, 0, 0, 0, 0, qJDD(1), t279, -g(2) * t247 + g(3) * t252, 0, 0, 0, 0, 0, 0, 0, t235, (t235 * t251 - t307) * pkin(1) + t284, ((-qJDD(1) - t235) * t246 + (-qJD(1) - t237) * t322) * pkin(1) - t368, 0, (t279 + (t246 ^ 2 + t251 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t149, t116, t190, t150, t191, 0, t262 * t245 + (-t147 + t265) * t250 + t309, t262 * t250 + (-t265 - t281) * t245 + t329, pkin(1) * qJD(2) * t373 + t218 * t294 + t286, t147 * t221 - g(2) * (-t230 * pkin(2) - t228 * pkin(7)) - g(3) * (-t228 * pkin(2) + t208) + t218 * t295 + (t370 * qJD(2) + t279) * pkin(1), t39, t14, t83, t40, t84, 0, t104 * t234 + t176 * t143 + t192 * t68 + t49 * t236 + t288, -t105 * t234 + t176 * t145 - t192 * t67 - t48 * t236 + t267, t104 * t67 - t105 * t68 - t48 * t143 - t49 * t145 + t271, -t290 * t105 + t77 * t48 + t26 * t104 + t76 * t49 + t101 * t192 + t146 * t176 - g(2) * (t291 - t357) - g(3) * (t272 - t360), t8, t1, t33, t9, t34, 0, -t131 * t263 + t43 * t225 + t7 * t226 + t96 * t92 + t289, -t131 * t23 - t44 * t225 - t6 * t226 - t274 * t96 + t268, t43 * t23 + t263 * t44 + t274 * t7 - t6 * t92 + t275, t4 * t44 + t28 * t6 + t5 * t43 + t27 * t7 + t47 * t131 + t97 * t96 - g(2) * (t292 - t357) - g(3) * (t273 - t360); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t284 + t287, (-t315 + (-qJD(2) + t237) * t324) * pkin(1) - t368, 0, 0, t149, t116, t190, t150, t191, 0, t264 * t245 + (-t147 + t266) * t250 + t309, t264 * t250 + (-t266 - t281) * t245 + t329, pkin(7) * t294 - t350 * t373 + t286, -g(3) * t208 + (-t147 + t281) * pkin(2) + (t295 + t211) * pkin(7) - t370 * t350, t39, t14, t83, t40, t84, 0, t126 * t234 + t270 * t143 - t220 * t68 + t344 * t236 + t288, -t127 * t234 + t270 * t145 + t220 * t67 - t345 * t236 + t267, t126 * t67 - t127 * t68 - t345 * t143 - t344 * t145 + t271, -g(2) * t291 - g(3) * t272 - t101 * t220 + t26 * t126 - t127 * t290 + t270 * t146 + t344 * t76 + t345 * t77, t8, t1, t33, t9, t34, 0, -t142 * t263 + t56 * t225 + t351 * t226 + t282 * t92 + t289, -t142 * t23 - t57 * t225 - t352 * t226 - t274 * t282 + t268, t56 * t23 + t263 * t57 + t274 * t351 - t352 * t92 + t275, -g(2) * t292 - g(3) * t273 + t47 * t142 + t27 * t351 + t28 * t352 + t282 * t97 + t4 * t57 + t5 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, t326 * t233, t333, t310, t330, qJDD(3), t269 * t245 - t366, g(1) * t245 + t269 * t250, 0, 0, t343, t71, t52, -t343, -t278, t234, -t79 * t236 + (-t143 * t337 + t234 * t249 - t236 * t319) * pkin(3) + t256, t80 * t236 + (-t145 * t337 - t234 * t244 - t236 * t318) * pkin(3) + t261, (t77 + t79) * t145 + (-t76 + t80) * t143 + (-t244 * t68 + t249 * t67 + (-t143 * t249 + t145 * t244) * qJD(4)) * pkin(3), -t76 * t79 - t77 * t80 + (-t366 - t244 * t290 + t249 * t26 + (-t244 * t76 + t249 * t77) * qJD(4) + (-t146 * t237 - t368) * t245) * pkin(3), -t376, t22, t18, t376, t19, t225, -t115 * t92 + t156 * t225 + t346 * t226 + t257, t115 * t274 - t157 * t225 - t347 * t226 + t259, t156 * t23 + t157 * t263 + t274 * t346 - t347 * t92 + t277, t4 * t157 + t5 * t156 - t97 * t115 - g(1) * t302 + t347 * t28 + t346 * t27 + t368 * (-t245 * pkin(3) - pkin(4) * t227); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, t71, t52, -t343, -t278, t234, t77 * t236 + t256, t76 * t236 + t261, 0, 0, -t376, t22, t18, t376, t19, t225, -t29 * t226 + (-t145 * t92 + t225 * t248 - t226 * t317) * pkin(4) + t257, t30 * t226 + (t145 * t274 - t225 * t243 - t226 * t316) * pkin(4) + t259, -t29 * t274 + t30 * t92 + (t23 * t248 + t263 * t243 + (-t243 * t274 - t248 * t92) * qJD(5)) * pkin(4) + t277, -t27 * t29 - t28 * t30 + (-t145 * t97 + t243 * t4 + t248 * t5 + (-t243 * t27 + t248 * t28) * qJD(5) + t377) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t376, t22, t18, t376, t19, t225, t28 * t226 + t257, t27 * t226 + t259, 0, 0;];
tau_reg = t2;
