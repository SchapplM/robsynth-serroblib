% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:21
% EndTime: 2019-03-09 10:06:34
% DurationCPUTime: 5.88s
% Computational Cost: add. (4198->530), mult. (8421->612), div. (0->0), fcn. (4952->6), ass. (0->267)
t195 = sin(qJ(2));
t306 = qJD(1) * t195;
t380 = -t306 - qJD(4);
t197 = cos(qJ(4));
t194 = sin(qJ(4));
t302 = qJD(2) * t194;
t198 = cos(qJ(2));
t305 = qJD(1) * t198;
t109 = t197 * t305 + t302;
t239 = t109 * t380;
t289 = qJD(1) * qJD(2);
t270 = t195 * t289;
t287 = t198 * qJDD(1);
t373 = -t270 + t287;
t45 = qJD(4) * t109 - t197 * qJDD(2) + t194 * t373;
t369 = t45 + t239;
t379 = t45 - t239;
t271 = t194 * t305;
t300 = qJD(2) * t197;
t111 = -t271 + t300;
t331 = t111 * t380;
t46 = -qJD(4) * t271 + t194 * qJDD(2) + (qJD(2) * qJD(4) + t373) * t197;
t378 = -t46 + t331;
t375 = t46 + t331;
t377 = t194 * t375 - t369 * t197;
t357 = pkin(3) + pkin(7);
t269 = t198 * t289;
t288 = t195 * qJDD(1);
t225 = t269 + t288;
t108 = qJDD(4) + t225;
t376 = t108 * qJ(5) - qJD(5) * t380;
t358 = pkin(2) + pkin(8);
t316 = qJ(6) - t358;
t168 = pkin(7) * t306;
t290 = pkin(3) * t306 + qJD(3) + t168;
t372 = t46 * qJ(6) + t109 * qJD(6);
t297 = qJD(4) * t197;
t86 = t194 * t108;
t371 = t297 * t380 - t86;
t370 = 0.2e1 * t376;
t298 = qJD(4) * t194;
t119 = t380 * t298;
t87 = t197 * t108;
t368 = t87 + t119;
t304 = qJD(2) * t109;
t367 = t87 - t304;
t177 = t195 * qJ(3);
t267 = -pkin(1) - t177;
t224 = -t198 * t358 + t267;
t69 = t224 * qJD(1);
t72 = -qJD(2) * t358 + t290;
t31 = -t194 * t69 + t197 * t72;
t314 = qJD(5) - t31;
t107 = t111 ^ 2;
t141 = t380 ^ 2;
t366 = -t141 - t107;
t365 = -t197 * qJD(5) + t290;
t196 = sin(qJ(1));
t199 = cos(qJ(1));
t253 = g(1) * t199 + g(2) * t196;
t99 = t108 * pkin(4);
t364 = t99 - qJDD(5);
t187 = g(3) * t195;
t319 = t198 * t199;
t324 = t196 * t198;
t363 = -g(1) * t319 - g(2) * t324 - t187;
t170 = pkin(7) * t305;
t117 = pkin(3) * t305 + t170;
t191 = qJD(2) * qJ(3);
t89 = t191 + t117;
t236 = t111 * qJ(5) - t89;
t356 = pkin(4) + pkin(5);
t21 = -t109 * t356 + qJD(6) + t236;
t149 = pkin(2) * t270;
t336 = qJ(3) * t198;
t248 = pkin(8) * t195 - t336;
t292 = t195 * qJD(3);
t214 = qJD(2) * t248 - t292;
t30 = qJD(1) * t214 + qJDD(1) * t224 + t149;
t147 = pkin(7) * t269;
t165 = pkin(7) * t288;
t268 = qJDD(3) + t147 + t165;
t50 = pkin(3) * t225 - qJDD(2) * t358 + t268;
t264 = -t194 * t30 + t197 * t50 - t69 * t297 - t72 * t298;
t322 = t197 * t198;
t321 = t197 * t199;
t91 = t194 * t196 - t195 * t321;
t325 = t196 * t197;
t93 = t194 * t199 + t195 * t325;
t218 = g(1) * t91 - g(2) * t93 + g(3) * t322 + t264;
t212 = t218 + t364;
t343 = qJ(6) * t45;
t362 = (qJD(6) + t21) * t111 + t212 - t343;
t323 = t197 * qJ(5);
t361 = t194 * t356 - t323;
t335 = qJ(5) * t194;
t230 = t197 * t356 + t335;
t360 = -0.2e1 * pkin(1);
t359 = t109 ^ 2;
t355 = pkin(5) * t108;
t202 = qJD(2) ^ 2;
t354 = pkin(7) * t202;
t353 = g(1) * t196;
t188 = g(3) * t198;
t166 = pkin(7) * t287;
t189 = qJDD(2) * qJ(3);
t190 = qJD(2) * qJD(3);
t70 = pkin(7) * t270 - t166 - t189 - t190;
t51 = pkin(3) * t373 - t70;
t7 = t46 * pkin(4) + t45 * qJ(5) - t111 * qJD(5) + t51;
t350 = t197 * t7;
t183 = t198 * pkin(2);
t32 = t194 * t72 + t197 * t69;
t349 = t230 * t380 - t365;
t173 = pkin(2) * t306;
t78 = qJD(1) * t248 + t173;
t348 = t194 * t117 + t197 * t78;
t260 = qJD(4) * t316;
t67 = t194 * t78;
t263 = -t197 * t117 + t67;
t291 = t197 * qJD(6);
t330 = t194 * t195;
t347 = t194 * t260 - t291 - (-qJ(6) * t330 - t198 * t356) * qJD(1) - t263;
t274 = t197 * t306;
t293 = t194 * qJD(6);
t36 = qJ(5) * t305 + t348;
t346 = qJ(6) * t274 + t197 * t260 + t293 - t36;
t251 = pkin(4) * t197 + t335;
t345 = -t251 * t380 + t365;
t344 = qJ(5) * t46;
t138 = t380 * qJ(5);
t25 = -t138 + t32;
t342 = t380 * t25;
t341 = t380 * t32;
t340 = t197 * t45;
t339 = t51 * t197;
t308 = t183 + t177;
t277 = -t198 * pkin(8) - t308;
t102 = -pkin(1) + t277;
t132 = t357 * t195;
t338 = t197 * t102 + t194 * t132;
t337 = pkin(7) * qJDD(2);
t334 = qJDD(2) * pkin(2);
t333 = t109 * qJ(5);
t332 = t111 * t109;
t329 = t194 * t198;
t328 = t195 * t196;
t327 = t195 * t197;
t326 = t195 * t199;
t320 = t198 * qJ(6);
t203 = qJD(1) ^ 2;
t318 = t198 * t203;
t317 = t358 * t108;
t19 = t111 * qJ(6) + t31;
t315 = qJD(5) - t19;
t143 = qJ(3) * t324;
t284 = pkin(4) * t329;
t312 = t196 * t284 + t143;
t145 = qJ(3) * t319;
t311 = t199 * t284 + t145;
t310 = g(1) * t324 - g(2) * t319;
t133 = t357 * t198;
t192 = t195 ^ 2;
t193 = t198 ^ 2;
t307 = t192 - t193;
t303 = qJD(2) * t111;
t301 = qJD(2) * t195;
t299 = qJD(2) * t198;
t296 = qJD(4) * t198;
t295 = qJD(4) * t358;
t294 = qJD(5) * t194;
t285 = -t194 * t50 - t197 * t30 - t72 * t297;
t283 = -t274 * t380 - t371;
t48 = t195 * qJ(5) + t338;
t282 = t21 * t297;
t281 = t195 * t318;
t280 = t196 * t322;
t279 = t197 * t319;
t278 = g(1) * t279 + g(2) * t280 + g(3) * t327;
t276 = -g(1) * t326 - g(2) * t328 + t188;
t275 = t194 * t306;
t273 = t194 * t301;
t272 = t197 * t295;
t92 = t194 * t326 + t325;
t266 = -t91 * pkin(4) + qJ(5) * t92;
t94 = -t194 * t328 + t321;
t265 = t93 * pkin(4) - qJ(5) * t94;
t262 = -t194 * t102 + t197 * t132;
t261 = -qJD(2) * pkin(2) + qJD(3);
t258 = t199 * pkin(1) + pkin(2) * t319 + t196 * pkin(7) + qJ(3) * t326;
t257 = -t165 - t276;
t256 = t166 + t363;
t255 = -g(1) * t93 - g(2) * t91;
t254 = -g(1) * t94 - g(2) * t92;
t20 = qJ(6) * t109 + t32;
t252 = -g(2) * t199 + t353;
t250 = -pkin(4) * t194 + t323;
t249 = pkin(5) * t194 - t323;
t22 = pkin(4) * t380 + t314;
t247 = t194 * t22 + t197 * t25;
t123 = t168 + t261;
t126 = -t170 - t191;
t245 = t123 * t198 + t126 * t195;
t244 = g(3) * (pkin(4) * t330 - t277);
t184 = t199 * pkin(7);
t243 = t199 * pkin(3) + t94 * pkin(4) + t93 * qJ(5) + t184;
t118 = t357 * t299;
t172 = pkin(2) * t301;
t64 = t172 + t214;
t242 = -t102 * t297 + t118 * t197 - t132 * t298 - t194 * t64;
t240 = t380 * t194;
t238 = qJDD(2) * t195 + t198 * t202;
t6 = -t264 - t364;
t237 = t267 - t183;
t235 = t283 + t303;
t4 = -pkin(5) * t46 + qJDD(6) - t7;
t233 = t197 * t4 - t21 * t298;
t90 = t237 * qJD(1);
t229 = t90 * t306 + qJDD(3) - t257;
t228 = -t298 * t69 - t285;
t227 = -qJ(3) * t299 - t292;
t226 = -t102 * t298 + t194 * t118 + t132 * t297 + t197 * t64;
t223 = t275 * t380 + t119 + t367;
t222 = t230 * t198;
t5 = t228 + t376;
t125 = -pkin(1) - t308;
t221 = t337 + (-qJD(1) * t125 - t90) * qJD(2);
t34 = pkin(4) * t109 - t236;
t220 = -t34 * t380 - t317;
t11 = qJ(5) * t299 + t195 * qJD(5) + t226;
t47 = qJD(1) * t227 + qJDD(1) * t237 + t149;
t81 = t172 + t227;
t219 = qJD(1) * t81 + qJDD(1) * t125 + t354 + t47;
t217 = -t198 * t253 - t187;
t216 = t196 * pkin(3) + t92 * pkin(4) + pkin(8) * t319 + qJ(5) * t91 + t258;
t1 = -qJD(6) * t111 + t343 - t355 + t6;
t14 = t356 * t380 + t315;
t17 = -t138 + t20;
t3 = t5 + t372;
t215 = -t1 * t197 + t3 * t194 + t276 + (t274 + t297) * t17 + (t275 + t298) * t14;
t213 = -t108 + t332;
t211 = g(1) * t92 - g(2) * t94 - g(3) * t329 - t228;
t77 = t268 - t334;
t210 = qJD(2) * t245 + t77 * t195 - t70 * t198;
t209 = -t295 * t380 + t217;
t206 = t111 * t34 - t212;
t205 = -t31 * t380 + t211;
t124 = qJ(3) - t250;
t122 = t316 * t197;
t121 = t316 * t194;
t116 = t357 * t301;
t114 = -qJ(3) * t305 + t173;
t98 = -qJ(3) - t361;
t63 = t198 * t251 + t133;
t53 = pkin(4) * t111 + t333;
t52 = -t222 - t133;
t49 = -pkin(4) * t195 - t262;
t38 = -pkin(4) * t305 + t263;
t37 = t197 * t320 + t48;
t35 = -t111 * t356 - t333;
t33 = t194 * t320 - t195 * t356 - t262;
t26 = (qJD(4) * t250 + t294) * t198 + (-t251 - t357) * t301;
t18 = (qJD(4) * t361 - t294) * t198 + (t230 + t357) * t301;
t13 = -pkin(4) * t299 - t242;
t9 = t198 * t291 + (-t194 * t296 - t195 * t300) * qJ(6) + t11;
t8 = -qJ(6) * t273 + (qJ(6) * t297 - qJD(2) * t356 + t293) * t198 - t242;
t2 = [qJDD(1), t252, t253, qJDD(1) * t192 + 0.2e1 * t195 * t269, 0.2e1 * t195 * t287 - 0.2e1 * t289 * t307, t238, qJDD(2) * t198 - t195 * t202, 0, 0.2e1 * pkin(1) * t373 - t238 * pkin(7) + t310 (t289 * t360 - t337) * t198 + (qJDD(1) * t360 - t252 + t354) * t195 (t192 + t193) * qJDD(1) * pkin(7) + t210 - t253, t195 * t221 + t198 * t219 - t310, t221 * t198 + (-t219 + t252) * t195, pkin(7) * t210 - g(1) * t184 - g(2) * t258 + t47 * t125 - t237 * t353 + t90 * t81, t45 * t329 + (-t197 * t296 + t273) * t111 (-t109 * t194 + t111 * t197) * t301 + (t194 * t46 + t340 + (t109 * t197 + t111 * t194) * qJD(4)) * t198 (-t302 * t380 - t45) * t195 + (t303 + t371) * t198 (-t300 * t380 - t46) * t195 + (-t304 - t368) * t198, t108 * t195 - t299 * t380, -t242 * t380 + t262 * t108 - t116 * t109 + t133 * t46 + (-t300 * t89 + t264) * t195 + (qJD(2) * t31 - t298 * t89 + t339) * t198 + t254, t226 * t380 - t338 * t108 - t116 * t111 - t133 * t45 + ((qJD(2) * t89 + qJD(4) * t69) * t194 + t285) * t195 + (-qJD(2) * t32 - t51 * t194 - t297 * t89) * t198 - t255, -t108 * t49 + t109 * t26 + t13 * t380 + t46 * t63 + (-t300 * t34 - t6) * t195 + (-qJD(2) * t22 - t298 * t34 + t350) * t198 + t254, -t109 * t11 + t111 * t13 - t45 * t49 - t46 * t48 + t247 * t301 + (-t194 * t6 - t197 * t5 + (t194 * t25 - t197 * t22) * qJD(4)) * t198 + t310, t108 * t48 - t11 * t380 - t111 * t26 + t45 * t63 + (-t302 * t34 + t5) * t195 + (qJD(2) * t25 + t194 * t7 + t297 * t34) * t198 + t255, -g(1) * t243 - g(2) * t216 + t25 * t11 + t22 * t13 - t224 * t353 + t34 * t26 + t5 * t48 + t6 * t49 + t7 * t63, -t108 * t33 - t109 * t18 + t380 * t8 - t46 * t52 + (t21 * t300 - t1) * t195 + (-qJD(2) * t14 - t233) * t198 + t254, t108 * t37 + t111 * t18 - t380 * t9 - t45 * t52 + (t21 * t302 + t3) * t195 + (qJD(2) * t17 - t194 * t4 - t282) * t198 + t255, t109 * t9 - t111 * t8 + t33 * t45 + t37 * t46 + (-t14 * t194 - t17 * t197) * t301 + (t1 * t194 + t197 * t3 + (t14 * t197 - t17 * t194) * qJD(4)) * t198 - t310, t3 * t37 + t17 * t9 + t1 * t33 + t14 * t8 + t4 * t52 + t21 * t18 - g(1) * (pkin(5) * t94 + t243) - g(2) * (pkin(5) * t92 - qJ(6) * t319 + t216) - (t198 * t316 + t267) * t353; 0, 0, 0, -t281, t307 * t203, t288, t287, qJDD(2), pkin(1) * t195 * t203 + t257, pkin(1) * t318 - t256 (-pkin(2) * t195 + t336) * qJDD(1) + ((-t126 - t191) * t195 + (-t123 + t261) * t198) * qJD(1), -t114 * t305 + t229 - 0.2e1 * t334, 0.2e1 * t189 + 0.2e1 * t190 + (t114 * t195 + t198 * t90) * qJD(1) + t256, -t70 * qJ(3) - t126 * qJD(3) - t77 * pkin(2) - t90 * t114 - g(1) * (-pkin(2) * t326 + t145) - g(2) * (-pkin(2) * t328 + t143) - g(3) * t308 - t245 * qJD(1) * pkin(7), t111 * t240 - t340, t379 * t194 + t378 * t197 (-t111 * t198 + t330 * t380) * qJD(1) + t368, t109 * t305 - t283, t380 * t305, -t31 * t305 + qJ(3) * t46 + t290 * t109 - t317 * t197 + (t51 + t209) * t194 - (t67 + (-t117 + t89) * t197) * t380, t32 * t305 - qJ(3) * t45 + t339 - (t272 + t348) * t380 + t290 * t111 + (t380 * t89 + t317) * t194 - t278, t22 * t305 + t124 * t46 - t380 * t38 + t345 * t109 + t220 * t197 + (t7 + t209) * t194, t109 * t36 - t111 * t38 + (-t25 * t306 - t358 * t45 + t6 + (t109 * t358 - t25) * qJD(4)) * t197 + (-t22 * t306 + t358 * t46 - t5 + (-t111 * t358 - t22) * qJD(4)) * t194 - t276, -t25 * t305 + t124 * t45 - t350 - (-t36 - t272) * t380 - t345 * t111 + t220 * t194 + t278, t7 * t124 - t25 * t36 - t22 * t38 - g(1) * (-qJ(5) * t279 + t311) - g(2) * (-qJ(5) * t280 + t312) - t244 + t345 * t34 + (g(3) * t323 + t253 * t358) * t195 - (qJD(4) * t247 + t5 * t194 - t6 * t197) * t358, -t282 + t108 * t122 - t46 * t98 + t347 * t380 - t349 * t109 + (t14 * t198 - t21 * t327) * qJD(1) + (-t4 + t217) * t194, t108 * t121 - t45 * t98 - t346 * t380 + t349 * t111 + (-t17 * t198 - t21 * t330) * qJD(1) + t233 + t278, t109 * t346 - t111 * t347 + t121 * t46 - t122 * t45 + t215, t3 * t121 - t1 * t122 + t4 * t98 - g(1) * t311 - g(2) * t312 - t244 + t349 * t21 + (g(3) * qJ(6) - t249 * t253) * t198 + t346 * t17 + t347 * t14 + (-g(3) * t249 - t253 * t316) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, qJDD(2) + t281, -t192 * t203 - t202, qJD(2) * t126 + t147 + t229 - t334, 0, 0, 0, 0, 0, -t240 * t380 + t367, -t141 * t197 - t303 - t86, t223, -t377, t235, -qJD(2) * t34 + (-t6 - t342) * t197 + (-t22 * t380 + t5) * t194 + t276, t223, t235, t377, qJD(2) * t21 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t107 - t359, -t369, -t375, t108, -t111 * t89 + t218 - t341, t109 * t89 + t205, -t109 * t53 - t206 - t341 + t99, pkin(4) * t45 - t344 + (t25 - t32) * t111 + (t22 - t314) * t109, -t109 * t34 + t111 * t53 - t205 + t370, -t6 * pkin(4) - g(1) * t266 - g(2) * t265 + t5 * qJ(5) + t188 * t251 - t22 * t32 + t25 * t314 - t34 * t53, t109 * t35 - t380 * t20 + (pkin(5) + t356) * t108 + t362, t109 * t21 - t111 * t35 + t19 * t380 - t211 + t370 + t372, t344 - t356 * t45 + (-t17 + t20) * t111 + (-t14 + t315) * t109, t3 * qJ(5) - t1 * t356 - t14 * t20 - t21 * t35 - g(1) * (-pkin(5) * t91 + t266) - g(2) * (pkin(5) * t93 + t265) + t315 * t17 + g(3) * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t369, t366, t206 + t342, t213, t366, t369, t17 * t380 - t355 - t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378, -t379, -t107 - t359, -t109 * t17 + t111 * t14 - t363 + t4;];
tau_reg  = t2;
