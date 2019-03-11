% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:48
% EndTime: 2019-03-09 08:24:00
% DurationCPUTime: 5.36s
% Computational Cost: add. (3171->567), mult. (7162->682), div. (0->0), fcn. (4859->8), ass. (0->282)
t196 = cos(qJ(2));
t304 = qJD(1) * t196;
t151 = -qJD(6) + t304;
t370 = qJD(6) + t151;
t369 = pkin(4) + qJ(3);
t193 = sin(qJ(2));
t194 = sin(qJ(1));
t197 = cos(qJ(1));
t250 = g(1) * t197 + g(2) * t194;
t353 = t193 * t250;
t191 = cos(pkin(9));
t305 = qJD(1) * t193;
t275 = t191 * t305;
t190 = sin(pkin(9));
t303 = qJD(2) * t190;
t113 = t275 + t303;
t297 = t113 * qJD(4);
t170 = t191 * qJDD(2);
t293 = qJD(1) * qJD(2);
t265 = t196 * t293;
t290 = t193 * qJDD(1);
t218 = t265 + t290;
t76 = t190 * t218 - t170;
t349 = t76 * pkin(3);
t164 = pkin(7) * t290;
t89 = -qJDD(2) * pkin(2) + pkin(7) * t265 + qJDD(3) + t164;
t204 = -t297 + t89 + t349;
t341 = -pkin(5) - qJ(4);
t277 = t190 * t305;
t296 = t191 * qJD(2);
t111 = t277 - t296;
t364 = t76 * qJ(5) + t111 * qJD(5);
t291 = t190 * qJDD(2);
t77 = t191 * t218 + t291;
t4 = t341 * t77 + t204 + t364;
t368 = -t4 + t353;
t192 = sin(qJ(6));
t195 = cos(qJ(6));
t52 = t111 * t192 - t195 * t113;
t320 = t193 * t197;
t321 = t193 * t194;
t367 = g(1) * t320 + g(2) * t321;
t174 = t196 * qJDD(1);
t266 = t193 * t293;
t357 = -t266 + t174;
t115 = qJDD(6) - t357;
t233 = t151 ^ 2;
t366 = t195 * t115 - t192 * t233;
t365 = -2 * pkin(1);
t187 = g(3) * t196;
t269 = t89 + t187;
t342 = pkin(3) + qJ(5);
t363 = t111 * t342;
t362 = t191 * t342;
t292 = qJD(1) * qJD(4);
t294 = qJ(4) * qJDD(1);
t361 = t196 * (t292 + t294);
t298 = qJD(4) * t190;
t167 = pkin(7) * t304;
t276 = t190 * t304;
t311 = pkin(3) * t276 + t167;
t360 = -qJD(5) * t191 - t298 - t311;
t179 = t196 * qJ(5);
t324 = t191 * t193;
t359 = pkin(4) * t324 + t179;
t175 = t196 * qJD(5);
t358 = qJ(5) * t174 + qJD(1) * t175;
t302 = qJD(2) * t193;
t356 = -qJ(4) * t302 + t196 * qJD(4);
t189 = t196 ^ 2;
t199 = qJD(1) ^ 2;
t329 = t189 * t199;
t351 = t113 ^ 2;
t355 = -t351 - t329;
t107 = t111 ^ 2;
t354 = t351 + t107;
t54 = t111 * t195 + t113 * t192;
t11 = qJD(6) * t54 + t192 * t76 - t195 * t77;
t301 = qJD(2) * t196;
t273 = t190 * t301;
t310 = pkin(3) * t273 + pkin(7) * t301;
t352 = -(qJD(4) * t191 - qJD(5) * t190) * t193 + t310;
t217 = t191 * t290 + t291;
t40 = (t111 + t296) * t304 + t217;
t350 = pkin(4) + pkin(8);
t348 = t77 * pkin(4);
t347 = pkin(1) * t199;
t346 = pkin(7) * t113;
t345 = g(1) * t194;
t183 = t196 * pkin(2);
t177 = t193 * qJ(3);
t282 = -pkin(1) - t183;
t229 = t282 - t177;
t245 = pkin(2) * t193 - qJ(3) * t196;
t93 = qJD(2) * t245 - t193 * qJD(3);
t49 = qJD(1) * t93 + qJDD(1) * t229;
t84 = t357 * pkin(7) + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t19 = t190 * t49 + t191 * t84;
t116 = t190 * t195 + t191 * t192;
t220 = t196 * t116;
t95 = t116 * qJD(6);
t340 = -qJD(1) * t220 + t95;
t274 = t191 * t304;
t323 = t191 * t195;
t117 = t190 * t192 - t323;
t96 = t117 * qJD(6);
t339 = -t192 * t276 + t195 * t274 + t96;
t338 = qJ(4) * t77;
t337 = t151 * t52;
t336 = t191 * t93;
t335 = t54 * t151;
t104 = t229 * qJD(1);
t133 = qJD(2) * qJ(3) + t167;
t57 = t190 * t104 + t191 * t133;
t332 = qJ(5) * t190;
t222 = t191 * t341 + t332;
t208 = t196 * t222;
t334 = -qJD(1) * t208 + t360;
t234 = qJ(4) * t191 - t332;
t219 = t234 * t196;
t333 = -qJD(1) * t219 - t360;
t331 = t111 * t113;
t328 = t190 * qJ(4);
t327 = t190 * t193;
t326 = t190 * t196;
t120 = t245 * qJD(1);
t325 = t191 * t120;
t322 = t191 * t196;
t319 = t194 * t196;
t317 = t196 * qJ(4);
t316 = t196 * t197;
t315 = t197 * t190;
t56 = t191 * t104 - t190 * t133;
t42 = pkin(3) * t304 + qJD(4) - t56;
t216 = qJ(5) * t304 + t42;
t22 = pkin(4) * t113 + t216;
t314 = -qJD(4) - t22;
t102 = t190 * t120;
t161 = qJ(4) * t305;
t313 = t102 + t161;
t312 = t367 * t190;
t309 = pkin(3) * t327 + t193 * pkin(7);
t307 = t183 + t177;
t124 = -pkin(1) - t307;
t79 = pkin(7) * t322 + t190 * t124;
t308 = g(1) * t321 - g(2) * t320;
t127 = t369 * t190;
t128 = t369 * t191;
t188 = t193 ^ 2;
t306 = t188 - t189;
t300 = qJD(3) * t111;
t299 = qJD(3) * t113;
t287 = -qJ(4) * t266 - t19;
t286 = pkin(7) * t302;
t285 = t190 * t350;
t284 = qJD(3) * t276 + t191 * t367;
t283 = g(1) * t316 + g(2) * t319 + g(3) * t193;
t281 = -pkin(7) * t190 - pkin(3);
t13 = t204 - t338;
t7 = -t13 - t364;
t280 = t7 - t187;
t18 = -t190 * t84 + t191 * t49;
t227 = pkin(3) * t174 + qJDD(4) - t18;
t256 = t342 * t302;
t201 = -qJD(1) * t256 + t227 + t358;
t5 = t350 * t77 + t201;
t260 = qJDD(5) - t287;
t6 = pkin(5) * t266 - t350 * t76 + (qJDD(1) * t341 - t292) * t196 + t260;
t279 = -t192 * t5 + t195 * t6;
t278 = qJ(3) * t302;
t272 = t196 * t296;
t271 = qJ(4) * t304;
t270 = t13 + t187;
t268 = t341 * t196;
t267 = qJ(3) * t174;
t264 = pkin(2) + t328;
t262 = pkin(4) * t274 - t325;
t143 = pkin(7) * t326;
t78 = t124 * t191 - t143;
t259 = pkin(3) * t322 + t190 * t317 + t307;
t258 = t197 * pkin(1) + pkin(2) * t316 + t194 * pkin(7) + qJ(3) * t320;
t184 = t197 * pkin(7);
t98 = t190 * t319 + t191 * t197;
t99 = t191 * t319 - t315;
t257 = -t99 * pkin(3) - t98 * qJ(4) + t184;
t100 = -t194 * t191 + t196 * t315;
t255 = g(1) * t98 - g(2) * t100;
t101 = t190 * t194 + t191 * t316;
t254 = g(1) * t99 - g(2) * t101;
t253 = pkin(4) * t272 + t175 - t336;
t252 = t281 * t193;
t251 = -pkin(2) - t362;
t249 = -g(2) * t197 + t345;
t248 = t192 * t6 + t195 * t5;
t247 = qJD(2) * pkin(2) - pkin(7) * t305 - qJD(3);
t65 = -t79 + t317;
t82 = t190 * t93;
t246 = t82 - t356;
t182 = t196 * pkin(3);
t69 = t182 - t78;
t244 = pkin(7) * t111 - t190 * t247;
t242 = t190 * t290 - t170;
t15 = t113 * t350 + t216;
t16 = qJD(1) * t268 - t111 * t350 + qJD(5) + t57;
t2 = t15 * t195 + t16 * t192;
t241 = t15 * t192 - t16 * t195;
t32 = t143 + t182 + (pkin(8) * t193 - t124) * t191 + t359;
t33 = -t193 * t285 + t268 + t79;
t240 = t192 * t33 + t195 * t32;
t239 = -t192 * t32 + t195 * t33;
t238 = t192 * t99 + t195 * t98;
t237 = t192 * t98 - t195 * t99;
t236 = -qJ(3) * t76 - t300;
t235 = qJ(3) * t77 + t299;
t63 = -t191 * t286 + t82;
t198 = qJD(2) ^ 2;
t231 = qJDD(2) * t196 - t193 * t198;
t71 = -pkin(7) * t275 + t102;
t230 = (-qJ(5) + t281) * t193;
t228 = t191 * pkin(3) + t264;
t47 = t271 - t57;
t225 = -t76 * pkin(4) + t260;
t105 = pkin(8) * t190 + t127;
t205 = (-pkin(7) * t191 + pkin(5)) * t193 - t196 * t285;
t224 = -qJD(1) * t205 + qJD(3) * t191 - qJD(6) * t105 - t313;
t106 = pkin(8) * t191 + t128;
t207 = pkin(8) * t322 + t230;
t223 = -qJD(1) * t207 + qJD(3) * t190 + qJD(6) * t106 - t262;
t10 = -qJD(6) * t52 + t192 * t77 + t195 * t76;
t221 = -pkin(4) * t326 - pkin(7) * t324;
t215 = t191 * t267 - t312;
t213 = t113 * qJ(4) + t247;
t212 = -t357 - t331;
t211 = t101 * pkin(3) + qJ(4) * t100 + t258;
t210 = t229 * t345;
t209 = t190 * t267 + t284;
t206 = -t192 * t115 - t195 * t233;
t38 = (-t113 + t303) * t304 + t242;
t203 = -g(1) * t100 - g(2) * t98 - g(3) * t327 + t227;
t202 = -t353 + t269;
t200 = t202 - t338 + t349;
t149 = qJ(3) * t316;
t144 = qJ(3) * t319;
t97 = t264 + t362;
t88 = t116 * t193;
t87 = t192 * t327 - t193 * t323;
t85 = -qJ(4) * t324 + t309;
t83 = t190 * t341 + t251;
t75 = -t191 * t271 + t311;
t70 = pkin(7) * t277 + t325;
t64 = t193 * t234 - t309;
t62 = t190 * t286 + t336;
t61 = qJD(1) * t252 - t325;
t60 = -t161 - t71;
t59 = (-qJ(4) * t301 - qJD(4) * t193) * t191 + t310;
t55 = -pkin(4) * t327 - t65;
t51 = t193 * t222 + t309;
t50 = qJD(2) * t252 - t336;
t48 = t69 + t359;
t46 = t100 * t195 + t101 * t192;
t45 = -t100 * t192 + t101 * t195;
t44 = qJD(1) * t221 + t313;
t41 = -t63 + t356;
t39 = (-t111 + t296) * t304 + t217;
t37 = t111 * pkin(3) - t213;
t35 = qJD(2) * t220 - t193 * t96;
t34 = t192 * t273 + t193 * t95 - t195 * t272;
t31 = qJD(1) * t230 + t262;
t30 = qJD(2) * t219 - t352;
t29 = qJD(2) * t221 + t246;
t26 = qJD(2) * t230 + t253;
t25 = -t111 * pkin(4) + qJD(5) - t47;
t24 = t213 - t363;
t23 = qJD(2) * t208 + t352;
t21 = qJD(2) * t205 + t246;
t20 = qJD(2) * t207 + t253;
t17 = t113 * t341 - t247 + t363;
t14 = -pkin(3) * t266 + t227;
t12 = t287 + t361;
t9 = t225 - t361;
t8 = t201 + t348;
t1 = [qJDD(1), t249, t250, qJDD(1) * t188 + 0.2e1 * t193 * t265, 0.2e1 * t174 * t193 - 0.2e1 * t293 * t306, qJDD(2) * t193 + t196 * t198, t231, 0 (-pkin(7) * qJDD(2) + t293 * t365) * t193 + (0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t198 + t249) * t196, -pkin(7) * t231 + t218 * t365 - t308 (pkin(7) * t76 + t190 * t89 + (qJD(1) * t78 + t56) * qJD(2)) * t193 + (-qJD(1) * t62 + qJD(2) * t244 - qJDD(1) * t78 - t18) * t196 + t254 (pkin(7) * t77 + t89 * t191 + (-qJD(1) * t79 - t57) * qJD(2)) * t193 + (qJD(1) * t63 + qJDD(1) * t79 + t19 + (-t191 * t247 + t346) * qJD(2)) * t196 - t255, -t63 * t111 - t62 * t113 - t79 * t76 - t78 * t77 + (-t18 * t191 - t19 * t190) * t193 + (-t190 * t57 - t191 * t56) * t301 + t308, t19 * t79 + t57 * t63 + t18 * t78 + t56 * t62 - g(1) * t184 - g(2) * t258 - t210 + (t193 * t89 - t247 * t301) * pkin(7), t41 * t111 + t50 * t113 + t65 * t76 + t69 * t77 + (t12 * t190 + t14 * t191) * t193 + (t190 * t47 + t191 * t42) * t301 + t308, -t111 * t59 - t76 * t85 + (-t13 * t190 + (qJD(1) * t69 + t42) * qJD(2)) * t193 + (-qJD(1) * t50 - qJDD(1) * t69 - t303 * t37 - t14) * t196 - t254, -t113 * t59 - t77 * t85 + (-t13 * t191 + (-qJD(1) * t65 - t47) * qJD(2)) * t193 + (qJD(1) * t41 + qJDD(1) * t65 - t296 * t37 + t12) * t196 + t255, -g(1) * t257 - g(2) * t211 + t12 * t65 + t13 * t85 + t14 * t69 + t37 * t59 + t47 * t41 + t42 * t50 - t210, t113 * t30 + t64 * t77 + (t191 * t7 + (qJD(1) * t55 + t25) * qJD(2)) * t193 + (-qJD(1) * t29 - qJDD(1) * t55 + t24 * t296 - t9) * t196 + t255, t29 * t111 - t26 * t113 - t48 * t77 + t55 * t76 + (t190 * t9 - t191 * t8) * t193 + (t190 * t25 - t191 * t22) * t301 - t308, -t111 * t30 - t64 * t76 + (-t190 * t7 + (-qJD(1) * t48 - t22) * qJD(2)) * t193 + (qJD(1) * t26 + qJDD(1) * t48 - t24 * t303 + t8) * t196 + t254, t8 * t48 + t22 * t26 + t7 * t64 + t24 * t30 + t9 * t55 + t25 * t29 - g(1) * (-t99 * qJ(5) + t257) - g(2) * (pkin(4) * t320 + qJ(5) * t101 + t211) - (-t193 * t369 + t282) * t345, t10 * t88 + t35 * t54, -t10 * t87 - t11 * t88 - t34 * t54 - t35 * t52, -t10 * t196 + t115 * t88 - t151 * t35 + t302 * t54, t11 * t196 - t115 * t87 + t151 * t34 - t302 * t52, -t115 * t196 - t151 * t302 -(-t192 * t20 + t195 * t21) * t151 + t239 * t115 - t279 * t196 - t241 * t302 + t23 * t52 + t51 * t11 + t4 * t87 + t17 * t34 + g(1) * t238 - g(2) * t46 + (t151 * t240 + t196 * t2) * qJD(6) (t192 * t21 + t195 * t20) * t151 - t240 * t115 + t248 * t196 - t2 * t302 + t23 * t54 + t51 * t10 + t4 * t88 + t17 * t35 - g(1) * t237 - g(2) * t45 + (t151 * t239 - t196 * t241) * qJD(6); 0, 0, 0, -t193 * t199 * t196, t306 * t199, t290, t174, qJDD(2), -t164 - t187 + (t250 + t347) * t193 (-pkin(7) * qJDD(1) + t347) * t196 + t283, -pkin(2) * t76 - t269 * t191 + ((-qJ(3) * t303 - t56) * t193 + (-t244 + t70) * t196) * qJD(1) + t209, -pkin(2) * t77 + t269 * t190 + ((-qJ(3) * t296 + t57) * t193 + (-t346 - t71 + (qJD(3) + t247) * t191) * t196) * qJD(1) + t215, t111 * t71 + t113 * t70 + (t304 * t56 + t19 + t236) * t191 + (t304 * t57 - t18 + t235) * t190 - t283, -t89 * pkin(2) - t57 * t71 - t56 * t70 + t247 * t167 - g(1) * (-pkin(2) * t320 + t149) - g(2) * (-pkin(2) * t321 + t144) - g(3) * t307 + (-t190 * t56 + t191 * t57) * qJD(3) + (-t18 * t190 + t19 * t191) * qJ(3), -t111 * t60 - t113 * t61 + (-t304 * t42 - t12 + t236) * t191 + (-t304 * t47 + t14 + t235) * t190 - t283, t228 * t76 + t270 * t191 + (t75 + t298) * t111 + (-t193 * t42 + t196 * t61 + (t196 * t37 + t278) * t190) * qJD(1) - t209, t113 * t75 + t228 * t77 + (-t270 + t297) * t190 + (t193 * t47 - t196 * t60 + (t278 + (-qJD(3) + t37) * t196) * t191) * qJD(1) - t215, -t37 * t75 - t47 * t60 - t42 * t61 - g(1) * t149 - g(2) * t144 - g(3) * t259 + (-qJ(3) * t12 - qJD(3) * t47) * t191 + (qJ(3) * t14 + qJD(3) * t42 - qJD(4) * t37) * t190 + (-t13 + t353) * t228, -t128 * t174 + t77 * t97 + t280 * t190 + t333 * t113 + ((qJD(2) * t128 - t25) * t193 + (t44 + (-qJD(3) - t24) * t191) * t196) * qJD(1) + t312, -t111 * t44 + t113 * t31 - t127 * t77 + t128 * t76 + (t22 * t304 + t300 - t9) * t191 + (-t25 * t304 - t299 - t8) * t190 + t283, t127 * t174 - t76 * t97 + t280 * t191 - t333 * t111 + ((t190 * t24 - t31) * t196 + (-qJD(2) * t127 + t22) * t193) * qJD(1) + t284, t8 * t127 + t7 * t97 + t9 * t128 - t22 * t31 - t25 * t44 - g(1) * (pkin(4) * t316 + t149) - g(2) * (pkin(4) * t319 + t144) - g(3) * (t179 * t191 + t259) + t333 * t24 + (t190 * t22 + t191 * t25) * qJD(3) + (-g(3) * pkin(4) + t250 * (-t251 + t328)) * t193, t10 * t117 + t340 * t54, t10 * t116 - t117 * t11 - t339 * t54 - t340 * t52, t117 * t115 - t151 * t340 - t305 * t54, t116 * t115 + t151 * t339 + t305 * t52, t151 * t305 (-t105 * t192 + t106 * t195) * t115 + t83 * t11 + t334 * t52 - g(3) * t220 + t339 * t17 + (t192 * t223 - t195 * t224) * t151 + t241 * t305 + t368 * t116 -(t105 * t195 + t106 * t192) * t115 + t83 * t10 + t334 * t54 + t340 * t17 + (t192 * t224 + t195 * t223) * t151 + t2 * t305 + (t187 - t368) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, -t354, t57 * t111 + t56 * t113 + t202, -t354, -t38, -t40, -t111 * t47 + (-qJD(4) - t42) * t113 + t200, -t40, t354, t38, t111 * t25 + t113 * t314 + t200 + t364, 0, 0, 0, 0, 0, t11 - t335, t10 + t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t212, t355, t37 * t113 + (-pkin(3) * t302 - t196 * t47) * qJD(1) + t203, t355, -t39, -t212, t348 - t113 * t24 + (t196 * t25 - t256) * qJD(1) + t203 + t358, 0, 0, 0, 0, 0, t113 * t52 + t206, t113 * t54 - t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t357 + t331 (t113 + t303) * t304 + t242, -t107 - t329, -g(3) * t324 - g(1) * t101 - g(2) * t99 + t111 * t24 + (qJD(1) * t314 - t294) * t196 + t225, 0, 0, 0, 0, 0, -t111 * t52 + t366, -t111 * t54 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t52 ^ 2 + t54 ^ 2, t10 - t337, -t11 - t335, t115, -g(1) * t45 + g(2) * t237 + g(3) * t87 - t17 * t54 - t370 * t2 + t279, g(1) * t46 + g(2) * t238 + g(3) * t88 + t17 * t52 + t370 * t241 - t248;];
tau_reg  = t1;
