% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:26:04
% EndTime: 2019-05-05 10:26:20
% DurationCPUTime: 7.40s
% Computational Cost: add. (46120->510), mult. (96788->750), div. (0->0), fcn. (74468->14), ass. (0->309)
t287 = sin(qJ(6));
t289 = sin(qJ(4));
t294 = cos(qJ(4));
t295 = cos(qJ(3));
t290 = sin(qJ(3));
t335 = qJD(2) * t290;
t250 = qJD(2) * t294 * t295 - t289 * t335;
t251 = (t289 * t295 + t290 * t294) * qJD(2);
t288 = sin(qJ(5));
t293 = cos(qJ(5));
t228 = t250 * t288 + t251 * t293;
t332 = qJD(2) * qJD(3);
t320 = t295 * t332;
t331 = t290 * qJDD(2);
t257 = t320 + t331;
t321 = t290 * t332;
t330 = t295 * qJDD(2);
t305 = -t321 + t330;
t212 = t250 * qJD(4) + t294 * t257 + t289 * t305;
t313 = t289 * t257 - t294 * t305;
t306 = -qJD(4) * t251 - t313;
t314 = t288 * t212 - t293 * t306;
t156 = -qJD(5) * t228 - t314;
t155 = qJDD(6) - t156;
t281 = qJD(3) + qJD(4);
t277 = qJD(5) + t281;
t292 = cos(qJ(6));
t207 = t228 * t287 - t292 * t277;
t209 = t228 * t292 + t277 * t287;
t175 = t209 * t207;
t370 = t155 - t175;
t377 = t287 * t370;
t226 = -t293 * t250 + t251 * t288;
t191 = t228 * t226;
t329 = qJDD(3) + qJDD(4);
t276 = qJDD(5) + t329;
t369 = -t191 + t276;
t376 = t288 * t369;
t234 = t250 * t251;
t368 = t234 + t329;
t375 = t289 * t368;
t374 = t292 * t370;
t373 = t293 * t369;
t372 = t294 * t368;
t284 = sin(pkin(6));
t285 = cos(pkin(6));
t359 = sin(pkin(12));
t360 = cos(pkin(12));
t304 = g(1) * t359 - g(2) * t360;
t336 = -g(3) + qJDD(1);
t371 = t284 * t336 + t285 * t304;
t185 = pkin(5) * t226 - pkin(11) * t228;
t367 = t277 ^ 2;
t261 = -g(1) * t360 - g(2) * t359;
t291 = sin(qJ(2));
t296 = cos(qJ(2));
t223 = t296 * t261 + t291 * t371;
t298 = qJD(2) ^ 2;
t215 = -t298 * pkin(2) + qJDD(2) * pkin(8) + t223;
t240 = -t284 * t304 + t285 * t336;
t192 = t290 * t215 - t295 * t240;
t337 = t295 * t298;
t271 = t290 * t337;
t262 = qJDD(3) + t271;
t172 = (-t257 + t320) * pkin(9) + t262 * pkin(3) - t192;
t325 = pkin(9) * t335;
t265 = qJD(3) * pkin(3) - t325;
t340 = t290 * t240;
t173 = t340 + (-t265 - t325) * qJD(3) + (-pkin(3) * t337 + pkin(9) * qJDD(2) + t215) * t295;
t132 = t289 * t172 + t294 * t173;
t239 = pkin(4) * t281 - pkin(10) * t251;
t248 = t250 ^ 2;
t111 = -t248 * pkin(4) + pkin(10) * t306 - t281 * t239 + t132;
t131 = -t294 * t172 + t289 * t173;
t343 = t281 * t250;
t198 = -t212 + t343;
t300 = pkin(4) * t368 + t198 * pkin(10) - t131;
t72 = t293 * t111 + t288 * t300;
t51 = -pkin(5) * t367 + t276 * pkin(11) - t226 * t185 + t72;
t310 = t291 * t261 - t296 * t371;
t214 = -qJDD(2) * pkin(2) - t298 * pkin(8) + t310;
t366 = t295 ^ 2;
t279 = t366 * t298;
t188 = -t305 * pkin(3) - pkin(9) * t279 + t265 * t335 + t214;
t147 = -t306 * pkin(4) - t248 * pkin(10) + t251 * t239 + t188;
t157 = -t226 * qJD(5) + t293 * t212 + t288 * t306;
t346 = t277 * t226;
t308 = t157 - t346;
t80 = (t228 * t277 - t156) * pkin(5) + t147 - t308 * pkin(11);
t39 = t287 * t51 - t292 * t80;
t40 = t287 * t80 + t292 * t51;
t22 = t287 * t39 + t292 * t40;
t221 = qJD(6) + t226;
t315 = t287 * t157 - t292 * t276;
t118 = (qJD(6) - t221) * t209 + t315;
t205 = t207 ^ 2;
t206 = t209 ^ 2;
t220 = t221 ^ 2;
t224 = t226 ^ 2;
t225 = t228 ^ 2;
t249 = t251 ^ 2;
t280 = t281 ^ 2;
t365 = pkin(5) * t288;
t71 = t111 * t288 - t293 * t300;
t50 = -t276 * pkin(5) - pkin(11) * t367 + t185 * t228 + t71;
t364 = -pkin(5) * t50 + pkin(11) * t22;
t47 = t287 * t50;
t42 = t288 * t72 - t293 * t71;
t363 = t289 * t42;
t83 = -t131 * t294 + t132 * t289;
t362 = t290 * t83;
t48 = t292 * t50;
t361 = t294 * t42;
t128 = t155 + t175;
t358 = t128 * t287;
t357 = t128 * t292;
t356 = t147 * t288;
t355 = t147 * t293;
t183 = t191 + t276;
t354 = t183 * t288;
t353 = t183 * t293;
t352 = t188 * t289;
t351 = t188 * t294;
t350 = t221 * t287;
t349 = t221 * t292;
t231 = -t234 + t329;
t348 = t231 * t289;
t347 = t231 * t294;
t345 = t277 * t288;
t344 = t277 * t293;
t342 = t281 * t289;
t341 = t281 * t294;
t339 = t290 * t262;
t263 = qJDD(3) - t271;
t338 = t295 * t263;
t333 = qJD(6) + t221;
t14 = t22 * t288 - t293 * t50;
t328 = pkin(4) * t14 + t364;
t309 = -t292 * t157 - t287 * t276;
t123 = t207 * t333 + t309;
t169 = -t206 - t220;
t94 = -t169 * t287 - t357;
t327 = pkin(5) * t123 + pkin(11) * t94 + t47;
t120 = -t209 * t333 - t315;
t160 = -t220 - t205;
t91 = t160 * t292 - t377;
t326 = pkin(5) * t120 + pkin(11) * t91 - t48;
t324 = t288 * t175;
t323 = t293 * t175;
t322 = -pkin(5) * t293 - pkin(4);
t43 = t288 * t71 + t293 * t72;
t154 = t205 + t206;
t134 = -qJD(6) * t207 - t309;
t180 = t221 * t207;
t122 = t134 + t180;
t76 = -t118 * t292 + t122 * t287;
t318 = pkin(5) * t154 + pkin(11) * t76 + t22;
t64 = t123 * t293 + t288 * t94;
t317 = pkin(4) * t64 + t327;
t61 = t120 * t293 + t288 * t91;
t316 = pkin(4) * t61 + t326;
t84 = t131 * t289 + t294 * t132;
t193 = t295 * t215 + t340;
t152 = t192 * t290 + t295 * t193;
t55 = t154 * t293 + t288 * t76;
t312 = pkin(4) * t55 + t318;
t181 = -t367 - t224;
t149 = t181 * t288 + t373;
t311 = pkin(4) * t149 - t71;
t21 = t287 * t40 - t292 * t39;
t307 = t212 + t343;
t258 = -0.2e1 * t321 + t330;
t303 = (-qJD(5) + t277) * t228 - t314;
t302 = (-qJD(4) + t281) * t251 - t313;
t210 = -t225 - t367;
t163 = t210 * t293 - t354;
t299 = pkin(4) * t163 - t72;
t297 = qJD(3) ^ 2;
t282 = t290 ^ 2;
t278 = t282 * t298;
t269 = -t279 - t297;
t268 = -t278 - t297;
t260 = t278 + t279;
t259 = (t282 + t366) * qJDD(2);
t256 = 0.2e1 * t320 + t331;
t242 = -t249 + t280;
t241 = t248 - t280;
t237 = -t249 - t280;
t236 = -t268 * t290 - t338;
t235 = t269 * t295 - t339;
t233 = t249 - t248;
t229 = -t280 - t248;
t217 = -t225 + t367;
t216 = t224 - t367;
t213 = -t248 - t249;
t201 = -t237 * t289 - t347;
t200 = t237 * t294 - t348;
t194 = (qJD(4) + t281) * t251 + t313;
t189 = -t224 + t225;
t187 = t229 * t294 - t375;
t186 = t229 * t289 + t372;
t179 = -t206 + t220;
t178 = t205 - t220;
t177 = (-t226 * t293 + t228 * t288) * t277;
t176 = (-t226 * t288 - t228 * t293) * t277;
t174 = t206 - t205;
t170 = -t224 - t225;
t168 = t216 * t293 - t354;
t167 = -t217 * t288 + t373;
t166 = t216 * t288 + t353;
t165 = t217 * t293 + t376;
t164 = -t210 * t288 - t353;
t161 = -t200 * t290 + t201 * t295;
t159 = -t198 * t289 + t294 * t302;
t158 = t198 * t294 + t289 * t302;
t151 = -t186 * t290 + t187 * t295;
t150 = t181 * t293 - t376;
t146 = (-t207 * t292 + t209 * t287) * t221;
t145 = (-t207 * t287 - t209 * t292) * t221;
t143 = -t157 - t346;
t139 = (qJD(5) + t277) * t228 + t314;
t138 = t157 * t293 - t228 * t345;
t137 = t157 * t288 + t228 * t344;
t136 = -t156 * t288 + t226 * t344;
t135 = t156 * t293 + t226 * t345;
t133 = -qJD(6) * t209 - t315;
t126 = -t163 * t289 + t164 * t294;
t125 = t163 * t294 + t164 * t289;
t124 = -t158 * t290 + t159 * t295;
t121 = t134 - t180;
t115 = t134 * t292 - t209 * t350;
t114 = t134 * t287 + t209 * t349;
t113 = -t133 * t287 + t207 * t349;
t112 = t133 * t292 + t207 * t350;
t109 = -pkin(10) * t163 + t355;
t108 = t146 * t293 + t155 * t288;
t107 = t146 * t288 - t155 * t293;
t106 = t178 * t292 - t358;
t105 = -t179 * t287 + t374;
t104 = t178 * t287 + t357;
t103 = t179 * t292 + t377;
t102 = -t149 * t289 + t150 * t294;
t101 = t149 * t294 + t150 * t289;
t100 = -pkin(10) * t149 + t356;
t99 = -t139 * t293 - t288 * t308;
t98 = -t143 * t288 + t293 * t303;
t97 = -t139 * t288 + t293 * t308;
t96 = t143 * t293 + t288 * t303;
t95 = pkin(4) * t96;
t93 = t169 * t292 - t358;
t90 = t160 * t287 + t374;
t88 = t115 * t293 + t324;
t87 = t113 * t293 - t324;
t86 = t115 * t288 - t323;
t85 = t113 * t288 + t323;
t82 = -pkin(4) * t308 + pkin(10) * t164 + t356;
t81 = -pkin(4) * t139 + pkin(10) * t150 - t355;
t78 = -t125 * t290 + t126 * t295;
t77 = t120 * t292 - t121 * t287;
t75 = t120 * t287 + t121 * t292;
t74 = -t118 * t287 - t122 * t292;
t69 = t106 * t293 - t118 * t288;
t68 = t105 * t293 + t122 * t288;
t67 = t106 * t288 + t118 * t293;
t66 = t105 * t288 - t122 * t293;
t65 = -t123 * t288 + t293 * t94;
t62 = -t120 * t288 + t293 * t91;
t59 = -t101 * t290 + t102 * t295;
t58 = t174 * t288 + t293 * t77;
t57 = -t174 * t293 + t288 * t77;
t56 = -t154 * t288 + t293 * t76;
t53 = -t289 * t96 + t294 * t98;
t52 = t289 * t98 + t294 * t96;
t46 = t295 * t84 - t362;
t45 = -pkin(11) * t93 + t48;
t44 = -pkin(11) * t90 + t47;
t41 = pkin(4) * t42;
t36 = -t289 * t64 + t294 * t65;
t35 = t289 * t65 + t294 * t64;
t34 = -t289 * t61 + t294 * t62;
t33 = t289 * t62 + t294 * t61;
t32 = -pkin(4) * t147 + pkin(10) * t43;
t31 = -t289 * t55 + t294 * t56;
t30 = t289 * t56 + t294 * t55;
t29 = -t290 * t52 + t295 * t53;
t28 = -pkin(10) * t96 - t42;
t27 = -pkin(5) * t93 + t40;
t26 = -pkin(5) * t90 + t39;
t25 = -pkin(4) * t170 + pkin(10) * t98 + t43;
t24 = t294 * t43 - t363;
t23 = t289 * t43 + t361;
t19 = -t290 * t35 + t295 * t36;
t18 = -t290 * t33 + t295 * t34;
t17 = -t290 * t30 + t295 * t31;
t16 = -pkin(11) * t74 - t21;
t15 = t22 * t293 + t288 * t50;
t12 = -pkin(10) * t64 - t27 * t288 + t293 * t45;
t11 = -pkin(10) * t61 - t26 * t288 + t293 * t44;
t10 = -pkin(4) * t93 + pkin(10) * t65 + t27 * t293 + t288 * t45;
t9 = -pkin(4) * t90 + pkin(10) * t62 + t26 * t293 + t288 * t44;
t8 = -pkin(10) * t55 + t16 * t293 + t365 * t74;
t7 = -t23 * t290 + t24 * t295;
t6 = pkin(10) * t56 + t288 * t16 + t322 * t74;
t5 = -t14 * t289 + t15 * t294;
t4 = t14 * t294 + t15 * t289;
t3 = -pkin(10) * t14 + (-pkin(11) * t293 + t365) * t21;
t2 = pkin(10) * t15 + (-pkin(11) * t288 + t322) * t21;
t1 = -t290 * t4 + t295 * t5;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t336, 0, 0, 0, 0, 0, 0, (qJDD(2) * t296 - t291 * t298) * t284, (-qJDD(2) * t291 - t296 * t298) * t284, 0, t285 * t240 + (t223 * t291 - t296 * t310) * t284, 0, 0, 0, 0, 0, 0, t285 * (t262 * t295 + t269 * t290) + (t235 * t291 + t258 * t296) * t284, t285 * (-t263 * t290 + t268 * t295) + (t236 * t291 - t256 * t296) * t284, (t259 * t291 + t260 * t296) * t284, t285 * (-t192 * t295 + t193 * t290) + (t152 * t291 - t214 * t296) * t284, 0, 0, 0, 0, 0, 0, t285 * (t186 * t295 + t187 * t290) + (t151 * t291 - t194 * t296) * t284, t285 * (t200 * t295 + t201 * t290) + (t161 * t291 - t296 * t307) * t284, t285 * (t158 * t295 + t159 * t290) + (t124 * t291 - t213 * t296) * t284, t285 * (t290 * t84 + t295 * t83) + (-t188 * t296 + t291 * t46) * t284, 0, 0, 0, 0, 0, 0, t285 * (t101 * t295 + t102 * t290) + (-t139 * t296 + t291 * t59) * t284, t285 * (t125 * t295 + t126 * t290) + (t291 * t78 - t296 * t308) * t284, t285 * (t290 * t53 + t295 * t52) + (-t170 * t296 + t29 * t291) * t284, t285 * (t23 * t295 + t24 * t290) + (-t147 * t296 + t291 * t7) * t284, 0, 0, 0, 0, 0, 0, t285 * (t290 * t34 + t295 * t33) + (t18 * t291 - t296 * t90) * t284, t285 * (t290 * t36 + t295 * t35) + (t19 * t291 - t296 * t93) * t284, t285 * (t290 * t31 + t295 * t30) + (t17 * t291 - t296 * t74) * t284, t285 * (t290 * t5 + t295 * t4) + (t1 * t291 - t21 * t296) * t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t310, -t223, 0, 0, (t257 + t320) * t290, t256 * t295 + t258 * t290, t339 + t295 * (-t278 + t297), t258 * t295, t290 * (t279 - t297) + t338, 0, pkin(2) * t258 + pkin(8) * t235 - t214 * t295, -pkin(2) * t256 + pkin(8) * t236 + t214 * t290, pkin(2) * t260 + pkin(8) * t259 + t152, -pkin(2) * t214 + pkin(8) * t152, t290 * (t212 * t294 - t251 * t342) + t295 * (t212 * t289 + t251 * t341), t290 * (-t194 * t294 - t289 * t307) + t295 * (-t194 * t289 + t294 * t307), t290 * (-t242 * t289 + t372) + t295 * (t242 * t294 + t375), t290 * (-t250 * t341 - t289 * t306) + t295 * (-t250 * t342 + t294 * t306), t290 * (t241 * t294 - t348) + t295 * (t241 * t289 + t347), (t290 * (t250 * t294 + t251 * t289) + t295 * (t250 * t289 - t251 * t294)) * t281, t290 * (-pkin(9) * t186 + t352) + t295 * (-pkin(3) * t194 + pkin(9) * t187 - t351) - pkin(2) * t194 + pkin(8) * t151, t290 * (-pkin(9) * t200 + t351) + t295 * (-pkin(3) * t307 + pkin(9) * t201 + t352) - pkin(2) * t307 + pkin(8) * t161, t290 * (-pkin(9) * t158 - t83) + t295 * (-pkin(3) * t213 + pkin(9) * t159 + t84) - pkin(2) * t213 + pkin(8) * t124, -pkin(9) * t362 + t295 * (-pkin(3) * t188 + pkin(9) * t84) - pkin(2) * t188 + pkin(8) * t46, t290 * (-t137 * t289 + t138 * t294) + t295 * (t137 * t294 + t138 * t289), t290 * (-t289 * t97 + t294 * t99) + t295 * (t289 * t99 + t294 * t97), t290 * (-t165 * t289 + t167 * t294) + t295 * (t165 * t294 + t167 * t289), t290 * (-t135 * t289 + t136 * t294) + t295 * (t135 * t294 + t136 * t289), t290 * (-t166 * t289 + t168 * t294) + t295 * (t166 * t294 + t168 * t289), t290 * (-t176 * t289 + t177 * t294) + t295 * (t176 * t294 + t177 * t289), t290 * (-pkin(9) * t101 + t100 * t294 - t289 * t81) + t295 * (-pkin(3) * t139 + pkin(9) * t102 + t100 * t289 + t294 * t81) - pkin(2) * t139 + pkin(8) * t59, t290 * (-pkin(9) * t125 + t109 * t294 - t289 * t82) + t295 * (-pkin(3) * t308 + pkin(9) * t126 + t109 * t289 + t294 * t82) - pkin(2) * t308 + pkin(8) * t78, t290 * (-pkin(9) * t52 - t25 * t289 + t28 * t294) + t295 * (-pkin(3) * t170 + pkin(9) * t53 + t25 * t294 + t28 * t289) - pkin(2) * t170 + pkin(8) * t29, t290 * (-pkin(9) * t23 - pkin(10) * t361 - t289 * t32) + t295 * (-pkin(3) * t147 + pkin(9) * t24 - pkin(10) * t363 + t294 * t32) - pkin(2) * t147 + pkin(8) * t7, t290 * (-t289 * t86 + t294 * t88) + t295 * (t289 * t88 + t294 * t86), t290 * (-t289 * t57 + t294 * t58) + t295 * (t289 * t58 + t294 * t57), t290 * (-t289 * t66 + t294 * t68) + t295 * (t289 * t68 + t294 * t66), t290 * (-t289 * t85 + t294 * t87) + t295 * (t289 * t87 + t294 * t85), t290 * (-t289 * t67 + t294 * t69) + t295 * (t289 * t69 + t294 * t67), t290 * (-t107 * t289 + t108 * t294) + t295 * (t107 * t294 + t108 * t289), t290 * (-pkin(9) * t33 + t11 * t294 - t289 * t9) + t295 * (-pkin(3) * t90 + pkin(9) * t34 + t11 * t289 + t294 * t9) - pkin(2) * t90 + pkin(8) * t18, t290 * (-pkin(9) * t35 - t10 * t289 + t12 * t294) + t295 * (-pkin(3) * t93 + pkin(9) * t36 + t10 * t294 + t12 * t289) - pkin(2) * t93 + pkin(8) * t19, t290 * (-pkin(9) * t30 - t289 * t6 + t294 * t8) + t295 * (-pkin(3) * t74 + pkin(9) * t31 + t289 * t8 + t294 * t6) - pkin(2) * t74 + pkin(8) * t17, t290 * (-pkin(9) * t4 - t2 * t289 + t294 * t3) + t295 * (-pkin(3) * t21 + pkin(9) * t5 + t2 * t294 + t289 * t3) - pkin(2) * t21 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, -t279 + t278, t331, t271, t330, qJDD(3), -t192, -t193, 0, 0, -t234, t233, -t198, t234, t302, t329, pkin(3) * t186 - t131, pkin(3) * t200 - t132, pkin(3) * t158, pkin(3) * t83, t191, t189, -t143, -t191, t303, t276, pkin(3) * t101 + t311, pkin(3) * t125 + t299, pkin(3) * t52 + t95, pkin(3) * t23 + t41, t114, t75, t103, t112, t104, t145, pkin(3) * t33 + t316, pkin(3) * t35 + t317, pkin(3) * t30 + t312, pkin(3) * t4 + t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t233, -t198, t234, t302, t329, -t131, -t132, 0, 0, t191, t189, -t143, -t191, t303, t276, t311, t299, t95, t41, t114, t75, t103, t112, t104, t145, t316, t317, t312, t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t189, -t143, -t191, t303, t276, -t71, -t72, 0, 0, t114, t75, t103, t112, t104, t145, t326, t327, t318, t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t174, t122, -t175, -t118, t155, -t39, -t40, 0, 0;];
tauJ_reg  = t13;
