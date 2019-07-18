% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:01:26
% EndTime: 2019-05-05 07:01:38
% DurationCPUTime: 7.09s
% Computational Cost: add. (41987->498), mult. (90657->749), div. (0->0), fcn. (69082->14), ass. (0->303)
t378 = 2 * qJD(5);
t289 = sin(pkin(12));
t295 = sin(qJ(4));
t299 = cos(qJ(4));
t300 = cos(qJ(3));
t296 = sin(qJ(3));
t335 = qJD(2) * t296;
t256 = qJD(2) * t299 * t300 - t295 * t335;
t257 = (t295 * t300 + t296 * t299) * qJD(2);
t291 = cos(pkin(12));
t230 = -t291 * t256 + t257 * t289;
t232 = t289 * t256 + t291 * t257;
t192 = t232 * t230;
t285 = qJDD(3) + qJDD(4);
t369 = -t192 + t285;
t377 = t289 * t369;
t376 = t291 * t369;
t294 = sin(qJ(6));
t330 = qJD(2) * qJD(3);
t323 = t300 * t330;
t329 = t296 * qJDD(2);
t261 = t323 + t329;
t281 = t300 * qJDD(2);
t324 = t296 * t330;
t311 = t281 - t324;
t315 = t261 * t295 - t299 * t311;
t211 = -qJD(4) * t257 - t315;
t212 = t256 * qJD(4) + t299 * t261 + t295 * t311;
t316 = -t291 * t211 + t212 * t289;
t171 = qJDD(6) + t316;
t286 = qJD(3) + qJD(4);
t298 = cos(qJ(6));
t208 = t232 * t294 - t298 * t286;
t210 = t232 * t298 + t286 * t294;
t177 = t210 * t208;
t370 = t171 - t177;
t375 = t294 * t370;
t238 = t256 * t257;
t368 = t238 + t285;
t374 = t295 * t368;
t373 = t298 * t370;
t372 = t299 * t368;
t290 = sin(pkin(6));
t292 = cos(pkin(6));
t358 = sin(pkin(11));
t359 = cos(pkin(11));
t307 = g(1) * t358 - g(2) * t359;
t336 = -g(3) + qJDD(1);
t371 = t290 * t336 + t292 * t307;
t347 = t232 * t286;
t151 = t316 + t347;
t187 = pkin(5) * t230 - pkin(10) * t232;
t366 = t286 ^ 2;
t265 = -g(1) * t359 - g(2) * t358;
t297 = sin(qJ(2));
t301 = cos(qJ(2));
t229 = t301 * t265 + t297 * t371;
t303 = qJD(2) ^ 2;
t216 = -t303 * pkin(2) + qJDD(2) * pkin(8) + t229;
t244 = -t290 * t307 + t292 * t336;
t194 = t216 * t296 - t300 * t244;
t337 = t300 * t303;
t274 = t296 * t337;
t266 = qJDD(3) + t274;
t172 = (-t261 + t323) * pkin(9) + t266 * pkin(3) - t194;
t269 = qJD(3) * pkin(3) - pkin(9) * t335;
t340 = t296 * t244;
t175 = t340 + t311 * pkin(9) - qJD(3) * t269 + (-pkin(3) * t337 + t216) * t300;
t131 = t172 * t295 + t175 * t299;
t243 = pkin(4) * t286 - qJ(5) * t257;
t254 = t256 ^ 2;
t102 = -pkin(4) * t254 + qJ(5) * t211 - t243 * t286 + t131;
t130 = -t299 * t172 + t175 * t295;
t248 = t286 * t256;
t201 = -t248 + t212;
t304 = t368 * pkin(4) - t201 * qJ(5) - t130;
t55 = -0.2e1 * qJD(5) * t230 + t291 * t102 + t289 * t304;
t51 = -pkin(5) * t366 + pkin(10) * t285 - t187 * t230 + t55;
t312 = t265 * t297 - t301 * t371;
t215 = -qJDD(2) * pkin(2) - t303 * pkin(8) + t312;
t365 = t300 ^ 2;
t283 = t365 * t303;
t191 = -t311 * pkin(3) - pkin(9) * t283 + t269 * t335 + t215;
t137 = -t211 * pkin(4) - t254 * qJ(5) + t243 * t257 + qJDD(5) + t191;
t174 = t211 * t289 + t212 * t291;
t220 = t286 * t230;
t155 = t174 - t220;
t81 = pkin(5) * t151 - pkin(10) * t155 + t137;
t37 = t294 * t51 - t298 * t81;
t38 = t294 * t81 + t298 * t51;
t22 = t294 * t37 + t298 * t38;
t367 = t248 + t212;
t319 = t102 * t289 - t291 * t304;
t54 = t232 * t378 + t319;
t225 = qJD(6) + t230;
t317 = t174 * t294 - t298 * t285;
t123 = (qJD(6) - t225) * t210 + t317;
t205 = t208 ^ 2;
t206 = t210 ^ 2;
t224 = t225 ^ 2;
t226 = t230 ^ 2;
t227 = t232 ^ 2;
t255 = t257 ^ 2;
t364 = pkin(5) * t289;
t50 = -t285 * pkin(5) - t366 * pkin(10) + (t378 + t187) * t232 + t319;
t47 = t294 * t50;
t31 = t289 * t55 - t291 * t54;
t363 = t295 * t31;
t83 = -t130 * t299 + t131 * t295;
t362 = t296 * t83;
t48 = t298 * t50;
t361 = t299 * t31;
t133 = t171 + t177;
t357 = t133 * t294;
t356 = t133 * t298;
t355 = t137 * t289;
t354 = t137 * t291;
t185 = t192 + t285;
t353 = t185 * t289;
t352 = t185 * t291;
t351 = t191 * t295;
t350 = t191 * t299;
t349 = t225 * t294;
t348 = t225 * t298;
t235 = -t238 + t285;
t346 = t235 * t295;
t345 = t235 * t299;
t344 = t286 * t289;
t343 = t286 * t291;
t342 = t286 * t295;
t341 = t286 * t299;
t339 = t296 * t266;
t267 = qJDD(3) - t274;
t338 = t300 * t267;
t331 = qJD(6) + t225;
t13 = t22 * t289 - t291 * t50;
t328 = pkin(4) * t13 - pkin(5) * t50 + pkin(10) * t22;
t327 = t289 * t177;
t326 = t291 * t177;
t325 = -pkin(5) * t291 - pkin(4);
t32 = t289 * t54 + t291 * t55;
t310 = -t174 * t298 - t285 * t294;
t128 = t208 * t331 + t310;
t168 = -t206 - t224;
t97 = -t168 * t294 - t356;
t67 = t128 * t291 + t289 * t97;
t321 = pkin(4) * t67 + pkin(5) * t128 + pkin(10) * t97 + t47;
t125 = -t210 * t331 - t317;
t160 = -t224 - t205;
t92 = t160 * t298 - t375;
t62 = t125 * t291 + t289 * t92;
t320 = pkin(4) * t62 + pkin(5) * t125 + pkin(10) * t92 - t48;
t84 = t130 * t295 + t299 * t131;
t195 = t300 * t216 + t340;
t148 = t194 * t296 + t300 * t195;
t214 = -t227 - t366;
t162 = t214 * t291 - t353;
t314 = pkin(4) * t162 - t55;
t150 = t205 + t206;
t142 = -qJD(6) * t208 - t310;
t182 = t225 * t208;
t127 = t142 + t182;
t77 = -t123 * t298 + t127 * t294;
t57 = t150 * t291 + t289 * t77;
t313 = pkin(4) * t57 + pkin(5) * t150 + pkin(10) * t77 + t22;
t262 = t281 - 0.2e1 * t324;
t21 = t294 * t38 - t298 * t37;
t183 = -t366 - t226;
t139 = t183 * t289 + t376;
t309 = pkin(4) * t139 - t54;
t308 = -t316 + t347;
t306 = (-qJD(4) + t286) * t257 - t315;
t302 = qJD(3) ^ 2;
t287 = t296 ^ 2;
t282 = t287 * t303;
t272 = -t283 - t302;
t271 = -t282 - t302;
t264 = t282 + t283;
t263 = (t287 + t365) * qJDD(2);
t260 = 0.2e1 * t323 + t329;
t246 = -t255 + t366;
t245 = t254 - t366;
t241 = -t255 - t366;
t240 = -t271 * t296 - t338;
t239 = t272 * t300 - t339;
t237 = t255 - t254;
t233 = -t366 - t254;
t218 = -t227 + t366;
t217 = t226 - t366;
t213 = -t254 - t255;
t203 = -t241 * t295 - t345;
t202 = t241 * t299 - t346;
t196 = (qJD(4) + t286) * t257 + t315;
t190 = t233 * t299 - t374;
t189 = t233 * t295 + t372;
t188 = t227 - t226;
t181 = -t206 + t224;
t180 = t205 - t224;
t179 = (-t230 * t291 + t232 * t289) * t286;
t178 = (-t230 * t289 - t232 * t291) * t286;
t176 = t206 - t205;
t169 = -t226 - t227;
t167 = t217 * t291 - t353;
t166 = -t218 * t289 + t376;
t165 = t217 * t289 + t352;
t164 = t218 * t291 + t377;
t163 = -t214 * t289 - t352;
t159 = -t202 * t296 + t203 * t300;
t158 = t201 * t295 + t299 * t306;
t157 = -t201 * t299 + t295 * t306;
t156 = t174 + t220;
t147 = t174 * t291 - t232 * t344;
t146 = t174 * t289 + t232 * t343;
t145 = t230 * t343 + t289 * t316;
t144 = t230 * t344 - t291 * t316;
t143 = -t189 * t296 + t190 * t300;
t141 = -qJD(6) * t210 - t317;
t140 = t183 * t291 - t377;
t136 = (-t208 * t298 + t210 * t294) * t225;
t135 = (-t208 * t294 - t210 * t298) * t225;
t126 = t142 - t182;
t120 = t142 * t298 - t210 * t349;
t119 = t142 * t294 + t210 * t348;
t118 = -t141 * t294 + t208 * t348;
t117 = t141 * t298 + t208 * t349;
t116 = -t162 * t295 + t163 * t299;
t115 = t162 * t299 + t163 * t295;
t114 = -t157 * t296 + t158 * t300;
t113 = t136 * t291 + t171 * t289;
t112 = t136 * t289 - t171 * t291;
t111 = t156 * t289 + t291 * t308;
t110 = -t151 * t291 - t155 * t289;
t109 = -t156 * t291 + t289 * t308;
t108 = -t151 * t289 + t155 * t291;
t107 = pkin(4) * t109;
t106 = t180 * t298 - t357;
t105 = -t181 * t294 + t373;
t104 = t180 * t294 + t356;
t103 = t181 * t298 + t375;
t98 = -qJ(5) * t162 + t354;
t96 = t168 * t298 - t357;
t94 = -t139 * t295 + t140 * t299;
t93 = t139 * t299 + t140 * t295;
t91 = t160 * t294 + t373;
t89 = -qJ(5) * t139 + t355;
t88 = t120 * t291 + t327;
t87 = t118 * t291 - t327;
t86 = t120 * t289 - t326;
t85 = t118 * t289 + t326;
t82 = -pkin(4) * t155 + qJ(5) * t163 + t355;
t79 = -pkin(4) * t151 + qJ(5) * t140 - t354;
t78 = t125 * t298 - t126 * t294;
t76 = t125 * t294 + t126 * t298;
t75 = -t123 * t294 - t127 * t298;
t73 = t106 * t291 - t123 * t289;
t72 = t105 * t291 + t127 * t289;
t71 = t106 * t289 + t123 * t291;
t70 = t105 * t289 - t127 * t291;
t69 = -t115 * t296 + t116 * t300;
t68 = -t128 * t289 + t291 * t97;
t65 = -t109 * t295 + t111 * t299;
t64 = t109 * t299 + t111 * t295;
t63 = -t125 * t289 + t291 * t92;
t60 = t176 * t289 + t291 * t78;
t59 = -t176 * t291 + t289 * t78;
t58 = -t150 * t289 + t291 * t77;
t52 = -t296 * t93 + t300 * t94;
t46 = t300 * t84 - t362;
t45 = -pkin(10) * t96 + t48;
t44 = -pkin(10) * t91 + t47;
t43 = -t295 * t67 + t299 * t68;
t42 = t295 * t68 + t299 * t67;
t41 = -t296 * t64 + t300 * t65;
t40 = -t295 * t62 + t299 * t63;
t39 = t295 * t63 + t299 * t62;
t34 = -t295 * t57 + t299 * t58;
t33 = t295 * t58 + t299 * t57;
t30 = pkin(4) * t31;
t29 = -pkin(4) * t137 + qJ(5) * t32;
t28 = -pkin(5) * t96 + t38;
t27 = -pkin(5) * t91 + t37;
t26 = -qJ(5) * t109 - t31;
t25 = -pkin(4) * t169 + qJ(5) * t111 + t32;
t24 = -t296 * t42 + t300 * t43;
t23 = -t296 * t39 + t300 * t40;
t19 = -t296 * t33 + t300 * t34;
t18 = t299 * t32 - t363;
t17 = t295 * t32 + t361;
t16 = -pkin(10) * t75 - t21;
t15 = -qJ(5) * t67 - t28 * t289 + t291 * t45;
t14 = t291 * t22 + t289 * t50;
t11 = -qJ(5) * t62 - t27 * t289 + t291 * t44;
t10 = -pkin(4) * t96 + qJ(5) * t68 + t28 * t291 + t289 * t45;
t9 = -pkin(4) * t91 + qJ(5) * t63 + t27 * t291 + t289 * t44;
t8 = -qJ(5) * t57 + t16 * t291 + t364 * t75;
t7 = qJ(5) * t58 + t16 * t289 + t325 * t75;
t6 = -t17 * t296 + t18 * t300;
t5 = -t13 * t295 + t14 * t299;
t4 = t13 * t299 + t14 * t295;
t3 = -qJ(5) * t13 + (-pkin(10) * t291 + t364) * t21;
t2 = qJ(5) * t14 + (-pkin(10) * t289 + t325) * t21;
t1 = -t296 * t4 + t300 * t5;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t336, 0, 0, 0, 0, 0, 0, (qJDD(2) * t301 - t297 * t303) * t290, (-qJDD(2) * t297 - t301 * t303) * t290, 0, t244 * t292 + (t229 * t297 - t301 * t312) * t290, 0, 0, 0, 0, 0, 0, t292 * (t266 * t300 + t272 * t296) + (t239 * t297 + t262 * t301) * t290, t292 * (-t267 * t296 + t271 * t300) + (t240 * t297 - t260 * t301) * t290, (t263 * t297 + t264 * t301) * t290, t292 * (-t194 * t300 + t195 * t296) + (t148 * t297 - t215 * t301) * t290, 0, 0, 0, 0, 0, 0, t292 * (t189 * t300 + t190 * t296) + (t143 * t297 - t196 * t301) * t290, t292 * (t202 * t300 + t203 * t296) + (t159 * t297 - t301 * t367) * t290, t292 * (t157 * t300 + t158 * t296) + (t114 * t297 - t213 * t301) * t290, t292 * (t296 * t84 + t300 * t83) + (-t191 * t301 + t297 * t46) * t290, 0, 0, 0, 0, 0, 0, t292 * (t296 * t94 + t300 * t93) + (-t151 * t301 + t297 * t52) * t290, t292 * (t115 * t300 + t116 * t296) + (-t155 * t301 + t297 * t69) * t290, t292 * (t296 * t65 + t300 * t64) + (-t169 * t301 + t297 * t41) * t290, t292 * (t17 * t300 + t18 * t296) + (-t137 * t301 + t297 * t6) * t290, 0, 0, 0, 0, 0, 0, t292 * (t296 * t40 + t300 * t39) + (t23 * t297 - t301 * t91) * t290, t292 * (t296 * t43 + t300 * t42) + (t24 * t297 - t301 * t96) * t290, t292 * (t296 * t34 + t300 * t33) + (t19 * t297 - t301 * t75) * t290, t292 * (t296 * t5 + t300 * t4) + (t1 * t297 - t21 * t301) * t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t312, -t229, 0, 0, (t261 + t323) * t296, t260 * t300 + t262 * t296, t339 + t300 * (-t282 + t302), t262 * t300, t296 * (t283 - t302) + t338, 0, pkin(2) * t262 + pkin(8) * t239 - t215 * t300, -pkin(2) * t260 + pkin(8) * t240 + t215 * t296, pkin(2) * t264 + pkin(8) * t263 + t148, -pkin(2) * t215 + pkin(8) * t148, t296 * (t212 * t299 - t257 * t342) + t300 * (t212 * t295 + t257 * t341), t296 * (-t196 * t299 - t295 * t367) + t300 * (-t196 * t295 + t299 * t367), t296 * (-t246 * t295 + t372) + t300 * (t246 * t299 + t374), t296 * (-t211 * t295 - t256 * t341) + t300 * (t211 * t299 - t256 * t342), t296 * (t245 * t299 - t346) + t300 * (t245 * t295 + t345), (t296 * (t256 * t299 + t257 * t295) + t300 * (t256 * t295 - t257 * t299)) * t286, t296 * (-pkin(9) * t189 + t351) + t300 * (-pkin(3) * t196 + pkin(9) * t190 - t350) - pkin(2) * t196 + pkin(8) * t143, t296 * (-pkin(9) * t202 + t350) + t300 * (-pkin(3) * t367 + pkin(9) * t203 + t351) - pkin(2) * t367 + pkin(8) * t159, t296 * (-pkin(9) * t157 - t83) + t300 * (-pkin(3) * t213 + pkin(9) * t158 + t84) - pkin(2) * t213 + pkin(8) * t114, -pkin(9) * t362 + t300 * (-pkin(3) * t191 + pkin(9) * t84) - pkin(2) * t191 + pkin(8) * t46, t296 * (-t146 * t295 + t147 * t299) + t300 * (t146 * t299 + t147 * t295), t296 * (-t108 * t295 + t110 * t299) + t300 * (t108 * t299 + t110 * t295), t296 * (-t164 * t295 + t166 * t299) + t300 * (t164 * t299 + t166 * t295), t296 * (-t144 * t295 + t145 * t299) + t300 * (t144 * t299 + t145 * t295), t296 * (-t165 * t295 + t167 * t299) + t300 * (t165 * t299 + t167 * t295), t296 * (-t178 * t295 + t179 * t299) + t300 * (t178 * t299 + t179 * t295), t296 * (-pkin(9) * t93 - t295 * t79 + t299 * t89) + t300 * (-pkin(3) * t151 + pkin(9) * t94 + t295 * t89 + t299 * t79) - pkin(2) * t151 + pkin(8) * t52, t296 * (-pkin(9) * t115 - t295 * t82 + t299 * t98) + t300 * (-pkin(3) * t155 + pkin(9) * t116 + t295 * t98 + t299 * t82) - pkin(2) * t155 + pkin(8) * t69, t296 * (-pkin(9) * t64 - t25 * t295 + t26 * t299) + t300 * (-pkin(3) * t169 + pkin(9) * t65 + t25 * t299 + t26 * t295) - pkin(2) * t169 + pkin(8) * t41, t296 * (-pkin(9) * t17 - qJ(5) * t361 - t29 * t295) + t300 * (-pkin(3) * t137 + pkin(9) * t18 - qJ(5) * t363 + t29 * t299) - pkin(2) * t137 + pkin(8) * t6, t296 * (-t295 * t86 + t299 * t88) + t300 * (t295 * t88 + t299 * t86), t296 * (-t295 * t59 + t299 * t60) + t300 * (t295 * t60 + t299 * t59), t296 * (-t295 * t70 + t299 * t72) + t300 * (t295 * t72 + t299 * t70), t296 * (-t295 * t85 + t299 * t87) + t300 * (t295 * t87 + t299 * t85), t296 * (-t295 * t71 + t299 * t73) + t300 * (t295 * t73 + t299 * t71), t296 * (-t112 * t295 + t113 * t299) + t300 * (t112 * t299 + t113 * t295), t296 * (-pkin(9) * t39 + t11 * t299 - t295 * t9) + t300 * (-pkin(3) * t91 + pkin(9) * t40 + t11 * t295 + t299 * t9) - pkin(2) * t91 + pkin(8) * t23, t296 * (-pkin(9) * t42 - t10 * t295 + t15 * t299) + t300 * (-pkin(3) * t96 + pkin(9) * t43 + t10 * t299 + t15 * t295) - pkin(2) * t96 + pkin(8) * t24, t296 * (-pkin(9) * t33 - t295 * t7 + t299 * t8) + t300 * (-pkin(3) * t75 + pkin(9) * t34 + t295 * t8 + t299 * t7) - pkin(2) * t75 + pkin(8) * t19, t296 * (-pkin(9) * t4 - t2 * t295 + t299 * t3) + t300 * (-pkin(3) * t21 + pkin(9) * t5 + t2 * t299 + t295 * t3) - pkin(2) * t21 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t274, -t283 + t282, t329, t274, t281, qJDD(3), -t194, -t195, 0, 0, -t238, t237, t201, t238, t306, t285, pkin(3) * t189 - t130, pkin(3) * t202 - t131, pkin(3) * t157, pkin(3) * t83, t192, t188, t156, -t192, t308, t285, pkin(3) * t93 + t309, pkin(3) * t115 + t314, pkin(3) * t64 + t107, pkin(3) * t17 + t30, t119, t76, t103, t117, t104, t135, pkin(3) * t39 + t320, pkin(3) * t42 + t321, pkin(3) * t33 + t313, pkin(3) * t4 + t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t237, t201, t238, t306, t285, -t130, -t131, 0, 0, t192, t188, t156, -t192, t308, t285, t309, t314, t107, t30, t119, t76, t103, t117, t104, t135, t320, t321, t313, t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t155, t169, t137, 0, 0, 0, 0, 0, 0, t91, t96, t75, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t176, t127, -t177, -t123, t171, -t37, -t38, 0, 0;];
tauJ_reg  = t12;