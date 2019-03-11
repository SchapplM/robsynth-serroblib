% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:44:00
% EndTime: 2019-03-09 08:44:13
% DurationCPUTime: 5.32s
% Computational Cost: add. (5503->521), mult. (11391->623), div. (0->0), fcn. (7364->10), ass. (0->260)
t346 = cos(qJ(5));
t276 = qJD(5) * t346;
t200 = sin(qJ(2));
t300 = qJD(1) * t200;
t368 = t346 * t300 + t276;
t199 = sin(qJ(5));
t294 = qJD(5) * t199;
t367 = -t199 * t300 - t294;
t202 = cos(qJ(2));
t291 = qJD(1) * qJD(2);
t273 = t202 * t291;
t290 = t200 * qJDD(1);
t227 = t273 + t290;
t123 = qJDD(5) + t227;
t189 = pkin(9) + qJ(5);
t179 = cos(t189);
t201 = sin(qJ(1));
t203 = cos(qJ(1));
t256 = g(1) * t203 + g(2) * t201;
t238 = t256 * t202;
t340 = g(3) * t200;
t361 = t238 + t340;
t195 = sin(pkin(9));
t338 = pkin(2) + qJ(4);
t336 = -pkin(8) - t338;
t133 = t336 * t195;
t196 = cos(pkin(9));
t134 = t336 * t196;
t79 = t346 * t133 + t199 * t134;
t366 = -t79 * t123 - t361 * t179;
t160 = qJD(5) + t300;
t231 = -t199 * t195 + t346 * t196;
t325 = -t195 * t368 + t367 * t196;
t263 = t231 * t123 + t160 * t325;
t299 = qJD(1) * t202;
t278 = t196 * t299;
t293 = t195 * qJD(2);
t120 = -t278 - t293;
t279 = t195 * t299;
t298 = qJD(2) * t196;
t121 = -t279 + t298;
t69 = -t346 * t120 + t121 * t199;
t365 = -qJD(2) * t69 + t263;
t364 = t69 ^ 2;
t234 = -t199 * t120 - t346 * t121;
t349 = t234 ^ 2;
t348 = pkin(3) + pkin(7);
t363 = t160 * t69;
t324 = t367 * t195 + t196 * t368;
t362 = t160 * t234;
t102 = t231 * t202;
t171 = pkin(7) * t300;
t360 = qJD(3) + t171;
t232 = -t346 * t195 - t199 * t196;
t233 = -t199 * t133 + t346 * t134;
t314 = t195 * t200;
t245 = pkin(4) * t202 - pkin(8) * t314;
t175 = pkin(2) * t300;
t322 = qJ(3) * t202;
t251 = qJ(4) * t200 - t322;
t106 = qJD(1) * t251 + t175;
t172 = pkin(7) * t299;
t131 = pkin(3) * t299 + t172;
t62 = -t106 * t195 + t196 * t131;
t44 = qJD(1) * t245 + t62;
t63 = t196 * t106 + t195 * t131;
t53 = pkin(8) * t196 * t300 + t63;
t359 = -qJD(4) * t232 - qJD(5) * t233 + t199 * t44 + t346 * t53;
t358 = -qJD(4) * t231 - qJD(5) * t79 + t199 * t53 - t346 * t44;
t167 = pkin(4) * t196 + pkin(3);
t305 = t167 * t300 + t360;
t228 = t123 * t232 - t324 * t160;
t357 = qJD(2) * t234 + t228;
t192 = qJD(2) * qJ(3);
t356 = qJD(4) + t192;
t289 = t202 * qJDD(1);
t147 = t195 * t289;
t266 = qJDD(2) * t196 - t147;
t274 = t200 * t291;
t218 = t195 * t274 + t266;
t303 = t195 * qJDD(2) + t196 * t289;
t226 = t196 * t274 - t303;
t25 = -t120 * t276 + t121 * t294 - t199 * t226 - t346 * t218;
t26 = -qJD(5) * t234 + t199 * t218 - t346 * t226;
t355 = t26 * pkin(5) + t25 * qJ(6) + qJD(6) * t234;
t354 = -t231 * t25 - t234 * t325;
t110 = t131 + t356;
t296 = qJD(2) * t202;
t353 = t200 * (-t110 + t356) + t338 * t296;
t292 = pkin(3) * t300 + t360;
t101 = -t338 * qJD(2) + t292;
t181 = t200 * qJ(3);
t268 = -pkin(1) - t181;
t224 = -t338 * t202 + t268;
t92 = t224 * qJD(1);
t48 = t196 * t101 - t195 * t92;
t49 = t195 * t101 + t196 * t92;
t253 = t195 * t48 - t196 * t49;
t188 = g(3) * t202;
t310 = t200 * t203;
t311 = t200 * t201;
t285 = -g(1) * t310 - g(2) * t311 + t188;
t159 = pkin(2) * t274;
t295 = qJD(3) * t200;
t209 = qJD(2) * t251 - qJD(4) * t202 - t295;
t41 = qJD(1) * t209 + qJDD(1) * t224 + t159;
t158 = pkin(7) * t273;
t168 = pkin(7) * t290;
t270 = qJDD(3) + t158 + t168;
t64 = t227 * pkin(3) - qJD(2) * qJD(4) - t338 * qJDD(2) + t270;
t21 = -t195 * t41 + t196 * t64;
t329 = t21 * t196;
t352 = t253 * t300 - t285 - t329;
t185 = t202 * pkin(2);
t302 = t185 + t181;
t321 = qJ(4) * t202;
t254 = t302 + t321;
t119 = -pkin(1) - t254;
t143 = t348 * t200;
t127 = t196 * t143;
t55 = pkin(4) * t200 + t127 + (pkin(8) * t202 - t119) * t195;
t313 = t196 * t202;
t75 = t196 * t119 + t195 * t143;
t61 = -pkin(8) * t313 + t75;
t240 = t199 * t55 + t346 * t61;
t132 = t348 * t296;
t297 = qJD(2) * t200;
t174 = pkin(2) * t297;
t80 = t174 + t209;
t51 = t196 * t132 - t195 * t80;
t35 = qJD(2) * t245 + t51;
t281 = t196 * t297;
t52 = t195 * t132 + t196 * t80;
t40 = pkin(8) * t281 + t52;
t351 = -qJD(5) * t240 - t199 * t40 + t346 * t35;
t12 = pkin(4) * t227 - pkin(8) * t218 + t21;
t22 = t195 * t64 + t196 * t41;
t17 = pkin(8) * t226 + t22;
t30 = pkin(4) * t300 - pkin(8) * t121 + t48;
t34 = pkin(8) * t120 + t49;
t239 = t199 * t12 + t346 * t17 + t30 * t276 - t294 * t34;
t320 = qJ(6) * t123;
t1 = qJD(6) * t160 + t239 + t320;
t269 = -t346 * t12 + t199 * t17 + t34 * t276 + t30 * t294;
t345 = pkin(5) * t123;
t2 = qJDD(6) + t269 - t345;
t10 = -t199 * t34 + t346 * t30;
t304 = qJD(6) - t10;
t6 = -t160 * pkin(5) + t304;
t11 = t199 * t30 + t346 * t34;
t7 = t160 * qJ(6) + t11;
t350 = -t1 * t232 - t2 * t231 + t324 * t7 - t325 * t6 + t285;
t344 = g(1) * t201;
t341 = g(2) * t203;
t180 = t195 * pkin(4);
t339 = t234 * t69;
t335 = t324 * pkin(5) - t325 * qJ(6) - qJD(6) * t231 + t305;
t333 = -qJ(6) * t299 - t359;
t332 = pkin(5) * t299 - t358;
t331 = t11 * t160;
t328 = t22 * t195;
t323 = pkin(7) * qJDD(2);
t317 = qJDD(2) * pkin(2);
t193 = t200 ^ 2;
t205 = qJD(1) ^ 2;
t315 = t193 * t205;
t309 = t200 * t205;
t308 = t201 * t202;
t307 = t202 * t203;
t165 = qJ(3) + t180;
t144 = t348 * t202;
t194 = t202 ^ 2;
t301 = t193 - t194;
t190 = qJDD(2) * qJ(3);
t288 = pkin(4) * t314;
t287 = t202 * t309;
t169 = pkin(7) * t289;
t191 = qJD(2) * qJD(3);
t286 = t169 + t190 + t191;
t111 = pkin(4) * t313 + t144;
t284 = t348 * qJD(2);
t277 = g(3) * t302;
t275 = t303 * pkin(4);
t272 = t195 * t290;
t271 = t196 * t290;
t267 = -qJD(2) * pkin(2) + qJD(3);
t265 = t203 * pkin(1) + pkin(2) * t307 + t201 * pkin(7) + qJ(3) * t310;
t264 = -t168 - t285;
t261 = qJD(1) * t284;
t178 = sin(t189);
t97 = t178 * t201 - t179 * t310;
t99 = t178 * t203 + t179 * t311;
t260 = g(1) * t99 + g(2) * t97;
t100 = -t178 * t311 + t179 * t203;
t98 = t178 * t310 + t179 * t201;
t259 = -g(1) * t100 - g(2) * t98;
t204 = qJD(2) ^ 2;
t258 = pkin(7) * t204 + t341;
t154 = qJ(3) * t308;
t157 = qJ(3) * t307;
t257 = -g(1) * t157 - g(2) * t154;
t255 = -t341 + t344;
t136 = t171 + t267;
t142 = -t172 - t192;
t250 = t136 * t202 + t142 * t200;
t249 = pkin(3) * t289 + qJDD(4) + t286;
t248 = t196 * t266;
t247 = (-pkin(7) - t167) * qJD(2);
t244 = t268 - t185;
t242 = -t199 * t61 + t346 * t55;
t118 = t244 * qJD(1);
t237 = t118 * t300 + qJDD(3) - t264;
t236 = -0.2e1 * pkin(1) * t291 - t323;
t235 = t199 * t35 + t55 * t276 - t294 * t61 + t346 * t40;
t95 = t200 * t247;
t230 = -qJ(3) * t296 - t295;
t76 = -pkin(4) * t120 + t110;
t229 = qJD(1) * t247;
t67 = -t200 * t261 + t249;
t225 = -t202 * t67 + t256;
t223 = 0.2e1 * qJDD(1) * pkin(1) - t258;
t222 = pkin(5) * t178 - qJ(6) * t179 + t180;
t103 = t232 * t202;
t137 = -pkin(1) - t302;
t220 = t323 + (-qJD(1) * t137 - t118) * qJD(2);
t219 = g(1) * t97 - g(2) * t99 + t179 * t188 - t269;
t217 = t25 + t363;
t215 = t67 - t361;
t108 = t174 + t230;
t66 = qJD(1) * t230 + qJDD(1) * t244 + t159;
t214 = qJD(1) * t108 + qJDD(1) * t137 + t258 + t66;
t212 = -t238 + t249;
t20 = pkin(5) * t69 + qJ(6) * t234 + t76;
t211 = -t20 * t234 + qJDD(6) - t219;
t105 = t270 - t317;
t93 = pkin(7) * t274 - t286;
t210 = qJD(2) * t250 + t105 * t200 - t93 * t202;
t208 = t123 * t233 - t178 * t361;
t207 = -g(1) * t98 + g(2) * t100 + t178 * t188 + t239;
t38 = t200 * t229 + t249 + t275;
t206 = t26 - t362;
t197 = -pkin(8) - qJ(4);
t186 = t203 * pkin(7);
t163 = g(1) * t308;
t130 = t200 * t284;
t128 = -qJ(3) * t299 + t175;
t74 = -t119 * t195 + t127;
t65 = -pkin(5) * t232 - qJ(6) * t231 + t165;
t58 = -t199 * t200 * t293 - qJD(5) * t103 + t346 * t281;
t57 = -qJD(5) * t102 - t232 * t297;
t39 = pkin(5) * t102 - qJ(6) * t103 + t111;
t29 = -pkin(5) * t234 + qJ(6) * t69;
t24 = -t200 * pkin(5) - t242;
t23 = qJ(6) * t200 + t240;
t16 = -pkin(5) * t58 - qJ(6) * t57 - qJD(6) * t103 + t95;
t13 = -t25 + t363;
t5 = -pkin(5) * t296 - t351;
t4 = qJ(6) * t296 + qJD(6) * t200 + t235;
t3 = t38 + t355;
t8 = [qJDD(1), t255, t256, qJDD(1) * t193 + 0.2e1 * t200 * t273, 0.2e1 * t200 * t289 - 0.2e1 * t291 * t301, qJDD(2) * t200 + t202 * t204, qJDD(2) * t202 - t200 * t204, 0, t200 * t236 + t202 * t223 + t163, t236 * t202 + (-t223 - t344) * t200 (t193 + t194) * qJDD(1) * pkin(7) + t210 - t256, t200 * t220 + t202 * t214 - t163, t220 * t202 + (-t214 + t344) * t200, pkin(7) * t210 - g(1) * t186 - g(2) * t265 + t118 * t108 + t66 * t137 - t244 * t344, t130 * t120 + t144 * t303 + (qJD(1) * t74 + t48) * t296 - t225 * t196 + (-t110 * t298 + qJDD(1) * t74 + t21 + t255 * t195 + (-t144 * t298 + t51) * qJD(1)) * t200, -t130 * t121 + t144 * t266 + (-qJD(1) * t75 - t49) * t296 + t225 * t195 + (t110 * t293 - t75 * qJDD(1) - t22 + t255 * t196 + (t144 * t293 - t52) * qJD(1)) * t200, t52 * t120 - t75 * t303 - t51 * t121 - t74 * t266 + t163 + ((-t195 * t74 + t196 * t75) * qJD(1) - t253) * t297 + (t195 * t21 - t196 * t22 - t341) * t202, t22 * t75 + t49 * t52 + t21 * t74 + t48 * t51 + t67 * t144 - t110 * t130 - g(1) * (pkin(3) * t203 + t186) - g(2) * (qJ(4) * t307 + t265) + (-g(1) * (t244 - t321) - g(2) * pkin(3)) * t201, -t103 * t25 - t234 * t57, t102 * t25 - t103 * t26 - t234 * t58 - t57 * t69, t103 * t123 + t160 * t57 - t200 * t25 - t234 * t296, -t102 * t123 + t160 * t58 - t200 * t26 - t296 * t69, t123 * t200 + t160 * t296, t10 * t296 + t38 * t102 + t111 * t26 + t242 * t123 + t160 * t351 - t269 * t200 - t76 * t58 + t95 * t69 + t259, t38 * t103 - t11 * t296 - t111 * t25 - t123 * t240 - t160 * t235 - t200 * t239 - t234 * t95 + t76 * t57 + t260, t102 * t3 - t123 * t24 + t16 * t69 - t160 * t5 - t2 * t200 - t20 * t58 + t26 * t39 - t296 * t6 + t259, -g(2) * t307 - t1 * t102 + t103 * t2 - t23 * t26 - t234 * t5 - t24 * t25 - t4 * t69 + t57 * t6 + t58 * t7 + t163, t1 * t200 - t103 * t3 + t123 * t23 + t16 * t234 + t160 * t4 - t20 * t57 + t25 * t39 + t296 * t7 - t260, t1 * t23 + t7 * t4 + t3 * t39 + t20 * t16 + t2 * t24 + t6 * t5 - g(1) * (pkin(5) * t100 + qJ(6) * t99 + t167 * t203 + t186) - g(2) * (pkin(5) * t98 + qJ(6) * t97 - t197 * t307 + t203 * t288 + t265) + (-g(1) * (t197 * t202 + t244 - t288) - g(2) * t167) * t201; 0, 0, 0, -t287, t301 * t205, t290, t289, qJDD(2), pkin(1) * t309 + t264, t340 - t169 + (pkin(1) * t205 + t256) * t202 (-pkin(2) * t200 + t322) * qJDD(1) + ((-t142 - t192) * t200 + (-t136 + t267) * t202) * qJD(1), -t128 * t299 + t237 - 0.2e1 * t317, t169 + 0.2e1 * t190 + 0.2e1 * t191 + (qJD(1) * t128 - g(3)) * t200 + (qJD(1) * t118 - t256) * t202, -t93 * qJ(3) - t142 * qJD(3) - t105 * pkin(2) - t118 * t128 - g(1) * (-pkin(2) * t310 + t157) - g(2) * (-pkin(2) * t311 + t154) - t277 - t250 * qJD(1) * pkin(7), -t338 * t271 + qJ(3) * t303 - t292 * t120 + t215 * t195 + (-t196 * t353 - t62 * t200 - t48 * t202) * qJD(1), t338 * t272 - qJ(3) * t147 + t292 * t121 + (t215 + t190) * t196 + (t195 * t353 + t63 * t200 + t49 * t202) * qJD(1), t338 * t248 - t63 * t120 + (t196 * qJD(4) + t62) * t121 + (-qJD(4) * t120 + t303 * t338 - t22) * t195 + t352, t67 * qJ(3) - t49 * t63 - t48 * t62 - g(3) * t254 + t292 * t110 + (-t49 * t195 - t48 * t196) * qJD(4) + t257 + (t256 * t200 - t328 - t329) * t338, t354, -t231 * t26 - t232 * t25 + t234 * t324 - t325 * t69, t234 * t299 + t263, t299 * t69 + t228, -t160 * t299, -t10 * t299 + t160 * t358 + t165 * t26 - t232 * t38 + t305 * t69 + t324 * t76 + t208, t11 * t299 + t160 * t359 - t165 * t25 + t231 * t38 - t234 * t305 + t325 * t76 + t366, -t160 * t332 + t20 * t324 - t232 * t3 + t26 * t65 + t299 * t6 + t335 * t69 + t208, t233 * t25 - t234 * t332 - t26 * t79 - t333 * t69 - t350, t333 * t160 - t325 * t20 - t231 * t3 + t335 * t234 + t25 * t65 - t7 * t299 - t366, t1 * t79 + t3 * t65 - t2 * t233 - t277 + t333 * t7 + t332 * t6 + t335 * t20 + (g(3) * t197 - t222 * t256) * t202 + (-g(3) * t222 + t256 * (pkin(2) - t197)) * t200 + t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, qJDD(2) + t287, -t204 - t315, qJD(2) * t142 + t158 + t237 - t317, t271 - t195 * t315 + (t120 + t278) * qJD(2), -t272 - t196 * t315 + (-t121 - t279) * qJD(2), -t195 * t303 - t248 + (t120 * t196 + t121 * t195) * t300, -qJD(2) * t110 + t328 - t352, 0, 0, 0, 0, 0, t365, t357, t365, t232 * t26 - t324 * t69 - t354, -t357, -qJD(2) * t20 + t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t121 - t298) * t300 + t303 (t120 + t293) * t300 + t266, -t120 ^ 2 - t121 ^ 2, -t120 * t49 + t121 * t48 + (-g(3) - t261) * t200 + t212, 0, 0, 0, 0, 0, t206, -t217, t206, -t349 - t364, t217, t275 + t7 * t69 + t6 * t234 + (-g(3) + t229) * t200 + t212 + t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, t349 - t364, t13, -t26 - t362, t123, t234 * t76 + t219 + t331, t10 * t160 + t69 * t76 - t207, -t29 * t69 - t211 + t331 + 0.2e1 * t345, pkin(5) * t25 - qJ(6) * t26 - (-t11 + t7) * t234 + (t6 - t304) * t69, 0.2e1 * t320 - t20 * t69 - t29 * t234 + (0.2e1 * qJD(6) - t10) * t160 + t207, t1 * qJ(6) - t2 * pkin(5) - t20 * t29 - t6 * t11 - g(1) * (-pkin(5) * t97 + qJ(6) * t98) - g(2) * (pkin(5) * t99 - qJ(6) * t100) + t304 * t7 - (-pkin(5) * t179 - qJ(6) * t178) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123 - t339, t13, -t160 ^ 2 - t349, -t160 * t7 + t211 - t345;];
tau_reg  = t8;
