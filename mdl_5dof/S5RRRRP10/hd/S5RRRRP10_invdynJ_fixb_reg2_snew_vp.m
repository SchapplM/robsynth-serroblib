% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:24
% EndTime: 2019-12-31 22:12:39
% DurationCPUTime: 5.49s
% Computational Cost: add. (25493->417), mult. (55146->581), div. (0->0), fcn. (43163->10), ass. (0->285)
t229 = cos(pkin(5));
t223 = t229 * qJD(1) + qJD(2);
t231 = sin(qJ(3));
t234 = cos(qJ(3));
t228 = sin(pkin(5));
t232 = sin(qJ(2));
t288 = qJD(1) * t232;
t272 = t228 * t288;
t200 = t231 * t223 + t234 * t272;
t230 = sin(qJ(4));
t233 = cos(qJ(4));
t235 = cos(qJ(2));
t287 = qJD(1) * t235;
t221 = t228 * t287;
t282 = t221 - qJD(3);
t181 = t230 * t200 + t233 * t282;
t183 = t233 * t200 - t230 * t282;
t153 = t183 * t181;
t281 = qJDD(1) * t228;
t206 = qJD(2) * t221 + t232 * t281;
t222 = t229 * qJDD(1) + qJDD(2);
t269 = t231 * t206 - t234 * t222;
t167 = -t200 * qJD(3) - t269;
t166 = qJDD(4) - t167;
t345 = -t153 + t166;
t352 = pkin(4) * t345;
t198 = -t234 * t223 + t231 * t272;
t258 = -t234 * t206 - t231 * t222;
t168 = -t198 * qJD(3) - t258;
t280 = qJDD(1) * t235;
t252 = qJD(2) * t288 - t280;
t245 = t252 * t228;
t243 = qJDD(3) + t245;
t129 = -t181 * qJD(4) + t233 * t168 + t230 * t243;
t194 = qJD(4) + t198;
t163 = t194 * t181;
t112 = t163 + t129;
t351 = qJ(5) * t112;
t301 = t230 * t345;
t295 = t233 * t345;
t180 = t183 ^ 2;
t192 = t194 ^ 2;
t150 = -t180 - t192;
t179 = t181 ^ 2;
t271 = t230 * t168 - t233 * t243;
t128 = -t183 * qJD(4) - t271;
t156 = t194 * pkin(4) - t183 * qJ(5);
t236 = qJD(1) ^ 2;
t328 = sin(qJ(1));
t329 = cos(qJ(1));
t249 = t329 * g(1) + t328 * g(2);
t202 = -t236 * pkin(1) + pkin(7) * t281 - t249;
t264 = -t235 * pkin(2) - t232 * pkin(8);
t289 = qJD(1) * t228;
t205 = t264 * t289;
t220 = t223 ^ 2;
t248 = t328 * g(1) - t329 * g(2);
t303 = t228 * t236;
t240 = qJDD(1) * pkin(1) + pkin(7) * t303 + t248;
t239 = t229 * t240;
t304 = t228 * t232;
t238 = -g(3) * t304 + t232 * t239;
t141 = t222 * pkin(8) - t220 * pkin(2) + (t205 * t289 + t202) * t235 + t238;
t323 = t229 * g(3);
t237 = -t206 * pkin(8) - t323 + (-t223 * pkin(8) * t287 + (-t280 + (qJD(2) + t223) * t288) * pkin(2) - t240) * t228;
t100 = t234 * t141 + t231 * t237;
t172 = t198 * pkin(3) - t200 * pkin(9);
t279 = t282 ^ 2;
t80 = -t279 * pkin(3) + t243 * pkin(9) - t198 * t172 + t100;
t270 = t232 * t202 - t235 * t239;
t322 = t235 * g(3);
t140 = -t222 * pkin(2) - t220 * pkin(8) + (t205 * t288 + t322) * t228 + t270;
t189 = t282 * t198;
t147 = t168 + t189;
t268 = t282 * t200;
t83 = -t147 * pkin(9) + (-t167 - t268) * pkin(3) + t140;
t48 = t230 * t83 + t233 * t80;
t253 = t128 * qJ(5) - 0.2e1 * qJD(5) * t181 - t194 * t156 + t48;
t350 = -t253 + (t150 + t179) * pkin(4);
t225 = t228 ^ 2;
t348 = t225 * t235;
t306 = t200 * t198;
t241 = t243 - t306;
t347 = t231 * t241;
t346 = t234 * t241;
t344 = -t163 + t129;
t170 = t228 * t322 + t270;
t171 = t235 * t202 + t238;
t343 = t232 * t170 + t235 * t171;
t109 = (qJD(4) - t194) * t183 + t271;
t196 = t198 ^ 2;
t197 = t200 ^ 2;
t131 = -t179 - t180;
t74 = -t109 * t233 + t230 * t112;
t56 = -t234 * t131 + t231 * t74;
t341 = pkin(2) * t56;
t108 = (qJD(4) + t194) * t183 + t271;
t136 = -t192 - t179;
t89 = t233 * t136 - t301;
t60 = -t234 * t108 + t231 * t89;
t340 = pkin(2) * t60;
t118 = t153 + t166;
t296 = t233 * t118;
t93 = -t230 * t150 - t296;
t63 = t231 * t93 - t234 * t344;
t339 = pkin(2) * t63;
t88 = t230 * t136 + t295;
t338 = pkin(3) * t88;
t302 = t230 * t118;
t92 = t233 * t150 - t302;
t337 = pkin(3) * t92;
t285 = qJD(5) * t183;
t176 = -0.2e1 * t285;
t47 = t230 * t80 - t233 * t83;
t247 = -t351 - t47 + t352;
t34 = t176 + t247;
t336 = pkin(4) * t34;
t335 = pkin(8) * t56;
t334 = pkin(8) * t60;
t333 = pkin(8) * t63;
t72 = -t109 * t230 - t233 * t112;
t332 = pkin(9) * t72;
t331 = pkin(9) * t88;
t330 = pkin(9) * t92;
t327 = pkin(3) * t231;
t326 = pkin(3) * t234;
t325 = pkin(4) * t112;
t324 = pkin(7) * t228;
t57 = t231 * t131 + t234 * t74;
t321 = pkin(1) * (-t228 * t56 + (t232 * t57 - t235 * t72) * t229) + (t232 * t72 + t235 * t57) * t324;
t61 = t231 * t108 + t234 * t89;
t320 = pkin(1) * (-t228 * t60 + (t232 * t61 - t235 * t88) * t229) + (t232 * t88 + t235 * t61) * t324;
t64 = t231 * t344 + t234 * t93;
t319 = pkin(1) * (-t228 * t63 + (t232 * t64 - t235 * t92) * t229) + (t232 * t92 + t235 * t64) * t324;
t318 = -pkin(2) * t72 + pkin(8) * t57;
t317 = -pkin(2) * t88 + pkin(8) * t61;
t316 = -pkin(2) * t92 + pkin(8) * t64;
t315 = t230 * t34;
t99 = t231 * t141 - t234 * t237;
t79 = -t243 * pkin(3) - t279 * pkin(9) + t200 * t172 + t99;
t314 = t230 * t79;
t313 = t233 * t34;
t312 = t233 * t79;
t311 = -pkin(3) * t131 + pkin(9) * t74;
t310 = -pkin(3) * t108 + pkin(9) * t89;
t309 = -pkin(3) * t344 + pkin(9) * t93;
t308 = t194 * t230;
t307 = t194 * t233;
t305 = t225 * t236;
t300 = t231 * t140;
t160 = -t306 - t243;
t299 = t231 * t160;
t215 = t235 * t232 * t305;
t203 = t215 + t222;
t297 = t232 * t203;
t294 = t234 * t140;
t293 = t234 * t160;
t204 = -t215 + t222;
t291 = t235 * t204;
t284 = 0.2e1 * qJD(3) - t221;
t226 = t232 ^ 2;
t278 = t226 * t305;
t227 = t235 ^ 2;
t277 = t227 * t305;
t276 = t231 * t153;
t275 = t234 * t153;
t274 = t235 * t306;
t210 = t223 * t221;
t273 = t210 + t206;
t28 = t230 * t47 + t233 * t48;
t66 = t234 * t100 + t231 * t99;
t267 = -pkin(3) * t79 + pkin(9) * t28;
t266 = t310 - t312;
t265 = t309 + t314;
t263 = t231 * t189;
t262 = t231 * t268;
t261 = t234 * t189;
t260 = t234 * t268;
t27 = t230 * t48 - t233 * t47;
t259 = t231 * t100 - t234 * t99;
t257 = t28 + t311;
t256 = qJD(1) * t223 - t229 * t236;
t255 = -pkin(1) + t264;
t144 = t221 * t200 + t269;
t49 = -t128 * pkin(4) - t179 * qJ(5) + t183 * t156 + qJDD(5) + t79;
t45 = -qJ(5) * t150 + t49;
t75 = -pkin(4) * t344 - qJ(5) * t118;
t251 = t230 * t45 + t233 * t75 + t309;
t26 = -qJ(5) * t109 + (-t131 - t179) * pkin(4) + t253;
t177 = 0.2e1 * t285;
t29 = t177 - t247 + t351;
t250 = t230 * t29 + t233 * t26 + t311;
t37 = -pkin(4) * t108 + qJ(5) * t136 - t49;
t246 = -qJ(5) * t301 + t233 * t37 + t310;
t244 = t247 + t352;
t36 = -t179 * pkin(4) + t253;
t15 = t233 * t36 - t315;
t19 = -pkin(4) * t49 + qJ(5) * t36;
t242 = -pkin(3) * t49 + pkin(9) * t15 - qJ(5) * t315 + t233 * t19;
t209 = t223 * t272;
t208 = (t226 - t227) * t305;
t207 = -t220 - t277;
t193 = -t278 - t220;
t190 = t228 * t240 + t323;
t188 = -t209 - t245;
t187 = t209 - t245;
t186 = -t210 + t206;
t185 = -t197 + t279;
t184 = t196 - t279;
t174 = -t197 - t279;
t173 = t197 - t196;
t169 = -t279 - t196;
t158 = -t180 + t192;
t157 = t179 - t192;
t155 = t196 + t197;
t154 = t263 + t260;
t151 = t180 - t179;
t149 = t284 * t198 + t258;
t148 = t168 - t189;
t145 = -t284 * t200 - t269;
t143 = t231 * t168 - t260;
t142 = t234 * t167 - t263;
t135 = t231 * t184 - t293;
t134 = t234 * t185 + t347;
t133 = -t231 * t174 + t293;
t132 = t234 * t174 + t299;
t126 = t234 * t169 - t347;
t125 = t231 * t169 + t346;
t121 = (-t181 * t233 + t183 * t230) * t194;
t120 = (-t181 * t230 - t183 * t233) * t194;
t116 = -t144 * t234 + t231 * t148;
t114 = t231 * t145 + t234 * t147;
t105 = t233 * t129 - t183 * t308;
t104 = t230 * t129 + t183 * t307;
t103 = -t230 * t128 + t181 * t307;
t102 = -t233 * t128 - t181 * t308;
t101 = t231 * t121 - t234 * t166;
t97 = t233 * t157 - t302;
t96 = -t230 * t158 + t295;
t95 = t230 * t157 + t296;
t94 = t233 * t158 + t301;
t85 = t231 * t105 - t275;
t84 = t231 * t103 + t275;
t78 = pkin(2) * t149 + pkin(8) * t133 + t300;
t76 = pkin(2) * t145 + pkin(8) * t126 - t294;
t73 = -t233 * t108 - t230 * t344;
t71 = -t230 * t108 + t233 * t344;
t68 = t234 * t109 + t231 * t97;
t67 = -t234 * t112 + t231 * t96;
t58 = -t234 * t151 + t231 * t73;
t54 = -pkin(2) * t140 + pkin(8) * t66;
t53 = t312 - t330;
t52 = t229 * t101 + (t232 * (t234 * t121 + t231 * t166) - t235 * t120) * t228;
t51 = t314 - t331;
t50 = -pkin(3) * t72 + t325;
t44 = pkin(2) * t155 + pkin(8) * t116 + t66;
t41 = t229 * t85 + (t232 * (t234 * t105 + t276) - t235 * t104) * t228;
t40 = t229 * t84 + (t232 * (t234 * t103 - t276) + t235 * t102) * t228;
t39 = t48 - t337;
t38 = t47 - t338;
t33 = t229 * t68 + (t232 * (-t231 * t109 + t234 * t97) - t235 * t95) * t228;
t32 = t229 * t67 + (t232 * (t231 * t112 + t234 * t96) - t235 * t94) * t228;
t25 = -t337 - t350;
t24 = -t230 * t75 + t233 * t45 - t330;
t23 = -qJ(5) * t295 - t230 * t37 - t331;
t22 = t229 * t58 + (t232 * (t231 * t151 + t234 * t73) - t235 * t71) * t228;
t21 = t177 - t244 - t338;
t18 = t231 * t79 + t234 * t28;
t17 = t231 * t28 - t234 * t79;
t16 = -t27 - t332;
t14 = t230 * t36 + t313;
t13 = t231 * t53 + t234 * t39 + t316;
t12 = t231 * t51 + t234 * t38 + t317;
t11 = t234 * t15 + t231 * t49;
t10 = t231 * t15 - t234 * t49;
t9 = -t230 * t26 + t233 * t29 - t332;
t8 = -pkin(3) * t14 - t336;
t7 = t231 * t16 - t72 * t326 + t318;
t6 = t231 * t24 + t234 * t25 + t316;
t5 = t234 * t21 + t231 * t23 + t317;
t4 = t231 * t9 + t234 * t50 + t318;
t3 = -pkin(9) * t14 - qJ(5) * t313 - t230 * t19;
t2 = pkin(8) * t18 + (-pkin(9) * t231 - pkin(2) - t326) * t27;
t1 = -pkin(2) * t14 + pkin(8) * t11 + t231 * t3 + t234 * t8;
t20 = [0, 0, 0, 0, 0, qJDD(1), t248, t249, 0, 0, (t206 * t228 + t256 * t348) * t232, t229 * t208 + (t232 * t188 + t235 * t273) * t228, t229 * t186 + (t297 + t235 * (t220 - t278)) * t228, (-t232 * t256 - t252) * t348, t229 * t187 + (t232 * (-t220 + t277) + t291) * t228, t229 * t222, (-t170 + pkin(1) * (t203 * t235 + t207 * t232)) * t229 + (t235 * t190 + pkin(1) * t188 + pkin(7) * (t235 * t207 - t297)) * t228, -t190 * t304 - t229 * t171 + pkin(1) * (-t228 * t273 + (t235 * t193 - t232 * t204) * t229) + (-t232 * t193 - t291) * t324, pkin(1) * ((-t235 * t186 + t232 * t187) * t229 - (-t226 - t227) * t225 * t303) + (t232 * t186 + t235 * t187) * t324 + t343 * t228, pkin(1) * (t228 * t190 + (-t170 * t235 + t171 * t232) * t229) + t343 * t324, t229 * t143 + (t232 * (t234 * t168 + t262) - t274) * t228, t229 * t114 + (t232 * (t234 * t145 - t231 * t147) - t235 * t173) * t228, t229 * t134 + (t232 * (-t231 * t185 + t346) - t235 * t148) * t228, t229 * t142 + (t232 * (-t231 * t167 - t261) + t274) * t228, t229 * t135 + (t232 * (t234 * t184 + t299) + t235 * t144) * t228, t229 * t154 + (t232 * (t261 - t262) - t243 * t235) * t228, (t76 + pkin(1) * (t126 * t232 + t145 * t235)) * t229 + (t232 * (-pkin(8) * t125 + t300) + t235 * (-pkin(2) * t125 + t99) - pkin(1) * t125 + pkin(7) * (t235 * t126 - t232 * t145)) * t228, (t78 + pkin(1) * (t133 * t232 + t149 * t235)) * t229 + (t232 * (-pkin(8) * t132 + t294) + t235 * (-pkin(2) * t132 + t100) - pkin(1) * t132 + pkin(7) * (t235 * t133 - t232 * t149)) * t228, (t44 + pkin(1) * (t116 * t232 + t155 * t235)) * t229 + (-t232 * t259 + pkin(7) * (t235 * t116 - t232 * t155) + t255 * (-t144 * t231 - t234 * t148)) * t228, (t54 + pkin(1) * (-t140 * t235 + t232 * t66)) * t229 + (pkin(7) * (t232 * t140 + t235 * t66) + t255 * t259) * t228, t41, t22, t32, t40, t33, t52, t229 * t12 + (t232 * (-t231 * t38 + t234 * t51 - t334) + t235 * (-t266 - t340)) * t228 + t320, t229 * t13 + (t232 * (-t231 * t39 + t234 * t53 - t333) + t235 * (-t265 - t339)) * t228 + t319, t229 * t7 + (t232 * (t234 * t16 + t72 * t327 - t335) + t235 * (-t257 - t341)) * t228 + t321, (t2 + pkin(1) * (t18 * t232 - t235 * t27)) * t229 + (t232 * (-pkin(8) * t17 + (-pkin(9) * t234 + t327) * t27) + t235 * (-pkin(2) * t17 - t267) - pkin(1) * t17 + pkin(7) * (t235 * t18 + t232 * t27)) * t228, t41, t22, t32, t40, t33, t52, t229 * t5 + (t232 * (-t231 * t21 + t234 * t23 - t334) + t235 * (-t246 - t340)) * t228 + t320, t229 * t6 + (t232 * (-t231 * t25 + t234 * t24 - t333) + t235 * (-t251 - t339)) * t228 + t319, t229 * t4 + (t232 * (-t231 * t50 + t234 * t9 - t335) + t235 * (-t250 - t341)) * t228 + t321, (t1 + pkin(1) * (t11 * t232 - t14 * t235)) * t229 + (t232 * (-pkin(8) * t10 - t231 * t8 + t234 * t3) + t235 * (-pkin(2) * t10 - t242) - pkin(1) * t10 + pkin(7) * (t235 * t11 + t232 * t14)) * t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t208, t186, t215, t187, t222, -t170, -t171, 0, 0, t143, t114, t134, t142, t135, t154, t76, t78, t44, t54, t85, t58, t67, t84, t68, t101, t12, t13, t7, t2, t85, t58, t67, t84, t68, t101, t5, t6, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, t173, t148, -t306, -t144, t243, -t99, -t100, 0, 0, t104, t71, t94, -t102, t95, t120, t266, t265, t257, t267, t104, t71, t94, -t102, t95, t120, t246, t251, t250, t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t151, t112, -t153, -t109, t166, -t47, -t48, 0, 0, t153, t151, t112, -t153, -t109, t166, t176 + t244, t350, -t325, t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t344, t131, t49;];
tauJ_reg = t20;
