% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:19
% EndTime: 2019-12-31 21:17:33
% DurationCPUTime: 5.47s
% Computational Cost: add. (30563->445), mult. (64082->636), div. (0->0), fcn. (45892->10), ass. (0->280)
t265 = sin(pkin(9));
t268 = sin(qJ(3));
t269 = sin(qJ(2));
t272 = cos(qJ(3));
t273 = cos(qJ(2));
t241 = (t273 * t268 + t269 * t272) * qJD(1);
t262 = qJD(2) + qJD(3);
t266 = cos(pkin(9));
t224 = t265 * t241 - t266 * t262;
t226 = t266 * t241 + t265 * t262;
t198 = t226 * t224;
t257 = t269 * qJDD(1);
t309 = qJD(1) * qJD(2);
t298 = t273 * t309;
t246 = t257 + t298;
t258 = t273 * qJDD(1);
t299 = t269 * t309;
t247 = t258 - t299;
t290 = t268 * t246 - t272 * t247;
t202 = t241 * qJD(3) + t290;
t349 = -t198 + t202;
t357 = t265 * t349;
t356 = t266 * t349;
t267 = sin(qJ(5));
t271 = cos(qJ(5));
t191 = t271 * t224 + t267 * t226;
t193 = -t267 * t224 + t271 * t226;
t148 = t193 * t191;
t201 = qJDD(5) + t202;
t351 = -t148 + t201;
t355 = t267 * t351;
t312 = qJD(1) * t269;
t239 = -t272 * t273 * qJD(1) + t268 * t312;
t218 = t241 * t239;
t308 = qJDD(2) + qJDD(3);
t347 = -t218 + t308;
t354 = t268 * t347;
t353 = t271 * t351;
t352 = t272 * t347;
t286 = t272 * t246 + t268 * t247;
t203 = -t239 * qJD(3) + t286;
t233 = t262 * t239;
t183 = t203 - t233;
t264 = t273 ^ 2;
t275 = qJD(1) ^ 2;
t270 = sin(qJ(1));
t342 = cos(qJ(1));
t297 = t270 * g(1) - t342 * g(2);
t282 = qJDD(1) * pkin(1) + t297;
t283 = qJD(2) * pkin(2) - pkin(7) * t312;
t205 = t247 * pkin(2) + (pkin(7) * t264 + pkin(6)) * t275 - t283 * t312 + t282;
t120 = -t183 * qJ(4) + (t262 * t241 + t202) * pkin(3) - t205;
t284 = t342 * g(1) + t270 * g(2);
t334 = qJDD(1) * pkin(6);
t278 = -t275 * pkin(1) - t284 + t334;
t230 = -t269 * g(3) + t273 * t278;
t260 = t264 * t275;
t199 = -pkin(2) * t260 + t247 * pkin(7) - qJD(2) * t283 + t230;
t316 = t272 * t199;
t277 = t269 * t278;
t318 = t269 * t275;
t340 = t246 * pkin(7);
t346 = qJDD(2) * pkin(2) - t277 + (pkin(2) * t318 + pkin(7) * t309 - g(3)) * t273 - t340;
t151 = t346 * t268 + t316;
t215 = t239 * pkin(3) - t241 * qJ(4);
t343 = t262 ^ 2;
t131 = -t343 * pkin(3) + t308 * qJ(4) - t239 * t215 + t151;
t75 = 0.2e1 * qJD(4) * t226 - t266 * t120 + t265 * t131;
t76 = -0.2e1 * qJD(4) * t224 + t265 * t120 + t266 * t131;
t45 = t265 * t75 + t266 * t76;
t194 = t266 * t203 + t265 * t308;
t292 = t265 * t203 - t266 * t308;
t119 = -t191 * qJD(5) + t271 * t194 - t267 * t292;
t235 = qJD(5) + t239;
t176 = t235 * t191;
t350 = -t176 + t119;
t209 = t239 * t224;
t163 = -t209 - t194;
t348 = -t209 + t194;
t293 = t267 * t194 + t271 * t292;
t100 = (qJD(5) - t235) * t193 + t293;
t188 = t191 ^ 2;
t189 = t193 ^ 2;
t345 = t224 ^ 2;
t223 = t226 ^ 2;
t234 = t235 ^ 2;
t344 = t239 ^ 2;
t238 = t241 ^ 2;
t341 = pkin(3) * t268;
t54 = pkin(4) * t349 + t163 * pkin(8) - t75;
t206 = t239 * pkin(4) - t226 * pkin(8);
t64 = -t345 * pkin(4) - t292 * pkin(8) - t239 * t206 + t76;
t25 = t267 * t64 - t271 * t54;
t26 = t267 * t54 + t271 * t64;
t15 = -t271 * t25 + t267 * t26;
t339 = t265 * t15;
t338 = t266 * t15;
t150 = t268 * t199 - t272 * t346;
t130 = -t308 * pkin(3) - t343 * qJ(4) + t241 * t215 + qJDD(4) + t150;
t80 = pkin(4) * t292 - t345 * pkin(8) + t226 * t206 + t130;
t337 = t267 * t80;
t336 = t271 * t80;
t335 = -pkin(3) * t130 + qJ(4) * t45;
t333 = t235 * t267;
t332 = t235 * t271;
t331 = t239 * t226;
t330 = t239 * t265;
t329 = t239 * t266;
t328 = t262 * t268;
t327 = t262 * t272;
t124 = t265 * t130;
t165 = t198 + t202;
t326 = t265 * t165;
t125 = t266 * t130;
t325 = t266 * t165;
t133 = t148 + t201;
t324 = t267 * t133;
t323 = t268 * t202;
t322 = t268 * t205;
t213 = t218 + t308;
t321 = t268 * t213;
t109 = -t272 * t150 + t268 * t151;
t320 = t269 * t109;
t253 = t273 * t318;
t249 = qJDD(2) + t253;
t319 = t269 * t249;
t317 = t271 * t133;
t315 = t272 * t205;
t314 = t272 * t213;
t313 = t273 * (qJDD(2) - t253);
t311 = qJD(3) + t262;
t307 = t268 * t148;
t306 = t268 * t198;
t305 = t272 * t148;
t304 = t272 * t198;
t195 = -t223 - t344;
t137 = -t265 * t195 - t325;
t303 = -pkin(3) * t348 + qJ(4) * t137 + t124;
t178 = -t344 - t345;
t127 = t266 * t178 - t357;
t158 = t292 + t331;
t302 = -pkin(3) * t158 + qJ(4) * t127 - t125;
t301 = -pkin(3) * t272 - pkin(2);
t16 = t267 * t25 + t271 * t26;
t122 = -t188 - t189;
t103 = t176 + t119;
t68 = -t100 * t271 + t267 * t103;
t10 = -pkin(4) * t122 + pkin(8) * t68 + t16;
t66 = -t100 * t267 - t271 * t103;
t12 = -pkin(8) * t66 - t15;
t31 = -t265 * t66 + t266 * t68;
t296 = -pkin(3) * t122 + qJ(4) * t31 + t266 * t10 + t265 * t12;
t146 = -t234 - t188;
t84 = t271 * t146 - t355;
t99 = (qJD(5) + t235) * t193 + t293;
t36 = -pkin(4) * t99 + pkin(8) * t84 - t336;
t83 = t267 * t146 + t353;
t47 = -pkin(8) * t83 + t337;
t50 = -t265 * t83 + t266 * t84;
t295 = -pkin(3) * t99 + qJ(4) * t50 + t265 * t47 + t266 * t36;
t167 = -t189 - t234;
t92 = -t267 * t167 - t317;
t40 = -pkin(4) * t350 + pkin(8) * t92 + t337;
t91 = t271 * t167 - t324;
t52 = -pkin(8) * t91 + t336;
t59 = -t265 * t91 + t266 * t92;
t294 = -pkin(3) * t350 + qJ(4) * t59 + t265 * t52 + t266 * t40;
t110 = t268 * t150 + t272 * t151;
t229 = t273 * g(3) + t277;
t291 = t269 * t229 + t273 * t230;
t159 = t292 - t331;
t115 = -t159 * t266 - t265 * t163;
t171 = -t223 - t345;
t289 = -pkin(3) * t171 + qJ(4) * t115 + t45;
t288 = t265 * t76 - t266 * t75;
t14 = -pkin(4) * t80 + pkin(8) * t16;
t8 = t266 * t16 - t339;
t285 = -pkin(3) * t80 - pkin(8) * t339 + qJ(4) * t8 + t266 * t14;
t281 = (-qJD(3) + t262) * t241 - t290;
t274 = qJD(2) ^ 2;
t263 = t269 ^ 2;
t259 = t263 * t275;
t248 = t258 - 0.2e1 * t299;
t245 = t257 + 0.2e1 * t298;
t242 = t275 * pkin(6) + t282;
t232 = -t238 + t343;
t231 = t344 - t343;
t228 = -t238 - t343;
t216 = t238 - t344;
t211 = -t344 - t343;
t208 = -t223 + t344;
t207 = -t344 + t345;
t204 = -t344 - t238;
t200 = t272 * t202;
t197 = -t223 + t345;
t186 = -t268 * t228 - t314;
t185 = t272 * t228 - t321;
t184 = t203 + t233;
t182 = -t311 * t239 + t286;
t179 = t311 * t241 + t290;
t175 = t272 * t211 - t354;
t174 = t268 * t211 + t352;
t173 = -t189 + t234;
t172 = t188 - t234;
t169 = (-t224 * t266 + t226 * t265) * t239;
t168 = (-t224 * t265 - t226 * t266) * t239;
t155 = t266 * t194 - t226 * t330;
t154 = t265 * t194 + t226 * t329;
t153 = t224 * t329 + t265 * t292;
t152 = t224 * t330 - t266 * t292;
t147 = t189 - t188;
t145 = t268 * t184 + t272 * t281;
t144 = -t272 * t184 + t268 * t281;
t143 = t266 * t207 - t326;
t142 = -t265 * t208 + t356;
t141 = t265 * t207 + t325;
t140 = t266 * t208 + t357;
t139 = (-t191 * t271 + t193 * t267) * t235;
t138 = (-t191 * t267 - t193 * t271) * t235;
t136 = t266 * t195 - t326;
t126 = t265 * t178 + t356;
t118 = -t193 * qJD(5) - t293;
t114 = -t266 * t158 - t265 * t348;
t113 = -t159 * t265 + t266 * t163;
t112 = -t265 * t158 + t266 * t348;
t108 = t271 * t172 - t324;
t107 = -t267 * t173 + t353;
t106 = t267 * t172 + t317;
t105 = t271 * t173 + t355;
t96 = t271 * t119 - t193 * t333;
t95 = t267 * t119 + t193 * t332;
t94 = -t267 * t118 + t191 * t332;
t93 = t271 * t118 + t191 * t333;
t90 = t272 * t137 + t268 * t348;
t89 = t268 * t137 - t272 * t348;
t88 = t272 * t127 + t268 * t158;
t87 = t268 * t127 - t272 * t158;
t86 = t272 * t115 + t268 * t171;
t85 = t268 * t115 - t272 * t171;
t82 = -t265 * t138 + t266 * t139;
t81 = t266 * t138 + t265 * t139;
t78 = -qJ(4) * t136 + t125;
t77 = -qJ(4) * t126 + t124;
t72 = -t265 * t106 + t266 * t108;
t71 = -t265 * t105 + t266 * t107;
t70 = t266 * t106 + t265 * t108;
t69 = t266 * t105 + t265 * t107;
t67 = -t267 * t350 - t271 * t99;
t65 = -t267 * t99 + t271 * t350;
t63 = -t265 * t95 + t266 * t96;
t62 = -t265 * t93 + t266 * t94;
t61 = t265 * t96 + t266 * t95;
t60 = t265 * t94 + t266 * t93;
t58 = t265 * t92 + t266 * t91;
t56 = -pkin(3) * t136 + t76;
t55 = -pkin(3) * t126 + t75;
t49 = t265 * t84 + t266 * t83;
t42 = t268 * t350 + t272 * t59;
t41 = t268 * t59 - t272 * t350;
t38 = t268 * t99 + t272 * t50;
t37 = t268 * t50 - t272 * t99;
t33 = -t272 * t130 + t268 * t45;
t32 = -qJ(4) * t113 - t288;
t30 = -t265 * t65 + t266 * t67;
t29 = t265 * t68 + t266 * t66;
t28 = t265 * t67 + t266 * t65;
t23 = t268 * t122 + t272 * t31;
t22 = -t272 * t122 + t268 * t31;
t21 = -pkin(3) * t29 - pkin(4) * t66;
t20 = -pkin(3) * t58 - pkin(4) * t91 + t26;
t19 = -qJ(4) * t58 - t265 * t40 + t266 * t52;
t18 = -pkin(3) * t49 - pkin(4) * t83 + t25;
t17 = -qJ(4) * t49 - t265 * t36 + t266 * t47;
t7 = t265 * t16 + t338;
t5 = t268 * t80 + t272 * t8;
t4 = t268 * t8 - t272 * t80;
t3 = -pkin(3) * t7 - pkin(4) * t15;
t2 = -qJ(4) * t29 - t265 * t10 + t266 * t12;
t1 = -pkin(8) * t338 - qJ(4) * t7 - t265 * t14;
t6 = [0, 0, 0, 0, 0, qJDD(1), t297, t284, 0, 0, (t246 + t298) * t269, t273 * t245 + t269 * t248, t319 + t273 * (-t259 + t274), (t247 - t299) * t273, t269 * (t260 - t274) + t313, 0, t273 * t242 + pkin(1) * t248 + pkin(6) * (t273 * (-t260 - t274) - t319), -t269 * t242 - pkin(1) * t245 + pkin(6) * (-t313 - t269 * (-t259 - t274)), pkin(1) * (t259 + t260) + (t263 + t264) * t334 + t291, pkin(1) * t242 + pkin(6) * t291, t269 * (t272 * t203 - t241 * t328) + t273 * (t268 * t203 + t241 * t327), t269 * (-t272 * t179 - t268 * t183) + t273 * (-t268 * t179 + t272 * t183), t269 * (-t268 * t232 + t352) + t273 * (t272 * t232 + t354), t269 * (t239 * t327 + t323) + t273 * (t239 * t328 - t200), t269 * (t272 * t231 - t321) + t273 * (t268 * t231 + t314), (t269 * (-t239 * t272 + t241 * t268) + t273 * (-t239 * t268 - t241 * t272)) * t262, t269 * (-pkin(7) * t174 - t322) + t273 * (-pkin(2) * t179 + pkin(7) * t175 + t315) - pkin(1) * t179 + pkin(6) * (-t269 * t174 + t273 * t175), t269 * (-pkin(7) * t185 - t315) + t273 * (-pkin(2) * t182 + pkin(7) * t186 - t322) - pkin(1) * t182 + pkin(6) * (-t269 * t185 + t273 * t186), t269 * (-pkin(7) * t144 - t109) + t273 * (-pkin(2) * t204 + pkin(7) * t145 + t110) - pkin(1) * t204 + pkin(6) * (-t269 * t144 + t273 * t145), -pkin(7) * t320 + t273 * (pkin(2) * t205 + pkin(7) * t110) + pkin(1) * t205 + pkin(6) * (t273 * t110 - t320), t269 * (t272 * t155 + t306) + t273 * (t268 * t155 - t304), t269 * (t272 * t114 - t268 * t197) + t273 * (t268 * t114 + t272 * t197), t269 * (t272 * t142 - t268 * t163) + t273 * (t268 * t142 + t272 * t163), t269 * (t272 * t153 - t306) + t273 * (t268 * t153 + t304), t269 * (t272 * t143 - t268 * t159) + t273 * (t268 * t143 + t272 * t159), t269 * (t272 * t169 + t323) + t273 * (t268 * t169 - t200), t269 * (-pkin(7) * t87 - t268 * t55 + t272 * t77) + t273 * (-pkin(2) * t126 + pkin(7) * t88 + t268 * t77 + t272 * t55) - pkin(1) * t126 + pkin(6) * (-t269 * t87 + t273 * t88), t269 * (-pkin(7) * t89 - t268 * t56 + t272 * t78) + t273 * (-pkin(2) * t136 + pkin(7) * t90 + t268 * t78 + t272 * t56) - pkin(1) * t136 + pkin(6) * (-t269 * t89 + t273 * t90), t269 * (-pkin(7) * t85 + t272 * t32) + t273 * (pkin(7) * t86 + t268 * t32) + pkin(6) * (-t269 * t85 + t273 * t86) + (t269 * t341 + t273 * t301 - pkin(1)) * t113, (t269 * (-qJ(4) * t272 + t341) + t273 * (-qJ(4) * t268 + t301) - pkin(1)) * t288 + (pkin(6) + pkin(7)) * (-t269 * t33 + t273 * (t268 * t130 + t272 * t45)), t269 * (t272 * t63 + t307) + t273 * (t268 * t63 - t305), t269 * (t268 * t147 + t272 * t30) + t273 * (-t272 * t147 + t268 * t30), t269 * (t268 * t103 + t272 * t71) + t273 * (-t272 * t103 + t268 * t71), t269 * (t272 * t62 - t307) + t273 * (t268 * t62 + t305), t269 * (-t268 * t100 + t272 * t72) + t273 * (t272 * t100 + t268 * t72), t269 * (t268 * t201 + t272 * t82) + t273 * (-t272 * t201 + t268 * t82), t269 * (-pkin(7) * t37 + t272 * t17 - t268 * t18) + t273 * (-pkin(2) * t49 + pkin(7) * t38 + t268 * t17 + t272 * t18) - pkin(1) * t49 + pkin(6) * (-t269 * t37 + t273 * t38), t269 * (-pkin(7) * t41 + t272 * t19 - t268 * t20) + t273 * (-pkin(2) * t58 + pkin(7) * t42 + t268 * t19 + t272 * t20) - pkin(1) * t58 + pkin(6) * (-t269 * t41 + t273 * t42), t269 * (-pkin(7) * t22 + t272 * t2 - t268 * t21) + t273 * (-pkin(2) * t29 + pkin(7) * t23 + t268 * t2 + t272 * t21) - pkin(1) * t29 + pkin(6) * (-t269 * t22 + t273 * t23), t269 * (-pkin(7) * t4 + t272 * t1 - t268 * t3) + t273 * (-pkin(2) * t7 + pkin(7) * t5 + t268 * t1 + t272 * t3) - pkin(1) * t7 + pkin(6) * (-t269 * t4 + t273 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t259 - t260, t257, t253, t258, qJDD(2), -t229, -t230, 0, 0, t218, t216, t184, -t218, t281, t308, pkin(2) * t174 - t150, -t316 - t268 * (pkin(7) * t298 - t229 - t340) + (-t268 * t249 + t185) * pkin(2), pkin(2) * t144, pkin(2) * t109, t154, t112, t140, t152, t141, t168, pkin(2) * t87 + t302, pkin(2) * t89 + t303, pkin(2) * t85 + t289, pkin(2) * t33 + t335, t61, t28, t69, t60, t70, t81, pkin(2) * t37 + t295, pkin(2) * t41 + t294, pkin(2) * t22 + t296, pkin(2) * t4 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, t216, t184, -t218, t281, t308, -t150, -t151, 0, 0, t154, t112, t140, t152, t141, t168, t302, t303, t289, t335, t61, t28, t69, t60, t70, t81, t295, t294, t296, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t348, t171, t130, 0, 0, 0, 0, 0, 0, t99, t350, t122, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t147, t103, -t148, -t100, t201, -t25, -t26, 0, 0;];
tauJ_reg = t6;
