% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:14:59
% EndTime: 2019-05-05 18:15:13
% DurationCPUTime: 6.59s
% Computational Cost: add. (37297->501), mult. (83104->733), div. (0->0), fcn. (59760->12), ass. (0->298)
t353 = 2 * qJD(4);
t267 = sin(pkin(11));
t269 = cos(pkin(11));
t277 = cos(qJ(3));
t274 = sin(qJ(3));
t311 = qJD(1) * t274;
t233 = qJD(1) * t269 * t277 - t267 * t311;
t234 = (t267 * t277 + t269 * t274) * qJD(1);
t219 = t233 * t234;
t342 = qJDD(3) + t219;
t352 = t267 * t342;
t351 = t269 * t342;
t272 = sin(qJ(6));
t273 = sin(qJ(5));
t276 = cos(qJ(5));
t210 = t233 * t273 + t234 * t276;
t303 = qJD(1) * qJD(3);
t293 = t277 * t303;
t302 = t274 * qJDD(1);
t242 = t293 + t302;
t294 = t274 * t303;
t301 = t277 * qJDD(1);
t285 = -t294 + t301;
t220 = -t242 * t267 + t269 * t285;
t221 = t269 * t242 + t267 * t285;
t289 = -t276 * t220 + t273 * t221;
t158 = -qJD(5) * t210 - t289;
t157 = qJDD(6) - t158;
t263 = qJD(3) + qJD(5);
t275 = cos(qJ(6));
t189 = t210 * t272 - t275 * t263;
t191 = t210 * t275 + t263 * t272;
t164 = t191 * t189;
t344 = t157 - t164;
t350 = t272 * t344;
t208 = -t276 * t233 + t234 * t273;
t180 = t210 * t208;
t262 = qJDD(3) + qJDD(5);
t343 = -t180 + t262;
t349 = t273 * t343;
t348 = t275 * t344;
t347 = t276 * t343;
t279 = qJD(1) ^ 2;
t338 = sin(qJ(1));
t339 = cos(qJ(1));
t284 = t339 * g(1) + t338 * g(2);
t240 = -t279 * pkin(1) - t284;
t268 = sin(pkin(10));
t270 = cos(pkin(10));
t283 = t338 * g(1) - t339 * g(2);
t281 = qJDD(1) * pkin(1) + t283;
t312 = t270 * t240 + t268 * t281;
t212 = -pkin(2) * t279 + qJDD(1) * pkin(7) + t312;
t265 = -g(3) + qJDD(2);
t197 = t274 * t212 - t277 * t265;
t314 = t277 * t279;
t252 = t274 * t314;
t247 = qJDD(3) + t252;
t177 = (-t242 + t293) * qJ(4) + t247 * pkin(3) - t197;
t295 = qJ(4) * t311;
t246 = qJD(3) * pkin(3) - t295;
t317 = t274 * t265;
t179 = t317 + (-t246 - t295) * qJD(3) + (-pkin(3) * t314 + qJ(4) * qJDD(1) + t212) * t277;
t127 = -t269 * t177 + t267 * t179 + t234 * t353;
t310 = qJD(3) * t233;
t287 = -t221 + t310;
t346 = t287 * pkin(8) - t127;
t176 = pkin(5) * t208 - pkin(9) * t210;
t341 = t263 ^ 2;
t280 = pkin(4) * t342 + t346;
t128 = t267 * t177 + t269 * t179 + t233 * t353;
t224 = qJD(3) * pkin(4) - pkin(8) * t234;
t231 = t233 ^ 2;
t119 = -pkin(4) * t231 + pkin(8) * t220 - qJD(3) * t224 + t128;
t316 = t276 * t119;
t68 = t273 * t280 + t316;
t52 = -t341 * pkin(5) + t262 * pkin(9) - t208 * t176 + t68;
t288 = -t268 * t240 + t270 * t281;
t211 = -qJDD(1) * pkin(2) - t279 * pkin(7) - t288;
t340 = t277 ^ 2;
t261 = t340 * t279;
t183 = -t285 * pkin(3) - qJ(4) * t261 + t246 * t311 + qJDD(4) + t211;
t134 = -t220 * pkin(4) - t231 * pkin(8) + t234 * t224 + t183;
t159 = -qJD(5) * t208 + t220 * t273 + t221 * t276;
t203 = t263 * t208;
t143 = t159 - t203;
t76 = (t210 * t263 - t158) * pkin(5) - t143 * pkin(9) + t134;
t38 = t272 * t52 - t275 * t76;
t39 = t272 * t76 + t275 * t52;
t21 = t272 * t38 + t275 * t39;
t205 = qJD(6) + t208;
t290 = t272 * t159 - t275 * t262;
t113 = (qJD(6) - t205) * t191 + t290;
t187 = t189 ^ 2;
t188 = t191 ^ 2;
t204 = t205 ^ 2;
t206 = t208 ^ 2;
t207 = t210 ^ 2;
t232 = t234 ^ 2;
t337 = pkin(5) * t273;
t67 = t119 * t273 - t276 * t280;
t51 = -t262 * pkin(5) - t341 * pkin(9) + t176 * t210 + t67;
t336 = -pkin(5) * t51 + pkin(9) * t21;
t40 = t273 * t68 - t276 * t67;
t335 = t267 * t40;
t334 = t269 * t40;
t48 = t272 * t51;
t78 = -t127 * t269 + t128 * t267;
t333 = t274 * t78;
t49 = t275 * t51;
t122 = t157 + t164;
t332 = t122 * t272;
t331 = t122 * t275;
t330 = t134 * t273;
t329 = t134 * t276;
t173 = t180 + t262;
t328 = t173 * t273;
t327 = t173 * t276;
t326 = t183 * t267;
t325 = t183 * t269;
t324 = t205 * t272;
t323 = t205 * t275;
t215 = qJDD(3) - t219;
t322 = t215 * t267;
t321 = t215 * t269;
t320 = t263 * t273;
t319 = t263 * t276;
t318 = t274 * t247;
t248 = qJDD(3) - t252;
t315 = t277 * t248;
t309 = qJD(3) * t234;
t308 = qJD(3) * t267;
t307 = qJD(3) * t269;
t304 = qJD(6) + t205;
t286 = -t275 * t159 - t272 * t262;
t118 = t304 * t189 + t286;
t154 = -t188 - t204;
t91 = -t154 * t272 - t331;
t300 = pkin(5) * t118 + pkin(9) * t91 + t48;
t114 = -t304 * t191 - t290;
t147 = -t204 - t187;
t86 = t147 * t275 - t350;
t299 = pkin(5) * t114 + pkin(9) * t86 - t49;
t298 = t273 * t164;
t297 = t276 * t164;
t296 = -pkin(5) * t276 - pkin(4);
t41 = t273 * t67 + t276 * t68;
t146 = t187 + t188;
t130 = -qJD(6) * t189 - t286;
t169 = t205 * t189;
t117 = t130 + t169;
t73 = -t113 * t275 + t117 * t272;
t292 = pkin(5) * t146 + pkin(9) * t73 + t21;
t79 = t127 * t267 + t269 * t128;
t198 = t277 * t212 + t317;
t162 = t274 * t197 + t277 * t198;
t20 = t272 * t39 - t275 * t38;
t194 = t220 + t309;
t243 = -0.2e1 * t294 + t301;
t282 = (-qJD(5) + t263) * t210 - t289;
t278 = qJD(3) ^ 2;
t264 = t274 ^ 2;
t260 = t264 * t279;
t251 = -t261 - t278;
t250 = -t260 - t278;
t245 = t260 + t261;
t244 = (t264 + t340) * qJDD(1);
t241 = 0.2e1 * t293 + t302;
t227 = -t232 - t278;
t226 = -t232 + t278;
t225 = t231 - t278;
t223 = -t250 * t274 - t315;
t222 = t251 * t277 - t318;
t213 = -t278 - t231;
t201 = -t207 + t341;
t200 = t206 - t341;
t199 = -t207 - t341;
t195 = t221 + t310;
t193 = -t220 + t309;
t192 = -t231 - t232;
t185 = -t227 * t267 - t321;
t184 = t227 * t269 - t322;
t182 = t213 * t269 - t352;
t181 = t213 * t267 + t351;
t178 = t207 - t206;
t171 = -t341 - t206;
t168 = -t188 + t204;
t167 = t187 - t204;
t166 = (-t208 * t276 + t210 * t273) * t263;
t165 = (-t208 * t273 - t210 * t276) * t263;
t163 = t188 - t187;
t161 = t194 * t269 - t267 * t287;
t160 = t194 * t267 + t269 * t287;
t156 = -t206 - t207;
t155 = -t184 * t274 + t185 * t277;
t153 = t200 * t276 - t328;
t152 = -t201 * t273 + t347;
t151 = t200 * t273 + t327;
t150 = t201 * t276 + t349;
t149 = -t199 * t273 - t327;
t148 = t199 * t276 - t328;
t144 = t159 + t203;
t139 = (qJD(5) + t263) * t210 + t289;
t138 = t159 * t276 - t210 * t320;
t137 = t159 * t273 + t210 * t319;
t136 = -t158 * t273 + t208 * t319;
t135 = t158 * t276 + t208 * t320;
t133 = -t181 * t274 + t182 * t277;
t132 = t171 * t276 - t349;
t131 = t171 * t273 + t347;
t129 = -qJD(6) * t191 - t290;
t126 = (-t189 * t275 + t191 * t272) * t205;
t125 = (-t189 * t272 - t191 * t275) * t205;
t120 = -t160 * t274 + t161 * t277;
t116 = t130 - t169;
t110 = t130 * t275 - t191 * t324;
t109 = t130 * t272 + t191 * t323;
t108 = -t129 * t272 + t189 * t323;
t107 = t129 * t275 + t189 * t324;
t105 = -t148 * t267 + t149 * t269;
t104 = t148 * t269 + t149 * t267;
t103 = -pkin(8) * t148 + t329;
t102 = t126 * t276 + t157 * t273;
t101 = t126 * t273 - t157 * t276;
t100 = t167 * t275 - t332;
t99 = -t168 * t272 + t348;
t98 = t167 * t272 + t331;
t97 = t168 * t275 + t350;
t96 = t144 * t273 + t276 * t282;
t95 = -t139 * t276 - t143 * t273;
t94 = -t144 * t276 + t273 * t282;
t93 = -t139 * t273 + t143 * t276;
t92 = -pkin(8) * t131 + t330;
t90 = t154 * t275 - t332;
t88 = -t131 * t267 + t132 * t269;
t87 = t131 * t269 + t132 * t267;
t85 = t147 * t272 + t348;
t83 = t110 * t276 + t298;
t82 = t108 * t276 - t298;
t81 = t110 * t273 - t297;
t80 = t108 * t273 + t297;
t77 = -pkin(4) * t143 + pkin(8) * t149 + t330;
t74 = -pkin(4) * t139 + pkin(8) * t132 - t329;
t72 = t114 * t275 - t116 * t272;
t71 = -t113 * t272 - t117 * t275;
t70 = t114 * t272 + t116 * t275;
t65 = -t104 * t274 + t105 * t277;
t64 = t100 * t276 - t113 * t273;
t63 = t117 * t273 + t276 * t99;
t62 = t100 * t273 + t113 * t276;
t61 = -t117 * t276 + t273 * t99;
t60 = -t118 * t273 + t276 * t91;
t59 = t118 * t276 + t273 * t91;
t58 = -t114 * t273 + t276 * t86;
t57 = t114 * t276 + t273 * t86;
t56 = -t267 * t94 + t269 * t96;
t55 = t267 * t96 + t269 * t94;
t54 = t163 * t273 + t276 * t72;
t53 = -t163 * t276 + t273 * t72;
t47 = -t146 * t273 + t276 * t73;
t46 = t146 * t276 + t273 * t73;
t45 = -t274 * t87 + t277 * t88;
t44 = t277 * t79 - t333;
t43 = -pkin(9) * t90 + t49;
t42 = -pkin(9) * t85 + t48;
t35 = -pkin(4) * t134 + pkin(8) * t41;
t34 = -t267 * t59 + t269 * t60;
t33 = t267 * t60 + t269 * t59;
t32 = -t267 * t57 + t269 * t58;
t31 = t267 * t58 + t269 * t57;
t30 = -t274 * t55 + t277 * t56;
t29 = -t267 * t46 + t269 * t47;
t28 = t267 * t47 + t269 * t46;
t27 = -pkin(8) * t94 - t40;
t26 = -pkin(5) * t90 + t39;
t25 = -pkin(5) * t85 + t38;
t24 = -pkin(4) * t156 + pkin(8) * t96 + t41;
t23 = t269 * t41 - t335;
t22 = t267 * t41 + t334;
t18 = -t274 * t33 + t277 * t34;
t17 = -t274 * t31 + t277 * t32;
t16 = -t274 * t28 + t277 * t29;
t15 = -pkin(9) * t71 - t20;
t14 = t21 * t276 + t273 * t51;
t13 = t21 * t273 - t276 * t51;
t12 = -pkin(8) * t59 - t26 * t273 + t276 * t43;
t11 = -pkin(8) * t57 - t25 * t273 + t276 * t42;
t10 = -pkin(4) * t90 + pkin(8) * t60 + t26 * t276 + t273 * t43;
t9 = -pkin(4) * t85 + pkin(8) * t58 + t25 * t276 + t273 * t42;
t8 = -pkin(8) * t46 + t15 * t276 + t71 * t337;
t7 = -t22 * t274 + t23 * t277;
t6 = pkin(8) * t47 + t273 * t15 + t296 * t71;
t5 = -t13 * t267 + t14 * t269;
t4 = t13 * t269 + t14 * t267;
t3 = -pkin(8) * t13 + (-pkin(9) * t276 + t337) * t20;
t2 = pkin(8) * t14 + (-pkin(9) * t273 + t296) * t20;
t1 = -t274 * t4 + t277 * t5;
t19 = [0, 0, 0, 0, 0, qJDD(1), t283, t284, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t270 - t268 * t279) + t288, pkin(1) * (-qJDD(1) * t268 - t270 * t279) - t312, 0, pkin(1) * (t268 * t312 + t270 * t288), (t242 + t293) * t274, t241 * t277 + t243 * t274, t318 + t277 * (-t260 + t278), t243 * t277, t274 * (t261 - t278) + t315, 0, -t277 * t211 + pkin(2) * t243 + pkin(7) * t222 + pkin(1) * (t222 * t268 + t243 * t270), t274 * t211 - pkin(2) * t241 + pkin(7) * t223 + pkin(1) * (t223 * t268 - t241 * t270), pkin(2) * t245 + pkin(7) * t244 + pkin(1) * (t244 * t268 + t245 * t270) + t162, -pkin(2) * t211 + pkin(7) * t162 + pkin(1) * (t162 * t268 - t211 * t270), t274 * (t221 * t269 - t234 * t308) + t277 * (t221 * t267 + t234 * t307), t274 * (-t193 * t269 - t195 * t267) + t277 * (-t193 * t267 + t195 * t269), t274 * (-t226 * t267 + t351) + t277 * (t226 * t269 + t352), t274 * (-t220 * t267 - t233 * t307) + t277 * (t220 * t269 - t233 * t308), t274 * (t225 * t269 - t322) + t277 * (t225 * t267 + t321), (t274 * (t233 * t269 + t234 * t267) + t277 * (t233 * t267 - t234 * t269)) * qJD(3), t274 * (-qJ(4) * t181 + t326) + t277 * (-pkin(3) * t193 + qJ(4) * t182 - t325) - pkin(2) * t193 + pkin(7) * t133 + pkin(1) * (t133 * t268 - t193 * t270), t274 * (-qJ(4) * t184 + t325) + t277 * (-pkin(3) * t195 + qJ(4) * t185 + t326) - pkin(2) * t195 + pkin(7) * t155 + pkin(1) * (t155 * t268 - t195 * t270), t274 * (-qJ(4) * t160 - t78) + t277 * (-pkin(3) * t192 + qJ(4) * t161 + t79) - pkin(2) * t192 + pkin(7) * t120 + pkin(1) * (t120 * t268 - t192 * t270), -qJ(4) * t333 + t277 * (-pkin(3) * t183 + qJ(4) * t79) - pkin(2) * t183 + pkin(7) * t44 + pkin(1) * (-t183 * t270 + t268 * t44), t274 * (-t137 * t267 + t138 * t269) + t277 * (t137 * t269 + t138 * t267), t274 * (-t267 * t93 + t269 * t95) + t277 * (t267 * t95 + t269 * t93), t274 * (-t150 * t267 + t152 * t269) + t277 * (t150 * t269 + t152 * t267), t274 * (-t135 * t267 + t136 * t269) + t277 * (t135 * t269 + t136 * t267), t274 * (-t151 * t267 + t153 * t269) + t277 * (t151 * t269 + t153 * t267), t274 * (-t165 * t267 + t166 * t269) + t277 * (t165 * t269 + t166 * t267), t274 * (-qJ(4) * t87 - t267 * t74 + t269 * t92) + t277 * (-pkin(3) * t139 + qJ(4) * t88 + t267 * t92 + t269 * t74) - pkin(2) * t139 + pkin(7) * t45 + pkin(1) * (-t139 * t270 + t268 * t45), t274 * (-qJ(4) * t104 + t103 * t269 - t267 * t77) + t277 * (-pkin(3) * t143 + qJ(4) * t105 + t103 * t267 + t269 * t77) - pkin(2) * t143 + pkin(7) * t65 + pkin(1) * (-t143 * t270 + t268 * t65), t274 * (-qJ(4) * t55 - t24 * t267 + t269 * t27) + t277 * (-pkin(3) * t156 + qJ(4) * t56 + t24 * t269 + t267 * t27) - pkin(2) * t156 + pkin(7) * t30 + pkin(1) * (-t156 * t270 + t268 * t30), t274 * (-pkin(8) * t334 - qJ(4) * t22 - t267 * t35) + t277 * (-pkin(3) * t134 - pkin(8) * t335 + qJ(4) * t23 + t269 * t35) - pkin(2) * t134 + pkin(7) * t7 + pkin(1) * (-t134 * t270 + t268 * t7), t274 * (-t267 * t81 + t269 * t83) + t277 * (t267 * t83 + t269 * t81), t274 * (-t267 * t53 + t269 * t54) + t277 * (t267 * t54 + t269 * t53), t274 * (-t267 * t61 + t269 * t63) + t277 * (t267 * t63 + t269 * t61), t274 * (-t267 * t80 + t269 * t82) + t277 * (t267 * t82 + t269 * t80), t274 * (-t267 * t62 + t269 * t64) + t277 * (t267 * t64 + t269 * t62), t274 * (-t101 * t267 + t102 * t269) + t277 * (t101 * t269 + t102 * t267), t274 * (-qJ(4) * t31 + t11 * t269 - t267 * t9) + t277 * (-pkin(3) * t85 + qJ(4) * t32 + t11 * t267 + t269 * t9) - pkin(2) * t85 + pkin(7) * t17 + pkin(1) * (t17 * t268 - t270 * t85), t274 * (-qJ(4) * t33 - t10 * t267 + t12 * t269) + t277 * (-pkin(3) * t90 + qJ(4) * t34 + t10 * t269 + t12 * t267) - pkin(2) * t90 + pkin(7) * t18 + pkin(1) * (t18 * t268 - t270 * t90), t274 * (-qJ(4) * t28 - t267 * t6 + t269 * t8) + t277 * (-pkin(3) * t71 + qJ(4) * t29 + t267 * t8 + t269 * t6) - pkin(2) * t71 + pkin(7) * t16 + pkin(1) * (t16 * t268 - t270 * t71), t274 * (-qJ(4) * t4 - t2 * t267 + t269 * t3) + t277 * (-pkin(3) * t20 + qJ(4) * t5 + t2 * t269 + t267 * t3) - pkin(2) * t20 + pkin(7) * t1 + pkin(1) * (t1 * t268 - t20 * t270); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, 0, 0, 0, 0, 0, 0, t247 * t277 + t251 * t274, -t248 * t274 + t250 * t277, 0, -t197 * t277 + t198 * t274, 0, 0, 0, 0, 0, 0, t181 * t277 + t182 * t274, t184 * t277 + t185 * t274, t160 * t277 + t161 * t274, t274 * t79 + t277 * t78, 0, 0, 0, 0, 0, 0, t274 * t88 + t277 * t87, t104 * t277 + t105 * t274, t274 * t56 + t277 * t55, t22 * t277 + t23 * t274, 0, 0, 0, 0, 0, 0, t274 * t32 + t277 * t31, t274 * t34 + t277 * t33, t274 * t29 + t277 * t28, t274 * t5 + t277 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252, -t261 + t260, t302, t252, t301, qJDD(3), -t197, -t198, 0, 0, -t219, t232 - t231, -t287, t219, t194, qJDD(3), pkin(3) * t181 - t127, pkin(3) * t184 - t128, pkin(3) * t160, pkin(3) * t78, t180, t178, t144, -t180, t282, t262, pkin(3) * t87 + pkin(4) * t131 - t67, pkin(3) * t104 - t316 - t273 * t346 + (-t273 * t342 + t148) * pkin(4), pkin(3) * t55 + pkin(4) * t94, pkin(3) * t22 + pkin(4) * t40, t109, t70, t97, t107, t98, t125, pkin(3) * t31 + pkin(4) * t57 + t299, pkin(3) * t33 + pkin(4) * t59 + t300, pkin(3) * t28 + pkin(4) * t46 + t292, pkin(3) * t4 + pkin(4) * t13 + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t195, t192, t183, 0, 0, 0, 0, 0, 0, t139, t143, t156, t134, 0, 0, 0, 0, 0, 0, t85, t90, t71, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t178, t144, -t180, t282, t262, -t67, -t68, 0, 0, t109, t70, t97, t107, t98, t125, t299, t300, t292, t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t163, t117, -t164, -t113, t157, -t38, -t39, 0, 0;];
tauJ_reg  = t19;
