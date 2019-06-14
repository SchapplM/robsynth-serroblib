% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:10:54
% EndTime: 2019-05-05 15:11:07
% DurationCPUTime: 6.66s
% Computational Cost: add. (32807->464), mult. (75359->692), div. (0->0), fcn. (57196->12), ass. (0->278)
t259 = sin(pkin(11));
t270 = t259 ^ 2;
t261 = cos(pkin(11));
t272 = t261 ^ 2;
t329 = qJD(1) ^ 2;
t238 = (t270 + t272) * t329;
t265 = sin(qJ(4));
t268 = cos(qJ(4));
t257 = -g(3) + qJDD(2);
t327 = sin(qJ(1));
t328 = cos(qJ(1));
t281 = t328 * g(1) + t327 * g(2);
t236 = -t329 * pkin(1) - t281;
t260 = sin(pkin(10));
t262 = cos(pkin(10));
t280 = t327 * g(1) - t328 * g(2);
t277 = qJDD(1) * pkin(1) + t280;
t306 = t262 * t236 + t260 * t277;
t339 = -t329 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t306;
t184 = t259 * t257 + t261 * t339;
t297 = t261 * qJDD(1);
t302 = t272 * t329;
t276 = -pkin(3) * t302 + pkin(7) * t297 + t184;
t245 = t261 * t257;
t334 = t245 + (pkin(3) * t329 * t261 - pkin(7) * qJDD(1) - t339) * t259;
t145 = t265 * t276 - t268 * t334;
t282 = t259 * t268 + t261 * t265;
t228 = t282 * qJDD(1);
t229 = (-t259 * t265 + t261 * t268) * qJD(1);
t304 = t229 * qJD(4);
t216 = t228 + t304;
t346 = -t145 + (-t216 + t304) * pkin(8);
t263 = sin(qJ(6));
t230 = t282 * qJD(1);
t264 = sin(qJ(5));
t267 = cos(qJ(5));
t205 = t229 * t264 + t230 * t267;
t298 = t259 * qJDD(1);
t185 = -t265 * t298 + t268 * t297;
t303 = t230 * qJD(4);
t214 = t185 - t303;
t286 = -t267 * t214 + t216 * t264;
t159 = -qJD(5) * t205 - t286;
t158 = qJDD(6) - t159;
t255 = qJD(4) + qJD(5);
t266 = cos(qJ(6));
t188 = t205 * t263 - t266 * t255;
t190 = t205 * t266 + t255 * t263;
t162 = t190 * t188;
t337 = t158 - t162;
t345 = t263 * t337;
t203 = -t267 * t229 + t230 * t264;
t174 = t205 * t203;
t254 = qJDD(4) + qJDD(5);
t336 = -t174 + t254;
t344 = t264 * t336;
t217 = t229 * t230;
t333 = qJDD(4) + t217;
t343 = t265 * t333;
t342 = t266 * t337;
t341 = t267 * t336;
t340 = t268 * t333;
t285 = -t260 * t236 + t262 * t277;
t200 = -qJDD(1) * pkin(2) - t329 * qJ(3) + qJDD(3) - t285;
t338 = (pkin(1) * t260 + qJ(3)) * t238 + t200 - (pkin(1) * t262 + pkin(2)) * qJDD(1);
t172 = pkin(5) * t203 - pkin(9) * t205;
t330 = t255 ^ 2;
t274 = pkin(4) * t333 + t346;
t146 = t265 * t334 + t268 * t276;
t218 = qJD(4) * pkin(4) - pkin(8) * t230;
t226 = t229 ^ 2;
t119 = -t226 * pkin(4) + t214 * pkin(8) - qJD(4) * t218 + t146;
t307 = t267 * t119;
t76 = t264 * t274 + t307;
t60 = -t330 * pkin(5) + t254 * pkin(9) - t203 * t172 + t76;
t160 = -qJD(5) * t203 + t214 * t264 + t216 * t267;
t196 = t255 * t203;
t138 = t160 - t196;
t182 = -pkin(3) * t297 + t200 + (-t270 * t329 - t302) * pkin(7);
t144 = -t214 * pkin(4) - t226 * pkin(8) + t218 * t230 + t182;
t73 = (t205 * t255 - t159) * pkin(5) - t138 * pkin(9) + t144;
t38 = t263 * t60 - t266 * t73;
t39 = t263 * t73 + t266 * t60;
t21 = t263 * t38 + t266 * t39;
t199 = qJD(6) + t203;
t287 = t160 * t263 - t266 * t254;
t112 = (qJD(6) - t199) * t190 + t287;
t186 = t188 ^ 2;
t187 = t190 ^ 2;
t198 = t199 ^ 2;
t201 = t203 ^ 2;
t202 = t205 ^ 2;
t227 = t230 ^ 2;
t326 = pkin(5) * t264;
t75 = t119 * t264 - t267 * t274;
t59 = -t254 * pkin(5) - t330 * pkin(9) + t172 * t205 + t75;
t325 = -pkin(5) * t59 + pkin(9) * t21;
t56 = t263 * t59;
t42 = t264 * t76 - t267 * t75;
t324 = t265 * t42;
t57 = t266 * t59;
t323 = t268 * t42;
t101 = -t145 * t268 + t146 * t265;
t322 = t101 * t259;
t121 = t158 + t162;
t321 = t121 * t263;
t320 = t121 * t266;
t319 = t144 * t264;
t318 = t144 * t267;
t170 = t174 + t254;
t317 = t170 * t264;
t316 = t170 * t267;
t315 = t182 * t265;
t314 = t182 * t268;
t313 = t199 * t263;
t312 = t199 * t266;
t211 = qJDD(4) - t217;
t311 = t211 * t265;
t310 = t211 * t268;
t309 = t255 * t264;
t308 = t255 * t267;
t300 = qJD(6) + t199;
t283 = -t160 * t266 - t254 * t263;
t117 = t300 * t188 + t283;
t154 = -t187 - t198;
t87 = -t154 * t263 - t320;
t296 = pkin(5) * t117 + pkin(9) * t87 + t56;
t113 = -t300 * t190 - t287;
t147 = -t198 - t186;
t84 = t147 * t266 - t345;
t295 = pkin(5) * t113 + pkin(9) * t84 - t57;
t293 = t264 * t162;
t292 = t267 * t162;
t290 = -pkin(5) * t267 - pkin(4);
t43 = t264 * t75 + t267 * t76;
t142 = t186 + t187;
t126 = -qJD(6) * t188 - t283;
t167 = t199 * t188;
t116 = t126 + t167;
t70 = -t112 * t266 + t116 * t263;
t289 = pkin(5) * t142 + pkin(9) * t70 + t21;
t102 = t145 * t265 + t268 * t146;
t183 = t259 * t339 - t245;
t156 = t259 * t183 + t261 * t184;
t20 = t263 * t39 - t266 * t38;
t279 = (-qJD(5) + t255) * t205 - t286;
t269 = qJD(4) ^ 2;
t252 = t272 * qJDD(1);
t251 = t270 * qJDD(1);
t237 = t252 + t251;
t221 = -t227 - t269;
t220 = -t227 + t269;
t219 = t226 - t269;
t215 = t228 + 0.2e1 * t304;
t213 = -t185 + 0.2e1 * t303;
t209 = -t269 - t226;
t194 = -t202 + t330;
t193 = t201 - t330;
t192 = -t202 - t330;
t191 = -t226 - t227;
t180 = -t221 * t265 - t310;
t179 = t221 * t268 - t311;
t178 = t185 * t268 + t228 * t265;
t177 = t185 * t265 - t228 * t268;
t176 = t209 * t268 - t343;
t175 = t209 * t265 + t340;
t173 = t202 - t201;
t168 = -t330 - t201;
t166 = -t187 + t198;
t165 = t186 - t198;
t164 = (-t203 * t267 + t205 * t264) * t255;
t163 = (-t203 * t264 - t205 * t267) * t255;
t161 = t187 - t186;
t157 = -t201 - t202;
t155 = -t179 * t259 + t180 * t261;
t153 = t193 * t267 - t317;
t152 = -t194 * t264 + t341;
t151 = t193 * t264 + t316;
t150 = t194 * t267 + t344;
t149 = -t192 * t264 - t316;
t148 = t192 * t267 - t317;
t140 = -t177 * t259 + t178 * t261;
t139 = t160 + t196;
t134 = (qJD(5) + t255) * t205 + t286;
t133 = -t175 * t259 + t176 * t261;
t132 = t160 * t267 - t205 * t309;
t131 = t160 * t264 + t205 * t308;
t130 = -t159 * t264 + t203 * t308;
t129 = t159 * t267 + t203 * t309;
t128 = t168 * t267 - t344;
t127 = t168 * t264 + t341;
t125 = -qJD(6) * t190 - t287;
t124 = (-t188 * t266 + t190 * t263) * t199;
t123 = (-t188 * t263 - t190 * t266) * t199;
t115 = t126 - t167;
t109 = t126 * t266 - t190 * t313;
t108 = t126 * t263 + t190 * t312;
t107 = -t125 * t263 + t188 * t312;
t106 = t125 * t266 + t188 * t313;
t105 = -t148 * t265 + t149 * t268;
t104 = t148 * t268 + t149 * t265;
t103 = -pkin(8) * t148 + t318;
t100 = t124 * t267 + t158 * t264;
t99 = t124 * t264 - t158 * t267;
t98 = t165 * t266 - t321;
t97 = -t166 * t263 + t342;
t96 = t165 * t263 + t320;
t95 = t166 * t266 + t345;
t94 = t139 * t264 + t267 * t279;
t93 = -t134 * t267 - t138 * t264;
t92 = -t139 * t267 + t264 * t279;
t91 = -t134 * t264 + t138 * t267;
t90 = -pkin(8) * t127 + t319;
t89 = -t127 * t265 + t128 * t268;
t88 = t127 * t268 + t128 * t265;
t86 = t154 * t266 - t321;
t83 = t147 * t263 + t342;
t81 = t109 * t267 + t293;
t80 = t107 * t267 - t293;
t79 = t109 * t264 - t292;
t78 = t107 * t264 + t292;
t77 = -pkin(4) * t138 + pkin(8) * t149 + t319;
t71 = -pkin(4) * t134 + pkin(8) * t128 - t318;
t69 = t113 * t266 - t115 * t263;
t68 = -t112 * t263 - t116 * t266;
t67 = t113 * t263 + t115 * t266;
t65 = -t104 * t259 + t105 * t261;
t64 = -t112 * t264 + t267 * t98;
t63 = t116 * t264 + t267 * t97;
t62 = t112 * t267 + t264 * t98;
t61 = -t116 * t267 + t264 * t97;
t55 = t102 * t261 - t322;
t54 = -t117 * t264 + t267 * t87;
t53 = t117 * t267 + t264 * t87;
t52 = -t113 * t264 + t267 * t84;
t51 = t113 * t267 + t264 * t84;
t50 = t161 * t264 + t267 * t69;
t49 = -t161 * t267 + t264 * t69;
t48 = -t265 * t92 + t268 * t94;
t47 = t265 * t94 + t268 * t92;
t46 = -t142 * t264 + t267 * t70;
t45 = t267 * t142 + t264 * t70;
t44 = -t259 * t88 + t261 * t89;
t41 = -pkin(9) * t86 + t57;
t40 = -pkin(9) * t83 + t56;
t35 = -pkin(4) * t144 + pkin(8) * t43;
t34 = -pkin(8) * t92 - t42;
t33 = -t265 * t53 + t268 * t54;
t32 = t265 * t54 + t268 * t53;
t31 = -t265 * t51 + t268 * t52;
t30 = t265 * t52 + t268 * t51;
t29 = -t259 * t47 + t261 * t48;
t28 = -pkin(4) * t157 + pkin(8) * t94 + t43;
t27 = -t265 * t45 + t268 * t46;
t26 = t265 * t46 + t268 * t45;
t25 = -pkin(5) * t86 + t39;
t24 = -pkin(5) * t83 + t38;
t23 = t268 * t43 - t324;
t22 = t265 * t43 + t323;
t18 = -t259 * t32 + t261 * t33;
t17 = -t259 * t30 + t261 * t31;
t16 = -pkin(9) * t68 - t20;
t15 = t21 * t267 + t264 * t59;
t14 = t21 * t264 - t267 * t59;
t13 = -t259 * t26 + t261 * t27;
t12 = -pkin(8) * t53 - t25 * t264 + t267 * t41;
t11 = -pkin(8) * t51 - t24 * t264 + t267 * t40;
t10 = -pkin(4) * t86 + pkin(8) * t54 + t25 * t267 + t264 * t41;
t9 = -pkin(4) * t83 + pkin(8) * t52 + t24 * t267 + t264 * t40;
t8 = -t22 * t259 + t23 * t261;
t7 = -pkin(8) * t45 + t16 * t267 + t68 * t326;
t6 = pkin(8) * t46 + t16 * t264 + t290 * t68;
t5 = -t14 * t265 + t15 * t268;
t4 = t14 * t268 + t15 * t265;
t3 = -pkin(8) * t14 + (-pkin(9) * t267 + t326) * t20;
t2 = pkin(8) * t15 + (-pkin(9) * t264 + t290) * t20;
t1 = -t259 * t4 + t261 * t5;
t19 = [0, 0, 0, 0, 0, qJDD(1), t280, t281, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t262 - t329 * t260) + t285, pkin(1) * (-qJDD(1) * t260 - t329 * t262) - t306, 0, pkin(1) * (t260 * t306 + t262 * t285), t251, 0.2e1 * t259 * t297, 0, t252, 0, 0, -t338 * t261, t338 * t259, pkin(2) * t238 + qJ(3) * t237 + pkin(1) * (t237 * t260 + t238 * t262) + t156, -pkin(2) * t200 + qJ(3) * t156 + pkin(1) * (t156 * t260 - t200 * t262), t259 * (t216 * t268 - t265 * t303) + t261 * (t216 * t265 + t268 * t303), t259 * (-t213 * t268 - t215 * t265) + t261 * (-t213 * t265 + t215 * t268), t259 * (-t220 * t265 + t340) + t261 * (t220 * t268 + t343), t259 * (-t214 * t265 - t268 * t304) + t261 * (t214 * t268 - t265 * t304), t259 * (t219 * t268 - t311) + t261 * (t219 * t265 + t310), (t259 * (t229 * t268 + t230 * t265) + t261 * (t229 * t265 - t230 * t268)) * qJD(4), t259 * (-pkin(7) * t175 + t315) + t261 * (-pkin(3) * t213 + pkin(7) * t176 - t314) - pkin(2) * t213 + qJ(3) * t133 + pkin(1) * (t133 * t260 - t213 * t262), t259 * (-pkin(7) * t179 + t314) + t261 * (-pkin(3) * t215 + pkin(7) * t180 + t315) - pkin(2) * t215 + qJ(3) * t155 + pkin(1) * (t155 * t260 - t215 * t262), t259 * (-pkin(7) * t177 - t101) + t261 * (-pkin(3) * t191 + pkin(7) * t178 + t102) - pkin(2) * t191 + qJ(3) * t140 + pkin(1) * (t140 * t260 - t191 * t262), -pkin(7) * t322 + t261 * (-pkin(3) * t182 + pkin(7) * t102) - pkin(2) * t182 + qJ(3) * t55 + pkin(1) * (-t182 * t262 + t260 * t55), t259 * (-t131 * t265 + t132 * t268) + t261 * (t131 * t268 + t132 * t265), t259 * (-t265 * t91 + t268 * t93) + t261 * (t265 * t93 + t268 * t91), t259 * (-t150 * t265 + t152 * t268) + t261 * (t150 * t268 + t152 * t265), t259 * (-t129 * t265 + t130 * t268) + t261 * (t129 * t268 + t130 * t265), t259 * (-t151 * t265 + t153 * t268) + t261 * (t151 * t268 + t153 * t265), t259 * (-t163 * t265 + t164 * t268) + t261 * (t163 * t268 + t164 * t265), t259 * (-pkin(7) * t88 - t265 * t71 + t268 * t90) + t261 * (-pkin(3) * t134 + pkin(7) * t89 + t265 * t90 + t268 * t71) - pkin(2) * t134 + qJ(3) * t44 + pkin(1) * (-t134 * t262 + t260 * t44), t259 * (-pkin(7) * t104 + t103 * t268 - t265 * t77) + t261 * (-pkin(3) * t138 + pkin(7) * t105 + t103 * t265 + t268 * t77) - pkin(2) * t138 + qJ(3) * t65 + pkin(1) * (-t138 * t262 + t260 * t65), t259 * (-pkin(7) * t47 - t265 * t28 + t268 * t34) + t261 * (-pkin(3) * t157 + pkin(7) * t48 + t265 * t34 + t268 * t28) - pkin(2) * t157 + qJ(3) * t29 + pkin(1) * (-t157 * t262 + t260 * t29), t259 * (-pkin(7) * t22 - pkin(8) * t323 - t265 * t35) + t261 * (-pkin(3) * t144 + pkin(7) * t23 - pkin(8) * t324 + t268 * t35) - pkin(2) * t144 + qJ(3) * t8 + pkin(1) * (-t144 * t262 + t260 * t8), t259 * (-t265 * t79 + t268 * t81) + t261 * (t265 * t81 + t268 * t79), t259 * (-t265 * t49 + t268 * t50) + t261 * (t265 * t50 + t268 * t49), t259 * (-t265 * t61 + t268 * t63) + t261 * (t265 * t63 + t268 * t61), t259 * (-t265 * t78 + t268 * t80) + t261 * (t265 * t80 + t268 * t78), t259 * (-t265 * t62 + t268 * t64) + t261 * (t265 * t64 + t268 * t62), t259 * (t100 * t268 - t265 * t99) + t261 * (t100 * t265 + t268 * t99), t259 * (-pkin(7) * t30 + t11 * t268 - t265 * t9) + t261 * (-pkin(3) * t83 + pkin(7) * t31 + t11 * t265 + t268 * t9) - pkin(2) * t83 + qJ(3) * t17 + pkin(1) * (t17 * t260 - t262 * t83), t259 * (-pkin(7) * t32 - t10 * t265 + t12 * t268) + t261 * (-pkin(3) * t86 + pkin(7) * t33 + t10 * t268 + t12 * t265) - pkin(2) * t86 + qJ(3) * t18 + pkin(1) * (t18 * t260 - t262 * t86), t259 * (-pkin(7) * t26 - t265 * t6 + t268 * t7) + t261 * (-pkin(3) * t68 + pkin(7) * t27 + t265 * t7 + t268 * t6) - pkin(2) * t68 + qJ(3) * t13 + pkin(1) * (t13 * t260 - t262 * t68), t259 * (-pkin(7) * t4 - t2 * t265 + t268 * t3) + t261 * (-pkin(3) * t20 + pkin(7) * t5 + t2 * t268 + t265 * t3) - pkin(2) * t20 + qJ(3) * t1 + pkin(1) * (t1 * t260 - t20 * t262); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183 * t261 + t184 * t259, 0, 0, 0, 0, 0, 0, t175 * t261 + t176 * t259, t179 * t261 + t180 * t259, t177 * t261 + t178 * t259, t101 * t261 + t102 * t259, 0, 0, 0, 0, 0, 0, t259 * t89 + t261 * t88, t104 * t261 + t105 * t259, t259 * t48 + t261 * t47, t22 * t261 + t23 * t259, 0, 0, 0, 0, 0, 0, t259 * t31 + t261 * t30, t259 * t33 + t261 * t32, t259 * t27 + t26 * t261, t259 * t5 + t261 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297, t298, -t238, t200, 0, 0, 0, 0, 0, 0, t213, t215, t191, t182, 0, 0, 0, 0, 0, 0, t134, t138, t157, t144, 0, 0, 0, 0, 0, 0, t83, t86, t68, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t227 - t226, t228, t217, t185, qJDD(4), -t145, -t146, 0, 0, t174, t173, t139, -t174, t279, t254, pkin(4) * t127 - t75, -t307 - t264 * t346 + (-t264 * t333 + t148) * pkin(4), pkin(4) * t92, pkin(4) * t42, t108, t67, t95, t106, t96, t123, pkin(4) * t51 + t295, pkin(4) * t53 + t296, pkin(4) * t45 + t289, pkin(4) * t14 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t173, t139, -t174, t279, t254, -t75, -t76, 0, 0, t108, t67, t95, t106, t96, t123, t295, t296, t289, t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t161, t116, -t162, -t112, t158, -t38, -t39, 0, 0;];
tauJ_reg  = t19;
