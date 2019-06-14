% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:35:28
% EndTime: 2019-05-06 08:35:38
% DurationCPUTime: 3.49s
% Computational Cost: add. (14021->371), mult. (31107->471), div. (0->0), fcn. (17412->8), ass. (0->228)
t307 = 2 * qJD(5);
t185 = sin(pkin(9));
t186 = cos(pkin(9));
t192 = cos(qJ(2));
t250 = qJD(1) * t192;
t133 = -qJD(2) * t186 + t185 * t250;
t134 = qJD(2) * t185 + t186 * t250;
t104 = t133 * t134;
t243 = qJD(1) * qJD(2);
t170 = t192 * t243;
t189 = sin(qJ(2));
t241 = t189 * qJDD(1);
t145 = t170 + t241;
t296 = -t104 + t145;
t306 = t185 * t296;
t305 = t186 * t296;
t188 = sin(qJ(6));
t136 = qJDD(6) + t145;
t191 = cos(qJ(6));
t101 = t133 * t188 - t134 * t191;
t99 = -t133 * t191 - t134 * t188;
t78 = t101 * t99;
t299 = t136 - t78;
t304 = t188 * t299;
t303 = t191 * t299;
t195 = qJD(1) ^ 2;
t297 = t192 * t195;
t235 = t189 * t297;
t154 = qJDD(2) + t235;
t144 = 0.2e1 * t170 + t241;
t194 = qJD(2) ^ 2;
t181 = t189 ^ 2;
t262 = t181 * t195;
t158 = t194 + t262;
t155 = -qJDD(2) + t235;
t257 = t192 * t155;
t302 = pkin(1) * t144 - pkin(7) * (t158 * t189 + t257);
t182 = t192 ^ 2;
t261 = t182 * t195;
t159 = t194 + t261;
t258 = t189 * t154;
t301 = pkin(7) * (t159 * t192 + t258);
t300 = pkin(7) - qJ(4);
t298 = qJDD(1) * pkin(7);
t169 = t189 * t243;
t240 = t192 * qJDD(1);
t146 = -t169 + t240;
t118 = -qJDD(2) * t185 - t146 * t186;
t251 = qJD(1) * t189;
t233 = t133 * t251;
t94 = t118 - t233;
t252 = t181 + t182;
t151 = t252 * t195;
t255 = pkin(1) * t151 + t252 * t298;
t295 = pkin(2) * t158 - qJ(3) * t155;
t294 = -pkin(3) * t261 - t146 * qJ(4);
t190 = sin(qJ(1));
t193 = cos(qJ(1));
t219 = g(1) * t193 + g(2) * t190;
t127 = -pkin(1) * t195 - t219 + t298;
t260 = t189 * qJ(3);
t218 = -pkin(2) * t192 - t260;
t142 = t218 * qJD(1);
t221 = qJD(1) * t142 + t127;
t293 = qJ(4) * qJDD(1) + 0.2e1 * qJD(1) * qJD(4) - t221;
t253 = g(1) * t190 - g(2) * t193;
t126 = qJDD(1) * pkin(1) + t195 * pkin(7) + t253;
t205 = -pkin(2) * t169 + t126;
t202 = pkin(3) * t146 - qJ(4) * t261 + qJDD(4) + t205;
t153 = -qJD(2) * pkin(3) - qJ(4) * t251;
t289 = 2 * qJD(3);
t246 = t289 + t153;
t222 = t246 * t189;
t282 = qJ(3) + pkin(4);
t283 = pkin(2) + qJ(5);
t43 = t283 * t146 + t282 * t145 + (t222 + (-t189 * qJ(5) + t192 * t282) * qJD(2)) * qJD(1) + t202;
t284 = t192 * g(3);
t113 = t189 * t127 + t284;
t228 = qJDD(2) * pkin(2) + qJ(3) * t194 - qJDD(3);
t210 = pkin(3) * t154 + t145 * qJ(4) + t228;
t201 = t113 - t210;
t217 = pkin(4) * t189 + qJ(5) * t192;
t249 = qJD(2) * t192;
t234 = qJ(4) * t249;
t245 = -0.2e1 * qJD(4) + t142;
t60 = -t194 * pkin(4) - qJDD(2) * qJ(5) + (t234 + (-qJD(1) * t217 + t245) * t189) * qJD(1) + t201;
t220 = t134 * t307 - t185 * t60 + t186 * t43;
t28 = t133 * t307 + t185 * t43 + t186 * t60;
t292 = t257 - (-t194 + t261) * t189;
t84 = t189 * t221 - t228 + t284;
t288 = pkin(2) + pkin(3);
t290 = t192 * t288 + pkin(1) + t260;
t97 = t99 ^ 2;
t98 = t101 ^ 2;
t130 = t133 ^ 2;
t131 = t134 ^ 2;
t163 = qJD(6) + t251;
t161 = t163 ^ 2;
t287 = t146 * pkin(2);
t197 = pkin(5) * t296 - pkin(8) * t94 + t220;
t117 = -qJDD(2) * t186 + t146 * t185;
t119 = pkin(5) * t251 + pkin(8) * t134;
t22 = -pkin(5) * t130 + pkin(8) * t117 - t119 * t251 + t28;
t10 = t188 * t22 - t191 * t197;
t11 = t188 * t197 + t191 * t22;
t7 = -t10 * t191 + t11 * t188;
t286 = t185 * t7;
t285 = t186 * t7;
t176 = t189 * g(3);
t254 = -pkin(2) * t194 - t176;
t196 = t246 * qJD(2) + (qJD(1) * t245 + t127) * t192 + t254 + t294;
t59 = -qJ(5) * t194 + qJDD(2) * t282 - t217 * t297 + qJDD(5) + t196;
t280 = t185 * t59;
t96 = t104 + t145;
t279 = t185 * t96;
t278 = t186 * t59;
t277 = t186 * t96;
t34 = -pkin(5) * t117 - pkin(8) * t130 - t134 * t119 + t59;
t276 = t188 * t34;
t72 = t136 + t78;
t275 = t188 * t72;
t232 = t134 * t251;
t92 = t117 - t232;
t67 = t185 * t94 + t186 * t92;
t274 = t189 * t67;
t273 = t191 * t34;
t272 = t191 * t72;
t271 = qJ(3) * t151;
t270 = qJ(3) * t159;
t269 = qJ(3) * t192;
t268 = t144 * t192;
t264 = t163 * t188;
t263 = t163 * t191;
t147 = -0.2e1 * t169 + t240;
t259 = t189 * t147;
t256 = pkin(1) * t147 - t301;
t244 = qJD(6) + t163;
t242 = qJDD(2) * qJ(3);
t239 = pkin(3) + t283;
t238 = t189 * t78;
t236 = t189 * t104;
t231 = t185 * t251;
t230 = t186 * t251;
t8 = t10 * t188 + t11 * t191;
t229 = qJD(3) * t251;
t114 = t192 * t127 - t176;
t225 = t113 * t189 + t114 * t192;
t224 = -t117 * t191 + t188 * t118;
t223 = -t131 - t262;
t216 = t145 + t170;
t12 = t185 * t28 + t186 * t220;
t13 = -t185 * t220 + t186 * t28;
t213 = t188 * t117 + t191 * t118;
t103 = t259 + t268;
t70 = -qJD(6) * t99 + t213;
t207 = qJD(2) * t289 + t242 + t254;
t206 = (-qJD(6) + t163) * t101 - t224;
t204 = t189 * t282 + t192 * t239 + pkin(1);
t203 = (t144 + t216) * t260 + t302;
t200 = t202 + t287;
t82 = t192 * t221 + t207;
t199 = t205 + 0.2e1 * t229 + t287;
t198 = t145 * qJ(3) + t200;
t65 = t196 + t242;
t68 = (t189 * t245 + t234) * qJD(1) + t201;
t152 = (t181 - t182) * t195;
t128 = t189 * t145;
t121 = -t131 + t262;
t120 = t130 - t262;
t111 = t170 * t189 + t128;
t110 = t258 + t192 * (t194 - t262);
t109 = (t146 - t169) * t192;
t102 = -t262 - t130;
t93 = t118 + t233;
t91 = -t117 - t232;
t90 = t163 * t99;
t89 = -t130 - t131;
t86 = -t98 + t161;
t85 = t97 - t161;
t83 = -t98 - t161;
t80 = -t185 * t223 - t277;
t79 = t186 * t223 - t279;
t77 = t98 - t97;
t76 = -t161 - t97;
t75 = t102 * t186 - t306;
t74 = t102 * t185 + t305;
t69 = -qJD(6) * t101 - t224;
t66 = t185 * t92 - t186 * t94;
t64 = (t101 * t188 - t191 * t99) * t163;
t63 = (-t101 * t191 - t188 * t99) * t163;
t62 = (qJ(3) * t249 + t222) * qJD(1) + t198;
t61 = -t97 - t98;
t57 = t70 + t90;
t56 = t70 - t90;
t55 = -t244 * t99 + t213;
t52 = t101 * t244 + t224;
t51 = t191 * t85 - t275;
t50 = -t188 * t86 + t303;
t49 = t188 * t85 + t272;
t48 = t191 * t86 + t304;
t47 = -t101 * t264 + t191 * t70;
t46 = t101 * t263 + t188 * t70;
t45 = -t188 * t69 + t263 * t99;
t44 = t191 * t69 + t264 * t99;
t42 = -t188 * t83 - t272;
t41 = t191 * t83 - t275;
t37 = t191 * t76 - t304;
t36 = t188 * t76 + t303;
t33 = t188 * t57 + t191 * t206;
t32 = -t188 * t56 - t191 * t52;
t31 = t188 * t206 - t191 * t57;
t30 = -t188 * t52 + t191 * t56;
t25 = -t185 * t41 + t186 * t42;
t24 = t185 * t42 + t186 * t41;
t23 = -pkin(8) * t41 + t273;
t21 = -t185 * t36 + t186 * t37;
t20 = t185 * t37 + t186 * t36;
t19 = -pkin(8) * t36 + t276;
t17 = -pkin(5) * t55 + pkin(8) * t42 + t276;
t16 = -pkin(5) * t52 + pkin(8) * t37 - t273;
t15 = -t185 * t31 + t186 * t33;
t14 = t185 * t33 + t186 * t31;
t5 = -pkin(5) * t34 + pkin(8) * t8;
t4 = -pkin(8) * t31 - t7;
t3 = -pkin(5) * t61 + pkin(8) * t33 + t8;
t2 = t186 * t8 - t286;
t1 = t185 * t8 + t285;
t6 = [0, 0, 0, 0, 0, qJDD(1), t253, t219, 0, 0, t111, t103, t110, t109, -t292, 0, t126 * t192 + t256, -t189 * t126 - t302, t225 + t255, pkin(1) * t126 + pkin(7) * t225, t111, t110, -t103, 0, t292, t109, t192 * (pkin(2) * t147 + t199) + (t192 * t216 + t259) * qJ(3) + t256, (pkin(2) * t151 + t82) * t192 + (t84 + t271) * t189 + t255, pkin(2) * t268 + t189 * t199 + t203, pkin(7) * (t189 * t84 + t192 * t82) + (pkin(1) - t218) * (qJ(3) * t216 + t199), t109, t103, -t292, t111, t110, 0, t192 * (qJ(4) * t155 + t144 * t288) + (qJ(4) * t158 + qJD(1) * t222 + t200) * t189 + t203, -qJ(4) * t258 + t192 * (-qJ(3) * t170 - qJ(4) * t159 - t153 * t251 - t198 - 0.2e1 * t229) + t301 - t290 * t147, (-qJD(2) * t153 - t288 * t151 + t293 * t192 - t207 - t294) * t192 + (t210 - t271 + (-qJ(4) * t243 - g(3)) * t192 + t293 * t189) * t189 - t255, t290 * t62 + t300 * (t189 * t68 + t192 * t65), t236 + t192 * (-t118 * t186 - t134 * t231), t189 * (t131 - t130) + t192 * (t185 * t93 + t186 * t91), t189 * t94 + t192 * (t121 * t185 - t305), -t236 + t192 * (t117 * t185 + t133 * t230), t189 * t92 + t192 * (-t120 * t186 + t279), t128 + (-t133 * t186 + t134 * t185) * t189 * t250, t189 * (-qJ(4) * t75 + t220) + t192 * (-qJ(4) * t91 - t280) + pkin(7) * (t189 * t75 + t192 * t91) + t204 * t74, t189 * (-qJ(4) * t80 - t28) + t192 * (-qJ(4) * t93 - t278) + pkin(7) * (t189 * t80 + t192 * t93) + t204 * t79, -qJ(4) * t274 + t192 * (-qJ(4) * t89 + t12) + pkin(7) * (t192 * t89 + t274) + t204 * t66, t12 * t204 + t300 * (t13 * t189 + t192 * t59), t238 + t192 * (t185 * t46 - t186 * t47), t189 * t77 + t192 * (t185 * t30 - t186 * t32), t189 * t57 + t192 * (t185 * t48 - t186 * t50), -t238 + t192 * (t185 * t44 - t186 * t45), t189 * t206 + t192 * (t185 * t49 - t186 * t51), t189 * t136 + t192 * (t185 * t63 - t186 * t64), t189 * (pkin(5) * t36 - qJ(4) * t21 - t10) + t192 * (-qJ(4) * t52 + t185 * t16 - t186 * t19) + pkin(7) * (t189 * t21 + t192 * t52) + t204 * t20, t189 * (pkin(5) * t41 - qJ(4) * t25 - t11) + t192 * (-qJ(4) * t55 + t185 * t17 - t186 * t23) + pkin(7) * (t189 * t25 + t192 * t55) + t204 * t24, t189 * (pkin(5) * t31 - qJ(4) * t15) + t192 * (-qJ(4) * t61 + t185 * t3 - t186 * t4) + pkin(7) * (t15 * t189 + t192 * t61) + t204 * t14, t189 * (pkin(5) * t7 - qJ(4) * t2) + t192 * (pkin(8) * t285 - qJ(4) * t34 + t185 * t5) + pkin(7) * (t189 * t2 + t192 * t34) + t204 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t152, t241, t235, t240, qJDD(2), -t113, -t114, 0, 0, -t235, t241, -t152, qJDD(2), -t240, t235, pkin(2) * t154 - t270 - t84, (-pkin(2) * t189 + t269) * qJDD(1), t82 + t295, -pkin(2) * t84 + qJ(3) * t82, t235, t152, t240, -t235, t241, qJDD(2), pkin(3) * t158 + t295 + t65, -t154 * t288 + t270 + t68, (t189 * t288 - t269) * qJDD(1), qJ(3) * t65 - t288 * t68, -t118 * t185 + t134 * t230, t185 * t91 - t186 * t93, -t121 * t186 - t306, -t117 * t186 + t133 * t231, -t120 * t185 - t277, (-t133 * t185 - t134 * t186) * t251, -t239 * t75 + t282 * t91 + t278, -t239 * t80 + t282 * t93 - t280, -t239 * t67 + t282 * t89 - t13, -t13 * t239 + t282 * t59, -t185 * t47 - t186 * t46, -t185 * t32 - t186 * t30, -t185 * t50 - t186 * t48, -t185 * t45 - t186 * t44, -t185 * t51 - t186 * t49, -t185 * t64 - t186 * t63, -t186 * t16 - t185 * t19 - t21 * t239 + t282 * t52, -t186 * t17 - t185 * t23 - t239 * t25 + t282 * t55, -t15 * t239 - t185 * t4 - t186 * t3 + t282 * t61, pkin(8) * t286 - t186 * t5 - t2 * t239 + t282 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t241, -t158, t84, 0, 0, 0, 0, 0, 0, -t158, t154, -t241, t68, 0, 0, 0, 0, 0, 0, t75, t80, t67, t13, 0, 0, 0, 0, 0, 0, t21, t25, t15, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t147, -t151, t62, 0, 0, 0, 0, 0, 0, t74, t79, t66, t12, 0, 0, 0, 0, 0, 0, t20, t24, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t93, t89, t59, 0, 0, 0, 0, 0, 0, t52, t55, t61, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t77, t57, -t78, t206, t136, -t10, -t11, 0, 0;];
tauJ_reg  = t6;
