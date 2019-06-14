% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:44:55
% EndTime: 2019-05-05 14:45:04
% DurationCPUTime: 3.67s
% Computational Cost: add. (13637->360), mult. (30798->498), div. (0->0), fcn. (22054->10), ass. (0->228)
t204 = sin(pkin(10));
t199 = t204 ^ 2;
t206 = cos(pkin(10));
t200 = t206 ^ 2;
t289 = qJD(1) ^ 2;
t187 = (t199 + t200) * t289;
t210 = sin(qJ(4));
t212 = cos(qJ(4));
t226 = t204 * t212 + t206 * t210;
t181 = t226 * qJD(1);
t209 = sin(qJ(5));
t211 = cos(qJ(5));
t166 = -qJD(4) * t211 + t181 * t209;
t168 = qJD(4) * t209 + t181 * t211;
t139 = t168 * t166;
t243 = t206 * qJDD(1);
t244 = t204 * qJDD(1);
t229 = t210 * t244 - t212 * t243;
t247 = t181 * qJD(4);
t156 = -t229 - t247;
t149 = qJDD(5) - t156;
t293 = -t139 + t149;
t303 = pkin(5) * t293;
t178 = t226 * qJDD(1);
t251 = qJD(1) * t206;
t179 = qJD(1) * t204 * t210 - t212 * t251;
t248 = t179 * qJD(4);
t158 = t178 - t248;
t125 = -qJD(5) * t166 + qJDD(4) * t209 + t158 * t211;
t175 = qJD(5) + t179;
t144 = t175 * t166;
t107 = t125 + t144;
t302 = qJ(6) * t107;
t159 = t181 * t179;
t291 = qJDD(4) - t159;
t301 = t210 * t291;
t300 = t212 * t291;
t263 = t293 * t209;
t262 = t293 * t211;
t165 = t168 ^ 2;
t174 = t175 ^ 2;
t130 = -t165 - t174;
t164 = t166 ^ 2;
t231 = -qJDD(4) * t211 + t209 * t158;
t124 = -qJD(5) * t168 - t231;
t140 = pkin(5) * t175 - qJ(6) * t168;
t150 = pkin(4) * t179 - pkin(8) * t181;
t288 = qJD(4) ^ 2;
t205 = sin(pkin(9));
t276 = sin(qJ(1));
t277 = cos(qJ(1));
t222 = g(1) * t276 - g(2) * t277;
t220 = qJDD(1) * pkin(1) + t222;
t223 = g(1) * t277 + g(2) * t276;
t185 = -pkin(1) * t289 - t223;
t207 = cos(pkin(9));
t256 = t207 * t185;
t217 = qJDD(1) * qJ(3) + t205 * t220 + t256;
t253 = -g(3) + qJDD(2);
t287 = 2 * qJD(3);
t134 = t206 * (-pkin(2) * t289 + t217) + t204 * t253 + t251 * t287;
t246 = t200 * t289;
t127 = -pkin(3) * t246 + pkin(7) * t243 + t134;
t218 = -t205 * t222 - t256;
t232 = t206 * t253;
t233 = pkin(1) * t205 + qJ(3);
t296 = t233 + pkin(7);
t215 = t232 + (-t296 * qJDD(1) + (-(2 * qJD(3)) + (pkin(3) * t206 + pkin(2)) * qJD(1)) * qJD(1) + t218) * t204;
t90 = t212 * t127 + t210 * t215;
t73 = -pkin(4) * t288 + qJDD(4) * pkin(8) - t179 * t150 + t90;
t182 = t207 * t220;
t230 = -t205 * t185 + t182;
t147 = -qJDD(1) * pkin(2) - qJ(3) * t289 + qJDD(3) - t230;
t132 = -pkin(3) * t243 + t147 + (-t199 * t289 - t246) * pkin(7);
t77 = (-t158 + t248) * pkin(8) + (-t156 + t247) * pkin(4) + t132;
t46 = t209 * t77 + t211 * t73;
t224 = qJ(6) * t124 - 0.2e1 * qJD(6) * t166 - t175 * t140 + t46;
t299 = -t224 + (t130 + t164) * pkin(5);
t237 = pkin(1) * t207 + pkin(2);
t298 = -t237 * qJDD(1) + t187 * t233 + t147;
t294 = t125 - t144;
t104 = (qJD(5) - t175) * t168 + t231;
t176 = t179 ^ 2;
t177 = t181 ^ 2;
t126 = -t174 - t164;
t81 = t126 * t209 + t262;
t286 = pkin(4) * t81;
t112 = t139 + t149;
t265 = t112 * t209;
t86 = t130 * t211 - t265;
t285 = pkin(4) * t86;
t249 = qJD(6) * t168;
t161 = -0.2e1 * t249;
t45 = t209 * t73 - t211 * t77;
t221 = -t302 - t45 + t303;
t28 = t161 + t221;
t284 = pkin(5) * t28;
t119 = -t164 - t165;
t70 = -t104 * t211 + t107 * t209;
t54 = -t119 * t212 + t210 * t70;
t283 = pkin(7) * t54;
t103 = (qJD(5) + t175) * t168 + t231;
t82 = t126 * t211 - t263;
t59 = -t103 * t212 + t210 * t82;
t282 = pkin(7) * t59;
t264 = t112 * t211;
t87 = -t130 * t209 - t264;
t63 = t210 * t87 - t212 * t294;
t281 = pkin(7) * t63;
t68 = -t104 * t209 - t107 * t211;
t280 = pkin(8) * t68;
t279 = pkin(8) * t81;
t278 = pkin(8) * t86;
t275 = pkin(4) * t210;
t274 = pkin(5) * t107;
t89 = t127 * t210 - t212 * t215;
t56 = t210 * t90 - t212 * t89;
t273 = t204 * t56;
t272 = t209 * t28;
t72 = -qJDD(4) * pkin(4) - pkin(8) * t288 + t150 * t181 + t89;
t271 = t209 * t72;
t270 = t211 * t28;
t269 = t211 * t72;
t268 = -pkin(4) * t119 + pkin(8) * t70;
t267 = -pkin(4) * t103 + pkin(8) * t82;
t266 = -pkin(4) * t294 + pkin(8) * t87;
t261 = t132 * t210;
t260 = t132 * t212;
t153 = qJDD(4) + t159;
t259 = t153 * t212;
t258 = t175 * t209;
t257 = t175 * t211;
t255 = t210 * t153;
t55 = t119 * t210 + t212 * t70;
t27 = -t204 * t54 + t206 * t55;
t242 = pkin(1) * (t205 * t27 - t207 * t68) + qJ(3) * t27 - pkin(2) * t68;
t60 = t103 * t210 + t212 * t82;
t33 = -t204 * t59 + t206 * t60;
t241 = pkin(1) * (t205 * t33 - t207 * t81) + qJ(3) * t33 - pkin(2) * t81;
t64 = t210 * t294 + t212 * t87;
t38 = -t204 * t63 + t206 * t64;
t240 = pkin(1) * (t205 * t38 - t207 * t86) + qJ(3) * t38 - pkin(2) * t86;
t239 = t212 * t139;
t238 = t210 * t139;
t236 = -pkin(4) * t212 - pkin(3);
t235 = -pkin(3) * t81 + pkin(7) * t60;
t234 = -pkin(3) * t86 + pkin(7) * t64;
t20 = t209 * t45 + t211 * t46;
t57 = t210 * t89 + t212 * t90;
t133 = -t232 + ((-pkin(2) * qJD(1) + t287) * qJD(1) + t217) * t204;
t96 = t204 * t133 + t134 * t206;
t19 = t209 * t46 - t211 * t45;
t219 = t221 + t303;
t47 = -pkin(5) * t124 - qJ(6) * t164 + t140 * t168 + qJDD(6) + t72;
t196 = t200 * qJDD(1);
t195 = t199 * qJDD(1);
t186 = t196 + t195;
t171 = -t177 - t288;
t170 = -t177 + t288;
t169 = t176 - t288;
t162 = 0.2e1 * t249;
t157 = t178 - 0.2e1 * t248;
t155 = t229 + 0.2e1 * t247;
t151 = -t288 - t176;
t142 = -t165 + t174;
t141 = t164 - t174;
t138 = -t176 - t177;
t136 = t165 - t164;
t129 = -t171 * t210 - t259;
t128 = t171 * t212 - t255;
t118 = t178 * t210 - t212 * t229;
t117 = -t178 * t212 - t210 * t229;
t115 = t151 * t212 - t301;
t114 = t151 * t210 + t300;
t110 = (-t166 * t211 + t168 * t209) * t175;
t109 = (-t166 * t209 - t168 * t211) * t175;
t100 = t125 * t211 - t168 * t258;
t99 = t125 * t209 + t168 * t257;
t98 = -t124 * t209 + t166 * t257;
t97 = t124 * t211 + t166 * t258;
t95 = -t128 * t204 + t129 * t206;
t94 = t141 * t211 - t265;
t93 = -t142 * t209 + t262;
t92 = t141 * t209 + t264;
t91 = t142 * t211 + t263;
t83 = -t117 * t204 + t118 * t206;
t80 = -t114 * t204 + t115 * t206;
t74 = -pkin(5) * t294 - qJ(6) * t112;
t69 = -t103 * t211 - t209 * t294;
t67 = -t103 * t209 + t211 * t294;
t61 = t204 * (t110 * t212 + t149 * t210) + t206 * (t110 * t210 - t149 * t212);
t53 = pkin(7) * t55;
t52 = -pkin(4) * t68 + t274;
t51 = t204 * (t100 * t212 + t238) + t206 * (t100 * t210 - t239);
t50 = t204 * (t212 * t98 - t238) + t206 * (t210 * t98 + t239);
t49 = t269 - t278;
t48 = t271 - t279;
t43 = -qJ(6) * t130 + t47;
t42 = t204 * (-t104 * t210 + t212 * t94) + t206 * (t104 * t212 + t210 * t94);
t41 = t204 * (t107 * t210 + t212 * t93) + t206 * (-t107 * t212 + t210 * t93);
t40 = t46 - t285;
t39 = t45 - t286;
t37 = t204 * t64 + t206 * t63;
t35 = -pkin(5) * t103 + qJ(6) * t126 - t47;
t34 = -pkin(5) * t164 + t224;
t32 = t204 * t60 + t206 * t59;
t30 = t204 * (t136 * t210 + t212 * t69) + t206 * (-t136 * t212 + t210 * t69);
t29 = t206 * t57 - t273;
t26 = t204 * t55 + t206 * t54;
t23 = t162 - t221 + t302;
t21 = -qJ(6) * t104 + (-t119 - t164) * pkin(5) + t224;
t18 = -t285 - t299;
t17 = -t209 * t74 + t211 * t43 - t278;
t16 = -qJ(6) * t262 - t209 * t35 - t279;
t15 = t162 - t219 - t286;
t13 = -pkin(5) * t47 + qJ(6) * t34;
t12 = t20 * t212 + t210 * t72;
t11 = t20 * t210 - t212 * t72;
t10 = -t19 - t280;
t9 = t211 * t34 - t272;
t8 = t209 * t34 + t270;
t7 = t210 * t47 + t212 * t9;
t6 = t210 * t9 - t212 * t47;
t5 = -t209 * t21 + t211 * t23 - t280;
t4 = -pkin(4) * t8 - t284;
t2 = -pkin(8) * t8 - qJ(6) * t270 - t13 * t209;
t1 = -t204 * t6 + t206 * t7;
t3 = [0, 0, 0, 0, 0, qJDD(1), t222, t223, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t207 - t205 * t289) + t230, (-0.2e1 * qJDD(1) * t205 - t207 * t289) * pkin(1) + t218, 0, pkin(1) * (t205 ^ 2 * t220 + t207 * t182), t195, 0.2e1 * t204 * t243, 0, t196, 0, 0, -t298 * t206, t298 * t204, pkin(2) * t187 + qJ(3) * t186 + pkin(1) * (t186 * t205 + t187 * t207) + t96, -pkin(2) * t147 + qJ(3) * t96 + pkin(1) * (-t147 * t207 + t205 * t96), t204 * (t158 * t212 - t210 * t247) + t206 * (t158 * t210 + t212 * t247), t204 * (-t155 * t212 - t157 * t210) + t206 * (-t155 * t210 + t157 * t212), t204 * (-t170 * t210 + t300) + t206 * (t170 * t212 + t301), t204 * (-t156 * t210 + t212 * t248) + t206 * (t156 * t212 + t210 * t248), t204 * (t169 * t212 - t255) + t206 * (t169 * t210 + t259), (t204 * (-t179 * t212 + t181 * t210) + t206 * (-t179 * t210 - t181 * t212)) * qJD(4), t204 * (-pkin(7) * t114 + t261) + t206 * (-pkin(3) * t155 + pkin(7) * t115 - t260) - pkin(2) * t155 + qJ(3) * t80 + pkin(1) * (-t155 * t207 + t205 * t80), t204 * (-pkin(7) * t128 + t260) + t206 * (-pkin(3) * t157 + pkin(7) * t129 + t261) - pkin(2) * t157 + qJ(3) * t95 + pkin(1) * (-t157 * t207 + t205 * t95), t204 * (-pkin(7) * t117 - t56) + t206 * (-pkin(3) * t138 + pkin(7) * t118 + t57) - pkin(2) * t138 + qJ(3) * t83 + pkin(1) * (-t138 * t207 + t205 * t83), -pkin(7) * t273 + t206 * (-pkin(3) * t132 + pkin(7) * t57) - pkin(2) * t132 + qJ(3) * t29 + pkin(1) * (-t132 * t207 + t205 * t29), t51, t30, t41, t50, t42, t61, t204 * (-t210 * t39 + t212 * t48 - t282) + t206 * (t210 * t48 + t212 * t39 + t235) + t241, t204 * (-t210 * t40 + t212 * t49 - t281) + t206 * (t210 * t49 + t212 * t40 + t234) + t240, t204 * (t10 * t212 + t275 * t68 - t283) + t206 * (t210 * t10 + t236 * t68 + t53) + t242, (t204 * (-pkin(8) * t212 + t275) + t206 * (-pkin(8) * t210 + t236) - t237) * t19 + t296 * (-t11 * t204 + t12 * t206), t51, t30, t41, t50, t42, t61, t204 * (-t15 * t210 + t16 * t212 - t282) + t206 * (t15 * t212 + t16 * t210 + t235) + t241, t204 * (t17 * t212 - t18 * t210 - t281) + t206 * (t17 * t210 + t18 * t212 + t234) + t240, t204 * (-t210 * t52 + t212 * t5 - t283) + t206 * (-pkin(3) * t68 + t210 * t5 + t212 * t52 + t53) + t242, t204 * (-pkin(7) * t6 + t2 * t212 - t210 * t4) + t206 * (-pkin(3) * t8 + pkin(7) * t7 + t2 * t210 + t212 * t4) - pkin(2) * t8 + qJ(3) * t1 + pkin(1) * (t1 * t205 - t207 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * t206 + t134 * t204, 0, 0, 0, 0, 0, 0, t114 * t206 + t115 * t204, t128 * t206 + t129 * t204, t117 * t206 + t118 * t204, t204 * t57 + t206 * t56, 0, 0, 0, 0, 0, 0, t32, t37, t26, t11 * t206 + t12 * t204, 0, 0, 0, 0, 0, 0, t32, t37, t26, t204 * t7 + t206 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t244, -t187, t147, 0, 0, 0, 0, 0, 0, t155, t157, t138, t132, 0, 0, 0, 0, 0, 0, t81, t86, t68, t19, 0, 0, 0, 0, 0, 0, t81, t86, t68, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t177 - t176, t178, -t159, -t229, qJDD(4), -t89, -t90, 0, 0, t99, t67, t91, t97, t92, t109, t267 - t269, t266 + t271, t20 + t268, -pkin(4) * t72 + pkin(8) * t20, t99, t67, t91, t97, t92, t109, -qJ(6) * t263 + t211 * t35 + t267, t209 * t43 + t211 * t74 + t266, t209 * t23 + t21 * t211 + t268, -pkin(4) * t47 + pkin(8) * t9 - qJ(6) * t272 + t13 * t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t136, t107, -t139, -t104, t149, -t45, -t46, 0, 0, t139, t136, t107, -t139, -t104, t149, t161 + t219, t299, -t274, t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t294, t119, t47;];
tauJ_reg  = t3;
