% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:56:33
% EndTime: 2019-05-04 20:56:42
% DurationCPUTime: 4.17s
% Computational Cost: add. (22947->377), mult. (42398->576), div. (0->0), fcn. (34252->16), ass. (0->248)
t197 = sin(qJ(6));
t198 = sin(qJ(5));
t202 = cos(qJ(5));
t199 = sin(qJ(4));
t250 = qJD(3) * t199;
t158 = -t202 * qJD(4) + t198 * t250;
t160 = qJD(4) * t198 + t202 * t250;
t201 = cos(qJ(6));
t137 = t201 * t158 + t160 * t197;
t139 = -t158 * t197 + t160 * t201;
t112 = t139 * t137;
t245 = qJD(3) * qJD(4);
t180 = t199 * t245;
t203 = cos(qJ(4));
t243 = t203 * qJDD(3);
t164 = -t180 + t243;
t157 = -qJDD(5) + t164;
t154 = -qJDD(6) + t157;
t281 = -t112 - t154;
t287 = t197 * t281;
t286 = t201 * t281;
t190 = sin(pkin(7));
t193 = cos(pkin(7));
t189 = sin(pkin(13));
t192 = cos(pkin(13));
t191 = sin(pkin(6));
t194 = cos(pkin(6));
t266 = sin(pkin(12));
t267 = cos(pkin(12));
t213 = t266 * g(1) - t267 * g(2);
t251 = -g(3) + qJDD(1);
t211 = t191 * t251 + t194 * t213;
t214 = -t267 * g(1) - t266 * g(2);
t208 = -t189 * t214 + t192 * t211;
t210 = -t191 * t213 + t194 * t251 + qJDD(2);
t285 = t190 * t210 + t193 * t208;
t249 = qJD(3) * t203;
t177 = -qJD(5) + t249;
t171 = -qJD(6) + t177;
t123 = t137 * t171;
t240 = t203 * t245;
t244 = t199 * qJDD(3);
t163 = t240 + t244;
t218 = -t198 * qJDD(4) - t202 * t163;
t132 = -qJD(5) * t158 - t218;
t238 = -t202 * qJDD(4) + t198 * t163;
t216 = qJD(5) * t160 + t238;
t96 = -t137 * qJD(6) + t201 * t132 - t197 * t216;
t284 = t123 + t96;
t261 = t160 * t158;
t212 = -t157 - t261;
t283 = t198 * t212;
t282 = t202 * t212;
t150 = t158 * t177;
t117 = t132 - t150;
t200 = sin(qJ(3));
t204 = cos(qJ(3));
t125 = t189 * t211 + t192 * t214;
t235 = t200 * t125 - t285 * t204;
t92 = t204 * t125 + t285 * t200;
t280 = t200 * t92 - t204 * t235;
t239 = t197 * t132 + t201 * t216;
t74 = (qJD(6) + t171) * t139 + t239;
t113 = (qJD(5) + t177) * t160 + t238;
t135 = t137 ^ 2;
t136 = t139 ^ 2;
t279 = t158 ^ 2;
t156 = t160 ^ 2;
t170 = t171 ^ 2;
t175 = t177 ^ 2;
t278 = qJD(4) ^ 2;
t237 = -pkin(4) * t203 - pkin(10) * t199;
t161 = t237 * qJD(3);
t206 = -t190 * t208 + t193 * t210;
t205 = qJD(3) ^ 2;
t90 = -t205 * pkin(3) + qJDD(3) * pkin(9) + t92;
t61 = t199 * t206 + t203 * t90;
t58 = -t278 * pkin(4) + qJDD(4) * pkin(10) + t161 * t249 + t61;
t233 = -t164 + t180;
t234 = t163 + t240;
t89 = -qJDD(3) * pkin(3) - t205 * pkin(9) + t235;
t64 = t233 * pkin(4) - t234 * pkin(10) + t89;
t42 = t198 * t58 - t202 * t64;
t30 = t212 * pkin(5) - pkin(11) * t117 - t42;
t147 = -pkin(5) * t177 - pkin(11) * t160;
t43 = t198 * t64 + t202 * t58;
t34 = -t279 * pkin(5) - pkin(11) * t216 + t177 * t147 + t43;
t18 = t197 * t34 - t201 * t30;
t19 = t197 * t30 + t201 * t34;
t11 = -t18 * t201 + t19 * t197;
t277 = pkin(5) * t11;
t77 = -t123 + t96;
t49 = -t197 * t74 - t201 * t77;
t276 = pkin(5) * t49;
t275 = t11 * t198;
t274 = t11 * t202;
t110 = t203 * t206;
t57 = -t110 - qJDD(4) * pkin(4) - t278 * pkin(10) + (qJD(3) * t161 + t90) * t199;
t53 = pkin(5) * t216 - t279 * pkin(11) + t160 * t147 + t57;
t273 = t197 * t53;
t272 = t198 * t57;
t270 = t201 * t53;
t269 = t202 * t57;
t100 = -t112 + t154;
t265 = t100 * t197;
t264 = t100 * t201;
t127 = t157 - t261;
t263 = t127 * t198;
t262 = t127 * t202;
t260 = t171 * t197;
t259 = t171 * t201;
t258 = t177 * t198;
t257 = t177 * t202;
t256 = t189 * t191;
t255 = t191 * t192;
t254 = t192 * t193;
t176 = t199 * t205 * t203;
t168 = qJDD(4) + t176;
t253 = t199 * t168;
t169 = qJDD(4) - t176;
t252 = t203 * t169;
t248 = qJD(5) - t177;
t242 = t203 * t112;
t241 = t203 * t261;
t12 = t18 * t197 + t201 * t19;
t22 = t198 * t42 + t202 * t43;
t60 = t199 * t90 - t110;
t40 = t199 * t60 + t203 * t61;
t5 = t12 * t202 - t275;
t3 = t199 * t53 + t203 * t5;
t4 = t12 * t198 + t274;
t236 = t200 * t3 - t204 * t4;
t16 = t199 * t57 + t203 * t22;
t21 = t198 * t43 - t202 * t42;
t232 = t16 * t200 - t204 * t21;
t51 = t197 * t77 - t201 * t74;
t27 = -t198 * t49 + t202 * t51;
t97 = -t135 - t136;
t24 = t199 * t97 + t203 * t27;
t26 = t198 * t51 + t202 * t49;
t231 = t200 * t24 - t204 * t26;
t107 = -t170 - t135;
t67 = t107 * t197 + t286;
t68 = t107 * t201 - t287;
t47 = -t198 * t67 + t202 * t68;
t73 = (qJD(6) - t171) * t139 + t239;
t32 = t199 * t73 + t203 * t47;
t46 = t198 * t68 + t202 * t67;
t230 = t200 * t32 - t204 * t46;
t119 = -t136 - t170;
t81 = t119 * t201 + t265;
t82 = -t119 * t197 + t264;
t55 = -t198 * t81 + t202 * t82;
t38 = t199 * t284 + t203 * t55;
t54 = t198 * t82 + t202 * t81;
t229 = t200 * t38 - t204 * t54;
t228 = t200 * t40 - t204 * t89;
t126 = t156 + t279;
t94 = -t113 * t202 + t117 * t198;
t66 = -t126 * t199 + t203 * t94;
t93 = -t113 * t198 - t117 * t202;
t227 = t200 * t66 - t204 * t93;
t133 = -t175 - t279;
t103 = t133 * t198 + t282;
t104 = t133 * t202 - t283;
t114 = -t248 * t160 - t238;
t80 = t104 * t203 - t114 * t199;
t226 = -t103 * t204 + t200 * t80;
t141 = -t156 - t175;
t108 = t141 * t202 + t263;
t109 = -t141 * t198 + t262;
t118 = t248 * t158 + t218;
t84 = t109 * t203 - t118 * t199;
t225 = -t108 * t204 + t200 * t84;
t186 = t203 ^ 2;
t184 = t186 * t205;
t174 = -t184 - t278;
t145 = t174 * t203 - t253;
t165 = -0.2e1 * t180 + t243;
t224 = t145 * t200 + t165 * t204;
t185 = t199 ^ 2;
t182 = t185 * t205;
t173 = -t182 - t278;
t146 = -t173 * t199 - t252;
t162 = 0.2e1 * t240 + t244;
t223 = t146 * t200 - t162 * t204;
t166 = (t185 + t186) * qJDD(3);
t167 = t182 + t184;
t222 = t166 * t200 + t167 * t204;
t221 = -pkin(3) + t237;
t220 = qJDD(3) * t204 - t200 * t205;
t219 = -qJDD(3) * t200 - t204 * t205;
t217 = pkin(5) * t67 - t18;
t215 = pkin(5) * t81 - t19;
t152 = t220 * t190;
t151 = t219 * t190;
t149 = -t156 + t175;
t148 = -t175 + t279;
t144 = -t169 * t199 + t173 * t203;
t143 = t168 * t203 + t174 * t199;
t142 = t156 - t279;
t134 = t222 * t190;
t122 = -t136 + t170;
t121 = t135 - t170;
t116 = t132 + t150;
t111 = t136 - t135;
t106 = t193 * t144 + t190 * t223;
t105 = t193 * t143 + t190 * t224;
t99 = (t137 * t201 - t139 * t197) * t171;
t98 = (t137 * t197 + t139 * t201) * t171;
t95 = -qJD(6) * t139 - t239;
t88 = t121 * t201 + t265;
t87 = -t122 * t197 + t286;
t86 = t121 * t197 - t264;
t85 = t122 * t201 + t287;
t83 = t109 * t199 + t118 * t203;
t79 = t104 * t199 + t114 * t203;
t72 = t139 * t260 + t201 * t96;
t71 = -t139 * t259 + t197 * t96;
t70 = -t137 * t259 - t197 * t95;
t69 = -t137 * t260 + t201 * t95;
t65 = t126 * t203 + t199 * t94;
t52 = -t197 * t284 - t201 * t73;
t50 = -t197 * t73 + t201 * t284;
t48 = t190 * t280 + t193 * t206;
t45 = t190 * t225 + t193 * t83;
t44 = t190 * t226 + t193 * t79;
t39 = t199 * t61 - t203 * t60;
t37 = t199 * t55 - t203 * t284;
t36 = -pkin(11) * t81 + t270;
t35 = t190 * t227 + t193 * t65;
t33 = -pkin(11) * t67 + t273;
t31 = t199 * t47 - t203 * t73;
t28 = -pkin(5) * t284 + pkin(11) * t82 + t273;
t25 = -pkin(5) * t73 + pkin(11) * t68 - t270;
t23 = t199 * t27 - t203 * t97;
t20 = t190 * t228 + t193 * t39;
t15 = t199 * t22 - t203 * t57;
t14 = t229 * t190 + t193 * t37;
t13 = t230 * t190 + t193 * t31;
t10 = t231 * t190 + t193 * t23;
t9 = -pkin(11) * t49 - t11;
t8 = -pkin(5) * t53 + pkin(11) * t12;
t7 = -pkin(5) * t97 + pkin(11) * t51 + t12;
t6 = t193 * t15 + t232 * t190;
t2 = t199 * t5 - t203 * t53;
t1 = t236 * t190 + t193 * t2;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125 * t256 + t194 * t210 + t208 * t255, 0, 0, 0, 0, 0, 0, t194 * t152 + (t189 * t219 + t220 * t254) * t191, t194 * t151 + (-t189 * t220 + t219 * t254) * t191, 0, (t200 * t235 + t204 * t92) * t256 + (-t190 * t206 + t193 * t280) * t255 + t194 * t48, 0, 0, 0, 0, 0, 0, t194 * t105 + (t189 * (t145 * t204 - t165 * t200) + t192 * (-t190 * t143 + t193 * t224)) * t191, t194 * t106 + (t189 * (t146 * t204 + t162 * t200) + t192 * (-t190 * t144 + t193 * t223)) * t191, t194 * t134 + (t189 * (t166 * t204 - t167 * t200) + t222 * t254) * t191, t194 * t20 + (t189 * (t200 * t89 + t204 * t40) + t192 * (-t190 * t39 + t193 * t228)) * t191, 0, 0, 0, 0, 0, 0, t194 * t44 + (t189 * (t103 * t200 + t204 * t80) + t192 * (-t190 * t79 + t193 * t226)) * t191, t194 * t45 + (t189 * (t108 * t200 + t204 * t84) + t192 * (-t190 * t83 + t193 * t225)) * t191, t194 * t35 + (t189 * (t200 * t93 + t204 * t66) + t192 * (-t190 * t65 + t193 * t227)) * t191, t194 * t6 + (t189 * (t16 * t204 + t200 * t21) + t192 * (-t190 * t15 + t193 * t232)) * t191, 0, 0, 0, 0, 0, 0, t194 * t13 + (t189 * (t200 * t46 + t204 * t32) + t192 * (-t190 * t31 + t193 * t230)) * t191, t194 * t14 + (t189 * (t200 * t54 + t204 * t38) + t192 * (-t190 * t37 + t193 * t229)) * t191, t194 * t10 + (t189 * (t200 * t26 + t204 * t24) + t192 * (-t190 * t23 + t193 * t231)) * t191, t194 * t1 + (t189 * (t200 * t4 + t204 * t3) + t192 * (-t190 * t2 + t193 * t236)) * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, 0, 0, 0, 0, 0, 0, t152, t151, 0, t48, 0, 0, 0, 0, 0, 0, t105, t106, t134, t20, 0, 0, 0, 0, 0, 0, t44, t45, t35, t6, 0, 0, 0, 0, 0, 0, t13, t14, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t235, -t92, 0, 0, t234 * t199, t162 * t203 + t165 * t199, t253 + t203 * (-t182 + t278), -t233 * t203, t199 * (t184 - t278) + t252, 0, pkin(3) * t165 + pkin(9) * t145 - t203 * t89, -pkin(3) * t162 + pkin(9) * t146 + t199 * t89, pkin(3) * t167 + pkin(9) * t166 + t40, -pkin(3) * t89 + pkin(9) * t40, t199 * (t132 * t202 + t160 * t258) - t241, t199 * (t114 * t202 - t116 * t198) - t203 * t142, t199 * (-t149 * t198 + t282) - t203 * t117, t199 * (-t158 * t257 + t198 * t216) + t241, t199 * (t148 * t202 + t263) + t203 * t113, t203 * t157 + t199 * (t158 * t202 - t160 * t198) * t177, t199 * (-pkin(10) * t103 + t272) + t203 * (-pkin(4) * t103 + t42) - pkin(3) * t103 + pkin(9) * t80, t199 * (-pkin(10) * t108 + t269) + t203 * (-pkin(4) * t108 + t43) - pkin(3) * t108 + pkin(9) * t84, pkin(9) * t66 - t199 * t21 + t221 * t93, pkin(9) * t16 + t21 * t221, t199 * (-t198 * t71 + t202 * t72) - t242, t199 * (-t198 * t50 + t202 * t52) - t203 * t111, t199 * (-t198 * t85 + t202 * t87) - t203 * t77, t199 * (-t198 * t69 + t202 * t70) + t242, t199 * (-t198 * t86 + t202 * t88) + t203 * t74, t199 * (-t198 * t98 + t202 * t99) + t203 * t154, t199 * (-pkin(10) * t46 - t198 * t25 + t202 * t33) + t203 * (-pkin(4) * t46 - t217) - pkin(3) * t46 + pkin(9) * t32, t199 * (-pkin(10) * t54 - t198 * t28 + t202 * t36) + t203 * (-pkin(4) * t54 - t215) - pkin(3) * t54 + pkin(9) * t38, t199 * (-pkin(10) * t26 - t198 * t7 + t202 * t9) + t203 * (-pkin(4) * t26 - t276) - pkin(3) * t26 + pkin(9) * t24, t199 * (-pkin(10) * t4 - pkin(11) * t274 - t198 * t8) + t203 * (-pkin(4) * t4 - t277) - pkin(3) * t4 + pkin(9) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t182 - t184, t244, t176, t243, qJDD(4), -t60, -t61, 0, 0, t132 * t198 - t160 * t257, t114 * t198 + t116 * t202, t149 * t202 + t283, -t158 * t258 - t202 * t216, t148 * t198 - t262, (t158 * t198 + t160 * t202) * t177, pkin(4) * t114 + pkin(10) * t104 - t269, pkin(4) * t118 + pkin(10) * t109 + t272, pkin(4) * t126 + pkin(10) * t94 + t22, -pkin(4) * t57 + pkin(10) * t22, t198 * t72 + t202 * t71, t198 * t52 + t202 * t50, t198 * t87 + t202 * t85, t198 * t70 + t202 * t69, t198 * t88 + t202 * t86, t198 * t99 + t202 * t98, -pkin(4) * t73 + pkin(10) * t47 + t198 * t33 + t202 * t25, -pkin(4) * t284 + pkin(10) * t55 + t198 * t36 + t202 * t28, -pkin(4) * t97 + pkin(10) * t27 + t198 * t9 + t202 * t7, -pkin(4) * t53 + pkin(10) * t5 - pkin(11) * t275 + t202 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t142, t117, -t261, -t113, -t157, -t42, -t43, 0, 0, t112, t111, t77, -t112, -t74, -t154, t217, t215, t276, t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, t77, -t112, -t74, -t154, -t18, -t19, 0, 0;];
tauJ_reg  = t17;
