% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:23:40
% EndTime: 2019-05-04 22:23:50
% DurationCPUTime: 3.25s
% Computational Cost: add. (15239->355), mult. (30203->537), div. (0->0), fcn. (22189->14), ass. (0->218)
t194 = sin(pkin(6));
t232 = -g(3) + qJDD(1);
t197 = cos(pkin(6));
t243 = sin(pkin(10));
t244 = cos(pkin(10));
t211 = t243 * g(1) - t244 * g(2);
t255 = t197 * t211;
t259 = t194 * t232 + t255;
t205 = qJD(2) ^ 2;
t192 = sin(pkin(12));
t195 = cos(pkin(12));
t200 = sin(qJ(4));
t231 = qJD(2) * t200;
t158 = -t195 * qJD(4) + t192 * t231;
t160 = qJD(4) * t192 + t195 * t231;
t199 = sin(qJ(6));
t202 = cos(qJ(6));
t132 = t202 * t158 + t160 * t199;
t203 = cos(qJ(4));
t230 = qJD(2) * t203;
t179 = -qJD(6) + t230;
t120 = t132 * t179;
t228 = qJD(2) * qJD(4);
t223 = t203 * t228;
t227 = t200 * qJDD(2);
t164 = t223 + t227;
t144 = qJDD(4) * t192 + t164 * t195;
t218 = -t195 * qJDD(4) + t164 * t192;
t99 = -t132 * qJD(6) + t202 * t144 - t199 * t218;
t258 = t120 + t99;
t182 = t200 * t228;
t226 = t203 * qJDD(2);
t165 = -t182 + t226;
t237 = t158 * t160;
t213 = -t165 - t237;
t257 = t192 * t213;
t256 = t195 * t213;
t161 = -qJDD(6) + t165;
t134 = -t158 * t199 + t160 * t202;
t238 = t134 * t132;
t210 = -t161 - t238;
t254 = t199 * t210;
t253 = t202 * t210;
t149 = t158 * t230;
t124 = -t144 + t149;
t150 = t160 * t230;
t122 = -t218 - t150;
t219 = t199 * t144 + t202 * t218;
t80 = (qJD(6) + t179) * t134 + t219;
t130 = t132 ^ 2;
t131 = t134 ^ 2;
t252 = t158 ^ 2;
t157 = t160 ^ 2;
t177 = t179 ^ 2;
t251 = qJD(4) ^ 2;
t215 = t197 * t232 + qJDD(3);
t207 = -t194 * t211 + t215;
t206 = t200 * t207;
t217 = -pkin(4) * t203 - qJ(5) * t200;
t171 = -t244 * g(1) - t243 * g(2);
t201 = sin(qJ(2));
t204 = cos(qJ(2));
t129 = t204 * t171 + t201 * t259;
t118 = -t205 * pkin(2) + t129;
t193 = sin(pkin(11));
t196 = cos(pkin(11));
t125 = -t171 * t201 + t204 * t259;
t212 = qJDD(2) * pkin(2) + t125;
t93 = t196 * t118 + t193 * t212;
t91 = -pkin(3) * t205 + qJDD(2) * pkin(8) + t93;
t221 = t205 * t217 + t91;
t60 = -t251 * pkin(4) + qJDD(4) * qJ(5) + t221 * t203 + t206;
t216 = t164 + t223;
t220 = t193 * t118 - t196 * t212;
t90 = -qJDD(2) * pkin(3) - t205 * pkin(8) + t220;
t63 = -t216 * qJ(5) + (-t165 + t182) * pkin(4) + t90;
t41 = 0.2e1 * qJD(5) * t160 + t192 * t60 - t195 * t63;
t34 = t213 * pkin(5) + t124 * pkin(9) - t41;
t145 = -pkin(5) * t230 - pkin(9) * t160;
t42 = -0.2e1 * qJD(5) * t158 + t192 * t63 + t195 * t60;
t35 = -t252 * pkin(5) - t218 * pkin(9) + t145 * t230 + t42;
t17 = t199 * t35 - t202 * t34;
t18 = t199 * t34 + t202 * t35;
t9 = -t17 * t202 + t18 * t199;
t250 = t192 * t9;
t249 = t195 * t9;
t137 = t203 * t207;
t59 = -qJDD(4) * pkin(4) - t251 * qJ(5) + t221 * t200 + qJDD(5) - t137;
t248 = t192 * t59;
t247 = t195 * t59;
t52 = t218 * pkin(5) - t252 * pkin(9) + t160 * t145 + t59;
t246 = t199 * t52;
t245 = t202 * t52;
t100 = t161 - t238;
t242 = t100 * t199;
t241 = t100 * t202;
t126 = t165 - t237;
t240 = t126 * t192;
t239 = t126 * t195;
t236 = t179 * t199;
t235 = t179 * t202;
t178 = t200 * t205 * t203;
t172 = qJDD(4) + t178;
t234 = t200 * t172;
t173 = qJDD(4) - t178;
t233 = t203 * t173;
t225 = t203 * t238;
t224 = t203 * t237;
t10 = t17 * t199 + t202 * t18;
t23 = t192 * t41 + t195 * t42;
t77 = t200 * t91 - t137;
t78 = t203 * t91 + t206;
t46 = t200 * t77 + t203 * t78;
t22 = t192 * t42 - t195 * t41;
t214 = -pkin(3) + t217;
t189 = t203 ^ 2;
t188 = t200 ^ 2;
t186 = t189 * t205;
t185 = t188 * t205;
t176 = -t186 - t251;
t175 = -t185 - t251;
t170 = t185 + t186;
t169 = (t188 + t189) * qJDD(2);
t168 = -qJDD(2) * t193 - t196 * t205;
t167 = qJDD(2) * t196 - t193 * t205;
t166 = -0.2e1 * t182 + t226;
t163 = 0.2e1 * t223 + t227;
t155 = t203 * t165;
t148 = -t157 - t186;
t147 = -t157 + t186;
t146 = -t186 + t252;
t142 = -t175 * t200 - t233;
t141 = t176 * t203 - t234;
t140 = -t173 * t200 + t175 * t203;
t139 = t172 * t203 + t176 * t200;
t136 = t169 * t193 + t170 * t196;
t135 = -t186 - t252;
t123 = t144 + t149;
t121 = -t150 + t218;
t117 = -t157 - t252;
t116 = -t131 + t177;
t115 = t130 - t177;
t114 = t142 * t193 - t163 * t196;
t113 = t141 * t193 + t166 * t196;
t109 = -t131 - t177;
t108 = -t148 * t192 + t239;
t107 = t148 * t195 + t240;
t106 = t131 - t130;
t105 = -t177 - t130;
t104 = t135 * t195 - t257;
t103 = t135 * t192 + t256;
t98 = -qJD(6) * t134 - t219;
t97 = t122 * t195 - t124 * t192;
t96 = t122 * t192 + t124 * t195;
t95 = (t132 * t202 - t134 * t199) * t179;
t94 = (t132 * t199 + t134 * t202) * t179;
t89 = -t130 - t131;
t88 = t108 * t203 + t123 * t200;
t87 = t108 * t200 - t123 * t203;
t86 = t104 * t203 + t121 * t200;
t85 = t104 * t200 - t121 * t203;
t84 = -t120 + t99;
t79 = (qJD(6) - t179) * t134 + t219;
t75 = t115 * t202 + t242;
t74 = -t116 * t199 + t253;
t73 = t115 * t199 - t241;
t72 = t116 * t202 + t254;
t71 = t134 * t236 + t202 * t99;
t70 = -t134 * t235 + t199 * t99;
t69 = -t132 * t235 - t199 * t98;
t68 = -t132 * t236 + t202 * t98;
t67 = t117 * t200 + t203 * t97;
t66 = -t117 * t203 + t200 * t97;
t65 = -t109 * t199 + t241;
t64 = t109 * t202 + t242;
t57 = t105 * t202 - t254;
t56 = t105 * t199 + t253;
t55 = -t107 * t196 + t193 * t88;
t54 = t193 * t93 - t196 * t220;
t53 = -t103 * t196 + t193 * t86;
t51 = t193 * t67 - t196 * t96;
t50 = t199 * t84 - t202 * t80;
t49 = -t199 * t258 - t202 * t79;
t48 = -t199 * t80 - t202 * t84;
t47 = -t199 * t79 + t202 * t258;
t45 = t200 * t78 - t203 * t77;
t44 = -t192 * t64 + t195 * t65;
t43 = t192 * t65 + t195 * t64;
t39 = -t192 * t56 + t195 * t57;
t38 = t192 * t57 + t195 * t56;
t37 = t193 * t46 - t196 * t90;
t36 = -pkin(9) * t64 + t245;
t33 = -pkin(9) * t56 + t246;
t31 = t200 * t258 + t203 * t44;
t30 = t200 * t44 - t203 * t258;
t29 = t200 * t79 + t203 * t39;
t28 = t200 * t39 - t203 * t79;
t27 = -pkin(5) * t258 + pkin(9) * t65 + t246;
t26 = -pkin(5) * t79 + pkin(9) * t57 - t245;
t25 = -t192 * t48 + t195 * t50;
t24 = t192 * t50 + t195 * t48;
t21 = t200 * t89 + t203 * t25;
t20 = t200 * t25 - t203 * t89;
t19 = t193 * t31 - t196 * t43;
t15 = t200 * t59 + t203 * t23;
t14 = t200 * t23 - t203 * t59;
t13 = t193 * t29 - t196 * t38;
t12 = t193 * t21 - t196 * t24;
t11 = t15 * t193 - t196 * t22;
t8 = -pkin(5) * t52 + pkin(9) * t10;
t7 = -pkin(9) * t48 - t9;
t6 = -pkin(5) * t89 + pkin(9) * t50 + t10;
t5 = t10 * t195 - t250;
t4 = t10 * t192 + t249;
t3 = t200 * t52 + t203 * t5;
t2 = t200 * t5 - t203 * t52;
t1 = t193 * t3 - t196 * t4;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t232, 0, 0, 0, 0, 0, 0, (qJDD(2) * t204 - t201 * t205) * t194, (-qJDD(2) * t201 - t204 * t205) * t194, 0, t197 ^ 2 * t232 + (t204 * t125 + t201 * t129 - t255) * t194, 0, 0, 0, 0, 0, 0, (t167 * t204 + t168 * t201) * t194, (-t167 * t201 + t168 * t204) * t194, 0, t197 * t215 + (t201 * (t193 * t220 + t196 * t93) + t204 * t54 - t255) * t194, 0, 0, 0, 0, 0, 0, t197 * t139 + (t201 * (t141 * t196 - t166 * t193) + t204 * t113) * t194, t197 * t140 + (t201 * (t142 * t196 + t163 * t193) + t204 * t114) * t194, (t201 * (t169 * t196 - t170 * t193) + t204 * t136) * t194, t197 * t45 + (t201 * (t193 * t90 + t196 * t46) + t204 * t37) * t194, 0, 0, 0, 0, 0, 0, t197 * t85 + (t201 * (t103 * t193 + t196 * t86) + t204 * t53) * t194, t197 * t87 + (t201 * (t107 * t193 + t196 * t88) + t204 * t55) * t194, t197 * t66 + (t201 * (t193 * t96 + t196 * t67) + t204 * t51) * t194, t197 * t14 + (t201 * (t15 * t196 + t193 * t22) + t204 * t11) * t194, 0, 0, 0, 0, 0, 0, t197 * t28 + (t201 * (t193 * t38 + t196 * t29) + t204 * t13) * t194, t197 * t30 + (t201 * (t193 * t43 + t196 * t31) + t204 * t19) * t194, t197 * t20 + (t201 * (t193 * t24 + t196 * t21) + t204 * t12) * t194, t197 * t2 + (t201 * (t193 * t4 + t196 * t3) + t204 * t1) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t125, -t129, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t167 - t220, pkin(2) * t168 - t93, 0, pkin(2) * t54, t216 * t200, t163 * t203 + t166 * t200, t234 + t203 * (-t185 + t251), -t200 * t223 + t155, t200 * (t186 - t251) + t233, 0, pkin(2) * t113 + pkin(3) * t166 + pkin(8) * t141 - t203 * t90, pkin(2) * t114 - pkin(3) * t163 + pkin(8) * t142 + t200 * t90, pkin(2) * t136 + pkin(3) * t170 + pkin(8) * t169 + t46, pkin(2) * t37 - pkin(3) * t90 + pkin(8) * t46, t200 * (t144 * t195 + t192 * t150) - t224, t200 * (-t121 * t195 - t123 * t192) + t203 * (-t157 + t252), t200 * (-t147 * t192 + t256) + t203 * t124, t200 * (-t195 * t149 + t192 * t218) + t224, t200 * (t146 * t195 + t240) - t203 * t122, t155 + t200 * (t158 * t195 - t160 * t192) * t230, t200 * (-qJ(5) * t103 + t248) + t203 * (-pkin(4) * t103 + t41) - pkin(3) * t103 + pkin(8) * t86 + pkin(2) * t53, t200 * (-qJ(5) * t107 + t247) + t203 * (-pkin(4) * t107 + t42) - pkin(3) * t107 + pkin(8) * t88 + pkin(2) * t55, pkin(2) * t51 + pkin(8) * t67 - t200 * t22 + t214 * t96, pkin(2) * t11 + pkin(8) * t15 + t214 * t22, t200 * (-t192 * t70 + t195 * t71) - t225, t200 * (-t192 * t47 + t195 * t49) - t203 * t106, t200 * (-t192 * t72 + t195 * t74) - t203 * t84, t200 * (-t192 * t68 + t195 * t69) + t225, t200 * (-t192 * t73 + t195 * t75) + t203 * t80, t200 * (-t192 * t94 + t195 * t95) + t203 * t161, t200 * (-qJ(5) * t38 - t192 * t26 + t195 * t33) + t203 * (-pkin(4) * t38 - pkin(5) * t56 + t17) - pkin(3) * t38 + pkin(8) * t29 + pkin(2) * t13, t200 * (-qJ(5) * t43 - t192 * t27 + t195 * t36) + t203 * (-pkin(4) * t43 - pkin(5) * t64 + t18) - pkin(3) * t43 + pkin(8) * t31 + pkin(2) * t19, t200 * (-qJ(5) * t24 - t192 * t6 + t195 * t7) + t203 * (-pkin(4) * t24 - pkin(5) * t48) - pkin(3) * t24 + pkin(8) * t21 + pkin(2) * t12, t200 * (-pkin(9) * t249 - qJ(5) * t4 - t192 * t8) + t203 * (-pkin(4) * t4 - pkin(5) * t9) - pkin(3) * t4 + pkin(8) * t3 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, 0, 0, 0, 0, 0, 0, t139, t140, 0, t45, 0, 0, 0, 0, 0, 0, t85, t87, t66, t14, 0, 0, 0, 0, 0, 0, t28, t30, t20, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t185 - t186, t227, t178, t226, qJDD(4), -t77, -t78, 0, 0, t144 * t192 - t195 * t150, -t121 * t192 + t123 * t195, t147 * t195 + t257, -t192 * t149 - t195 * t218, t146 * t192 - t239, (t158 * t192 + t160 * t195) * t230, -pkin(4) * t121 + qJ(5) * t104 - t247, -pkin(4) * t123 + qJ(5) * t108 + t248, -pkin(4) * t117 + qJ(5) * t97 + t23, -pkin(4) * t59 + qJ(5) * t23, t192 * t71 + t195 * t70, t192 * t49 + t195 * t47, t192 * t74 + t195 * t72, t192 * t69 + t195 * t68, t192 * t75 + t195 * t73, t192 * t95 + t195 * t94, -pkin(4) * t79 + qJ(5) * t39 + t192 * t33 + t195 * t26, -pkin(4) * t258 + qJ(5) * t44 + t192 * t36 + t195 * t27, -pkin(4) * t89 + qJ(5) * t25 + t192 * t7 + t195 * t6, -pkin(4) * t52 - pkin(9) * t250 + qJ(5) * t5 + t195 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t123, t117, t59, 0, 0, 0, 0, 0, 0, t79, t258, t89, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t106, t84, -t238, -t80, -t161, -t17, -t18, 0, 0;];
tauJ_reg  = t16;
