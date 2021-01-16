% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:23
% EndTime: 2021-01-15 16:45:38
% DurationCPUTime: 3.05s
% Computational Cost: add. (2224->386), mult. (5251->537), div. (0->0), fcn. (4021->10), ass. (0->194)
t132 = cos(qJ(4));
t130 = sin(qJ(3));
t191 = qJD(3) * t130;
t129 = sin(qJ(4));
t234 = pkin(7) * t129;
t125 = sin(pkin(5));
t131 = sin(qJ(2));
t133 = cos(qJ(3));
t134 = cos(qJ(2));
t200 = t133 * t134;
t57 = (-t129 * t200 + t131 * t132) * t125;
t158 = pkin(3) * t130 - pkin(8) * t133;
t95 = t158 * qJD(3);
t254 = -qJD(1) * t57 + t132 * t95 + t191 * t234;
t188 = qJD(4) * t132;
t58 = (t129 * t131 + t132 * t200) * t125;
t99 = -pkin(3) * t133 - pkin(8) * t130 - pkin(2);
t253 = -qJD(1) * t58 + t129 * t95 + t99 * t188;
t194 = qJD(2) * t130;
t252 = qJD(4) * t194 - qJDD(3);
t182 = t130 * qJDD(2);
t193 = qJD(2) * t133;
t25 = ((qJD(4) + t193) * qJD(3) + t182) * t129 + t252 * t132;
t201 = t132 * t133;
t113 = pkin(7) * t201;
t153 = pkin(4) * t130 - qJ(5) * t201;
t187 = qJD(5) * t132;
t216 = qJ(5) * t130;
t251 = -t130 * t187 + t153 * qJD(3) + (-t113 + (-t99 + t216) * t129) * qJD(4) + t254;
t204 = t130 * t132;
t233 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t204 + (-qJD(5) * t130 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t133) * t129 + t253;
t112 = -qJD(4) + t193;
t186 = t132 * qJD(3);
t90 = t129 * t194 - t186;
t227 = t112 * t90;
t184 = qJD(2) * qJD(3);
t170 = t133 * t184;
t24 = -qJD(4) * t186 + (-t170 - t182) * t132 + t252 * t129;
t250 = -t24 + t227;
t192 = qJD(3) * t129;
t92 = t132 * t194 + t192;
t226 = t112 * t92;
t249 = t25 - t226;
t127 = cos(pkin(5));
t196 = qJD(1) * t127;
t197 = qJD(1) * t125;
t97 = qJD(2) * pkin(7) + t131 * t197;
t248 = -t130 * t97 + t133 * t196;
t135 = qJD(3) ^ 2;
t185 = qJD(1) * qJD(2);
t171 = t131 * t185;
t209 = t125 * t134;
t157 = -qJDD(1) * t209 + t125 * t171;
t124 = sin(pkin(9));
t126 = cos(pkin(9));
t206 = t127 * t134;
t71 = t124 * t131 - t126 * t206;
t73 = t124 * t206 + t126 * t131;
t159 = g(1) * t73 + g(2) * t71;
t247 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t135 + (-g(3) * t134 + t171) * t125 - t157 + t159;
t202 = t131 * t133;
t77 = t125 * t202 + t127 * t130;
t40 = -t129 * t77 - t132 * t209;
t242 = g(3) * t40;
t211 = t125 * t130;
t207 = t127 * t131;
t70 = t124 * t207 - t126 * t134;
t29 = -t124 * t211 + t133 * t70;
t72 = t124 * t134 + t126 * t207;
t31 = -t126 * t211 + t133 * t72;
t62 = t71 * t132;
t65 = t73 * t132;
t245 = -t242 - g(1) * (t129 * t29 + t65) - g(2) * (-t129 * t31 + t62);
t244 = t92 ^ 2;
t243 = pkin(4) * t90;
t150 = t129 * t209 - t132 * t77;
t241 = g(3) * t150;
t240 = g(3) * t57;
t239 = g(3) * t58;
t205 = t130 * t131;
t76 = t125 * t205 - t127 * t133;
t238 = g(3) * t76;
t120 = t133 * qJDD(2);
t87 = t130 * t184 + qJDD(4) - t120;
t237 = t87 * pkin(4);
t110 = t130 * t196;
t52 = t133 * t97 + t110;
t46 = qJD(3) * pkin(8) + t52;
t179 = t134 * t197;
t54 = t99 * qJD(2) - t179;
t17 = -t129 * t46 + t132 * t54;
t11 = -qJ(5) * t92 + t17;
t6 = -pkin(4) * t112 + t11;
t236 = t11 - t6;
t128 = -qJ(5) - pkin(8);
t94 = t158 * qJD(2);
t232 = t129 * t94 + t132 * t248;
t165 = qJD(4) * t128;
t175 = t129 * t193;
t231 = qJ(5) * t175 + t129 * t165 + t187 - t232;
t81 = t132 * t94;
t230 = -t153 * qJD(2) + t132 * t165 - t81 + (-qJD(5) + t248) * t129;
t228 = qJD(2) * pkin(2);
t18 = t129 * t54 + t132 * t46;
t12 = -qJ(5) * t90 + t18;
t225 = t12 * t112;
t224 = t129 * t24;
t221 = t129 * t71;
t220 = t129 * t73;
t219 = t132 * t92;
t66 = t70 * t129;
t63 = t72 * t129;
t217 = t129 * t99 + t113;
t215 = qJD(3) * t90;
t214 = qJD(3) * t92;
t212 = t124 * t127;
t210 = t125 * t133;
t208 = t126 * t127;
t203 = t130 * t134;
t199 = qJDD(1) - g(3);
t122 = t130 ^ 2;
t198 = -t133 ^ 2 + t122;
t195 = qJD(2) * t125;
t190 = qJD(3) * t133;
t189 = qJD(4) * t129;
t183 = qJDD(1) * t127;
t166 = -qJD(5) - t243;
t45 = -qJD(3) * pkin(3) - t248;
t23 = -t166 + t45;
t181 = t23 * t188;
t180 = pkin(4) * t129 + pkin(7);
t177 = t131 * t195;
t176 = t134 * t195;
t174 = t112 * t186;
t173 = t112 * t189;
t172 = t112 * t188;
t169 = t134 * t184;
t167 = t130 * t183;
t60 = qJDD(2) * pkin(7) + (qJDD(1) * t131 + t134 * t185) * t125;
t14 = qJDD(3) * pkin(8) + qJD(3) * t248 + t133 * t60 + t167;
t22 = qJD(2) * t95 + t99 * qJDD(2) + t157;
t164 = -t129 * t22 - t132 * t14 - t54 * t188 + t46 * t189;
t162 = t90 * t179;
t161 = t92 * t179;
t30 = t126 * t210 + t130 * t72;
t32 = t124 * t210 + t130 * t70;
t160 = g(1) * t32 - g(2) * t30;
t143 = qJD(3) * t110 + t130 * t60 - t133 * t183 + t97 * t190;
t15 = -qJDD(3) * pkin(3) + t143;
t5 = t25 * pkin(4) + qJDD(5) + t15;
t156 = -t160 - t5;
t117 = pkin(4) * t132 + pkin(3);
t155 = -t117 * t133 + t128 * t130;
t136 = qJD(2) ^ 2;
t154 = qJDD(2) * t134 - t131 * t136;
t151 = t25 * qJ(5) + t164;
t147 = t129 * t87 - t172;
t146 = t132 * t87 + t173;
t145 = -g(1) * t29 + g(2) * t31 + g(3) * t77;
t144 = -t160 + t238;
t142 = -pkin(8) * t87 - t112 * t45;
t75 = -t127 * t205 - t210;
t141 = -pkin(8) * qJD(4) * t112 - g(1) * (t124 * t75 + t126 * t203) - g(2) * (t124 * t203 - t126 * t75) + t15;
t21 = t132 * t22;
t139 = -t18 * qJD(4) - t129 * t14 + t21;
t98 = -t179 - t228;
t138 = -pkin(7) * qJDD(3) + (t179 + t98 - t228) * qJD(3);
t137 = t24 * qJ(5) + t139;
t101 = t128 * t132;
t100 = t128 * t129;
t96 = t180 * t130;
t89 = t132 * t99;
t86 = t90 ^ 2;
t78 = t127 * t202 - t211;
t67 = t70 * t132;
t64 = t72 * t132;
t61 = t129 * t238;
t56 = t73 * t133;
t55 = t71 * t133;
t53 = pkin(7) * t190 + (t129 * t190 + t130 * t188) * pkin(4);
t42 = -t129 * t216 + t217;
t39 = t77 * qJD(3) + t130 * t176;
t38 = -t76 * qJD(3) + t133 * t176;
t37 = -t124 * t78 + t126 * t200;
t35 = t124 * t200 + t126 * t78;
t27 = pkin(4) * t175 + t52;
t26 = -qJ(5) * t204 + t89 + (-pkin(4) - t234) * t133;
t9 = t40 * qJD(4) + t129 * t177 + t132 * t38;
t8 = t150 * qJD(4) - t129 * t38 + t132 * t177;
t4 = t112 * t9 + t150 * t87 - t24 * t76 + t39 * t92;
t3 = -t112 * t8 + t25 * t76 + t39 * t90 + t40 * t87;
t2 = -t90 * qJD(5) - t151;
t1 = -t92 * qJD(5) + t137 + t237;
t7 = [t199, 0, t154 * t125, (-qJDD(2) * t131 - t134 * t136) * t125, 0, 0, 0, 0, 0, -qJD(3) * t39 - qJDD(3) * t76 + (-t130 * t169 + t154 * t133) * t125, -qJD(3) * t38 - qJDD(3) * t77 + (-t154 * t130 - t133 * t169) * t125, 0, 0, 0, 0, 0, t3, t4, t3, t4, t150 * t25 + t24 * t40 - t8 * t92 - t9 * t90, t1 * t40 + t12 * t9 - t150 * t2 + t23 * t39 + t5 * t76 + t6 * t8 - g(3); 0, qJDD(2), t199 * t209 + t159, -t199 * t131 * t125 - g(1) * t70 + g(2) * t72, qJDD(2) * t122 + 0.2e1 * t130 * t170, 0.2e1 * t130 * t120 - 0.2e1 * t198 * t184, qJDD(3) * t130 + t133 * t135, qJDD(3) * t133 - t130 * t135, 0, t138 * t130 + t247 * t133, -t247 * t130 + t138 * t133, t133 * t92 * t186 + (-t132 * t24 - t92 * t189) * t130, (-t129 * t92 - t132 * t90) * t190 + (t224 - t132 * t25 + (t129 * t90 - t219) * qJD(4)) * t130, (t24 - t174) * t133 + (t146 + t214) * t130, (t112 * t192 + t25) * t133 + (-t147 - t215) * t130, -t112 * t191 - t133 * t87, t89 * t87 - g(1) * (-t132 * t56 - t66) - g(2) * (-t132 * t55 + t63) - t239 + (t189 * t99 - t254) * t112 + (t46 * t188 - t21 + (t172 + t215) * pkin(7) + (-pkin(7) * t87 + qJD(3) * t45 + qJD(4) * t54 + t14) * t129) * t133 + (pkin(7) * t25 + t17 * qJD(3) + t15 * t129 + t188 * t45 - t162) * t130, -t217 * t87 - g(1) * (t129 * t56 - t67) - g(2) * (t129 * t55 + t64) - t240 + t253 * t112 + (t45 * t186 + (-t173 + t214) * pkin(7) - t164) * t133 + (-t161 - t45 * t189 - t18 * qJD(3) + t15 * t132 + (-t24 - t174) * pkin(7)) * t130, g(1) * t66 - g(2) * t63 - t239 + t96 * t25 + t26 * t87 + t53 * t90 + (t132 * t159 + t192 * t23 - t1) * t133 - t251 * t112 + (qJD(3) * t6 + t129 * t5 - t162 + t181) * t130, g(1) * t67 - g(2) * t64 - t240 - t96 * t24 - t42 * t87 + t53 * t92 + (-t129 * t159 + t186 * t23 + t2) * t133 + t233 * t112 + (-qJD(3) * t12 + t132 * t5 - t189 * t23 - t161) * t130, t24 * t26 - t25 * t42 - t251 * t92 - t233 * t90 + (-t12 * t129 - t132 * t6) * t190 + (-g(3) * t209 - t1 * t132 - t129 * t2 + (-t12 * t132 + t129 * t6) * qJD(4) + t159) * t130, t2 * t42 + t1 * t26 + t5 * t96 - g(1) * (-pkin(4) * t66 - (pkin(2) * t126 + pkin(7) * t212) * t131 + (-pkin(2) * t212 + pkin(7) * t126) * t134 + t155 * t73) - g(2) * (pkin(4) * t63 - (pkin(2) * t124 - pkin(7) * t208) * t131 + (pkin(2) * t208 + pkin(7) * t124) * t134 + t155 * t71) - g(3) * (t180 * t131 + (pkin(2) - t155) * t134) * t125 + t251 * t6 + (-t130 * t179 + t53) * t23 + t233 * t12; 0, 0, 0, 0, -t130 * t136 * t133, t198 * t136, t182, t120, qJDD(3), qJD(3) * t52 - t98 * t194 - t143 + t144, -t167 + (-qJD(2) * t98 - t60) * t133 + t145, -t112 * t219 - t224, -t249 * t129 + t250 * t132, (t112 * t201 - t130 * t92) * qJD(2) + t147, (-t112 * t129 * t133 + t130 * t90) * qJD(2) + t146, t112 * t194, -t17 * t194 - pkin(3) * t25 + t81 * t112 - t52 * t90 + (-t112 * t248 + t142) * t129 + (-t141 + t238) * t132, pkin(3) * t24 - t112 * t232 + t129 * t141 + t132 * t142 + t18 * t194 - t52 * t92 - t61, -t6 * t194 + t100 * t87 - t117 * t25 - t27 * t90 - t230 * t112 + (-t23 * t193 + (t23 + t243) * qJD(4)) * t129 + (t144 - t5) * t132, t181 + t101 * t87 + t117 * t24 - t27 * t92 - t61 + t231 * t112 + (t12 * t130 - t201 * t23) * qJD(2) + (pkin(4) * qJD(4) * t92 - t156) * t129, t100 * t24 + t101 * t25 - t230 * t92 - t231 * t90 + (t112 * t6 + t2) * t132 + (-t1 + t225) * t129 - t145, -t2 * t101 + t1 * t100 - t5 * t117 - g(1) * (t117 * t32 + t128 * t29) - g(2) * (-t117 * t30 - t128 * t31) - g(3) * (-t117 * t76 - t128 * t77) + t230 * t6 + (pkin(4) * t189 - t27) * t23 + t231 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92 * t90, -t86 + t244, -t24 - t227, -t226 - t25, t87, -t18 * t112 - t45 * t92 - g(1) * (-t129 * t37 + t65) - g(2) * (-t129 * t35 + t62) - t242 + t139, -t17 * t112 + t45 * t90 - g(1) * (-t132 * t37 - t220) - g(2) * (-t132 * t35 - t221) - t241 + t164, 0.2e1 * t237 - t225 + (t166 - t23) * t92 + t137 + t245, -t11 * t112 - t244 * pkin(4) - g(1) * (t132 * t29 - t220) - g(2) * (-t132 * t31 - t221) - t241 + (qJD(5) + t23) * t90 + t151, pkin(4) * t24 + t236 * t90, -t236 * t12 + (-t23 * t92 + t1 + t245) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t250, -t86 - t244, t12 * t90 + t6 * t92 - t156 - t238;];
tau_reg = t7;
