% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:27
% EndTime: 2019-03-09 02:16:35
% DurationCPUTime: 3.27s
% Computational Cost: add. (4610->396), mult. (9120->477), div. (0->0), fcn. (6467->10), ass. (0->186)
t131 = sin(pkin(9));
t132 = cos(pkin(9));
t136 = sin(qJ(4));
t244 = cos(qJ(4));
t156 = -t244 * t131 - t136 * t132;
t210 = qJD(1) * t131;
t196 = t136 * t210;
t147 = qJD(4) * t196 + qJDD(1) * t156;
t197 = t244 * t132;
t185 = qJD(1) * t197;
t56 = qJD(4) * t185 - t147;
t53 = qJDD(5) + t56;
t83 = t156 * qJD(1);
t255 = qJD(5) - t83;
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t134 = -pkin(1) - qJ(3);
t102 = t134 * qJD(1) + qJD(2);
t191 = -pkin(7) * qJD(1) + t102;
t74 = t191 * t131;
t75 = t191 * t132;
t42 = t136 * t75 + t244 * t74;
t39 = qJD(4) * pkin(8) + t42;
t85 = t185 - t196;
t113 = qJD(1) * qJ(2) + qJD(3);
t96 = pkin(3) * t210 + t113;
t40 = -pkin(4) * t83 - pkin(8) * t85 + t96;
t17 = t135 * t40 + t138 * t39;
t8 = qJ(6) * t255 + t17;
t262 = t255 * t8;
t209 = qJD(4) * t135;
t64 = t138 * t85 + t209;
t194 = t64 * t255;
t155 = -t136 * t131 + t197;
t151 = t155 * qJDD(1);
t256 = t83 * qJD(4);
t141 = t151 + t256;
t261 = -qJD(4) * qJD(5) - t141;
t127 = pkin(9) + qJ(4);
t114 = sin(t127);
t115 = cos(t127);
t137 = sin(qJ(1));
t139 = cos(qJ(1));
t254 = g(1) * t137 - g(2) * t139;
t149 = g(3) * t114 - t254 * t115;
t207 = qJD(5) * t135;
t27 = -t135 * qJDD(4) + t138 * t261 + t207 * t85;
t195 = qJD(4) * t244;
t208 = qJD(4) * t136;
t86 = -t131 * t195 - t132 * t208;
t173 = -t155 * t27 + t86 * t64;
t188 = -qJD(5) * t156 + qJD(1);
t48 = t138 * t53;
t87 = -t131 * t208 + t132 * t195;
t260 = t156 * t48 + (t135 * t188 - t138 * t87) * t255 - t173;
t259 = -t136 * t74 + t244 * t75;
t243 = g(1) * t139;
t178 = g(2) * t137 + t243;
t128 = qJDD(1) * qJ(2);
t129 = qJD(1) * qJD(2);
t253 = t128 + t129;
t98 = qJDD(3) + t253;
t258 = t98 - t178;
t211 = t131 ^ 2 + t132 ^ 2;
t257 = t102 * t211;
t249 = -qJD(1) * qJD(3) + qJDD(1) * t134;
t92 = qJDD(2) + t249;
t190 = -pkin(7) * qJDD(1) + t92;
t68 = t190 * t131;
t69 = t190 * t132;
t163 = t136 * t69 + t244 * t68;
t14 = qJDD(4) * pkin(8) + qJD(4) * t259 + t163;
t206 = qJD(5) * t138;
t203 = t131 * qJDD(1);
t89 = pkin(3) * t203 + t98;
t23 = t56 * pkin(4) - pkin(8) * t141 + t89;
t192 = t135 * t14 - t138 * t23 + t39 * t206 + t40 * t207;
t246 = pkin(5) * t53;
t2 = qJDD(6) + t192 - t246;
t252 = -t2 + t262;
t38 = -qJD(4) * pkin(4) - t259;
t62 = -t138 * qJD(4) + t135 * t85;
t18 = t62 * pkin(5) - t64 * qJ(6) + t38;
t245 = pkin(8) * t53;
t251 = t18 * t255 - t245;
t241 = g(3) * t115;
t148 = -t114 * t254 - t241;
t248 = t64 ^ 2;
t247 = 0.2e1 * t129;
t240 = g(3) * t138;
t119 = t131 * pkin(3);
t239 = t17 * t255;
t238 = t62 * t83;
t237 = t64 * t62;
t236 = t64 * t85;
t235 = t85 * t62;
t234 = -pkin(7) + t134;
t176 = -t138 * qJDD(4) + t206 * t85;
t28 = t135 * t151 + (qJD(5) + t83) * t209 + t176;
t233 = -t135 * t28 - t62 * t206;
t55 = pkin(4) * t85 - pkin(8) * t83;
t232 = t135 * t55 + t138 * t259;
t107 = qJ(2) + t119;
t54 = -pkin(4) * t156 - pkin(8) * t155 + t107;
t93 = t234 * t131;
t94 = t234 * t132;
t59 = t136 * t94 + t244 * t93;
t231 = t135 * t54 + t138 * t59;
t230 = t86 * qJD(4) + qJDD(4) * t155;
t169 = pkin(5) * t135 - qJ(6) * t138;
t229 = -t135 * qJD(6) + t169 * t255 - t42;
t228 = pkin(8) * qJD(5);
t46 = t135 * t53;
t227 = t138 * t255;
t225 = t27 * t135;
t224 = t53 * qJ(6);
t216 = t139 * t138;
t223 = g(2) * t115 * t216 + t114 * t240;
t222 = pkin(1) * qJDD(1);
t221 = qJD(5) * t155;
t220 = t115 * t137;
t219 = t135 * t137;
t218 = t135 * t139;
t217 = t137 * t138;
t16 = -t135 * t39 + t138 * t40;
t214 = qJD(6) - t16;
t213 = t139 * pkin(1) + t137 * qJ(2);
t202 = g(1) * t220;
t201 = t255 * t228;
t199 = t155 * t207;
t198 = t155 * t206;
t193 = t211 * t92;
t189 = t135 * t255;
t187 = qJDD(2) - t222;
t159 = t136 * t68 + t74 * t195 + t75 * t208 - t244 * t69;
t15 = -qJDD(4) * pkin(4) + t159;
t3 = t28 * pkin(5) + t27 * qJ(6) - t64 * qJD(6) + t15;
t186 = -t3 - t201;
t78 = t114 * t219 - t216;
t80 = t114 * t218 + t217;
t184 = g(1) * t80 + g(2) * t78;
t79 = t114 * t217 + t218;
t81 = t114 * t216 - t219;
t183 = -g(1) * t81 - g(2) * t79;
t158 = t135 * t23 + t138 * t14 + t40 * t206 - t207 * t39;
t1 = qJD(6) * t255 + t158 + t224;
t182 = -t1 * t156 + t8 * t87;
t7 = -pkin(5) * t255 + t214;
t181 = -t156 * t2 + t7 * t87;
t180 = t155 * t3 + t18 * t86;
t179 = t255 * t7 + t1;
t175 = -t135 * t8 + t138 * t7;
t174 = t15 * t155 + t38 * t86;
t172 = t156 * t27 + t64 * t87;
t171 = t156 * t28 - t62 * t87;
t170 = -t155 * t53 - t255 * t86;
t168 = t206 * t255 - t83 * t227 + t46;
t167 = -qJD(4) * t87 + qJDD(4) * t156;
t166 = t48 + (t135 * t83 - t207) * t255;
t164 = pkin(5) * t138 + qJ(6) * t135 + pkin(4);
t161 = -t136 * t93 + t244 * t94;
t160 = t255 * t38 - t245;
t35 = qJD(3) * t156 + qJD(4) * t161;
t51 = pkin(4) * t87 - pkin(8) * t86 + qJD(2);
t157 = t135 * t51 + t138 * t35 + t54 * t206 - t207 * t59;
t153 = pkin(4) * t114 - pkin(8) * t115 + t119;
t150 = g(1) * t78 - g(2) * t80 + t135 * t241 - t192;
t145 = t253 + t258;
t144 = t18 * t64 + qJDD(6) - t150;
t143 = -g(1) * t79 + g(2) * t81 - t115 * t240 + t158;
t142 = t156 * t46 - t155 * t28 - t86 * t62 + (-t135 * t87 - t138 * t188) * t255;
t36 = qJD(3) * t155 + qJD(4) * t59;
t140 = qJD(1) ^ 2;
t133 = -pkin(7) - qJ(3);
t121 = t139 * qJ(2);
t30 = pkin(5) * t64 + qJ(6) * t62;
t29 = t155 * t169 - t161;
t22 = pkin(5) * t156 + t135 * t59 - t138 * t54;
t21 = -qJ(6) * t156 + t231;
t13 = -pkin(5) * t85 + t135 * t259 - t138 * t55;
t12 = qJ(6) * t85 + t232;
t9 = t255 * t62 - t27;
t6 = (pkin(5) * t86 + qJ(6) * t221) * t135 + (-qJ(6) * t86 - (-pkin(5) * qJD(5) + qJD(6)) * t155) * t138 + t36;
t5 = -t87 * pkin(5) + t231 * qJD(5) + t135 * t35 - t138 * t51;
t4 = qJ(6) * t87 - qJD(6) * t156 + t157;
t10 = [qJDD(1), t254, t178, qJDD(2) - t254 - 0.2e1 * t222, 0.2e1 * t128 + t247 - t178, -t187 * pkin(1) - g(1) * (-t137 * pkin(1) + t121) - g(2) * t213 + (t128 + t247) * qJ(2), t145 * t131, t145 * t132, t254 + t211 * (-t249 - t92) t98 * qJ(2) + t113 * qJD(2) - g(1) * (t134 * t137 + t121) - g(2) * (qJ(3) * t139 + t213) + t134 * t193 - qJD(3) * t257, t141 * t155 + t85 * t86, t141 * t156 - t155 * t56 + t86 * t83 - t85 * t87, t230, t167, 0, -qJD(2) * t83 - t36 * qJD(4) + qJDD(4) * t161 + t107 * t56 - t114 * t178 - t156 * t89 + t96 * t87, -g(2) * t220 + qJD(2) * t85 - t35 * qJD(4) - t59 * qJDD(4) + t107 * t141 - t115 * t243 + t155 * t89 + t96 * t86, t138 * t173 - t199 * t64 (-t135 * t64 - t138 * t62) * t86 - (-t225 + t138 * t28 + (-t135 * t62 + t138 * t64) * qJD(5)) * t155, -t138 * t170 - t199 * t255 + t172, t135 * t170 - t198 * t255 + t171, -t156 * t53 + t255 * t87, t192 * t156 + t16 * t87 + t36 * t62 - t161 * t28 + ((-qJD(5) * t59 + t51) * t255 + t54 * t53 + t38 * t221) * t138 + ((-qJD(5) * t54 - t35) * t255 - t59 * t53 + t174) * t135 + t183, t174 * t138 + t156 * t158 - t157 * t255 + t161 * t27 - t17 * t87 - t38 * t199 - t231 * t53 + t36 * t64 + t184, t135 * t180 + t18 * t198 - t22 * t53 - t255 * t5 + t29 * t28 + t6 * t62 - t181 + t183, -t21 * t28 - t22 * t27 - t4 * t62 + t5 * t64 + t175 * t86 + t178 * t115 - (t1 * t135 - t2 * t138 + (t135 * t7 + t138 * t8) * qJD(5)) * t155, -t138 * t180 + t18 * t199 + t21 * t53 + t255 * t4 + t29 * t27 - t6 * t64 + t182 - t184, t1 * t21 + t8 * t4 + t3 * t29 + t18 * t6 + t2 * t22 + t7 * t5 - g(1) * (t81 * pkin(5) + t80 * qJ(6) + t121) - g(2) * (t79 * pkin(5) + t78 * qJ(6) + t213) + (-g(1) * t153 + g(2) * t133) * t139 + (-g(1) * (-pkin(1) + t133) - g(2) * t153) * t137; 0, 0, 0, qJDD(1), -t140, -qJ(2) * t140 + t187 - t254, -t140 * t131, -t140 * t132, -t211 * qJDD(1), -t113 * qJD(1) + t193 - t254, 0, 0, 0, 0, 0, qJD(1) * t83 + t230, -qJD(1) * t85 + t167, 0, 0, 0, 0, 0, t142, t260, t142 (t188 * t64 + t171) * t138 + (t188 * t62 + t172) * t135, -t260 (t188 * t7 + t182) * t138 + (-t188 * t8 + t181) * t135 - t180 - t254; 0, 0, 0, 0, 0, 0, t203, t132 * qJDD(1), -t211 * t140, qJD(1) * t257 + t258, 0, 0, 0, 0, 0 (t185 + t85) * qJD(4) - t147, t151 + 0.2e1 * t256, 0, 0, 0, 0, 0, t166 - t235, -t227 * t255 - t236 - t46, -t189 * t255 - t235 + t48 (t27 + t238) * t138 + t135 * t194 + t233, t168 + t236, t179 * t135 + t252 * t138 - t18 * t85 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t83, -t83 ^ 2 + t85 ^ 2, t151 (-t185 + t85) * qJD(4) + t147, qJDD(4), t42 * qJD(4) - t96 * t85 + t149 - t159, -t96 * t83 - t148 - t163, t138 * t194 - t225 (-t27 + t238) * t138 - t64 * t189 + t233, t168 - t236, t166 + t235, -t255 * t85, -pkin(4) * t28 - t16 * t85 - t42 * t62 + (-t202 - t15 + (-t55 - t228) * t255) * t138 + (t255 * t259 + t160) * t135 + t223, pkin(4) * t27 + t232 * t255 + t17 * t85 - t42 * t64 + t160 * t138 + (-t149 + t15 + t201) * t135, t13 * t255 - t164 * t28 + t7 * t85 + t229 * t62 + (t186 - t202) * t138 + t251 * t135 + t223, t12 * t62 - t13 * t64 + ((qJD(5) * t64 - t28) * pkin(8) + t179) * t138 + ((qJD(5) * t62 - t27) * pkin(8) - t252) * t135 + t148, -t12 * t255 - t164 * t27 - t8 * t85 - t229 * t64 - t251 * t138 + (t149 + t186) * t135, -t8 * t12 - t7 * t13 + t229 * t18 + (qJD(5) * t175 + t1 * t138 + t2 * t135 + t148) * pkin(8) + (-t3 + t149) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, -t62 ^ 2 + t248, t9, t135 * t261 - t176 + t194, t53, -t38 * t64 + t150 + t239, t16 * t255 + t38 * t62 - t143, -t30 * t62 - t144 + t239 + 0.2e1 * t246, pkin(5) * t27 - t28 * qJ(6) + (-t17 + t8) * t64 + (t7 - t214) * t62, 0.2e1 * t224 - t18 * t62 + t30 * t64 + (0.2e1 * qJD(6) - t16) * t255 + t143, t1 * qJ(6) - t2 * pkin(5) - t18 * t30 - t7 * t17 - g(1) * (-pkin(5) * t78 + qJ(6) * t79) - g(2) * (pkin(5) * t80 - qJ(6) * t81) + t214 * t8 + t169 * t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237 - t53, t9, -t255 ^ 2 - t248, t144 - t246 - t262;];
tau_reg  = t10;
