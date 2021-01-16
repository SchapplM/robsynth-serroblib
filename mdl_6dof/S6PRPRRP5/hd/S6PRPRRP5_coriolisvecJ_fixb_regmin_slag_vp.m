% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:49
% EndTime: 2021-01-16 01:51:58
% DurationCPUTime: 2.09s
% Computational Cost: add. (2145->299), mult. (4804->417), div. (0->0), fcn. (3221->8), ass. (0->168)
t94 = cos(qJ(5));
t145 = t94 * qJD(4);
t95 = cos(qJ(4));
t152 = qJD(2) * t95;
t91 = sin(qJ(5));
t68 = t91 * t152 - t145;
t92 = sin(qJ(4));
t153 = qJD(2) * t92;
t83 = qJD(5) + t153;
t183 = t68 * t83;
t146 = qJD(5) * t95;
t131 = t91 * t146;
t133 = t92 * t145;
t105 = t131 + t133;
t41 = qJD(2) * t105 - qJD(5) * t145;
t203 = -t41 - t183;
t89 = sin(pkin(6));
t156 = qJD(1) * t89;
t93 = sin(qJ(2));
t172 = t91 * t93;
t114 = pkin(4) * t95 + pkin(9) * t92;
t64 = t114 * qJD(4) + qJD(3);
t96 = cos(qJ(2));
t202 = (-t92 * t172 + t94 * t96) * t156 - t94 * t64;
t97 = -pkin(2) - pkin(8);
t169 = t92 * t97;
t74 = pkin(4) * t92 - pkin(9) * t95 + qJ(3);
t161 = t94 * t169 + t91 * t74;
t135 = t94 * t152;
t151 = qJD(4) * t91;
t70 = t135 + t151;
t182 = t70 * t83;
t150 = qJD(4) * t92;
t134 = t91 * t150;
t42 = -qJD(2) * t134 + qJD(5) * t70;
t201 = t42 + t182;
t147 = qJD(5) * t94;
t168 = t93 * t94;
t200 = (t92 * t168 + t91 * t96) * t156 - t95 * t97 * t145 - t74 * t147 - t91 * t64;
t155 = qJD(1) * t92;
t128 = t96 * t156;
t110 = qJD(3) - t128;
t57 = t97 * qJD(2) + t110;
t90 = cos(pkin(6));
t37 = -t90 * t155 + t57 * t95;
t199 = t70 ^ 2;
t198 = pkin(5) * t68;
t175 = t90 * t95;
t82 = qJD(1) * t175;
t38 = t92 * t57 + t82;
t30 = qJD(4) * pkin(9) + t38;
t129 = t93 * t156;
t48 = t74 * qJD(2) + t129;
t16 = -t30 * t91 + t94 * t48;
t10 = -qJ(6) * t70 + t16;
t7 = pkin(5) * t83 + t10;
t197 = t10 - t7;
t171 = t91 * t97;
t127 = pkin(5) - t171;
t144 = t94 * qJD(6);
t148 = qJD(5) * t91;
t196 = -qJ(6) * t133 + t161 * qJD(5) - (qJ(6) * t148 + t127 * qJD(4) - t144) * t95 + t202;
t130 = t94 * t146;
t195 = qJ(6) * t130 - (-qJD(6) * t95 + (qJ(6) * qJD(4) - qJD(5) * t97) * t92) * t91 + t200;
t174 = t91 * t48;
t17 = t30 * t94 + t174;
t11 = -qJ(6) * t68 + t17;
t194 = t11 * t83;
t154 = qJD(2) * t89;
t138 = t93 * t154;
t118 = t95 * t138;
t162 = -qJD(4) * t82 - t57 * t150;
t22 = -qJD(1) * t118 - t162;
t14 = pkin(5) * t42 + t22;
t193 = t14 * t91;
t192 = t14 * t94;
t191 = t22 * t91;
t190 = t22 * t94;
t29 = -qJD(4) * pkin(4) - t37;
t189 = t29 * t91;
t188 = t29 * t94;
t187 = t41 * t91;
t186 = t41 * t92;
t185 = t42 * t92;
t181 = t70 * t91;
t143 = qJD(2) * qJ(3);
t73 = t129 + t143;
t180 = t73 * t96;
t179 = t83 * t91;
t178 = t83 * t94;
t177 = t89 * t96;
t99 = qJD(2) ^ 2;
t176 = t89 * t99;
t173 = t91 * t92;
t170 = t92 * t94;
t167 = t95 * t41;
t166 = -qJ(6) - pkin(9);
t122 = qJD(5) * t166;
t72 = t114 * qJD(2);
t124 = -t91 * t37 + t94 * t72;
t165 = (pkin(5) * t95 + qJ(6) * t170) * qJD(2) + t124 + t91 * qJD(6) - t94 * t122;
t136 = t91 * t153;
t163 = t94 * t37 + t91 * t72;
t164 = qJ(6) * t136 - t91 * t122 - t144 + t163;
t88 = t95 ^ 2;
t160 = t92 ^ 2 - t88;
t98 = qJD(4) ^ 2;
t159 = -t98 - t99;
t158 = qJ(6) * t95;
t157 = qJD(2) * pkin(2);
t149 = qJD(4) * t95;
t142 = qJD(2) * qJD(4);
t141 = t93 * t176;
t140 = t96 * t176;
t137 = t96 * t154;
t132 = t83 * t148;
t126 = t95 * t142;
t21 = t57 * t149 + (-qJD(4) * t90 + t138) * t155;
t39 = (t64 + t128) * qJD(2);
t125 = -t91 * t21 + t94 * t39;
t123 = -qJD(6) - t198;
t121 = -t48 * t147 + t30 * t148 - t94 * t21 - t91 * t39;
t120 = t83 + t153;
t119 = t95 * t129;
t117 = t92 * t138;
t116 = qJD(5) * t92 + qJD(2);
t115 = pkin(5) * t126;
t113 = -t73 + t129;
t112 = t11 * t94 - t7 * t91;
t111 = t11 * t91 + t7 * t94;
t109 = qJD(2) * t88 - t83 * t92;
t56 = -t92 * t177 + t175;
t34 = t89 * t168 - t56 * t91;
t35 = t89 * t172 + t56 * t94;
t55 = t95 * t177 + t90 * t92;
t108 = -pkin(9) * t149 + t29 * t92;
t107 = t42 * qJ(6) + t121;
t106 = t113 * qJD(2);
t104 = t113 - t143;
t102 = -t17 * qJD(5) + t125;
t65 = (qJD(3) + t128) * qJD(2);
t101 = t110 * qJD(2) - t97 * t98 + t65;
t100 = t41 * qJ(6) + t102;
t85 = -pkin(5) * t94 - pkin(4);
t77 = t166 * t94;
t76 = t166 * t91;
t67 = t110 - t157;
t66 = (pkin(5) * t91 - t97) * t95;
t63 = t68 ^ 2;
t62 = t94 * t74;
t47 = t70 * t119;
t46 = t68 * t119;
t43 = t97 * t150 + (t130 - t134) * pkin(5);
t33 = t56 * qJD(4) - t118;
t32 = -t55 * qJD(4) + t117;
t31 = -t91 * t158 + t161;
t25 = t127 * t92 - t94 * t158 + t62;
t23 = -pkin(5) * t136 + t38;
t19 = -t123 + t29;
t13 = t34 * qJD(5) + t91 * t137 + t32 * t94;
t12 = -t35 * qJD(5) + t94 * t137 - t32 * t91;
t6 = -t95 * t42 - t116 * t178 + (-t120 * t95 * t91 + t68 * t92) * qJD(4);
t5 = t167 + t116 * t179 + (-t95 * t178 + (t70 - t135) * t92) * qJD(4);
t4 = -qJD(6) * t68 - t107;
t3 = -t70 * qJD(6) + t100 + t115;
t2 = -t35 * t126 - t13 * t83 + t33 * t70 - t41 * t55;
t1 = t12 * t83 + t34 * t126 + t33 * t68 + t42 * t55;
t8 = [0, 0, -t141, -t140, t141, t140, (t65 * t93 + (t180 + (t67 - t128) * t93) * qJD(2)) * t89, 0, 0, 0, 0, 0, t92 * t140 + (-t33 + t118) * qJD(4), t95 * t140 + (-t32 - t117) * qJD(4), 0, 0, 0, 0, 0, t1, t2, t1, t2, -t12 * t70 - t13 * t68 + t34 * t41 - t35 * t42, t11 * t13 + t12 * t7 + t14 * t55 + t19 * t33 + t3 * t34 + t35 * t4; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t65 * qJ(3) + t73 * qJD(3) + (-t180 + (-t67 - t157) * t93) * t156, -0.2e1 * t92 * t126, 0.2e1 * t160 * t142, -t98 * t92, -t98 * t95, 0, t101 * t92 - t104 * t149, t101 * t95 + t104 * t150, -t105 * t70 - t94 * t167, (t68 * t94 + t181) * t150 + (t187 - t42 * t94 + (t68 * t91 - t70 * t94) * qJD(5)) * t95, -t83 * t131 - t186 + (t109 * t94 + t70 * t95) * qJD(4), -t83 * t130 - t185 + (-t109 * t91 - t68 * t95) * qJD(4), t120 * t149, t46 + (-t74 * t148 - t202) * t83 + ((t68 * t97 - t189) * qJD(4) + (-t174 + (-t83 * t97 - t30) * t94) * qJD(5) + t125) * t92 + (t29 * t147 + t191 - t97 * t42 + (-t83 * t171 + (-t91 * t169 + t62) * qJD(2) + t16) * qJD(4)) * t95, t47 + t200 * t83 + (t97 * t132 + (t70 * t97 - t188) * qJD(4) + t121) * t92 + (-t29 * t148 + t190 + t97 * t41 + (-t161 * qJD(2) - t17) * qJD(4)) * t95, t42 * t66 + t43 * t68 + t46 + (-t19 * t151 + t3) * t92 - t196 * t83 + (t19 * t147 + t193 + (qJD(2) * t25 + t7) * qJD(4)) * t95, -t41 * t66 + t43 * t70 + t47 + (-t19 * t145 - t4) * t92 + t195 * t83 + (-t19 * t148 + t192 + (-qJD(2) * t31 - t11) * qJD(4)) * t95, t25 * t41 - t31 * t42 + t196 * t70 + t195 * t68 + t111 * t150 + (-t112 * qJD(5) - t3 * t94 - t4 * t91) * t95, t14 * t66 + t3 * t25 + t4 * t31 - t196 * t7 + (t43 + t119) * t19 - t195 * t11; 0, 0, 0, 0, 0, -t99, t106, 0, 0, 0, 0, 0, t159 * t92, t159 * t95, 0, 0, 0, 0, 0, t6, t5, t6, t5, (t116 * t70 - t149 * t68 - t185) * t94 + (t116 * t68 + t149 * t70 - t186) * t91, -t111 * qJD(2) + (qJD(4) * t112 - t14) * t95 + (qJD(4) * t19 - qJD(5) * t111 - t3 * t91 + t4 * t94) * t92; 0, 0, 0, 0, 0, 0, 0, t95 * t99 * t92, -t160 * t99, 0, 0, 0, t38 * qJD(4) + t106 * t95 + t162, -t113 * t153, t70 * t178 - t187, -t201 * t91 + t203 * t94, t83 * t147 + (t83 * t170 + (-t70 + t151) * t95) * qJD(2), -t132 + (-t83 * t173 + (t68 + t145) * t95) * qJD(2), -t83 * t152, -pkin(4) * t42 - t190 - t124 * t83 - t38 * t68 + (-pkin(9) * t178 + t189) * qJD(5) + (t108 * t91 - t16 * t95) * qJD(2), pkin(4) * t41 + t191 + t163 * t83 - t38 * t70 + (pkin(9) * t179 + t188) * qJD(5) + (t108 * t94 + t17 * t95) * qJD(2), -t192 - t23 * t68 + t42 * t85 - t165 * t83 + (t19 + t198) * t148 + (t19 * t173 + (qJD(4) * t76 - t7) * t95) * qJD(2), t193 - t23 * t70 - t41 * t85 + t164 * t83 + (pkin(5) * t181 + t19 * t94) * qJD(5) + (t19 * t170 + (qJD(4) * t77 + t11) * t95) * qJD(2), t41 * t76 + t42 * t77 + t165 * t70 + t164 * t68 + (-t7 * t83 + t4) * t94 + (-t3 - t194) * t91, t14 * t85 + t3 * t76 - t4 * t77 - t165 * t7 + (pkin(5) * t148 - t23) * t19 - t164 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t68, -t63 + t199, -t41 + t183, t182 - t42, t126, t17 * t83 - t29 * t70 + t102, t16 * t83 + t29 * t68 + t121, 0.2e1 * t115 + t194 + (t123 - t19) * t70 + t100, -t199 * pkin(5) + t10 * t83 + (qJD(6) + t19) * t68 + t107, t41 * pkin(5) + t197 * t68, -t197 * t11 + (-t19 * t70 + t3) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t203, -t63 - t199, t11 * t68 + t7 * t70 + t14;];
tauc_reg = t8;
