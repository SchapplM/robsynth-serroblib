% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:18
% DurationCPUTime: 1.49s
% Computational Cost: add. (3433->255), mult. (5942->311), div. (0->0), fcn. (3723->6), ass. (0->164)
t133 = cos(qJ(3));
t218 = cos(qJ(4));
t169 = t218 * t133;
t130 = sin(qJ(4));
t131 = sin(qJ(3));
t193 = t130 * t131;
t144 = t169 - t193;
t221 = -pkin(8) - pkin(7);
t107 = t221 * t131;
t124 = t133 * pkin(8);
t108 = pkin(7) * t133 + t124;
t145 = t218 * t107 - t130 * t108;
t170 = qJD(3) * t221;
t134 = cos(qJ(2));
t203 = pkin(1) * qJD(1);
t176 = t134 * t203;
t192 = t130 * t133;
t94 = t131 * t170;
t207 = t145 * qJD(4) - t144 * t176 + t170 * t192 + t218 * t94;
t128 = t131 ^ 2;
t129 = t133 ^ 2;
t185 = t128 + t129;
t153 = qJD(3) * t169;
t162 = qJD(4) * t218;
t224 = -t133 * t162 - t153;
t127 = qJD(1) + qJD(2);
t92 = t218 * t131 + t192;
t81 = t92 * t127;
t126 = qJD(3) + qJD(4);
t223 = t207 * t126;
t222 = t81 ^ 2;
t148 = t126 * t193;
t188 = t224 * t127;
t47 = t127 * t148 + t188;
t64 = t126 * t92;
t48 = t64 * t127;
t132 = sin(qJ(2));
t202 = pkin(1) * qJD(2);
t173 = qJD(1) * t202;
t113 = t132 * t173;
t184 = qJD(3) * t131;
t168 = t127 * t184;
t85 = pkin(3) * t168 + t113;
t12 = pkin(4) * t48 + qJ(5) * t47 - qJD(5) * t81 + t85;
t172 = t127 * t193;
t79 = -t127 * t169 + t172;
t119 = -pkin(3) * t133 - pkin(2);
t83 = t119 * t127 - t176;
t35 = pkin(4) * t79 - qJ(5) * t81 + t83;
t63 = t148 + t224;
t220 = -t12 * t92 + t35 * t63;
t219 = -t12 * t144 + t35 * t64;
t217 = pkin(1) * t134;
t177 = t132 * t203;
t98 = pkin(7) * t127 + t177;
t163 = pkin(8) * t127 + t98;
t149 = qJD(3) * t163;
t156 = t134 * t173;
t54 = -t131 * t149 + t133 * t156;
t55 = -t131 * t156 - t133 * t149;
t160 = t130 * t54 - t218 * t55;
t75 = t163 * t133;
t174 = t218 * t75;
t74 = t163 * t131;
t70 = qJD(3) * pkin(3) - t74;
t41 = t130 * t70 + t174;
t11 = qJD(4) * t41 + t160;
t116 = pkin(1) * t132 + pkin(7);
t208 = -pkin(8) - t116;
t89 = t208 * t131;
t90 = t116 * t133 + t124;
t146 = -t130 * t90 + t218 * t89;
t216 = t11 * t146;
t215 = t11 * t145;
t214 = t11 * t92;
t213 = t35 * t79;
t212 = t35 * t81;
t211 = t81 * t79;
t210 = t83 * t79;
t209 = t83 * t81;
t72 = t130 * t107 + t218 * t108;
t206 = t72 * qJD(4) + t130 * t94 - t221 * t153 - t92 * t176;
t205 = -t144 * t85 + t83 * t64;
t204 = -t83 * t63 + t85 * t92;
t159 = qJD(3) * t208;
t175 = t134 * t202;
t139 = -t131 * t175 + t133 * t159;
t73 = t131 * t159 + t133 * t175;
t17 = t146 * qJD(4) + t130 * t139 + t218 * t73;
t201 = t126 * t17;
t57 = t130 * t89 + t218 * t90;
t18 = t57 * qJD(4) + t130 * t73 - t218 * t139;
t200 = t126 * t18;
t199 = t126 * t64;
t198 = t130 * t75;
t183 = qJD(3) * t133;
t99 = -pkin(2) * t127 - t176;
t196 = t131 * t113 + t99 * t183;
t43 = -t218 * t74 - t198;
t195 = pkin(3) * t162 + qJD(5) - t43;
t194 = t127 * t131;
t191 = t132 * t133;
t135 = qJD(3) ^ 2;
t190 = t135 * t131;
t123 = t135 * t133;
t40 = t218 * t70 - t198;
t189 = qJD(5) - t40;
t187 = t185 * t156;
t186 = t128 - t129;
t182 = qJD(3) * t134;
t181 = qJD(4) * t130;
t180 = -qJD(1) - t127;
t179 = -qJD(2) + t127;
t178 = pkin(3) * t194;
t121 = t132 * t202;
t120 = pkin(3) * t184;
t36 = t79 ^ 2 - t222;
t125 = t127 ^ 2;
t171 = t131 * t125 * t133;
t167 = t127 * t183;
t166 = t131 * t182;
t158 = -t130 * t55 - t70 * t162 + t75 * t181 - t218 * t54;
t165 = -t144 * t158 + t40 * t63 - t41 * t64 + t214;
t33 = -t126 * pkin(4) + t189;
t34 = t126 * qJ(5) + t41;
t122 = t126 * qJD(5);
t9 = t122 - t158;
t164 = t144 * t9 - t33 * t63 - t34 * t64 + t214;
t161 = t206 * t126;
t157 = t185 * qJD(2);
t155 = t131 * t167;
t26 = pkin(4) * t64 + qJ(5) * t63 - qJD(5) * t92 + t120;
t152 = -t26 + t177;
t42 = -t130 * t74 + t174;
t151 = pkin(3) * t181 - t42;
t49 = pkin(4) * t81 + qJ(5) * t79;
t150 = -t144 * t48 + t64 * t79;
t147 = t126 * t40 + t158;
t143 = t146 * t47 - t17 * t79 + t18 * t81 - t57 * t48;
t142 = t120 - t177;
t141 = -t127 * t99 - t156;
t140 = -t132 * t194 + t133 * t182;
t58 = -pkin(4) * t144 - qJ(5) * t92 + t119;
t2 = -t144 * t47 - t48 * t92 + t63 * t79 - t64 * t81;
t138 = t47 * t145 + t206 * t81 - t207 * t79 - t72 * t48;
t137 = -t160 + (-qJD(4) + t126) * t41;
t136 = t42 * t126 + (-t174 + (-pkin(3) * t126 - t70) * t130) * qJD(4) - t160;
t118 = -pkin(2) - t217;
t117 = -t218 * pkin(3) - pkin(4);
t114 = pkin(3) * t130 + qJ(5);
t102 = t119 - t217;
t97 = -0.2e1 * t155;
t96 = 0.2e1 * t155;
t95 = t121 + t120;
t86 = t99 * t184;
t78 = -0.2e1 * t186 * t127 * qJD(3);
t59 = t63 * t126;
t53 = t58 - t217;
t44 = t49 + t178;
t31 = -t81 * t126 + t48;
t30 = -t188 + (-t172 + t79) * t126;
t23 = t121 + t26;
t14 = -t47 * t92 - t63 * t81;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 * t121 - t113, t180 * t175, 0, 0, t96, t78, t123, t97, -t190, 0, t118 * t168 - t116 * t123 + t86 + (t180 * t191 - t166) * t202, t116 * t190 + t118 * t167 - t140 * t202 + t196, t127 * t157 * t217 + t187, ((qJD(1) * t118 + t99) * t132 + (qJD(1) * t116 + t98) * t134 * t185) * t202, t14, t2, -t59, t150, -t199, 0, t102 * t48 + t79 * t95 - t200 + t205, -t102 * t47 + t81 * t95 - t201 + t204, t143 + t165, t102 * t85 - t158 * t57 + t17 * t41 - t18 * t40 + t83 * t95 - t216, t14, -t59, -t2, 0, t199, t150, t23 * t79 + t48 * t53 - t200 + t219, t143 + t164, -t23 * t81 + t47 * t53 + t201 + t220, t12 * t53 + t17 * t34 + t18 * t33 + t23 * t35 + t57 * t9 - t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * t177 - t113, t179 * t176, 0, 0, t96, t78, t123, t97, -t190, 0, -pkin(2) * t168 - pkin(7) * t123 + t86 + (t179 * t191 + t166) * t203, -pkin(2) * t167 + pkin(7) * t190 + t140 * t203 + t196, -t185 * t127 * t176 + t187, ((-pkin(2) * qJD(2) - t99) * t132 + (pkin(7) * t157 - t185 * t98) * t134) * t203, t14, t2, -t59, t150, -t199, 0, t119 * t48 + t142 * t79 - t161 + t205, -t119 * t47 + t142 * t81 + t204 - t223, t138 + t165, t119 * t85 + t142 * t83 - t158 * t72 - t206 * t40 + t207 * t41 - t215, t14, -t59, -t2, 0, t199, t150, -t152 * t79 + t48 * t58 - t161 + t219, t138 + t164, t152 * t81 + t47 * t58 + t220 + t223, t12 * t58 - t152 * t35 + t206 * t33 + t207 * t34 + t72 * t9 - t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t186 * t125, 0, t171, 0, 0, t141 * t131, t141 * t133, 0, 0, t211, -t36, t30, -t211, 0, 0, -t79 * t178 + t136 - t209, t43 * t126 + t210 + (-t126 * t162 - t81 * t194) * pkin(3) + t158, (t41 - t42) * t81 + (-t40 + t43) * t79 + (t218 * t47 - t130 * t48 + (t130 * t81 - t218 * t79) * qJD(4)) * pkin(3), t40 * t42 - t41 * t43 + (-t83 * t194 - t218 * t11 - t158 * t130 + (-t130 * t40 + t218 * t41) * qJD(4)) * pkin(3), t211, t30, t36, 0, t31, -t211, -t44 * t79 + t136 - t212, -t114 * t48 - t117 * t47 + (t151 + t34) * t81 + (t33 - t195) * t79, t195 * t126 + t44 * t81 - t213 + t9, t11 * t117 + t114 * t9 + t151 * t33 + t195 * t34 - t35 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, -t36, t30, -t211, 0, 0, t137 - t209, t147 + t210, 0, 0, t211, t30, t36, 0, t31, -t211, -t49 * t79 + t137 - t212, pkin(4) * t47 - qJ(5) * t48 + (t34 - t41) * t81 + (t33 - t189) * t79, t49 * t81 + 0.2e1 * t122 - t147 - t213, -pkin(4) * t11 + qJ(5) * t9 + t189 * t34 - t33 * t41 - t35 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t30, -t126 ^ 2 - t222, -t34 * t126 + t11 + t212;];
tauc_reg = t1;
