% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:17
% EndTime: 2019-12-31 21:14:24
% DurationCPUTime: 2.07s
% Computational Cost: add. (3570->262), mult. (9431->372), div. (0->0), fcn. (6982->8), ass. (0->170)
t127 = sin(pkin(9));
t195 = cos(pkin(9));
t132 = cos(qJ(3));
t133 = cos(qJ(2));
t186 = qJD(1) * t133;
t175 = t132 * t186;
t129 = sin(qJ(3));
t130 = sin(qJ(2));
t187 = qJD(1) * t130;
t176 = t129 * t187;
t85 = -t175 + t176;
t87 = -t129 * t186 - t132 * t187;
t170 = t127 * t87 - t195 * t85;
t218 = qJD(5) - t170;
t228 = qJD(5) - t218;
t124 = qJD(2) + qJD(3);
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t147 = -t127 * t85 - t195 * t87;
t53 = -t131 * t124 + t128 * t147;
t227 = t218 * t53;
t162 = t131 * t218;
t99 = t129 * t133 + t132 * t130;
t221 = qJD(1) * t99;
t138 = t124 * t221;
t182 = qJD(1) * qJD(2);
t172 = t133 * t182;
t66 = qJD(3) * t175 - t124 * t176 + t132 * t172;
t40 = t127 * t66 + t195 * t138;
t200 = t128 * t40;
t226 = -t218 * t162 - t200;
t225 = -0.2e1 * t182;
t122 = pkin(2) * t187;
t215 = t87 * pkin(3);
t33 = pkin(4) * t147 - pkin(8) * t170 - t215;
t120 = t132 * pkin(2) + pkin(3);
t168 = t195 * t129;
t81 = pkin(2) * t168 + t127 * t120;
t76 = pkin(8) + t81;
t224 = (qJD(5) * t76 + t122 + t33) * t218;
t216 = pkin(6) + pkin(7);
t107 = t216 * t133;
t102 = qJD(1) * t107;
t88 = t129 * t102;
t106 = t216 * t130;
t100 = qJD(1) * t106;
t203 = qJD(2) * pkin(2);
t94 = -t100 + t203;
t163 = t132 * t94 - t88;
t82 = t87 * qJ(4);
t51 = t163 + t82;
t184 = qJD(5) * t128;
t222 = -t131 * t40 + t184 * t218;
t177 = qJD(2) * t216;
t159 = qJD(1) * t177;
t95 = t130 * t159;
t220 = (qJD(3) * t94 - t95) * t132;
t101 = t130 * t177;
t103 = t133 * t177;
t152 = t129 * t106 - t132 * t107;
t141 = t152 * qJD(3) + t129 * t101 - t132 * t103;
t98 = t129 * t130 - t132 * t133;
t70 = t124 * t98;
t139 = t70 * qJ(4) - t99 * qJD(4) + t141;
t185 = qJD(3) * t129;
t193 = t132 * t106;
t146 = -qJD(3) * t193 - t132 * t101 - t129 * t103 - t107 * t185;
t144 = t99 * qJD(3);
t71 = t99 * qJD(2) + t144;
t29 = -t71 * qJ(4) - t98 * qJD(4) + t146;
t11 = t127 * t139 + t195 * t29;
t143 = -t99 * qJ(4) - t129 * t107 - t193;
t60 = -t98 * qJ(4) - t152;
t36 = t127 * t143 + t195 * t60;
t92 = t132 * t102;
t153 = -t129 * t94 - t92;
t96 = t133 * t159;
t154 = t129 * t95 - t132 * t96;
t140 = t153 * qJD(3) + t154;
t137 = -t66 * qJ(4) + t87 * qJD(4) + t140;
t169 = -t102 * t185 - t129 * t96;
t18 = -qJ(4) * t138 - t85 * qJD(4) + t169 + t220;
t6 = t127 * t18 - t195 * t137;
t69 = -t127 * t98 + t195 * t99;
t158 = -t36 * t40 + t6 * t69;
t196 = t85 * qJ(4);
t52 = -t153 - t196;
t201 = t127 * t52;
t45 = t124 * pkin(3) + t51;
t23 = t195 * t45 - t201;
t21 = -t124 * pkin(4) - t23;
t121 = -t133 * pkin(2) - pkin(1);
t105 = t121 * qJD(1);
t72 = t85 * pkin(3) + qJD(4) + t105;
t28 = -pkin(4) * t170 - pkin(8) * t147 + t72;
t151 = t98 * pkin(3) + t121;
t68 = t127 * t99 + t195 * t98;
t34 = t68 * pkin(4) - t69 * pkin(8) + t151;
t43 = -t127 * t71 - t195 * t70;
t7 = t127 * t137 + t195 * t18;
t217 = -(qJD(5) * t34 + t11) * t218 - (qJD(5) * t28 + t7) * t68 + t21 * t43 + t158;
t214 = t21 * t170;
t213 = t21 * t69;
t212 = t34 * t40;
t55 = t128 * t124 + t131 * t147;
t211 = t55 * t147;
t210 = t218 * t147;
t209 = t147 * t53;
t208 = t87 * t85;
t48 = t195 * t52;
t24 = t127 * t45 + t48;
t164 = t129 * t100 - t92;
t148 = t164 + t196;
t204 = pkin(2) * qJD(3);
t205 = -t132 * t100 - t88;
t56 = t82 + t205;
t207 = -t127 * t56 + t195 * t148 + (t127 * t132 + t168) * t204;
t194 = t127 * t129;
t206 = -t127 * t148 - t195 * t56 + (t195 * t132 - t194) * t204;
t202 = t105 * t87;
t41 = -t127 * t138 + t195 * t66;
t199 = t128 * t41;
t198 = t128 * t170;
t183 = qJD(5) * t131;
t15 = t124 * t183 + t131 * t41 - t147 * t184;
t197 = t15 * t128;
t135 = qJD(1) ^ 2;
t192 = t133 * t135;
t134 = qJD(2) ^ 2;
t191 = t134 * t130;
t190 = t134 * t133;
t188 = t130 ^ 2 - t133 ^ 2;
t123 = t130 * t203;
t180 = t69 * t184;
t179 = t218 * t183;
t22 = t124 * pkin(8) + t24;
t155 = t128 * t22 - t131 * t28;
t178 = t147 * t155 + t21 * t184;
t174 = -pkin(2) * t124 - t94;
t173 = t71 * pkin(3) + t123;
t161 = pkin(1) * t225;
t9 = t128 * t28 + t131 * t22;
t160 = t6 * t128 + t147 * t9 + t21 * t183;
t157 = t147 * t24 + t170 * t23;
t156 = t218 * t43 + t40 * t69;
t150 = t198 * t218 - t222;
t149 = t105 * t85 - t169;
t80 = -pkin(2) * t194 + t195 * t120;
t142 = -t206 * t218 - t76 * t40 - t214;
t136 = pkin(3) * t138 + qJD(2) * t122;
t118 = -t195 * pkin(3) - pkin(4);
t117 = t127 * pkin(3) + pkin(8);
t75 = -pkin(4) - t80;
t57 = -t85 ^ 2 + t87 ^ 2;
t47 = -t87 * t124 - t138;
t46 = t85 * t124 + t66;
t42 = -t127 * t70 + t195 * t71;
t35 = t127 * t60 - t195 * t143;
t26 = t195 * t51 - t201;
t25 = t127 * t51 + t48;
t16 = t55 * qJD(5) + t199;
t14 = t42 * pkin(4) - t43 * pkin(8) + t173;
t13 = t40 * pkin(4) - t41 * pkin(8) + t136;
t12 = t131 * t13;
t10 = t127 * t29 - t195 * t139;
t4 = t55 * t162 + t197;
t3 = -t211 - t226;
t2 = t150 + t209;
t1 = (t15 - t227) * t131 + (-t218 * t55 - t16) * t128;
t5 = [0, 0, 0, 0.2e1 * t130 * t172, t188 * t225, t190, -t191, 0, -pkin(6) * t190 + t130 * t161, pkin(6) * t191 + t133 * t161, t66 * t99 + t87 * t70, -t99 * t138 - t66 * t98 + t70 * t85 + t87 * t71, -t70 * t124, -t71 * t124, 0, t85 * t123 + t105 * t71 + t141 * t124 + (t121 * t144 + (t130 * pkin(2) * t98 + t121 * t99) * qJD(2)) * qJD(1), t121 * t66 - t105 * t70 - t146 * t124 + (-t87 + t221) * t123, t10 * t147 + t11 * t170 - t23 * t43 - t24 * t42 + t35 * t41 - t7 * t68 + t158, -t23 * t10 + t24 * t11 + t136 * t151 + t72 * t173 + t6 * t35 + t7 * t36, -t55 * t180 + (t15 * t69 + t43 * t55) * t131, (-t128 * t55 - t131 * t53) * t43 + (-t197 - t131 * t16 + (t128 * t53 - t131 * t55) * qJD(5)) * t69, t156 * t131 + t15 * t68 - t180 * t218 + t55 * t42, -t156 * t128 - t16 * t68 - t69 * t179 - t53 * t42, t218 * t42 + t40 * t68, t10 * t53 + t12 * t68 + t35 * t16 - t155 * t42 + (t14 * t218 + t212 + (-t218 * t36 - t22 * t68 + t213) * qJD(5)) * t131 + t217 * t128, t10 * t55 + t35 * t15 - t9 * t42 + (-(-qJD(5) * t36 + t14) * t218 - t212 - (-qJD(5) * t22 + t13) * t68 - qJD(5) * t213) * t128 + t217 * t131; 0, 0, 0, -t130 * t192, t188 * t135, 0, 0, 0, t135 * pkin(1) * t130, pkin(1) * t192, -t208, t57, t46, t47, 0, -t85 * t122 + t202 - t164 * t124 + (t174 * t129 - t92) * qJD(3) + t154, t87 * t122 + t205 * t124 + (t174 * qJD(3) + t95) * t132 + t149, t147 * t207 + t170 * t206 - t81 * t40 - t80 * t41 + t157, t7 * t81 - t6 * t80 - t72 * (t122 - t215) + t206 * t24 - t207 * t23, t4, t1, t3, t2, -t210, t75 * t16 + t207 * t53 + (-t6 - t224) * t131 + t142 * t128 + t178, t128 * t224 + t142 * t131 + t75 * t15 + t207 * t55 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t57, t46, t47, 0, -t153 * t124 + t140 + t202, t163 * t124 + t149 - t220, -t25 * t147 - t26 * t170 + (-t127 * t40 - t195 * t41) * pkin(3) + t157, t23 * t25 - t24 * t26 + (t127 * t7 - t195 * t6 + t72 * t87) * pkin(3), t4, t1, t3, t2, -t210, t118 * t16 - t6 * t131 - (-t128 * t26 + t131 * t33) * t218 - t25 * t53 - t21 * t198 + (-t179 - t200) * t117 + t178, t118 * t15 + (t128 * t33 + t131 * t26) * t218 - t25 * t55 - t131 * t214 + t222 * t117 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147 ^ 2 - t170 ^ 2, t147 * t23 - t170 * t24 + t136, 0, 0, 0, 0, 0, t150 - t209, -t211 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t53, -t53 ^ 2 + t55 ^ 2, t15 + t227, -t228 * t55 - t199, t40, -t128 * t7 - t21 * t55 - t228 * t9 + t12, -t128 * t13 - t131 * t7 + t228 * t155 + t21 * t53;];
tauc_reg = t5;
