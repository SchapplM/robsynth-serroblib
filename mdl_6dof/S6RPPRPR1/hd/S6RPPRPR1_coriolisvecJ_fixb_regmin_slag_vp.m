% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:02
% EndTime: 2019-03-09 01:40:08
% DurationCPUTime: 1.81s
% Computational Cost: add. (2950->258), mult. (7339->373), div. (0->0), fcn. (5640->10), ass. (0->141)
t135 = cos(pkin(10));
t197 = cos(qJ(4));
t164 = qJD(1) * t197;
t121 = t135 * t164;
t132 = sin(pkin(10));
t138 = sin(qJ(4));
t176 = t138 * t132;
t165 = qJD(1) * t176;
t100 = -t121 + t165;
t97 = qJD(6) + t100;
t122 = sin(pkin(9)) * pkin(1) + qJ(3);
t117 = t122 * qJD(1);
t128 = t135 * qJD(2);
t190 = pkin(7) * qJD(1);
t87 = t128 + (-t117 - t190) * t132;
t92 = t132 * qJD(2) + t135 * t117;
t88 = t135 * t190 + t92;
t51 = t138 * t87 + t197 * t88;
t202 = qJD(4) * t51;
t201 = qJD(6) - t97;
t139 = cos(qJ(6));
t175 = t138 * t135;
t113 = t197 * t132 + t175;
t102 = t113 * qJD(1);
t131 = sin(pkin(11));
t134 = cos(pkin(11));
t83 = -t134 * qJD(4) + t102 * t131;
t200 = t139 * t83;
t137 = sin(qJ(6));
t85 = qJD(4) * t131 + t102 * t134;
t151 = t137 * t83 - t139 * t85;
t199 = t151 * t97;
t112 = t131 * t139 + t134 * t137;
t182 = t97 * t112;
t50 = -t138 * t88 + t197 * t87;
t195 = pkin(7) + t122;
t107 = t195 * t132;
t108 = t195 * t135;
t71 = t197 * t107 + t138 * t108;
t174 = t139 * t134;
t177 = t131 * t137;
t110 = -t174 + t177;
t183 = t97 * t110;
t106 = t113 * qJD(4);
t94 = qJD(1) * t106;
t198 = -t112 * t94 + t183 * t97;
t120 = qJD(4) * t121;
t93 = -qJD(4) * t165 + t120;
t21 = -t151 * qJD(6) + t112 * t93;
t99 = t100 ^ 2;
t196 = pkin(8) * t134;
t194 = pkin(8) + qJ(5);
t147 = t197 * t135 - t176;
t185 = t131 * t93;
t150 = -qJD(6) * t85 - t185;
t171 = qJD(6) * t139;
t191 = -t83 * t171 + t93 * t174;
t20 = t150 * t137 + t191;
t193 = -t106 * t151 - t147 * t20;
t105 = t147 * qJD(4);
t172 = qJD(6) * t113;
t178 = t113 * t134;
t29 = t112 * t105 + t171 * t178 - t172 * t177;
t69 = t112 * t113;
t192 = -t29 * t97 - t69 * t94;
t142 = t147 * qJD(3);
t33 = qJD(1) * t142 + (qJD(5) + t50) * qJD(4);
t44 = pkin(4) * t94 - qJ(5) * t93 - qJD(5) * t102;
t9 = t131 * t44 + t134 * t33;
t43 = qJD(4) * qJ(5) + t51;
t114 = -cos(pkin(9)) * pkin(1) - pkin(3) * t135 - pkin(2);
t98 = t114 * qJD(1) + qJD(3);
t59 = pkin(4) * t100 - qJ(5) * t102 + t98;
t16 = t131 * t59 + t134 * t43;
t73 = pkin(4) * t102 + qJ(5) * t100;
t26 = t131 * t73 + t134 * t50;
t52 = -t71 * qJD(4) + t142;
t60 = pkin(4) * t106 - qJ(5) * t105 - qJD(5) * t113;
t19 = t131 * t60 + t134 * t52;
t66 = -pkin(4) * t147 - qJ(5) * t113 + t114;
t72 = -t138 * t107 + t197 * t108;
t31 = t131 * t66 + t134 * t72;
t47 = t137 * t85 + t200;
t189 = t102 * t47;
t188 = t102 * t151;
t187 = t147 * t93;
t184 = t134 * t93;
t181 = t100 * t131;
t180 = t105 * t131;
t179 = t113 * t131;
t173 = t132 ^ 2 + t135 ^ 2;
t170 = t105 * qJD(4);
t169 = qJD(1) * qJD(3);
t8 = -t131 * t33 + t134 * t44;
t4 = pkin(5) * t94 - pkin(8) * t184 + t8;
t5 = -pkin(8) * t185 + t9;
t167 = -t137 * t5 + t139 * t4;
t15 = -t131 * t43 + t134 * t59;
t25 = -t131 * t50 + t134 * t73;
t18 = -t131 * t52 + t134 * t60;
t30 = -t131 * t72 + t134 * t66;
t163 = qJD(1) * t173;
t34 = t132 * qJD(3) * t164 + t169 * t175 + t202;
t162 = -t110 * t94 - t182 * t97;
t161 = -t131 * t8 + t134 * t9;
t160 = t137 * t4 + t139 * t5;
t28 = -t110 * t105 - t112 * t172;
t70 = t110 * t113;
t159 = -t28 * t97 + t70 * t94;
t10 = -pkin(8) * t83 + t16;
t6 = pkin(5) * t100 - pkin(8) * t85 + t15;
t2 = t10 * t139 + t137 * t6;
t158 = t10 * t137 - t139 * t6;
t157 = -t106 * t47 + t147 * t21;
t156 = -t131 * t15 + t134 * t16;
t155 = t131 * t85 - t134 * t83;
t154 = t132 * (-t117 * t132 + t128) - t135 * t92;
t17 = -pkin(5) * t147 - pkin(8) * t178 + t30;
t23 = -pkin(8) * t179 + t31;
t153 = -t137 * t23 + t139 * t17;
t152 = t137 * t17 + t139 * t23;
t149 = -t100 * t105 - t113 * t94;
t116 = t194 * t134;
t146 = pkin(5) * t102 + qJD(5) * t131 + qJD(6) * t116 + t100 * t196 + t25;
t115 = t194 * t131;
t145 = pkin(8) * t181 - qJD(5) * t134 + qJD(6) * t115 + t26;
t41 = -qJD(4) * pkin(4) + qJD(5) - t50;
t144 = t105 * t41 + t113 * t34 + t71 * t93;
t143 = t149 - t187;
t141 = -pkin(4) * t93 - qJ(5) * t94 + (-qJD(5) + t41) * t100;
t53 = t113 * qJD(3) + t72 * qJD(4);
t125 = -pkin(5) * t134 - pkin(4);
t96 = t106 * qJD(4);
t56 = pkin(5) * t179 + t71;
t36 = pkin(5) * t180 + t53;
t35 = -pkin(5) * t181 + t51;
t27 = t83 * pkin(5) + t41;
t24 = pkin(5) * t185 + t34;
t12 = -pkin(8) * t180 + t19;
t7 = pkin(5) * t106 - t105 * t196 + t18;
t1 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t163 (t122 * t163 - t154) * qJD(3), t102 * t105 + t113 * t93, -t102 * t106 + t149 + t187, t170, -t96, 0, -qJD(4) * t53 + t106 * t98 + t114 * t94, -qJD(4) * t52 + t105 * t98 + t114 * t93, t100 * t18 + t106 * t15 + t144 * t131 - t147 * t8 + t30 * t94 + t53 * t83, -t100 * t19 - t106 * t16 + t144 * t134 + t147 * t9 - t31 * t94 + t53 * t85, -t18 * t85 - t19 * t83 + (-t105 * t15 - t113 * t8 - t30 * t93) * t134 + (-t105 * t16 - t113 * t9 - t31 * t93) * t131, t15 * t18 + t16 * t19 + t30 * t8 + t31 * t9 + t34 * t71 + t41 * t53, -t151 * t28 - t20 * t70, t151 * t29 - t20 * t69 + t21 * t70 - t28 * t47, -t159 + t193, t157 + t192, t106 * t97 - t147 * t94 (-t12 * t137 + t139 * t7) * t97 + t153 * t94 - t167 * t147 - t158 * t106 + t36 * t47 + t56 * t21 + t24 * t69 + t27 * t29 + (t147 * t2 - t152 * t97) * qJD(6) -(t12 * t139 + t137 * t7) * t97 - t152 * t94 + t160 * t147 - t2 * t106 - t36 * t151 + t56 * t20 - t24 * t70 + t27 * t28 + (-t147 * t158 - t153 * t97) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t170, t106 * t83 + t143 * t131, t106 * t85 + t143 * t134, t155 * t105, t156 * t105 + t106 * t41 + t161 * t113 - t147 * t34, 0, 0, 0, 0, 0, -t157 + t192, t159 + t193; 0, 0, 0, 0, 0, 0, -t173 * qJD(1) ^ 2, t154 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t102 * qJD(4), t120 + (-t100 - t165) * qJD(4), -t102 * t83 - t99 * t131 + t134 * t94, -t102 * t85 - t131 * t94 - t134 * t99 (-t131 ^ 2 - t134 ^ 2) * t93 + t155 * t100, t156 * t100 - t102 * t41 + t131 * t9 + t134 * t8, 0, 0, 0, 0, 0, t162 - t189, t188 + t198; 0, 0, 0, 0, 0, 0, 0, 0, t102 * t100, t102 ^ 2 - t99, t120 + (t100 - t165) * qJD(4), 0, 0, -t102 * t98 + t202 - t34, t98 * t100 - t147 * t169, -t100 * t25 - t102 * t15 + t141 * t131 - t134 * t34 - t51 * t83, t100 * t26 + t102 * t16 + t131 * t34 + t141 * t134 - t51 * t85, t25 * t85 + t26 * t83 + (-qJD(5) * t83 - t100 * t15 + t9) * t134 + (qJD(5) * t85 - t100 * t16 - t8) * t131, -pkin(4) * t34 + t161 * qJ(5) + t156 * qJD(5) - t15 * t25 - t16 * t26 - t41 * t51, t112 * t20 + t151 * t183, -t110 * t20 - t112 * t21 + t151 * t182 + t183 * t47, t188 - t198, t162 + t189, -t97 * t102 (-t115 * t139 - t116 * t137) * t94 + t125 * t21 + t24 * t110 + t158 * t102 - t35 * t47 + (t137 * t145 - t139 * t146) * t97 + t182 * t27 -(-t115 * t137 + t116 * t139) * t94 + t125 * t20 + t24 * t112 + t2 * t102 + t35 * t151 + (t137 * t146 + t139 * t145) * t97 - t183 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t85 + t185, -t100 * t83 + t184, -t83 ^ 2 - t85 ^ 2, t15 * t85 + t16 * t83 + t34, 0, 0, 0, 0, 0, t21 - t199, -t97 * t200 + (-t85 * t97 + t150) * t137 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151 * t47, t151 ^ 2 - t47 ^ 2, t47 * t97 + t20, -t21 - t199, t94, t27 * t151 - t201 * t2 + t167, t201 * t158 + t27 * t47 - t160;];
tauc_reg  = t1;
