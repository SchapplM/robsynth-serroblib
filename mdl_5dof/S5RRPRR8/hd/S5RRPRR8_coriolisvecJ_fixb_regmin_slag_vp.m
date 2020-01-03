% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:20
% DurationCPUTime: 2.13s
% Computational Cost: add. (3139->242), mult. (8302->345), div. (0->0), fcn. (6437->8), ass. (0->156)
t140 = cos(qJ(5));
t175 = qJD(5) * t140;
t135 = sin(pkin(9));
t136 = cos(pkin(9));
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t116 = t135 * t142 + t136 * t139;
t107 = t116 * qJD(1);
t138 = sin(qJ(4));
t183 = t136 * t142;
t115 = -t135 * t139 + t183;
t105 = t115 * qJD(1);
t141 = cos(qJ(4));
t93 = t141 * t105;
t64 = -t138 * t107 + t93;
t223 = t140 * t64;
t230 = t175 - t223;
t132 = qJD(2) + qJD(4);
t186 = t64 * t132;
t177 = qJD(4) * t138;
t106 = t116 * qJD(2);
t96 = qJD(1) * t106;
t173 = qJD(1) * qJD(2);
t168 = t139 * t173;
t97 = -t135 * t168 + t173 * t183;
t32 = qJD(4) * t93 - t107 * t177 - t138 * t96 + t141 * t97;
t229 = t32 - t186;
t179 = -qJD(5) + t64;
t228 = qJD(5) + t179;
t137 = sin(qJ(5));
t154 = t138 * t105 + t141 * t107;
t176 = qJD(5) * t137;
t14 = t132 * t175 + t140 * t32 - t154 * t176;
t52 = t137 * t132 + t140 * t154;
t15 = qJD(5) * t52 + t137 * t32;
t49 = -t140 * t132 + t137 * t154;
t227 = -t137 * t15 + t14 * t140 - t230 * t49;
t12 = t14 * t137;
t226 = t230 * t52 + t12;
t33 = t154 * qJD(4) + t138 * t97 + t141 * t96;
t29 = t137 * t33;
t55 = t179 * t175;
t195 = t29 - t55;
t200 = t52 * t154;
t225 = t179 * t223 + t195 - t200;
t204 = t107 * pkin(7);
t198 = -qJ(3) - pkin(6);
t124 = t198 * t142;
t121 = qJD(1) * t124;
t110 = t135 * t121;
t123 = t198 * t139;
t120 = qJD(1) * t123;
t193 = qJD(2) * pkin(2);
t114 = t120 + t193;
t66 = t136 * t114 + t110;
t44 = qJD(2) * pkin(3) - t204 + t66;
t205 = t105 * pkin(7);
t184 = t136 * t121;
t67 = t135 * t114 - t184;
t47 = t67 + t205;
t20 = -t138 * t47 + t141 * t44;
t18 = -t132 * pkin(4) - t20;
t224 = t18 * t64;
t222 = t154 * t64;
t221 = t137 * t179;
t187 = t154 * t132;
t220 = -t33 + t187;
t218 = t154 ^ 2 - t64 ^ 2;
t35 = pkin(4) * t154 - t64 * pkin(8);
t162 = qJD(2) * t198;
t102 = t142 * qJD(3) + t139 * t162;
t84 = t102 * qJD(1);
t103 = -t139 * qJD(3) + t142 * t162;
t85 = t103 * qJD(1);
t48 = -t135 * t84 + t136 * t85;
t39 = -t97 * pkin(7) + t48;
t51 = t135 * t85 + t136 * t84;
t40 = -pkin(7) * t96 + t51;
t2 = (qJD(4) * t44 + t40) * t141 + t138 * t39 - t47 * t177;
t169 = -t142 * pkin(2) - pkin(1);
t156 = t169 * qJD(1);
t122 = qJD(3) + t156;
t72 = -t105 * pkin(3) + t122;
t217 = -t72 * t64 - t2;
t215 = -0.2e1 * t173;
t199 = t154 * t49;
t213 = t179 * t154;
t31 = t140 * t33;
t212 = -t176 * t179 - t31;
t21 = t138 * t44 + t141 * t47;
t19 = t132 * pkin(8) + t21;
t22 = -pkin(4) * t64 - pkin(8) * t154 + t72;
t155 = t137 * t19 - t140 * t22;
t211 = t154 * t155 + t18 * t176;
t3 = t21 * qJD(4) + t138 * t40 - t141 * t39;
t5 = t137 * t22 + t140 * t19;
t210 = t3 * t137 + t5 * t154 + t18 * t175;
t209 = -t154 * t72 - t3;
t153 = t141 * t115 - t138 * t116;
t69 = t138 * t115 + t141 * t116;
t87 = -t115 * pkin(3) + t169;
t26 = -pkin(4) * t153 - t69 * pkin(8) + t87;
t74 = t136 * t123 + t135 * t124;
t58 = -t116 * pkin(7) + t74;
t75 = t135 * t123 - t136 * t124;
t59 = t115 * pkin(7) + t75;
t28 = t138 * t58 + t141 * t59;
t109 = t115 * qJD(2);
t36 = t153 * qJD(4) - t138 * t106 + t141 * t109;
t27 = t138 * t59 - t141 * t58;
t56 = -t135 * t102 + t136 * t103;
t41 = -t109 * pkin(7) + t56;
t57 = t136 * t102 + t135 * t103;
t42 = -t106 * pkin(7) + t57;
t6 = -t27 * qJD(4) + t138 * t41 + t141 * t42;
t208 = (qJD(5) * t26 + t6) * t179 + (qJD(5) * t22 + t2) * t153 + t18 * t36 - t28 * t33 + t3 * t69;
t206 = pkin(2) * t135;
t203 = t18 * t69;
t202 = t26 * t33;
t201 = t33 * t69;
t129 = t136 * pkin(2) + pkin(3);
t150 = t141 * t129 - t138 * t206;
t70 = -t135 * t120 + t184;
t53 = t70 - t205;
t71 = t136 * t120 + t110;
t54 = t71 - t204;
t196 = -t150 * qJD(4) + t138 * t53 + t141 * t54;
t151 = t138 * t129 + t141 * t206;
t194 = t151 * qJD(4) - t138 * t54 + t141 * t53;
t191 = t137 * t52;
t144 = qJD(1) ^ 2;
t182 = t142 * t144;
t143 = qJD(2) ^ 2;
t181 = t143 * t139;
t180 = t143 * t142;
t178 = t139 ^ 2 - t142 ^ 2;
t174 = t139 * qJD(1);
t131 = t139 * t193;
t170 = t69 * t176;
t127 = pkin(2) * t168;
t73 = t96 * pkin(3) + t127;
t80 = t106 * pkin(3) + t131;
t79 = pkin(2) * t174 + t107 * pkin(3);
t159 = pkin(1) * t215;
t101 = pkin(8) + t151;
t158 = qJD(5) * t101 + t35 + t79;
t157 = -t179 * t36 + t201;
t152 = -t221 * t64 - t212;
t148 = -t101 * t33 - t179 * t196 - t224;
t100 = -pkin(4) - t150;
t37 = t69 * qJD(4) + t141 * t106 + t138 * t109;
t10 = t37 * pkin(4) - t36 * pkin(8) + t80;
t9 = t33 * pkin(4) - t32 * pkin(8) + t73;
t8 = t140 * t9;
t7 = t28 * qJD(4) + t138 * t42 - t141 * t41;
t1 = [0, 0, 0, 0.2e1 * t142 * t168, t178 * t215, t180, -t181, 0, -pkin(6) * t180 + t139 * t159, pkin(6) * t181 + t142 * t159, t57 * t105 - t67 * t106 - t56 * t107 - t66 * t109 + t51 * t115 - t48 * t116 - t74 * t97 - t75 * t96, t48 * t74 + t51 * t75 + t66 * t56 + t67 * t57 + (t122 + t156) * t131, t154 * t36 + t32 * t69, t153 * t32 - t154 * t37 + t36 * t64 - t201, t36 * t132, -t37 * t132, 0, -t7 * t132 - t153 * t73 + t87 * t33 + t72 * t37 - t64 * t80, -t6 * t132 + t154 * t80 + t87 * t32 + t72 * t36 + t73 * t69, -t52 * t170 + (t14 * t69 + t36 * t52) * t140, (-t140 * t49 - t191) * t36 + (-t12 - t140 * t15 + (t137 * t49 - t140 * t52) * qJD(5)) * t69, -t14 * t153 + t157 * t140 + t170 * t179 + t52 * t37, -t157 * t137 + t15 * t153 - t49 * t37 + t69 * t55, -t153 * t33 - t179 * t37, t27 * t15 - t155 * t37 + t7 * t49 - t8 * t153 + (-t10 * t179 + t202 + (t153 * t19 + t179 * t28 + t203) * qJD(5)) * t140 + t208 * t137, t27 * t14 - t5 * t37 + t7 * t52 + ((-qJD(5) * t28 + t10) * t179 - t202 + (-qJD(5) * t19 + t9) * t153 - qJD(5) * t203) * t137 + t208 * t140; 0, 0, 0, -t139 * t182, t178 * t144, 0, 0, 0, t144 * pkin(1) * t139, pkin(1) * t182, (t67 + t70) * t107 + (t66 - t71) * t105 + (-t135 * t96 - t136 * t97) * pkin(2), -t66 * t70 - t67 * t71 + (-t122 * t174 + t135 * t51 + t136 * t48) * pkin(2), -t222, t218, t229, t220, 0, -t194 * t132 + t64 * t79 + t209, t196 * t132 - t154 * t79 + t217, t226, t179 * t191 + t227, t225, t152 + t199, t213, t100 * t15 + t194 * t49 + (t158 * t179 - t3) * t140 + t148 * t137 + t211, t100 * t14 + t148 * t140 - t158 * t221 + t194 * t52 + t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 ^ 2 - t107 ^ 2, -t67 * t105 + t66 * t107 + t127, 0, 0, 0, 0, 0, t33 + t187, t32 + t186, 0, 0, 0, 0, 0, t152 - t199, -t140 * t179 ^ 2 - t200 - t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, t218, t229, t220, 0, t21 * t132 + t209, t20 * t132 + t217, t226, t221 * t52 + t227, t225, -t179 * t221 + t199 + t31, t213, -pkin(4) * t15 - t3 * t140 + (-t137 * t20 + t140 * t35) * t179 - t21 * t49 - t137 * t224 - t195 * pkin(8) + t211, -pkin(4) * t14 - (t137 * t35 + t140 * t20) * t179 - t21 * t52 - t18 * t223 + t212 * pkin(8) + t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t49, -t49 ^ 2 + t52 ^ 2, -t179 * t49 + t14, -t179 * t52 - t15, t33, -t137 * t2 - t18 * t52 - t228 * t5 + t8, -t137 * t9 - t140 * t2 + t228 * t155 + t18 * t49;];
tauc_reg = t1;
