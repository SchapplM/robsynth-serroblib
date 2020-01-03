% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:15
% EndTime: 2019-12-31 19:39:20
% DurationCPUTime: 1.61s
% Computational Cost: add. (2468->261), mult. (5952->363), div. (0->0), fcn. (3956->6), ass. (0->153)
t156 = sin(qJ(5));
t158 = cos(qJ(5));
t153 = sin(pkin(8));
t154 = cos(pkin(8));
t157 = sin(qJ(2));
t159 = cos(qJ(2));
t106 = t153 * t157 + t154 * t159;
t89 = t106 * qJD(1);
t186 = qJD(1) * t159;
t178 = t153 * t186;
t187 = qJD(1) * t157;
t92 = t154 * t187 - t178;
t46 = -t156 * t89 + t158 * t92;
t215 = t46 ^ 2;
t207 = t156 * t92 + t158 * t89;
t214 = t207 ^ 2;
t148 = qJD(2) - qJD(5);
t213 = t46 * t148;
t200 = t46 * t207;
t212 = t148 * t207;
t211 = t214 - t215;
t149 = qJD(2) * qJD(3);
t185 = qJD(2) * t157;
t198 = pkin(6) - qJ(4);
t85 = -t159 * qJD(4) - t198 * t185;
t61 = t85 * qJD(1) + t149;
t180 = qJD(1) * qJD(2);
t176 = t159 * t180;
t133 = pkin(6) * t176;
t181 = t157 * qJD(4);
t184 = qJD(2) * t159;
t69 = t133 + (-qJ(4) * t184 - t181) * qJD(1);
t174 = t153 * t61 - t154 * t69;
t98 = t106 * qJD(2);
t82 = qJD(1) * t98;
t18 = -pkin(7) * t82 - t174;
t26 = t153 * t69 + t154 * t61;
t177 = t157 * t180;
t128 = t154 * t177;
t81 = t153 * t176 - t128;
t19 = -pkin(7) * t81 + t26;
t203 = pkin(7) * t92;
t140 = pkin(6) * t186;
t116 = -qJ(4) * t186 + t140;
t150 = qJD(2) * qJ(3);
t102 = t116 + t150;
t139 = pkin(6) * t187;
t114 = qJ(4) * t187 - t139;
t160 = -pkin(2) - pkin(3);
t179 = t160 * qJD(2);
t75 = qJD(3) + t179 - t114;
t32 = -t102 * t153 + t154 * t75;
t21 = -qJD(2) * pkin(4) - t203 + t32;
t204 = pkin(7) * t89;
t33 = t154 * t102 + t153 * t75;
t22 = t33 - t204;
t6 = t156 * t21 + t158 * t22;
t2 = -t6 * qJD(5) - t156 * t19 + t158 * t18;
t103 = -qJD(1) * pkin(1) - pkin(2) * t186 - qJ(3) * t187;
t65 = pkin(3) * t186 + qJD(4) - t103;
t39 = pkin(4) * t89 + t65;
t210 = -t39 * t46 + t2;
t182 = qJD(5) * t158;
t183 = qJD(5) * t156;
t1 = t156 * t18 + t158 * t19 + t21 * t182 - t22 * t183;
t209 = t207 * t39 - t1;
t11 = t156 * t81 - t158 * t82 + t89 * t182 + t92 * t183;
t208 = t11 + t212;
t206 = t89 ^ 2;
t205 = t92 ^ 2;
t105 = -t153 * t156 + t154 * t158;
t50 = -t114 * t153 + t154 * t116;
t27 = t50 - t204;
t51 = t154 * t114 + t153 * t116;
t28 = t51 + t203;
t117 = -qJ(3) * t153 + t154 * t160;
t113 = -pkin(4) + t117;
t118 = qJ(3) * t154 + t153 * t160;
t52 = t113 * t158 - t118 * t156;
t202 = t105 * qJD(3) + t52 * qJD(5) - t156 * t27 - t158 * t28;
t108 = t153 * t158 + t154 * t156;
t53 = t113 * t156 + t118 * t158;
t201 = -t108 * qJD(3) - t53 * qJD(5) + t156 * t28 - t158 * t27;
t199 = t92 * t89;
t125 = t198 * t159;
t86 = qJD(2) * t125 - t181;
t35 = t153 * t86 + t154 * t85;
t197 = t148 * t108;
t196 = t148 * t105;
t195 = qJD(2) * pkin(2);
t162 = qJD(1) ^ 2;
t192 = t159 * t162;
t161 = qJD(2) ^ 2;
t191 = t161 * t157;
t144 = t161 * t159;
t124 = t198 * t157;
t58 = t153 * t124 + t154 * t125;
t143 = t157 * qJD(3);
t190 = qJ(3) * t176 + qJD(1) * t143;
t189 = qJ(3) * t184 + t143;
t151 = t157 ^ 2;
t188 = t159 ^ 2 - t151;
t121 = -t159 * pkin(2) - t157 * qJ(3) - pkin(1);
t34 = -t153 * t85 + t154 * t86;
t173 = t156 * t82 + t158 * t81;
t172 = -0.2e1 * pkin(1) * t180;
t171 = qJD(3) - t195;
t170 = qJD(3) * t153 + t50;
t169 = qJD(3) * t154 - t51;
t57 = t154 * t124 - t125 * t153;
t104 = t159 * pkin(3) - t121;
t168 = qJD(1) * t121 + t103;
t167 = t157 * t179;
t166 = t157 * t176;
t165 = t153 * t32 - t154 * t33;
t107 = -t153 * t159 + t154 * t157;
t36 = -pkin(7) * t107 + t57;
t37 = -pkin(7) * t106 + t58;
t9 = -t156 * t37 + t158 * t36;
t10 = t156 * t36 + t158 * t37;
t49 = -t106 * t156 + t107 * t158;
t137 = qJ(3) * t186;
t84 = t160 * t187 + t137;
t70 = pkin(2) * t177 - t190;
t87 = pkin(2) * t185 - t189;
t164 = -pkin(6) * t161 - qJD(1) * t87 - t70;
t64 = t167 + t189;
t56 = qJD(1) * t167 + t190;
t31 = t81 * pkin(4) + t56;
t12 = t46 * qJD(5) + t173;
t119 = -pkin(6) * t177 + t149;
t120 = t139 + t171;
t123 = t140 + t150;
t163 = t119 * t159 + (t120 * t159 + (-t123 + t140) * t157) * qJD(2);
t132 = t157 * t192;
t127 = -0.2e1 * t166;
t126 = 0.2e1 * t166;
t122 = t188 * t162;
t115 = pkin(2) * t187 - t137;
t112 = t188 * t180;
t97 = t153 * t184 - t154 * t185;
t54 = pkin(4) * t106 + t104;
t48 = t158 * t106 + t107 * t156;
t47 = -t92 * pkin(4) + t84;
t38 = t97 * pkin(4) + t64;
t24 = -pkin(7) * t97 + t35;
t23 = -pkin(7) * t98 + t34;
t17 = t49 * qJD(5) + t156 * t98 + t158 * t97;
t16 = t106 * t182 + t107 * t183 + t156 * t97 - t158 * t98;
t5 = -t156 * t22 + t158 * t21;
t4 = -t10 * qJD(5) - t156 * t24 + t158 * t23;
t3 = t9 * qJD(5) + t156 * t23 + t158 * t24;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0.2e1 * t112, t144, t127, -t191, 0, -pkin(6) * t144 + t157 * t172, pkin(6) * t191 + t159 * t172, 0, 0, t126, t144, -0.2e1 * t112, 0, t191, t127, t164 * t159 + t168 * t185, t163, t164 * t157 - t168 * t184, t163 * pkin(6) + t103 * t87 + t70 * t121, t107 * t82 + t92 * t98, -t106 * t82 - t107 * t81 - t89 * t98 - t92 * t97, -t98 * qJD(2), t106 * t81 + t89 * t97, t97 * qJD(2), 0, -qJD(2) * t34 + t104 * t81 + t106 * t56 + t64 * t89 + t65 * t97, qJD(2) * t35 + t104 * t82 + t107 * t56 + t64 * t92 + t65 * t98, -t106 * t26 + t107 * t174 - t32 * t98 - t33 * t97 - t34 * t92 - t35 * t89 - t57 * t82 - t58 * t81, t104 * t56 - t174 * t57 + t26 * t58 + t32 * t34 + t33 * t35 + t64 * t65, -t11 * t49 - t16 * t46, t11 * t48 - t12 * t49 + t16 * t207 - t17 * t46, t16 * t148, t12 * t48 + t17 * t207, t17 * t148, 0, t12 * t54 - t148 * t4 + t17 * t39 + t207 * t38 + t31 * t48, -t11 * t54 + t148 * t3 - t16 * t39 + t31 * t49 + t38 * t46, -t1 * t48 - t10 * t12 + t11 * t9 + t16 * t5 - t17 * t6 - t2 * t49 - t207 * t3 - t4 * t46, t1 * t10 + t2 * t9 + t3 * t6 + t31 * t54 + t38 * t39 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t122, 0, t132, 0, 0, t162 * pkin(1) * t157, pkin(1) * t192, 0, 0, -t132, 0, t122, 0, 0, t132, (-t103 * t157 + t115 * t159) * qJD(1), ((t123 - t150) * t157 + (-t120 + t171) * t159) * qJD(1), 0.2e1 * t149 + (t103 * t159 + t115 * t157) * qJD(1), t119 * qJ(3) + t123 * qJD(3) - t103 * t115 + (t123 * t157 + (-t120 - t195) * t159) * qJD(1) * pkin(6), -t199, -t205 + t206, 0, t199, -t128 + (t92 + t178) * qJD(2), 0, t170 * qJD(2) + t65 * t92 - t84 * t89 + t174, t169 * qJD(2) - t65 * t89 - t84 * t92 + t26, -t117 * t82 - t118 * t81 + (t170 - t33) * t92 + (-t169 + t32) * t89, -t165 * qJD(3) - t117 * t174 + t118 * t26 - t32 * t50 - t33 * t51 - t65 * t84, -t200, t211, t208, t200, t12 + t213, 0, -t201 * t148 - t207 * t47 - t210, t202 * t148 - t46 * t47 - t209, t11 * t52 - t12 * t53 + (-t201 - t6) * t46 + (-t202 + t5) * t207, t1 * t53 + t2 * t52 + t201 * t5 + t202 * t6 - t39 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, 0, -t151 * t162 - t161, -qJD(2) * t123 + t103 * t187 + t133, 0, 0, 0, 0, 0, 0, -t153 * t161 - t89 * t187, -t154 * t161 - t92 * t187, -t153 * t81 - t154 * t82 + (-t153 * t92 + t154 * t89) * qJD(2), t165 * qJD(2) + t26 * t153 - t154 * t174 - t65 * t187, 0, 0, 0, 0, 0, 0, -t197 * t148 - t187 * t207, -t196 * t148 - t46 * t187, t105 * t11 - t108 * t12 + t196 * t207 - t197 * t46, t1 * t108 + t105 * t2 - t39 * t187 - t196 * t6 + t197 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 + (-t92 + t178) * qJD(2), 0.2e1 * t89 * qJD(2), -t205 - t206, t32 * t92 + t33 * t89 + t56, 0, 0, 0, 0, 0, 0, t12 - t213, -t11 + t212, -t214 - t215, t207 * t6 + t46 * t5 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t211, -t208, -t200, -t173 + (-qJD(5) - t148) * t46, 0, -t6 * t148 + t210, -t148 * t5 + t209, 0, 0;];
tauc_reg = t7;
