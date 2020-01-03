% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:16
% DurationCPUTime: 2.35s
% Computational Cost: add. (4053->293), mult. (8376->399), div. (0->0), fcn. (5345->6), ass. (0->166)
t200 = sin(qJ(4));
t147 = t200 * qJD(4);
t208 = t200 * qJD(3) + t147;
t105 = sin(qJ(3));
t107 = cos(qJ(4));
t108 = cos(qJ(3));
t71 = t107 * t105 + t200 * t108;
t66 = t71 * qJD(1);
t205 = qJD(5) + t66;
t104 = sin(qJ(5));
t106 = cos(qJ(5));
t162 = qJD(3) + qJD(4);
t141 = t106 * t162;
t168 = qJD(1) * t108;
t151 = t107 * t168;
t153 = t200 * t105;
t68 = -qJD(1) * t153 + t151;
t54 = t104 * t68 - t141;
t207 = t205 * t54;
t109 = -pkin(1) - pkin(6);
t85 = t109 * qJD(1) + qJD(2);
t145 = pkin(7) * qJD(1) - t85;
t202 = t145 * t105;
t180 = t107 * t202;
t206 = t200 * t202;
t204 = t208 * qJD(1);
t177 = t107 * t108;
t123 = t162 * t177;
t51 = -t208 * t105 + t123;
t203 = t51 * t162;
t63 = -pkin(7) * t168 + t108 * t85;
t60 = qJD(3) * pkin(3) + t63;
t39 = t200 * t60 - t180;
t36 = t162 * pkin(8) + t39;
t169 = qJD(1) * t105;
t77 = pkin(3) * t169 + qJD(1) * qJ(2);
t40 = t66 * pkin(4) - t68 * pkin(8) + t77;
t119 = t104 * t36 - t106 * t40;
t166 = qJD(3) * t108;
t61 = t145 * t166;
t144 = qJD(4) * t60 - t61;
t58 = t202 * t147;
t16 = qJD(3) * t206 + t144 * t107 + t58;
t45 = t162 * t107 * t169 + t204 * t108;
t190 = t204 * t105;
t46 = qJD(1) * t123 - t190;
t100 = qJD(1) * qJD(2);
t163 = qJD(1) * qJD(3);
t146 = t108 * t163;
t74 = pkin(3) * t146 + t100;
t23 = t46 * pkin(4) + t45 * pkin(8) + t74;
t3 = -t119 * qJD(5) + t104 * t23 + t106 * t16;
t201 = 0.2e1 * t100;
t135 = -qJD(3) * t180 - t200 * t61;
t17 = t39 * qJD(4) + t135;
t192 = pkin(7) - t109;
t75 = t192 * t105;
t76 = t192 * t108;
t52 = t107 * t76 - t200 * t75;
t199 = t17 * t52;
t72 = -t153 + t177;
t198 = t17 * t72;
t1 = t3 * t106;
t197 = t46 * t71;
t56 = t104 * t162 + t106 * t68;
t196 = t56 * t54;
t195 = t205 * t68;
t194 = t68 * t66;
t193 = t72 * t46;
t41 = t200 * t63 - t180;
t191 = t39 - t41;
t189 = pkin(3) * qJD(4);
t188 = t104 * t45;
t187 = t104 * t46;
t186 = t104 * t54;
t185 = t104 * t66;
t184 = t106 * t46;
t183 = t106 * t56;
t182 = t106 * t66;
t181 = t107 * t60;
t165 = qJD(5) * t104;
t27 = -qJD(5) * t141 + t106 * t45 + t68 * t165;
t179 = t27 * t104;
t28 = t56 * qJD(5) - t188;
t178 = t28 * t106;
t92 = t105 * pkin(3) + qJ(2);
t110 = qJD(3) ^ 2;
t176 = t110 * t105;
t175 = t110 * t108;
t111 = qJD(1) ^ 2;
t174 = t111 * qJ(2);
t173 = t111 * t108;
t86 = pkin(3) * t166 + qJD(2);
t171 = t105 ^ 2 - t108 ^ 2;
t170 = -t110 - t111;
t167 = qJD(3) * t105;
t164 = qJD(5) * t106;
t14 = t104 * t40 + t106 * t36;
t161 = t119 * t182 - t14 * t185 + t1;
t160 = 0.2e1 * qJD(1);
t159 = t107 * t189;
t158 = pkin(3) * t168;
t156 = t72 * t165;
t155 = t72 * t164;
t154 = t105 * t173;
t38 = t206 + t181;
t35 = -t162 * pkin(4) - t38;
t31 = t35 * t165;
t32 = t35 * t164;
t150 = t119 * t68 + t31;
t149 = t77 * t66 - t58;
t143 = t106 * t205;
t142 = qJD(5) * t71 + qJD(1);
t140 = t17 * t104 + t14 * t68 + t32;
t139 = t105 * t146;
t47 = t68 * pkin(4) + t66 * pkin(8);
t136 = t192 * t167;
t4 = -qJD(5) * t14 - t104 * t16 + t106 * t23;
t134 = t119 * t51 - t4 * t71;
t133 = t14 * t51 + t3 * t71;
t50 = t162 * t71;
t132 = -t35 * t50 + t198;
t131 = -t72 * t27 - t50 * t56;
t130 = -t27 * t71 + t56 * t51;
t129 = t72 * t28 - t50 * t54;
t128 = -t28 * t71 - t54 * t51;
t127 = -t72 * t45 - t68 * t50;
t126 = -t205 * t50 + t193;
t125 = t51 * t66 + t197;
t122 = t104 * t14 - t106 * t119;
t121 = -t104 * t119 - t106 * t14;
t48 = t71 * pkin(4) - t72 * pkin(8) + t92;
t53 = -t107 * t75 - t200 * t76;
t25 = -t104 * t53 + t106 * t48;
t26 = t104 * t48 + t106 * t53;
t118 = pkin(3) * t147 - t41;
t117 = -t77 * t68 - t135;
t94 = t200 * pkin(3) + pkin(8);
t115 = -t159 * t205 + t35 * t66 - t46 * t94;
t114 = t16 * t71 - t38 * t50 + t39 * t51 - t198;
t113 = -t122 * qJD(5) - t4 * t104;
t112 = t113 + t1;
t98 = qJ(2) * t201;
t95 = -t107 * pkin(3) - pkin(4);
t70 = qJD(3) * t76;
t49 = t50 * t162;
t43 = t47 + t158;
t42 = t107 * t63 + t206;
t37 = -t66 ^ 2 + t68 ^ 2;
t34 = t190 + (-t151 + t68) * t162;
t33 = t66 * t162 - t45;
t30 = t53 * qJD(4) - t107 * t136 - t200 * t70;
t29 = -t52 * qJD(4) - t107 * t70 + t200 * t136;
t24 = t51 * pkin(4) + t50 * pkin(8) + t86;
t22 = t104 * t47 + t106 * t38;
t21 = -t104 * t38 + t106 * t47;
t19 = t104 * t43 + t106 * t42;
t18 = -t104 * t42 + t106 * t43;
t10 = t143 * t205 - t56 * t68 + t187;
t9 = -t104 * t205 ^ 2 + t54 * t68 + t184;
t8 = t186 * t205 - t178;
t7 = t56 * t143 - t179;
t6 = -t26 * qJD(5) - t104 * t29 + t106 * t24;
t5 = t25 * qJD(5) + t104 * t24 + t106 * t29;
t2 = (-t27 - t207) * t106 + (-t205 * t56 - t28) * t104;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t98, -0.2e1 * t139, 0.2e1 * t171 * t163, -t176, 0.2e1 * t139, -t175, 0, -t109 * t176 + (qJ(2) * t166 + qJD(2) * t105) * t160, -t109 * t175 + (-qJ(2) * t167 + qJD(2) * t108) * t160, 0, t98, t127, t45 * t71 + t50 * t66 - t68 * t51 - t193, -t49, t125, -t203, 0, -t162 * t30 + t92 * t46 + t77 * t51 + t86 * t66 + t74 * t71, -t162 * t29 - t92 * t45 - t77 * t50 + t86 * t68 + t74 * t72, -t29 * t66 + t30 * t68 - t52 * t45 - t53 * t46 - t114, t16 * t53 + t39 * t29 - t38 * t30 + t74 * t92 + t77 * t86 + t199, t106 * t131 - t156 * t56, (t104 * t56 + t106 * t54) * t50 + (t179 - t178 + (-t183 + t186) * qJD(5)) * t72, t106 * t126 - t156 * t205 + t130, t129 * t104 + t155 * t54, -t126 * t104 - t155 * t205 + t128, t205 * t51 + t197, t132 * t104 + t205 * t6 + t25 * t46 + t52 * t28 + t30 * t54 + t32 * t72 - t134, t106 * t132 - t205 * t5 - t26 * t46 - t52 * t27 + t30 * t56 - t31 * t72 - t133, t25 * t27 - t26 * t28 - t5 * t54 - t6 * t56 + t122 * t50 + (qJD(5) * t121 - t104 * t3 - t106 * t4) * t72, -t119 * t6 + t14 * t5 + t4 * t25 + t3 * t26 + t35 * t30 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t174, 0, 0, 0, 0, 0, 0, t170 * t105, t170 * t108, 0, -t174, 0, 0, 0, 0, 0, 0, -qJD(1) * t66 - t49, -qJD(1) * t68 - t203, -t125 - t127, -t77 * qJD(1) + t114, 0, 0, 0, 0, 0, 0, -t71 * t187 + (-t104 * t51 - t106 * t142) * t205 - t129, -t71 * t184 + (t142 * t104 - t106 * t51) * t205 - t131, (t142 * t56 + t128) * t106 + (t142 * t54 + t130) * t104, (t119 * t142 + t133) * t106 + (-t14 * t142 + t134) * t104 - t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t171 * t111, 0, -t154, 0, 0, -qJ(2) * t173, t105 * t174, 0, 0, t194, t37, t33, -t194, t34, 0, t41 * qJD(3) - t191 * qJD(4) + (-t147 * t162 - t168 * t66) * pkin(3) + t117, -t68 * t158 + t42 * qJD(4) + (t42 - t206) * qJD(3) + (-t162 * t189 - t144) * t107 + t149, t191 * t68 + (-t38 + t42) * t66 + (-t200 * t46 + t107 * t45 + (-t107 * t66 + t200 * t68) * qJD(4)) * pkin(3), t38 * t41 - t39 * t42 + (-t77 * t168 + t200 * t16 - t107 * t17 + (t107 * t39 - t200 * t38) * qJD(4)) * pkin(3), t7, t2, t10, t8, t9, -t195, -t18 * t205 + t95 * t28 + t118 * t54 + (-qJD(5) * t205 * t94 - t17) * t106 + t115 * t104 + t150, -t95 * t27 + (t165 * t94 + t19) * t205 + t118 * t56 + t115 * t106 + t140, t18 * t56 + t19 * t54 + (-t54 * t159 - t28 * t94 + (t56 * t94 + t119) * qJD(5)) * t106 + (t56 * t159 - t27 * t94 - t4 + (t54 * t94 - t14) * qJD(5)) * t104 + t161, t119 * t18 - t14 * t19 + t17 * t95 - t35 * t41 + (-t121 * t107 + t200 * t35) * t189 + t112 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t37, t33, -t194, t34, 0, t39 * qJD(3) + t117, t107 * t61 + (t38 - t181) * qJD(4) + (t38 - t206) * qJD(3) + t149, 0, 0, t7, t2, t10, t8, t9, -t195, t35 * t185 - pkin(4) * t28 - t17 * t106 - t21 * t205 - t39 * t54 + (-t164 * t205 - t187) * pkin(8) + t150, t35 * t182 + pkin(4) * t27 + t22 * t205 - t39 * t56 + (t165 * t205 - t184) * pkin(8) + t140, t21 * t56 + t22 * t54 + (-t179 - t178 + (t183 + t186) * qJD(5)) * pkin(8) + t113 + t161, -t17 * pkin(4) + pkin(8) * t112 + t119 * t21 - t14 * t22 - t35 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, -t54 ^ 2 + t56 ^ 2, -t27 + t207, -t196, t188 + (-qJD(5) + t205) * t56, t46, t14 * t205 - t35 * t56 + t4, -t119 * t205 + t35 * t54 - t3, 0, 0;];
tauc_reg = t11;
