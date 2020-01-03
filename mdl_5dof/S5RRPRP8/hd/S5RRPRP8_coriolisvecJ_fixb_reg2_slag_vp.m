% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:27
% EndTime: 2019-12-31 20:04:32
% DurationCPUTime: 1.69s
% Computational Cost: add. (1837->241), mult. (4399->295), div. (0->0), fcn. (2651->4), ass. (0->142)
t137 = qJD(2) - qJD(4);
t143 = sin(qJ(4));
t144 = sin(qJ(2));
t145 = cos(qJ(4));
t146 = cos(qJ(2));
t88 = t144 * t143 + t146 * t145;
t72 = t88 * qJD(1);
t190 = t72 * t137;
t172 = qJD(1) * qJD(2);
t167 = t146 * t172;
t168 = t144 * t172;
t35 = qJD(4) * t72 - t143 * t168 - t145 * t167;
t17 = t35 + t190;
t173 = qJD(4) * t146;
t174 = qJD(4) * t145;
t176 = qJD(2) * t146;
t201 = t144 * t174 + (-t173 + t176) * t143;
t178 = qJD(1) * t146;
t179 = qJD(1) * t144;
t74 = -t143 * t178 + t145 * t179;
t188 = t74 * t137;
t36 = t201 * qJD(1) - t145 * t168;
t19 = t36 + t188;
t199 = t74 ^ 2;
t71 = t72 ^ 2;
t25 = -t71 + t199;
t147 = -pkin(2) - pkin(3);
t100 = -t143 * qJ(3) + t145 * t147;
t198 = pkin(6) - pkin(7);
t197 = t74 * t72;
t189 = t74 * qJ(5);
t170 = t147 * qJD(2);
t128 = pkin(6) * t179;
t95 = pkin(7) * t179 - t128;
t62 = qJD(3) + t170 - t95;
t139 = qJD(2) * qJ(3);
t129 = pkin(6) * t178;
t97 = -pkin(7) * t178 + t129;
t76 = t139 + t97;
t27 = -t143 * t76 + t145 * t62;
t13 = t27 - t189;
t12 = -t137 * pkin(4) + t13;
t196 = t12 - t13;
t43 = t143 * t97 + t145 * t95;
t108 = t198 * t144;
t109 = t198 * t146;
t49 = t143 * t108 + t145 * t109;
t195 = qJD(2) * pkin(2);
t191 = t72 * qJ(5);
t28 = t143 * t62 + t145 * t76;
t14 = t28 - t191;
t194 = t14 * t137;
t193 = t27 * t137;
t192 = t28 * t137;
t149 = qJD(1) ^ 2;
t186 = t146 * t149;
t148 = qJD(2) ^ 2;
t185 = t148 * t144;
t133 = t148 * t146;
t132 = t144 * qJD(3);
t182 = qJ(3) * t167 + qJD(1) * t132;
t181 = qJ(3) * t176 + t132;
t140 = t144 ^ 2;
t180 = t146 ^ 2 - t140;
t177 = qJD(2) * t144;
t175 = qJD(4) * t143;
t101 = t145 * qJ(3) + t143 * t147;
t63 = t145 * qJD(3) + t100 * qJD(4);
t64 = -t143 * qJD(3) - t101 * qJD(4);
t171 = -t101 * t36 - t63 * t72 - t64 * t74;
t103 = -t146 * pkin(2) - t144 * qJ(3) - pkin(1);
t77 = -qJD(1) * pkin(1) - pkin(2) * t178 - qJ(3) * t179;
t166 = -t72 * pkin(4) - qJD(5);
t42 = -t143 * t95 + t145 * t97;
t138 = qJD(2) * qJD(3);
t96 = t198 * t177;
t65 = -qJD(1) * t96 + t138;
t121 = pkin(6) * t167;
t84 = -pkin(7) * t167 + t121;
t165 = -t143 * t84 - t145 * t65 - t62 * t174 + t76 * t175;
t164 = t143 * t65 - t145 * t84 + t76 * t174 + t62 * t175;
t48 = t145 * t108 - t143 * t109;
t163 = -0.2e1 * pkin(1) * t172;
t162 = qJD(3) - t195;
t161 = qJD(1) * t103 + t77;
t85 = t146 * pkin(3) - t103;
t160 = t137 ^ 2;
t59 = pkin(3) * t178 - t77;
t159 = t144 * t170;
t158 = t144 * t167;
t3 = -t36 * qJ(5) - t72 * qJD(5) - t165;
t124 = qJ(3) * t178;
t68 = t147 * t179 + t124;
t58 = pkin(2) * t168 - t182;
t70 = pkin(2) * t177 - t181;
t157 = -pkin(6) * t148 - qJD(1) * t70 - t58;
t156 = -t59 * t74 - t164;
t155 = t59 * t72 + t165;
t154 = t35 * qJ(5) - t164;
t98 = qJD(2) * t109;
t15 = t108 * t174 - t109 * t175 + t143 * t98 - t145 * t96;
t54 = t159 + t181;
t31 = -t166 + t59;
t152 = t31 * t72 - t3;
t47 = qJD(1) * t159 + t182;
t11 = t36 * pkin(4) + t47;
t16 = -t49 * qJD(4) + t143 * t96 + t145 * t98;
t102 = t128 + t162;
t105 = t129 + t139;
t99 = -pkin(6) * t168 + t138;
t150 = t99 * t146 + (t102 * t146 + (-t105 + t129) * t144) * qJD(2);
t117 = t144 * t186;
t107 = -0.2e1 * t158;
t106 = 0.2e1 * t158;
t104 = t180 * t149;
t94 = -pkin(4) + t100;
t93 = pkin(2) * t179 - t124;
t89 = -t146 * t143 + t144 * t145;
t87 = t180 * t172;
t53 = t64 * t137;
t52 = t63 * t137;
t46 = t88 * pkin(4) + t85;
t45 = t88 * qJD(2) - t144 * t175 - t145 * t173;
t44 = -t145 * t177 + t201;
t41 = t45 * t137;
t40 = t44 * t137;
t37 = -t74 * pkin(4) + t68;
t34 = -t145 * t160 - t74 * t179;
t33 = -t143 * t160 - t72 * t179;
t30 = -t88 * qJ(5) + t49;
t29 = -t89 * qJ(5) + t48;
t23 = t189 + t43;
t22 = t42 - t191;
t21 = t44 * pkin(4) + t54;
t8 = t36 * t88 + t72 * t44;
t7 = -t35 * t89 + t74 * t45;
t6 = -t45 * qJ(5) - t89 * qJD(5) + t16;
t5 = -t44 * qJ(5) - t88 * qJD(5) + t15;
t4 = -t74 * qJD(5) + t154;
t2 = -t19 * t143 + t17 * t145;
t1 = t35 * t88 - t89 * t36 - t74 * t44 - t45 * t72;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0.2e1 * t87, t133, t107, -t185, 0, -pkin(6) * t133 + t144 * t163, pkin(6) * t185 + t146 * t163, 0, 0, t106, t133, -0.2e1 * t87, 0, t185, t107, t157 * t146 + t161 * t177, t150, t157 * t144 - t161 * t176, t150 * pkin(6) + t58 * t103 + t77 * t70, t7, t1, -t41, t8, t40, 0, -t16 * t137 + t85 * t36 + t59 * t44 + t47 * t88 + t54 * t72, t15 * t137 - t85 * t35 + t59 * t45 + t47 * t89 + t54 * t74, -t15 * t72 - t16 * t74 + t164 * t89 + t165 * t88 - t27 * t45 - t28 * t44 + t48 * t35 - t49 * t36, t28 * t15 + t27 * t16 - t164 * t48 - t165 * t49 + t47 * t85 + t59 * t54, t7, t1, -t41, t8, t40, 0, t11 * t88 - t6 * t137 + t21 * t72 + t31 * t44 + t46 * t36, t11 * t89 + t5 * t137 + t21 * t74 + t31 * t45 - t46 * t35, -t12 * t45 - t14 * t44 + t29 * t35 - t3 * t88 - t30 * t36 - t4 * t89 - t5 * t72 - t6 * t74, t11 * t46 + t12 * t6 + t14 * t5 + t31 * t21 + t4 * t29 + t3 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t104, 0, t117, 0, 0, t149 * pkin(1) * t144, pkin(1) * t186, 0, 0, -t117, 0, t104, 0, 0, t117, (-t144 * t77 + t146 * t93) * qJD(1), ((t105 - t139) * t144 + (-t102 + t162) * t146) * qJD(1), 0.2e1 * t138 + (t144 * t93 + t146 * t77) * qJD(1), t99 * qJ(3) + t105 * qJD(3) - t77 * t93 + (t105 * t144 + (-t102 - t195) * t146) * qJD(1) * pkin(6), -t197, -t25, t17, t197, t19, 0, t42 * t137 - t68 * t72 - t156 - t53, -t43 * t137 - t68 * t74 - t155 + t52, t100 * t35 + (-t28 + t42) * t74 + (t27 + t43) * t72 + t171, -t164 * t100 - t165 * t101 - t59 * t68 + (-t43 + t63) * t28 + (-t42 + t64) * t27, -t197, -t25, t17, t197, t19, 0, t22 * t137 - t37 * t72 - t53 + (qJD(5) + t31) * t74 - t154, -t23 * t137 - t37 * t74 - t152 + t52, t94 * t35 + (-t14 + t22) * t74 + (t12 + t23) * t72 + t171, t3 * t101 - t31 * t37 + t4 * t94 + (-t23 + t63) * t14 + (-t22 + t64) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, 0, -t140 * t149 - t148, -t105 * qJD(2) + t77 * t179 + t121, 0, 0, 0, 0, 0, 0, t33, t34, t2, -t59 * t179 + (-t164 - t192) * t145 + (-t165 + t193) * t143, 0, 0, 0, 0, 0, 0, t33, t34, t2, -t31 * t179 + (t4 - t194) * t145 + (t137 * t12 + t3) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t25, -t17, -t197, -t19, 0, t156 - t192, t155 - t193, 0, 0, t197, t25, -t17, -t197, -t19, 0, -t194 + (t166 - t31) * t74 + t154, -t199 * pkin(4) - t13 * t137 + t152, t35 * pkin(4) - t196 * t72, t196 * t14 + (-t31 * t74 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 - t188, -t35 + t190, -t71 - t199, t12 * t74 + t14 * t72 + t11;];
tauc_reg = t9;
