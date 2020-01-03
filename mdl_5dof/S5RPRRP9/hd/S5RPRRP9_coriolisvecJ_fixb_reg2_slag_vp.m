% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:17
% DurationCPUTime: 1.76s
% Computational Cost: add. (3377->241), mult. (9232->293), div. (0->0), fcn. (7056->6), ass. (0->133)
t106 = sin(qJ(4));
t173 = cos(qJ(4));
t104 = sin(pkin(8));
t105 = cos(pkin(8));
t107 = sin(qJ(3));
t108 = cos(qJ(3));
t87 = t108 * t104 + t107 * t105;
t182 = t87 * qJD(1);
t151 = t108 * t105;
t152 = t107 * t104;
t130 = -t151 + t152;
t80 = t130 * qJD(1);
t126 = t106 * t80 - t173 * t182;
t142 = qJD(1) * t152;
t148 = qJD(3) * t108;
t94 = t105 * qJD(1) * t148;
t73 = qJD(3) * t142 - t94;
t74 = qJD(3) * t182;
t114 = t126 * qJD(4) + t106 * t73 - t173 * t74;
t103 = qJD(3) + qJD(4);
t183 = t103 * t126;
t191 = t114 - t183;
t141 = qJD(4) * t173;
t147 = qJD(4) * t106;
t125 = -t106 * t74 - t80 * t141 - t147 * t182 - t173 * t73;
t53 = -t106 * t182 - t173 * t80;
t154 = t53 * t103;
t190 = t125 - t154;
t170 = t53 ^ 2;
t189 = t126 ^ 2;
t132 = t189 - t170;
t98 = -t105 * pkin(2) - pkin(1);
t90 = t98 * qJD(1) + qJD(2);
t60 = t80 * pkin(3) + t90;
t17 = -pkin(4) * t53 + qJ(5) * t126 + t60;
t188 = t17 * t53;
t187 = t60 * t53;
t186 = t126 * t60;
t172 = t17 * t126;
t185 = t53 * t126;
t163 = pkin(6) + qJ(2);
t91 = t163 * t104;
t92 = t163 * t105;
t62 = -t107 * t91 + t108 * t92;
t184 = t130 * pkin(7) - t62;
t26 = -pkin(4) * t126 - t53 * qJ(5);
t123 = t87 * qJD(2);
t181 = t184 * qJD(3) - t123;
t179 = t182 ^ 2;
t178 = qJD(3) ^ 2;
t174 = t87 * pkin(7);
t155 = t107 * t92;
t61 = -t108 * t91 - t155;
t48 = t61 - t174;
t127 = t106 * t184 + t173 * t48;
t145 = qJD(1) * qJD(2);
t140 = t104 * t145;
t88 = qJD(1) * t91;
t96 = qJD(2) * t151;
t161 = qJD(1) * t96 - t88 * t148;
t89 = qJD(1) * t92;
t39 = (-qJD(3) * t89 - t140) * t107 + t161;
t33 = -t74 * pkin(7) + t39;
t118 = qJD(1) * t123;
t59 = -t107 * t88 + t108 * t89;
t40 = -t59 * qJD(3) - t118;
t34 = t73 * pkin(7) + t40;
t139 = t106 * t33 - t173 * t34;
t44 = -t80 * pkin(7) + t59;
t143 = t173 * t44;
t156 = t107 * t89;
t58 = -t108 * t88 - t156;
t43 = -pkin(7) * t182 + t58;
t42 = qJD(3) * pkin(3) + t43;
t16 = t106 * t42 + t143;
t3 = t16 * qJD(4) + t139;
t177 = t3 * t127;
t176 = t74 * pkin(3);
t175 = t182 * pkin(3);
t146 = t104 * qJD(2);
t160 = -t91 * t148 + t96;
t36 = -t107 * t146 + (-t155 - t174) * qJD(3) + t160;
t7 = t127 * qJD(4) + t181 * t106 + t173 * t36;
t166 = t7 * t103;
t23 = t106 * t48 - t173 * t184;
t8 = t23 * qJD(4) + t106 * t36 - t173 * t181;
t165 = t8 * t103;
t164 = t182 * t80;
t157 = t106 * t44;
t19 = t173 * t43 - t157;
t162 = -pkin(3) * t141 - qJD(5) + t19;
t120 = t130 * qJD(3);
t122 = t87 * qJD(3);
t57 = -t106 * t130 + t173 * t87;
t32 = t57 * qJD(4) - t106 * t120 + t173 * t122;
t159 = t103 * t32;
t15 = t173 * t42 - t157;
t150 = qJD(5) - t15;
t149 = t104 ^ 2 + t105 ^ 2;
t137 = -t106 * t34 - t42 * t141 + t44 * t147 - t173 * t33;
t136 = t149 * qJD(1) ^ 2;
t18 = t106 * t43 + t143;
t135 = pkin(3) * t147 - t18;
t100 = t103 * qJD(5);
t1 = t100 - t137;
t119 = t173 * t130;
t56 = t106 * t87 + t119;
t134 = -t114 * t56 - t32 * t53;
t133 = -t189 - t170;
t129 = 0.2e1 * t149 * t145;
t128 = t15 * t103 + t137;
t117 = t114 * t23 - t125 * t127 - t126 * t8 + t3 * t57 + t53 * t7;
t116 = t125 + t154;
t31 = t103 * t119 + t106 * t122 + t87 * t147;
t115 = t114 * t57 - t125 * t56 + t126 * t32 - t31 * t53;
t5 = -pkin(4) * t114 - qJ(5) * t125 + qJD(5) * t126 + t176;
t66 = t130 * pkin(3) + t98;
t112 = -t139 + (-qJD(4) + t103) * t16;
t111 = -t114 - t183;
t110 = t18 * t103 + (-t143 + (-pkin(3) * t103 - t42) * t106) * qJD(4) - t139;
t99 = -t173 * pkin(3) - pkin(4);
t97 = t106 * pkin(3) + qJ(5);
t77 = t80 ^ 2;
t46 = -qJD(3) * t62 - t123;
t45 = (-qJD(3) * t92 - t146) * t107 + t160;
t27 = t31 * t103;
t21 = t56 * pkin(4) - t57 * qJ(5) + t66;
t20 = t175 + t26;
t14 = t103 * qJ(5) + t16;
t13 = -t103 * pkin(4) + t150;
t9 = pkin(3) * t122 + t32 * pkin(4) + t31 * qJ(5) - t57 * qJD(5);
t6 = t125 * t57 + t126 * t31;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, qJ(2) * t129, -t120 * t182 - t73 * t87, -t87 * t74 + t73 * t130 + (t130 * t80 - t182 * t87) * qJD(3), -t130 * t178, t80 * t122 + t74 * t130, -t87 * t178, 0, t98 * t74 + (t90 * t87 + t46) * qJD(3), -t98 * t73 + (-t90 * t130 - t45) * qJD(3), -t45 * t80 - t62 * t74 - t39 * t130 - t46 * t182 + t61 * t73 - t40 * t87 + (t58 * t130 - t59 * t87) * qJD(3), t39 * t62 + t40 * t61 + t59 * t45 + t58 * t46, t6, t115, -t27, t134, -t159, 0, -t165 - t66 * t114 + t60 * t32 + (-t122 * t53 + t74 * t56) * pkin(3), -t166 + t66 * t125 - t60 * t31 + (-t122 * t126 + t74 * t57) * pkin(3), t137 * t56 + t15 * t31 - t16 * t32 + t117, -t15 * t8 + t16 * t7 - t137 * t23 - t177 + (t60 * t122 + t74 * t66) * pkin(3), t6, -t27, -t115, 0, t159, t134, -t114 * t21 + t17 * t32 + t5 * t56 - t53 * t9 - t165, -t1 * t56 - t13 * t31 - t14 * t32 + t117, -t125 * t21 + t126 * t9 + t17 * t31 - t5 * t57 + t166, t1 * t23 + t13 * t8 + t14 * t7 + t17 * t9 + t5 * t21 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -qJ(2) * t136, 0, 0, 0, 0, 0, 0, 0.2e1 * t74, t94 + (-t80 - t142) * qJD(3), -t77 - t179, t182 * t58 + t59 * t80, 0, 0, 0, 0, 0, 0, t111, t116, t133, -t126 * t15 - t16 * t53 + t176, 0, 0, 0, 0, 0, 0, t111, t133, -t116, t126 * t13 - t14 * t53 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t77 + t179, t94 + (t80 - t142) * qJD(3), -t164, 0, 0, -t182 * t90 - t118, t107 * t140 + t90 * t80 + (t58 + t156) * qJD(3) - t161, 0, 0, t185, t132, t190, -t185, t191, 0, t175 * t53 + t110 + t186, t19 * t103 - t187 + (-t103 * t141 + t126 * t182) * pkin(3) + t137, t15 * t53 - t16 * t126 + t18 * t126 - t19 * t53 + (-t173 * t125 + t106 * t114 + (-t106 * t126 + t173 * t53) * qJD(4)) * pkin(3), t15 * t18 - t16 * t19 + (-t173 * t3 - t106 * t137 - t60 * t182 + (-t106 * t15 + t173 * t16) * qJD(4)) * pkin(3), t185, t190, -t132, 0, -t191, -t185, t20 * t53 + t110 + t172, t125 * t99 + t97 * t114 + (-t13 - t162) * t53 + (-t135 - t14) * t126, -t162 * t103 - t126 * t20 + t1 + t188, t1 * t97 + t135 * t13 - t162 * t14 - t17 * t20 + t3 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t132, t190, -t185, t191, 0, t112 + t186, t128 - t187, 0, 0, t185, t190, -t132, 0, -t191, -t185, t26 * t53 + t112 + t172, -pkin(4) * t125 + t114 * qJ(5) - (t14 - t16) * t126 - (t13 - t150) * t53, -t126 * t26 + 0.2e1 * t100 - t128 + t188, -t3 * pkin(4) + t1 * qJ(5) - t13 * t16 + t150 * t14 - t17 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t190, -t103 ^ 2 - t189, -t14 * t103 - t172 + t3;];
tauc_reg = t2;
