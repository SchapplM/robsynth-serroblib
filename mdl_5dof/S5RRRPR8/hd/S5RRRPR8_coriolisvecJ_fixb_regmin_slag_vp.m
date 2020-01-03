% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:32
% DurationCPUTime: 1.84s
% Computational Cost: add. (2154->251), mult. (5370->334), div. (0->0), fcn. (3684->6), ass. (0->151)
t120 = sin(qJ(3));
t121 = sin(qJ(2));
t123 = cos(qJ(2));
t190 = cos(qJ(3));
t91 = t120 * t123 + t190 * t121;
t175 = qJD(1) * t91;
t198 = qJD(5) + t175;
t115 = qJD(2) + qJD(3);
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t155 = t190 * t123;
t145 = qJD(1) * t155;
t165 = qJD(1) * t121;
t154 = t120 * t165;
t78 = -t145 + t154;
t59 = t119 * t115 - t122 * t78;
t200 = t198 * t59;
t199 = qJD(5) - t198;
t161 = qJD(1) * qJD(2);
t197 = -0.2e1 * t161;
t153 = qJD(3) * t190;
t193 = -pkin(7) - pkin(6);
t99 = t193 * t123;
t95 = qJD(1) * t99;
t82 = t120 * t95;
t98 = t193 * t121;
t93 = qJD(1) * t98;
t56 = t190 * t93 + t82;
t176 = pkin(2) * t153 + qJD(4) - t56;
t164 = qJD(3) * t120;
t85 = t190 * t95;
t55 = t120 * t93 - t85;
t143 = pkin(2) * t164 - t55;
t182 = qJD(2) * pkin(2);
t86 = t93 + t182;
t53 = -t190 * t86 - t82;
t169 = qJD(4) + t53;
t195 = t175 ^ 2;
t194 = pkin(3) + pkin(8);
t192 = t78 * pkin(4);
t191 = t175 * pkin(4);
t54 = t120 * t86 - t85;
t45 = -t115 * qJ(4) - t54;
t29 = -t45 - t192;
t174 = t120 * t121;
t90 = -t155 + t174;
t189 = t29 * t90;
t111 = -t123 * pkin(2) - pkin(1);
t138 = -t91 * qJ(4) + t111;
t36 = t194 * t90 + t138;
t140 = t115 * t174;
t167 = t115 * t145;
t49 = qJD(1) * t140 - t167;
t188 = t36 * t49;
t97 = t111 * qJD(1);
t127 = -qJ(4) * t175 + t97;
t40 = t78 * pkin(3) + t127;
t187 = t40 * t175;
t186 = t49 * t90;
t185 = t198 * t78;
t184 = t175 * t78;
t183 = t97 * t175;
t47 = t122 * t49;
t181 = t122 * t198;
t162 = qJD(5) * t122;
t163 = qJD(5) * t119;
t58 = t115 * t91;
t50 = t58 * qJD(1);
t21 = -t115 * t163 + t119 * t50 + t78 * t162;
t180 = t21 * t122;
t156 = qJD(2) * t193;
t94 = t121 * t156;
t96 = t123 * t156;
t26 = -t120 * t96 - t98 * t153 - t99 * t164 - t190 * t94;
t179 = t26 * t115;
t133 = t120 * t98 - t190 * t99;
t27 = t133 * qJD(3) + t120 * t94 - t190 * t96;
t178 = t27 * t115;
t177 = t191 + t176;
t126 = qJD(1) ^ 2;
t173 = t123 * t126;
t125 = qJD(2) ^ 2;
t172 = t125 * t121;
t171 = t125 * t123;
t170 = t191 + t169;
t166 = t121 ^ 2 - t123 ^ 2;
t113 = t121 * t182;
t146 = qJD(1) * t156;
t87 = t121 * t146;
t88 = t123 * t146;
t150 = -t120 * t88 - t86 * t153 - t95 * t164 - t190 * t87;
t18 = -t115 * qJD(4) + t150;
t11 = -t50 * pkin(4) - t18;
t23 = -t194 * t115 + t170;
t25 = t194 * t78 + t127;
t7 = t119 * t23 + t122 * t25;
t159 = t11 * t122 - t7 * t78;
t158 = t198 * t163;
t157 = t90 * t162;
t152 = t121 * t161;
t151 = t198 * t29;
t20 = t120 * t87 - t95 * t153 + t86 * t164 - t190 * t88;
t149 = t119 * t198;
t148 = pkin(1) * t197;
t110 = -t190 * pkin(2) - pkin(3);
t107 = -pkin(8) + t110;
t51 = pkin(3) * t175 + t78 * qJ(4);
t42 = pkin(2) * t165 + t51;
t74 = t175 * pkin(8);
t147 = -qJD(5) * t107 + t42 + t74;
t144 = t192 + t143;
t62 = -t120 * t99 - t190 * t98;
t142 = t198 * t58 - t186;
t139 = t119 * t25 - t122 * t23;
t141 = t11 * t119 - t139 * t78 + (t122 * t175 + t162) * t29;
t61 = t122 * t115 + t119 * t78;
t137 = t20 + t187;
t136 = t97 * t78 + t150;
t135 = -t149 * t198 - t47;
t134 = t115 * t54 - t20;
t132 = pkin(2) * t152 + t49 * qJ(4) - qJD(4) * t175;
t57 = -qJD(2) * t155 - t123 * t153 + t140;
t131 = t57 * qJ(4) - t91 * qJD(4) + t113;
t130 = -t40 * t78 - t18;
t43 = t91 * pkin(4) + t62;
t129 = t11 * t90 + t29 * t58 + t43 * t49;
t128 = t119 * t49 - t181 * t198;
t30 = t167 + (-t154 + t78) * t115;
t108 = t120 * pkin(2) + qJ(4);
t52 = t90 * pkin(3) + t138;
t48 = t122 * t50;
t44 = -t90 * pkin(4) + t133;
t41 = -t115 * pkin(3) + t169;
t37 = -t78 ^ 2 + t195;
t35 = t49 * t91;
t34 = t54 - t192;
t32 = t51 + t74;
t22 = t61 * qJD(5) - t48;
t17 = t58 * pkin(3) + t131;
t16 = -t57 * pkin(4) + t27;
t15 = -t58 * pkin(4) - t26;
t14 = t50 * pkin(3) + t132;
t13 = -t49 * pkin(4) + t20;
t12 = t122 * t13;
t10 = t194 * t58 + t131;
t5 = t194 * t50 + t132;
t4 = -t59 * t78 + t128;
t3 = t61 * t78 + t135;
t2 = -t61 * t149 + t180;
t1 = (-t198 * t61 - t22) * t122 + (-t21 + t200) * t119;
t6 = [0, 0, 0, 0.2e1 * t123 * t152, t166 * t197, t171, -t172, 0, -pkin(6) * t171 + t121 * t148, pkin(6) * t172 + t123 * t148, -t175 * t57 - t35, -t175 * t58 - t91 * t50 + t57 * t78 + t186, -t57 * t115, -t58 * t115, 0, t111 * t50 - t178 + t97 * t58 + (qJD(1) * t90 + t78) * t113, -t111 * t49 + 0.2e1 * t175 * t113 - t97 * t57 + t179, -t133 * t50 + t175 * t27 + t18 * t90 + t20 * t91 + t26 * t78 - t41 * t57 + t45 * t58 - t49 * t62, -t14 * t90 - t17 * t78 - t40 * t58 - t50 * t52 + t178, -t14 * t91 - t17 * t175 + t40 * t57 + t49 * t52 - t179, -t133 * t18 + t14 * t52 + t17 * t40 + t20 * t62 + t26 * t45 + t27 * t41, t61 * t157 + (t21 * t90 + t58 * t61) * t119, (-t119 * t59 + t122 * t61) * t58 + (-t119 * t22 + t180 + (-t119 * t61 - t122 * t59) * qJD(5)) * t90, t119 * t142 + t157 * t198 + t21 * t91 - t61 * t57, t122 * t142 - t90 * t158 - t22 * t91 + t59 * t57, -t198 * t57 - t35, t12 * t91 + t15 * t59 + t44 * t22 + t139 * t57 + (-t10 * t198 - t5 * t91 + t188) * t119 + (t16 * t198 - t129) * t122 + ((-t119 * t43 - t122 * t36) * t198 - t7 * t91 + t119 * t189) * qJD(5), t15 * t61 + t44 * t21 + t7 * t57 + (-(qJD(5) * t43 + t10) * t198 + t188 - (qJD(5) * t23 + t5) * t91 + qJD(5) * t189) * t122 + (-(-qJD(5) * t36 + t16) * t198 - (-qJD(5) * t25 + t13) * t91 + t129) * t119; 0, 0, 0, -t121 * t173, t166 * t126, 0, 0, 0, t126 * pkin(1) * t121, pkin(1) * t173, t184, t37, t30, 0, 0, t55 * t115 - t183 + (-t115 * t164 - t78 * t165) * pkin(2) - t20, t56 * t115 + (-t115 * t153 - t165 * t175) * pkin(2) + t136, -t108 * t50 - t110 * t49 + (t143 - t45) * t175 + (t41 - t176) * t78, t143 * t115 + t42 * t78 + t137, t176 * t115 + t175 * t42 + t130, -t18 * t108 + t20 * t110 + t143 * t41 - t176 * t45 - t40 * t42, t2, t1, t3, t4, t185, -t107 * t47 + t108 * t22 + t177 * t59 + (t119 * t147 + t122 * t144) * t198 + t141, t108 * t21 + t177 * t61 + t147 * t181 + (t107 * t49 - t144 * t198 - t151) * t119 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t37, t30, 0, 0, t134 - t183, -t115 * t53 + t136, pkin(3) * t49 - qJ(4) * t50 + (-t45 - t54) * t175 + (t41 - t169) * t78, t51 * t78 - t134 + t187, t169 * t115 + t175 * t51 + t130, -t20 * pkin(3) - t18 * qJ(4) - t169 * t45 - t40 * t51 - t41 * t54, t2, t1, t3, t4, t185, qJ(4) * t22 - (-t119 * t32 + t122 * t34) * t198 + t170 * t59 - (-t158 - t47) * t194 + t141, qJ(4) * t21 + t170 * t61 + (qJD(5) * t194 + t32) * t181 + (-t194 * t49 + t198 * t34 - t151) * t119 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t184, -t115 ^ 2 - t195, t115 * t45 + t137, 0, 0, 0, 0, 0, -t115 * t59 + t135, -t115 * t61 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t59, -t59 ^ 2 + t61 ^ 2, t21 + t200, -t199 * t61 + t48, -t49, -t119 * t5 - t199 * t7 - t29 * t61 + t12, -t119 * t13 - t122 * t5 + t199 * t139 + t29 * t59;];
tauc_reg = t6;
