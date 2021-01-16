% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:58
% EndTime: 2021-01-15 20:30:06
% DurationCPUTime: 1.63s
% Computational Cost: add. (2792->287), mult. (7189->389), div. (0->0), fcn. (5009->6), ass. (0->150)
t109 = sin(qJ(4));
t111 = cos(qJ(4));
t108 = sin(pkin(8));
t110 = sin(qJ(2));
t112 = cos(qJ(2));
t160 = cos(pkin(8));
t89 = t108 * t112 + t110 * t160;
t157 = qJD(1) * t89;
t62 = qJD(2) * t109 + t111 * t157;
t166 = t109 * t62;
t146 = qJD(1) * qJD(2);
t189 = -0.2e1 * t146;
t150 = qJD(1) * t110;
t136 = t160 * t112;
t98 = qJD(1) * t136;
t76 = t108 * t150 - t98;
t104 = -pkin(2) * t112 - pkin(1);
t151 = qJD(1) * t104;
t93 = qJD(3) + t151;
t31 = pkin(3) * t76 - pkin(7) * t157 + t93;
t172 = -qJ(3) - pkin(6);
t96 = t172 * t112;
t92 = qJD(1) * t96;
t140 = t160 * t92;
t167 = qJD(2) * pkin(2);
t95 = t172 * t110;
t91 = qJD(1) * t95;
t85 = t91 + t167;
t52 = t108 * t85 - t140;
t47 = qJD(2) * pkin(7) + t52;
t15 = t109 * t31 + t111 * t47;
t139 = t110 * t146;
t97 = t108 * t139;
t119 = qJD(2) * t98 - t97;
t78 = t89 * qJD(2);
t72 = qJD(1) * t78;
t99 = pkin(2) * t139;
t30 = t72 * pkin(3) - pkin(7) * t119 + t99;
t23 = t111 * t30;
t137 = qJD(2) * t172;
t75 = t112 * qJD(3) + t110 * t137;
t69 = t75 * qJD(1);
t118 = -t110 * qJD(3) + t112 * t137;
t70 = t118 * qJD(1);
t26 = t108 * t70 + t160 * t69;
t116 = -qJD(4) * t15 - t109 * t26 + t23;
t147 = t111 * qJD(2);
t149 = qJD(4) * t109;
t123 = qJD(4) * t147 + t111 * t119 - t149 * t157;
t115 = -qJ(5) * t123 + t116;
t184 = pkin(4) * t72;
t1 = -qJD(5) * t62 + t115 + t184;
t60 = t109 * t157 - t147;
t10 = -qJ(5) * t60 + t15;
t74 = qJD(4) + t76;
t180 = t10 * t74;
t188 = t1 + t180;
t187 = t74 * t166;
t28 = qJD(4) * t62 + t109 * t119;
t186 = t62 ^ 2;
t14 = -t109 * t47 + t111 * t31;
t9 = -qJ(5) * t62 + t14;
t5 = pkin(4) * t74 + t9;
t185 = t5 - t9;
t183 = t60 * pkin(4);
t101 = pkin(2) * t108 + pkin(7);
t153 = qJ(5) + t101;
t131 = qJD(4) * t153;
t158 = qJ(5) * t111;
t144 = pkin(2) * t150;
t41 = pkin(3) * t157 + pkin(7) * t76 + t144;
t35 = t111 * t41;
t82 = t108 * t92;
t54 = t160 * t91 + t82;
t182 = -pkin(4) * t157 - t111 * t131 - t158 * t76 - t35 + (-qJD(5) + t54) * t109;
t181 = pkin(2) * t110;
t25 = t108 * t69 - t160 * t70;
t179 = t25 * t89;
t51 = t160 * t85 + t82;
t46 = -qJD(2) * pkin(3) - t51;
t120 = -t108 * t110 + t136;
t81 = t120 * qJD(2);
t178 = t46 * t81;
t177 = t60 * t74;
t176 = t60 * t76;
t175 = t60 * t157;
t174 = t62 * t74;
t173 = t62 * t157;
t159 = qJ(5) * t109;
t169 = t109 * t41 + t111 * t54;
t171 = -t111 * qJD(5) + t109 * t131 + t159 * t76 + t169;
t148 = qJD(4) * t111;
t170 = -t109 * t28 - t60 * t148;
t50 = -pkin(3) * t120 - pkin(7) * t89 + t104;
t58 = t108 * t95 - t160 * t96;
t55 = t111 * t58;
t168 = t109 * t50 + t55;
t165 = t109 * t72;
t163 = t109 * t76;
t138 = qJD(5) + t183;
t18 = t138 + t46;
t162 = t111 * t18;
t161 = t123 * t109;
t114 = qJD(1) ^ 2;
t156 = t112 * t114;
t113 = qJD(2) ^ 2;
t155 = t113 * t110;
t154 = t113 * t112;
t152 = t110 ^ 2 - t112 ^ 2;
t40 = t108 * t118 + t160 * t75;
t143 = t110 * t167;
t42 = pkin(3) * t78 - pkin(7) * t81 + t143;
t145 = t109 * t42 + t111 * t40 + t50 * t148;
t142 = t89 * t149;
t141 = t89 * t148;
t39 = t108 * t75 - t160 * t118;
t53 = t108 * t91 - t140;
t57 = -t108 * t96 - t160 * t95;
t135 = -t109 * t30 - t111 * t26 - t31 * t148 + t47 * t149;
t134 = 0.2e1 * t157;
t133 = t111 * t74;
t132 = pkin(1) * t189;
t102 = -pkin(2) * t160 - pkin(3);
t130 = qJD(4) * t101 * t74 + t25;
t122 = qJ(5) * t28 + t135;
t2 = -qJD(5) * t60 - t122;
t129 = -t5 * t74 + t2;
t13 = pkin(4) * t28 + t25;
t128 = t13 * t89 + t18 * t81;
t127 = -t58 * t72 + t179;
t126 = t72 * t89 + t74 * t81;
t125 = -qJ(5) * t81 - qJD(5) * t89;
t124 = t111 * t72 + (-t149 - t163) * t74;
t121 = -t101 * t72 + t46 * t74;
t94 = -t111 * pkin(4) + t102;
t87 = t153 * t111;
t86 = t153 * t109;
t59 = t60 ^ 2;
t45 = t111 * t50;
t37 = pkin(4) * t109 * t89 + t57;
t36 = t111 * t42;
t21 = -pkin(4) * t163 + t53;
t17 = (t109 * t81 + t141) * pkin(4) + t39;
t16 = -t159 * t89 + t168;
t11 = -pkin(4) * t120 - t109 * t58 - t158 * t89 + t45;
t7 = -t111 * t74 ^ 2 - t165 - t173;
t6 = t124 - t175;
t4 = -qJ(5) * t141 + (-qJD(4) * t58 + t125) * t109 + t145;
t3 = pkin(4) * t78 - t109 * t40 + t36 + t125 * t111 + (-t55 + (qJ(5) * t89 - t50) * t109) * qJD(4);
t8 = [0, 0, 0, 0.2e1 * t112 * t139, t152 * t189, t154, -t155, 0, -pkin(6) * t154 + t110 * t132, pkin(6) * t155 + t112 * t132, t104 * t72 + t78 * t93 + (-t39 + (-qJD(1) * t120 + t76) * t181) * qJD(2), -t104 * t97 + t93 * t81 + (t104 * t98 + t134 * t181 - t40) * qJD(2), t119 * t57 + t120 * t26 + t157 * t39 - t40 * t76 - t51 * t81 - t52 * t78 + t127, t25 * t57 + t26 * t58 - t39 * t51 + t40 * t52 + (t93 + t151) * t143, -t62 * t142 + (t123 * t89 + t62 * t81) * t111, (-t111 * t60 - t166) * t81 + (-t161 - t111 * t28 + (t109 * t60 - t111 * t62) * qJD(4)) * t89, t111 * t126 - t120 * t123 - t142 * t74 + t62 * t78, -t109 * t126 + t120 * t28 - t141 * t74 - t60 * t78, -t120 * t72 + t74 * t78, (-t148 * t58 + t36) * t74 + t45 * t72 - (-t148 * t47 + t23) * t120 + t14 * t78 + t39 * t60 + t57 * t28 + t46 * t141 + ((-qJD(4) * t50 - t40) * t74 - (-qJD(4) * t31 - t26) * t120 + t178 + t127) * t109, -(-t149 * t58 + t145) * t74 - t168 * t72 - t135 * t120 - t15 * t78 + t39 * t62 + t57 * t123 - t46 * t142 + (t178 + t179) * t111, -t1 * t120 + t109 * t128 + t11 * t72 + t141 * t18 + t17 * t60 + t28 * t37 + t3 * t74 + t5 * t78, -t10 * t78 + t111 * t128 + t120 * t2 + t123 * t37 - t142 * t18 - t16 * t72 + t17 * t62 - t4 * t74, -t11 * t123 - t16 * t28 - t3 * t62 - t4 * t60 + (-t10 * t109 - t111 * t5) * t81 + (-t1 * t111 - t109 * t2 + (-t10 * t111 + t109 * t5) * qJD(4)) * t89, t1 * t11 + t10 * t4 + t13 * t37 + t16 * t2 + t17 * t18 + t3 * t5; 0, 0, 0, -t110 * t156, t152 * t114, 0, 0, 0, t114 * pkin(1) * t110, pkin(1) * t156, qJD(2) * t53 - t144 * t76 - t157 * t93 - t25, t54 * qJD(2) - t144 * t157 + t93 * t76 - t26, (t52 - t53) * t157 + (t54 - t51) * t76 + (-t108 * t72 - t119 * t160) * pkin(2), t51 * t53 - t52 * t54 + (t108 * t26 - t150 * t93 - t160 * t25) * pkin(2), t133 * t62 + t161, (t123 - t176) * t111 - t187 + t170, t133 * t74 + t165 - t173, t124 + t175, -t74 * t157, t102 * t28 - t14 * t157 - t35 * t74 - t53 * t60 - t130 * t111 + (t54 * t74 + t121) * t109, t102 * t123 + t109 * t130 + t111 * t121 + t15 * t157 + t169 * t74 - t53 * t62, -t111 * t13 - t21 * t60 + t28 * t94 - t5 * t157 - t72 * t86 + t182 * t74 + (t18 * t76 + (t18 + t183) * qJD(4)) * t109, t76 * t162 + t10 * t157 + t109 * t13 - t21 * t62 + t123 * t94 - t72 * t87 + t171 * t74 + (pkin(4) * t166 + t162) * qJD(4), -t188 * t109 + t129 * t111 + t123 * t86 + t171 * t60 - t182 * t62 - t28 * t87, -t1 * t86 + t13 * t94 + t2 * t87 + t182 * t5 + (pkin(4) * t149 - t21) * t18 - t171 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 * qJD(2), -t97 + (t98 - t76) * qJD(2), -t157 ^ 2 - t76 ^ 2, t157 * t51 + t52 * t76 + t99, 0, 0, 0, 0, 0, t6, t7, t6, t7, (-t123 - t176) * t111 + t187 + t170, t129 * t109 + t188 * t111 - t157 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t59 + t186, t123 + t177, t174 - t28, t72, t15 * t74 - t46 * t62 + t116, t14 * t74 + t46 * t60 + t135, 0.2e1 * t184 + t180 + (-t138 - t18) * t62 + t115, -pkin(4) * t186 + t74 * t9 + (qJD(5) + t18) * t60 + t122, -pkin(4) * t123 - t185 * t60, t185 * t10 + (-t18 * t62 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 + t174, t123 - t177, -t59 - t186, t10 * t60 + t5 * t62 + t13;];
tauc_reg = t8;
