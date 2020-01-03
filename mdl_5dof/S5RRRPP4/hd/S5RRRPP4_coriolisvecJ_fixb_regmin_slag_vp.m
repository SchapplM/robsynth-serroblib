% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:45
% EndTime: 2019-12-31 20:55:50
% DurationCPUTime: 1.21s
% Computational Cost: add. (2875->221), mult. (7615->295), div. (0->0), fcn. (5297->6), ass. (0->142)
t114 = qJD(2) + qJD(3);
t117 = sin(pkin(8));
t118 = cos(pkin(8));
t121 = cos(qJ(3));
t122 = cos(qJ(2));
t162 = qJD(1) * t122;
t154 = t121 * t162;
t119 = sin(qJ(3));
t120 = sin(qJ(2));
t163 = qJD(1) * t120;
t155 = t119 * t163;
t75 = -t154 + t155;
t77 = -t119 * t162 - t121 * t163;
t50 = t117 * t77 - t118 * t75;
t187 = t114 * t50;
t138 = -t117 * t75 - t118 * t77;
t186 = t138 ^ 2;
t159 = qJD(1) * qJD(2);
t185 = -0.2e1 * t159;
t182 = pkin(6) + pkin(7);
t97 = t182 * t122;
t92 = qJD(1) * t97;
t82 = t121 * t92;
t96 = t182 * t120;
t90 = qJD(1) * t96;
t146 = t119 * t90 - t82;
t172 = qJ(4) * t75;
t132 = t146 + t172;
t158 = pkin(2) * t117 * t119;
t160 = qJD(3) * t121;
t78 = t119 * t92;
t173 = -t121 * t90 - t78;
t72 = t77 * qJ(4);
t43 = t72 + t173;
t174 = -qJD(3) * t158 - t117 * t132 + (pkin(2) * t160 - t43) * t118;
t171 = qJD(2) * pkin(2);
t84 = -t90 + t171;
t149 = t121 * t84 - t78;
t41 = t149 + t72;
t89 = t119 * t122 + t120 * t121;
t184 = qJD(1) * t89;
t157 = qJD(2) * t182;
t142 = qJD(1) * t157;
t85 = t120 * t142;
t183 = (qJD(3) * t84 - t85) * t121;
t137 = -t119 * t84 - t82;
t86 = t122 * t142;
t148 = t119 * t85 - t121 * t86;
t128 = t137 * qJD(3) + t148;
t150 = t122 * t159;
t53 = qJD(3) * t154 - t114 * t155 + t121 * t150;
t126 = -qJ(4) * t53 + qJD(4) * t77 + t128;
t161 = qJD(3) * t119;
t147 = -t119 * t86 - t92 * t161;
t59 = t114 * t89;
t54 = t59 * qJD(1);
t13 = -qJ(4) * t54 - qJD(4) * t75 + t147 + t183;
t2 = t117 * t13 - t118 * t126;
t110 = -pkin(2) * t122 - pkin(1);
t95 = t110 * qJD(1);
t60 = t75 * pkin(3) + qJD(4) + t95;
t21 = -pkin(4) * t50 - qJ(5) * t138 + t60;
t134 = -t138 * t21 - t2;
t181 = pkin(3) * t77;
t130 = -qJ(4) * t89 - t119 * t97 - t121 * t96;
t136 = t119 * t96 - t121 * t97;
t88 = t119 * t120 - t121 * t122;
t46 = -qJ(4) * t88 - t136;
t28 = t117 * t46 - t118 * t130;
t180 = t2 * t28;
t42 = -t137 - t172;
t38 = t118 * t42;
t18 = t117 * t41 + t38;
t179 = t18 * t138;
t178 = t77 * t75;
t177 = t95 * t77;
t3 = t117 * t126 + t118 * t13;
t176 = -qJD(5) - t174;
t35 = pkin(3) * t114 + t41;
t17 = t117 * t35 + t38;
t169 = t118 * t119;
t175 = t117 * t43 - t118 * t132 - (t117 * t121 + t169) * qJD(3) * pkin(2);
t170 = t117 * t42;
t109 = pkin(2) * t121 + pkin(3);
t71 = pkin(2) * t169 + t117 * t109;
t124 = qJD(1) ^ 2;
t168 = t122 * t124;
t123 = qJD(2) ^ 2;
t167 = t123 * t120;
t166 = t123 * t122;
t19 = t118 * t41 - t170;
t165 = qJD(5) - t19;
t164 = t120 ^ 2 - t122 ^ 2;
t112 = t120 * t171;
t111 = pkin(2) * t163;
t156 = t175 * t138;
t153 = -pkin(2) * t114 - t84;
t152 = pkin(3) * t54 + qJD(2) * t111;
t151 = pkin(3) * t59 + t112;
t30 = t117 * t53 + t118 * t54;
t144 = pkin(1) * t185;
t143 = t21 * t50 + t3;
t16 = t118 * t35 - t170;
t14 = -pkin(4) * t114 + qJD(5) - t16;
t15 = qJ(5) * t114 + t17;
t141 = t138 * t15 - t14 * t50;
t140 = t138 * t17 + t16 * t50;
t139 = -t50 ^ 2 - t186;
t31 = -t117 * t54 + t118 * t53;
t135 = pkin(3) * t88 + t110;
t133 = t95 * t75 - t147;
t70 = t109 * t118 - t158;
t91 = t120 * t157;
t93 = t122 * t157;
t131 = -t119 * t93 - t121 * t91 - t96 * t160 - t97 * t161;
t26 = pkin(4) * t138 - qJ(5) * t50 - t181;
t29 = t117 * t130 + t118 * t46;
t57 = -t117 * t88 + t118 * t89;
t127 = t136 * qJD(3) + t119 * t91 - t121 * t93;
t58 = t114 * t88;
t125 = qJ(4) * t58 - qJD(4) * t89 + t127;
t22 = -qJ(4) * t59 - qJD(4) * t88 + t131;
t6 = t117 * t22 - t118 * t125;
t7 = t117 * t125 + t118 * t22;
t129 = t138 * t6 + t2 * t57 + t28 * t31 - t29 * t30 + t50 * t7;
t5 = pkin(4) * t30 - qJ(5) * t31 - qJD(5) * t138 + t152;
t113 = t114 * qJD(5);
t107 = -pkin(3) * t118 - pkin(4);
t106 = pkin(3) * t117 + qJ(5);
t65 = -pkin(4) - t70;
t64 = qJ(5) + t71;
t56 = t117 * t89 + t118 * t88;
t44 = -t75 ^ 2 + t77 ^ 2;
t37 = (-t184 - t77) * t114;
t36 = t75 * t114 + t53;
t33 = -t117 * t59 - t118 * t58;
t32 = -t117 * t58 + t118 * t59;
t27 = pkin(4) * t56 - qJ(5) * t57 + t135;
t25 = t111 + t26;
t8 = pkin(4) * t32 - qJ(5) * t33 - qJD(5) * t57 + t151;
t1 = t113 + t3;
t4 = [0, 0, 0, 0.2e1 * t120 * t150, t164 * t185, t166, -t167, 0, -pkin(6) * t166 + t120 * t144, pkin(6) * t167 + t122 * t144, t53 * t89 + t58 * t77, -t53 * t88 - t54 * t89 + t58 * t75 + t59 * t77, -t58 * t114, -t59 * t114, 0, t110 * t54 + t95 * t59 + t127 * t114 + (qJD(1) * t88 + t75) * t112, t110 * t53 - t95 * t58 - t131 * t114 + (-t77 + t184) * t112, -t16 * t33 - t17 * t32 - t3 * t56 + t129, t152 * t135 + t60 * t151 - t16 * t6 + t17 * t7 + t3 * t29 + t180, -t114 * t6 + t21 * t32 + t27 * t30 + t5 * t56 - t50 * t8, -t1 * t56 + t14 * t33 - t15 * t32 + t129, t114 * t7 - t138 * t8 - t21 * t33 - t27 * t31 - t5 * t57, t1 * t29 + t14 * t6 + t15 * t7 + t21 * t8 + t27 * t5 + t180; 0, 0, 0, -t120 * t168, t164 * t124, 0, 0, 0, t124 * pkin(1) * t120, pkin(1) * t168, -t178, t44, t36, t37, 0, -t75 * t111 + t177 - t146 * t114 + (t153 * t119 - t82) * qJD(3) + t148, t77 * t111 + t173 * t114 + (t153 * qJD(3) + t85) * t121 + t133, t174 * t50 - t30 * t71 - t31 * t70 + t140 - t156, t3 * t71 - t2 * t70 - t60 * (t111 - t181) + t174 * t17 + t175 * t16, t175 * t114 + t25 * t50 + t134, -t176 * t50 - t30 * t64 + t31 * t65 + t141 - t156, -t176 * t114 + t138 * t25 + t113 + t143, t1 * t64 - t175 * t14 - t176 * t15 + t2 * t65 - t21 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t44, t36, t37, 0, -t137 * t114 + t128 + t177, t149 * t114 + t133 - t183, -t179 - t19 * t50 + (-t117 * t30 - t118 * t31) * pkin(3) + t140, t16 * t18 - t17 * t19 + (t117 * t3 - t118 * t2 + t60 * t77) * pkin(3), t114 * t18 + t26 * t50 + t134, -t106 * t30 + t107 * t31 + t165 * t50 + t141 - t179, -t114 * t19 + t138 * t26 + 0.2e1 * t113 + t143, t1 * t106 + t107 * t2 - t14 * t18 + t15 * t165 - t21 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t138 * t16 - t17 * t50 + t152, t114 * t138 + t30, t139, -t31 - t187, -t138 * t14 - t15 * t50 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138 * t50, t31 - t187, -t114 ^ 2 - t186, -t114 * t15 - t134;];
tauc_reg = t4;
