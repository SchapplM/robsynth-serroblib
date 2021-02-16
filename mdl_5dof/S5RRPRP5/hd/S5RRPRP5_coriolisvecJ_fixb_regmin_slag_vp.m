% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP5
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
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:18
% EndTime: 2021-01-15 20:18:25
% DurationCPUTime: 1.35s
% Computational Cost: add. (2525->220), mult. (6658->290), div. (0->0), fcn. (4796->6), ass. (0->136)
t128 = sin(qJ(2));
t194 = 0.2e1 * t128;
t127 = sin(qJ(4));
t125 = sin(pkin(8));
t126 = cos(pkin(8));
t129 = cos(qJ(2));
t103 = t125 * t129 + t126 * t128;
t157 = qJD(1) * t103;
t165 = t126 * t129;
t102 = t125 * t128 - t165;
t158 = qJD(1) * t102;
t181 = cos(qJ(4));
t140 = t127 * t158 - t157 * t181;
t153 = qJD(1) * qJD(2);
t149 = t129 * t153;
t150 = t128 * t153;
t160 = t125 * t149 + t126 * t150;
t113 = t125 * t150;
t84 = t126 * t149 - t113;
t133 = t140 * qJD(4) - t127 * t84 - t181 * t160;
t122 = qJD(2) + qJD(4);
t168 = t140 * t122;
t193 = t133 - t168;
t151 = qJD(4) * t181;
t154 = qJD(4) * t127;
t139 = -t127 * t160 - t151 * t158 - t154 * t157 + t181 * t84;
t49 = -t127 * t157 - t158 * t181;
t167 = t49 * t122;
t7 = t139 - t167;
t175 = t49 ^ 2;
t191 = t140 ^ 2;
t192 = -t175 + t191;
t190 = pkin(2) * t194;
t118 = -t129 * pkin(2) - pkin(1);
t156 = qJD(1) * t118;
t110 = qJD(3) + t156;
t58 = pkin(3) * t158 + t110;
t12 = -pkin(4) * t49 + qJ(5) * t140 + t58;
t189 = t12 * t49;
t188 = t58 * t49;
t178 = t12 * t140;
t173 = t140 * t49;
t187 = t140 * t58;
t21 = -pkin(4) * t140 - t49 * qJ(5);
t186 = -0.2e1 * t153;
t117 = t126 * pkin(2) + pkin(3);
t180 = pkin(2) * t125;
t152 = t127 * t180;
t183 = t158 * pkin(7);
t172 = -qJ(3) - pkin(6);
t111 = t172 * t128;
t107 = qJD(1) * t111;
t112 = t172 * t129;
t108 = qJD(1) * t112;
t166 = t126 * t108;
t56 = -t125 * t107 + t166;
t40 = t56 + t183;
t182 = t157 * pkin(7);
t97 = t125 * t108;
t57 = t126 * t107 + t97;
t41 = t57 - t182;
t185 = qJD(4) * t152 - t117 * t151 + t127 * t40 + t181 * t41;
t60 = t126 * t111 + t125 * t112;
t44 = -t103 * pkin(7) + t60;
t61 = t125 * t111 - t126 * t112;
t45 = -t102 * pkin(7) + t61;
t141 = -t127 * t45 + t181 * t44;
t137 = t102 * qJD(2);
t147 = qJD(2) * t172;
t90 = t129 * qJD(3) + t128 * t147;
t91 = -t128 * qJD(3) + t129 * t147;
t42 = -t125 * t90 + t126 * t91;
t32 = pkin(7) * t137 + t42;
t138 = t103 * qJD(2);
t43 = t125 * t91 + t126 * t90;
t33 = -pkin(7) * t138 + t43;
t4 = t141 * qJD(4) + t127 * t32 + t181 * t33;
t177 = t4 * t122;
t18 = t127 * t44 + t181 * t45;
t5 = t18 * qJD(4) + t127 * t33 - t181 * t32;
t174 = t5 * t122;
t136 = t127 * t117 + t181 * t180;
t171 = -t136 * qJD(4) + t127 * t41 - t181 * t40;
t170 = -qJD(5) + t185;
t71 = t90 * qJD(1);
t72 = t91 * qJD(1);
t39 = t125 * t72 + t126 * t71;
t169 = qJD(2) * pkin(2);
t101 = t107 + t169;
t53 = t125 * t101 - t166;
t131 = qJD(1) ^ 2;
t164 = t129 * t131;
t130 = qJD(2) ^ 2;
t163 = t130 * t128;
t162 = t130 * t129;
t52 = t126 * t101 + t97;
t36 = qJD(2) * pkin(3) - t182 + t52;
t37 = t53 - t183;
t10 = -t127 * t37 + t181 * t36;
t161 = qJD(5) - t10;
t159 = t128 ^ 2 - t129 ^ 2;
t155 = qJD(1) * t128;
t120 = t128 * t169;
t119 = pkin(2) * t155;
t65 = pkin(3) * t157 + t119;
t38 = -t125 * t71 + t126 * t72;
t28 = -t84 * pkin(7) + t38;
t29 = -t160 * pkin(7) + t39;
t146 = -t127 * t28 - t36 * t151 + t37 * t154 - t181 * t29;
t2 = t127 * t29 + t37 * t151 + t36 * t154 - t181 * t28;
t145 = pkin(1) * t186;
t121 = t122 * qJD(5);
t1 = t121 - t146;
t116 = pkin(2) * t150;
t59 = t160 * pkin(3) + t116;
t74 = t102 * pkin(3) + t118;
t143 = t10 * t122 + t146;
t11 = t127 * t36 + t181 * t37;
t142 = t11 * t122 - t2;
t55 = -t127 * t102 + t181 * t103;
t135 = t171 * t122 - t2;
t134 = t139 + t167;
t66 = pkin(3) * t138 + t120;
t3 = -pkin(4) * t133 - qJ(5) * t139 + qJD(5) * t140 + t59;
t132 = -t133 - t168;
t88 = -t181 * t117 - pkin(4) + t152;
t87 = qJ(5) + t136;
t54 = t181 * t102 + t127 * t103;
t23 = t55 * qJD(4) - t127 * t137 + t181 * t138;
t22 = t102 * t151 + t103 * t154 + t127 * t138 + t181 * t137;
t16 = t54 * pkin(4) - t55 * qJ(5) + t74;
t15 = t21 + t65;
t9 = t122 * qJ(5) + t11;
t8 = -t122 * pkin(4) + t161;
t6 = t23 * pkin(4) + t22 * qJ(5) - t55 * qJD(5) + t66;
t13 = [0, 0, 0, t149 * t194, t159 * t186, t162, -t163, 0, -pkin(6) * t162 + t128 * t145, pkin(6) * t163 + t129 * t145, t118 * t160 + (t110 * t103 + t158 * t190 + t42) * qJD(2), t118 * t84 + (-t110 * t102 + t157 * t190 - t43) * qJD(2), -t43 * t158 - t61 * t160 - t39 * t102 - t42 * t157 - t60 * t84 - t38 * t103 + (t52 * t102 - t53 * t103) * qJD(2), t38 * t60 + t39 * t61 + t52 * t42 + t53 * t43 + (t110 + t156) * t120, t139 * t55 + t140 * t22, t133 * t55 - t139 * t54 + t140 * t23 - t22 * t49, -t22 * t122, -t23 * t122, 0, -t133 * t74 + t58 * t23 - t49 * t66 + t59 * t54 - t174, t139 * t74 - t140 * t66 - t58 * t22 + t59 * t55 - t177, t12 * t23 - t133 * t16 + t3 * t54 - t49 * t6 - t174, -t1 * t54 + t133 * t18 - t139 * t141 - t140 * t5 + t2 * t55 - t8 * t22 - t9 * t23 + t4 * t49, t12 * t22 - t139 * t16 + t140 * t6 - t3 * t55 + t177, t1 * t18 + t12 * t6 - t141 * t2 + t3 * t16 + t9 * t4 + t8 * t5; 0, 0, 0, -t128 * t164, t159 * t131, 0, 0, 0, t131 * pkin(1) * t128, pkin(1) * t164, -t56 * qJD(2) - t110 * t157 - t119 * t158 + t38, t57 * qJD(2) + t110 * t158 - t119 * t157 - t39, (t53 + t56) * t157 + (t57 - t52) * t158 + (-t125 * t160 - t126 * t84) * pkin(2), -t52 * t56 - t53 * t57 + (-t110 * t155 + t125 * t39 + t126 * t38) * pkin(2), t173, t192, t7, t193, 0, t49 * t65 + t135 + t187, t185 * t122 + t140 * t65 + t146 - t188, t15 * t49 + t135 + t178, t139 * t88 + t87 * t133 + (-t170 - t8) * t49 + (t171 - t9) * t140, -t170 * t122 - t140 * t15 + t1 + t189, t1 * t87 - t12 * t15 - t170 * t9 - t171 * t8 + t2 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t157 + t160, -t113 + (qJD(1) * t165 - t158) * qJD(2), -t157 ^ 2 - t158 ^ 2, t157 * t52 + t158 * t53 + t116, 0, 0, 0, 0, 0, t132, t134, t132, -t175 - t191, -t134, t140 * t8 - t9 * t49 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t192, t7, t193, 0, t142 + t187, t143 - t188, t21 * t49 + t142 + t178, -pkin(4) * t139 + t133 * qJ(5) - (-t11 + t9) * t140 - (t8 - t161) * t49, -t140 * t21 + 0.2e1 * t121 - t143 + t189, -t2 * pkin(4) + t1 * qJ(5) - t8 * t11 - t12 * t21 + t161 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t7, -t122 ^ 2 - t191, -t9 * t122 - t178 + t2;];
tauc_reg = t13;
