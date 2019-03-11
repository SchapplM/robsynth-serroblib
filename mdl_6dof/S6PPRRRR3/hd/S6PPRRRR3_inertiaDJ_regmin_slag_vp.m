% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:41
% EndTime: 2019-03-08 19:11:48
% DurationCPUTime: 2.48s
% Computational Cost: add. (2315->271), mult. (7521->522), div. (0->0), fcn. (8128->16), ass. (0->158)
t103 = sin(qJ(5));
t193 = -0.4e1 * t103;
t102 = sin(qJ(6));
t106 = cos(qJ(6));
t92 = t106 ^ 2;
t186 = t102 ^ 2 - t92;
t139 = t186 * qJD(6);
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t101 = cos(pkin(6));
t100 = cos(pkin(7));
t97 = sin(pkin(6));
t187 = t97 * cos(pkin(14));
t157 = t100 * t187;
t96 = sin(pkin(7));
t125 = t101 * t96 + t157;
t188 = sin(pkin(14)) * t97;
t49 = t105 * t125 + t109 * t188;
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t154 = t105 * t188;
t48 = t109 * t125 - t154;
t70 = t100 * t101 - t187 * t96;
t95 = sin(pkin(8));
t99 = cos(pkin(8));
t130 = t48 * t99 + t70 * t95;
t192 = -t49 * t104 + t108 * t130;
t107 = cos(qJ(5));
t179 = t108 * t95;
t158 = pkin(10) * t179;
t65 = t158 + (pkin(3) * t104 + pkin(11)) * t99;
t66 = (-pkin(4) * t108 - pkin(11) * t104 - pkin(3)) * t95;
t128 = t103 * t66 + t107 * t65;
t67 = (pkin(4) * t104 - pkin(11) * t108) * t95 * qJD(4);
t170 = qJD(4) * t104;
t148 = t95 * t170;
t169 = qJD(4) * t108;
t149 = t99 * t169;
t68 = -pkin(3) * t149 + pkin(10) * t148;
t27 = -qJD(5) * t128 + t103 * t68 + t107 * t67;
t191 = 0.2e1 * t95;
t190 = pkin(11) * t95;
t189 = t102 * pkin(11);
t91 = t103 ^ 2;
t185 = -t107 ^ 2 + t91;
t184 = t104 * t95;
t183 = t104 * t99;
t182 = t105 * t96;
t181 = t107 * t95;
t46 = t49 * qJD(3);
t180 = t108 * t46;
t178 = t109 * t96;
t147 = t95 * t169;
t72 = t103 * t184 - t107 * t99;
t52 = -qJD(5) * t72 + t107 * t147;
t73 = t103 * t99 + t104 * t181;
t54 = t102 * t73 + t106 * t179;
t29 = -qJD(6) * t54 + t102 * t148 + t106 * t52;
t177 = t29 * t102;
t176 = t29 * t106;
t174 = t104 * t109;
t173 = t105 * t108;
t172 = t106 * t107;
t171 = qJD(3) * t109;
t168 = qJD(5) * t102;
t167 = qJD(5) * t103;
t166 = qJD(5) * t106;
t165 = qJD(5) * t107;
t164 = qJD(5) * t108;
t163 = qJD(6) * t102;
t162 = qJD(6) * t106;
t161 = qJD(6) * t107;
t160 = -0.2e1 * pkin(4) * qJD(5);
t159 = -0.2e1 * pkin(5) * qJD(6);
t156 = t107 * t189;
t155 = pkin(11) * t172;
t152 = t102 * t179;
t151 = t104 * t182;
t89 = t95 ^ 2;
t150 = t89 * t169;
t146 = qJD(3) * t182;
t145 = t96 * t171;
t144 = t102 * t161;
t143 = t106 * t161;
t142 = t102 * t162;
t141 = t103 * t165;
t140 = t106 * t165;
t138 = t185 * qJD(5);
t137 = 0.2e1 * t141;
t136 = t89 * t146;
t135 = t95 * t146;
t134 = t104 * t150;
t133 = t102 * t140;
t132 = -pkin(5) * t107 - pkin(12) * t103;
t131 = pkin(5) * t103 - pkin(12) * t107;
t23 = t104 * t130 + t49 * t108;
t36 = -t48 * t95 + t70 * t99;
t14 = t103 * t36 + t107 * t23;
t10 = -t102 * t192 + t106 * t14;
t9 = -t102 * t14 - t106 * t192;
t124 = t100 * t95 + t178 * t99;
t51 = t104 * t124 + t173 * t96;
t71 = t100 * t99 - t178 * t95;
t38 = t103 * t71 + t107 * t51;
t50 = -t108 * t124 + t151;
t25 = t102 * t50 + t106 * t38;
t64 = pkin(10) * t184 + (-pkin(3) * t108 - pkin(4)) * t99;
t39 = t72 * pkin(5) - t73 * pkin(12) + t64;
t41 = -pkin(12) * t179 + t128;
t19 = t102 * t39 + t106 * t41;
t55 = t106 * t73 - t152;
t129 = -t102 * t55 - t106 * t54;
t13 = t103 * t23 - t107 * t36;
t127 = -t103 * t65 + t107 * t66;
t85 = -pkin(4) + t132;
t63 = t102 * t85 + t155;
t45 = qJD(3) * t154 - t101 * t145 - t157 * t171;
t12 = qJD(4) * t192 - t45 * t108 - t46 * t183;
t3 = qJD(5) * t14 + t12 * t103 - t181 * t46;
t123 = t102 * t3 + t13 * t162;
t122 = -t106 * t3 + t13 * t163;
t34 = -t149 * t178 - t100 * t147 - t108 * t145 + (qJD(3) * t99 + qJD(4)) * t151;
t17 = qJD(5) * t38 - t103 * t34 - t107 * t135;
t37 = t103 * t51 - t107 * t71;
t121 = t102 * t17 + t162 * t37;
t120 = -t106 * t17 + t163 * t37;
t21 = -pkin(5) * t148 - t27;
t40 = pkin(5) * t179 - t127;
t119 = t21 * t102 + t162 * t40;
t118 = -t21 * t106 + t163 * t40;
t53 = qJD(5) * t73 + t103 * t147;
t117 = t102 * t53 + t162 * t72;
t116 = -t106 * t53 + t163 * t72;
t115 = t131 * t102;
t26 = -t103 * t67 + t107 * t68 - t165 * t66 + t167 * t65;
t114 = t103 * t164 + t107 * t170;
t113 = t103 * t170 - t107 * t164;
t112 = t103 * t166 + t144;
t69 = (pkin(3) * t183 + t158) * qJD(4);
t111 = pkin(12) * t148 - t26;
t110 = t53 * pkin(5) - t52 * pkin(12) + t69;
t62 = t106 * t85 - t156;
t43 = -t63 * qJD(6) + (t103 * t189 + t106 * t131) * qJD(5);
t42 = pkin(11) * t112 - qJD(5) * t115 - t162 * t85;
t35 = t100 * t148 + ((t174 * t99 + t173) * qJD(4) + (t173 * t99 + t174) * qJD(3)) * t96;
t30 = -qJD(6) * t152 + t102 * t52 - t106 * t148 + t162 * t73;
t24 = -t102 * t38 + t106 * t50;
t18 = -t102 * t41 + t106 * t39;
t16 = -t103 * t135 + t107 * t34 - t165 * t71 + t167 * t51;
t11 = qJD(4) * t23 - t45 * t104 + t180 * t99;
t8 = -qJD(6) * t25 + t102 * t16 + t106 * t35;
t7 = -t102 * t35 + t106 * t16 - t162 * t50 + t163 * t38;
t6 = -qJD(6) * t19 - t102 * t111 + t106 * t110;
t5 = -t102 * t110 - t106 * t111 - t162 * t39 + t163 * t41;
t4 = t46 * t95 * t103 - qJD(5) * t13 + t12 * t107;
t2 = qJD(6) * t9 + t11 * t102 + t4 * t106;
t1 = -qJD(6) * t10 - t4 * t102 + t11 * t106;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t46, t45, 0, 0, 0, 0, 0, -t11 * t99 + t148 * t36 - t180 * t89, t104 * t46 * t89 - t12 * t99 + t147 * t36, 0, 0, 0, 0, 0, t11 * t72 - t192 * t53 + (t108 * t3 - t13 * t170) * t95, t11 * t73 - t192 * t52 + (t108 * t4 - t14 * t170) * t95, 0, 0, 0, 0, 0, t1 * t72 + t13 * t30 + t3 * t54 + t53 * t9, -t10 * t53 + t13 * t29 - t2 * t72 + t3 * t55; 0, 0, 0, -t146, -t145, 0, 0, 0, 0, 0, -t108 * t136 + t148 * t71 - t35 * t99, t104 * t136 + t147 * t71 + t34 * t99, 0, 0, 0, 0, 0, t35 * t72 + t50 * t53 + (t108 * t17 - t170 * t37) * t95, t35 * t73 + t50 * t52 + (-t108 * t16 - t170 * t38) * t95, 0, 0, 0, 0, 0, t17 * t54 + t24 * t53 + t30 * t37 + t72 * t8, t17 * t55 - t25 * t53 + t29 * t37 + t7 * t72; 0, 0, 0, 0, 0, 0.2e1 * t134, 0.2e1 * (-t104 ^ 2 + t108 ^ 2) * t89 * qJD(4), 0.2e1 * t99 * t147, -0.2e1 * t99 * t148, 0, -0.2e1 * pkin(3) * t170 * t89 - 0.2e1 * t69 * t99, -0.2e1 * pkin(3) * t150 + 0.2e1 * t68 * t99, 0.2e1 * t73 * t52, -0.2e1 * t52 * t72 - 0.2e1 * t53 * t73 (-t108 * t52 + t170 * t73) * t191 (t108 * t53 - t170 * t72) * t191, -0.2e1 * t134, 0.2e1 * t64 * t53 + 0.2e1 * t69 * t72 + 0.2e1 * (-t108 * t27 + t127 * t170) * t95, 0.2e1 * t64 * t52 + 0.2e1 * t69 * t73 + 0.2e1 * (-t26 * t108 - t128 * t170) * t95, 0.2e1 * t55 * t29, -0.2e1 * t29 * t54 - 0.2e1 * t30 * t55, 0.2e1 * t29 * t72 + 0.2e1 * t53 * t55, -0.2e1 * t30 * t72 - 0.2e1 * t53 * t54, 0.2e1 * t72 * t53, 0.2e1 * t18 * t53 + 0.2e1 * t21 * t54 + 0.2e1 * t30 * t40 + 0.2e1 * t6 * t72, -0.2e1 * t19 * t53 + 0.2e1 * t21 * t55 + 0.2e1 * t29 * t40 + 0.2e1 * t5 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t107 * t11 - t167 * t192, t103 * t11 - t165 * t192, 0, 0, 0, 0, 0 (t13 * t168 - t1) * t107 + (qJD(5) * t9 + t123) * t103 (t13 * t166 + t2) * t107 + (-qJD(5) * t10 - t122) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, 0, 0, 0, 0, -t107 * t35 + t167 * t50, t103 * t35 + t165 * t50, 0, 0, 0, 0, 0 (t168 * t37 - t8) * t107 + (qJD(5) * t24 + t121) * t103 (t166 * t37 - t7) * t107 + (-qJD(5) * t25 - t120) * t103; 0, 0, 0, 0, 0, 0, 0, t147, -t148, 0, -t69, t68, t103 * t52 + t165 * t73, -t103 * t53 + t52 * t107 + (-t103 * t73 - t107 * t72) * qJD(5), t113 * t95, t114 * t95, 0, -pkin(4) * t53 - t69 * t107 - t113 * t190 + t167 * t64, -pkin(4) * t52 + t69 * t103 - t114 * t190 + t165 * t64, t55 * t140 + (-t163 * t55 + t176) * t103, t129 * t165 + (-t177 - t106 * t30 + (t102 * t54 - t106 * t55) * qJD(6)) * t103 (t166 * t72 - t29) * t107 + (qJD(5) * t55 - t116) * t103 (-t168 * t72 + t30) * t107 + (-qJD(5) * t54 - t117) * t103, -t107 * t53 + t167 * t72, t43 * t72 + t62 * t53 + (-t6 + (pkin(11) * t54 + t102 * t40) * qJD(5)) * t107 + (pkin(11) * t30 + qJD(5) * t18 + t119) * t103, t42 * t72 - t63 * t53 + (-t5 + (pkin(11) * t55 + t106 * t40) * qJD(5)) * t107 + (pkin(11) * t29 - qJD(5) * t19 - t118) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -0.2e1 * t138, 0, 0, 0, t103 * t160, t107 * t160, 0.2e1 * t141 * t92 - 0.2e1 * t142 * t91, t133 * t193 + 0.2e1 * t139 * t91, 0.2e1 * t103 * t144 + 0.2e1 * t166 * t185, -0.2e1 * t102 * t138 + 0.2e1 * t103 * t143, -0.2e1 * t141, 0.2e1 * t62 * t167 - 0.2e1 * t43 * t107 + 0.2e1 * (t102 * t137 + t162 * t91) * pkin(11), -0.2e1 * t63 * t167 - 0.2e1 * t42 * t107 + 0.2e1 * (t106 * t137 - t163 * t91) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, t122, t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, 0, 0, 0, t120, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t53, t148, t27, t26, t162 * t55 + t177, qJD(6) * t129 - t102 * t30 + t176, t117, -t116, 0, -pkin(5) * t30 - pkin(12) * t117 + t118, -pkin(5) * t29 + pkin(12) * t116 + t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t167, 0, -pkin(11) * t165, pkin(11) * t167, -t103 * t139 + t133, t142 * t193 - t165 * t186, t102 * t167 - t143, t112, 0 (pkin(12) * t172 + (-pkin(5) * t106 + t189) * t103) * qJD(6) + (t102 * t132 - t155) * qJD(5) (pkin(11) * t103 * t106 + t115) * qJD(6) + (t106 * t132 + t156) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t142, -0.2e1 * t139, 0, 0, 0, t102 * t159, t106 * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t30, t53, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * t163 + t140, -t102 * t165 - t103 * t162, t167, t43, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, -t163, 0, -pkin(12) * t162, pkin(12) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
