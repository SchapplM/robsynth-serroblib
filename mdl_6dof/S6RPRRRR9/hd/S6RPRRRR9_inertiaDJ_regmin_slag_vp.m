% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x34]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:19
% EndTime: 2019-03-09 07:25:26
% DurationCPUTime: 2.32s
% Computational Cost: add. (2875->253), mult. (6822->450), div. (0->0), fcn. (6228->8), ass. (0->155)
t109 = sin(qJ(4));
t110 = sin(qJ(3));
t157 = t110 * qJD(3);
t142 = t109 * t157;
t113 = cos(qJ(4));
t114 = cos(qJ(3));
t163 = qJD(4) * t114;
t143 = t113 * t163;
t68 = t142 - t143;
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t78 = t108 * t109 - t112 * t113;
t60 = t78 * t110;
t79 = t108 * t113 + t109 * t112;
t59 = t79 * t114;
t139 = t113 * t157;
t145 = t109 * t163;
t66 = -t139 - t145;
t104 = t110 ^ 2;
t106 = t114 ^ 2;
t127 = (t104 - t106) * qJD(3);
t105 = t113 ^ 2;
t168 = t109 ^ 2 - t105;
t128 = t168 * qJD(4);
t185 = qJD(4) + qJD(5);
t184 = 2 * qJD(2);
t183 = pkin(8) + pkin(9);
t182 = t110 * pkin(5);
t181 = t112 * pkin(4);
t115 = -pkin(1) - pkin(7);
t172 = t109 * t115;
t133 = pkin(4) - t172;
t170 = t113 * t114;
t123 = pkin(3) * t110 - pkin(8) * t114;
t85 = qJ(2) + t123;
t75 = t113 * t85;
t46 = -pkin(9) * t170 + t110 * t133 + t75;
t173 = t109 * t114;
t74 = t109 * t85;
t171 = t110 * t113;
t94 = t115 * t171;
t178 = t94 + t74;
t54 = -pkin(9) * t173 + t178;
t48 = t112 * t54;
t180 = t108 * t46 + t48;
t92 = t183 * t109;
t93 = t183 * t113;
t179 = -t108 * t92 + t112 * t93;
t111 = cos(qJ(6));
t20 = -pkin(10) * t59 + t180;
t177 = t111 * t20;
t53 = t185 * t79;
t176 = t53 * t110;
t107 = sin(qJ(6));
t175 = t107 * t108;
t174 = t108 * t111;
t169 = t114 * t115;
t166 = t104 + t106;
t165 = qJD(4) * t109;
t164 = qJD(4) * t113;
t162 = qJD(4) * t115;
t161 = qJD(5) * t108;
t160 = qJD(5) * t112;
t159 = qJD(6) * t107;
t158 = qJD(6) * t111;
t102 = t114 * qJD(3);
t156 = -0.2e1 * pkin(3) * qJD(4);
t155 = pkin(5) * t102;
t154 = pkin(4) * t165;
t153 = pkin(4) * t161;
t152 = pkin(4) * t160;
t151 = pkin(5) * t159;
t150 = pkin(5) * t158;
t149 = t78 * t157;
t148 = t79 * t102;
t101 = -pkin(4) * t113 - pkin(3);
t29 = -t185 * t59 + t149;
t124 = pkin(3) * t114 + pkin(8) * t110;
t76 = qJD(3) * t124 + qJD(2);
t63 = t113 * t76;
t24 = t63 + (-t94 + (pkin(9) * t114 - t85) * t109) * qJD(4) + (pkin(9) * t171 + t114 * t133) * qJD(3);
t136 = t115 * t102;
t144 = t109 * t162;
t34 = -t109 * t76 + t110 * t144 - t113 * t136 - t164 * t85;
t26 = pkin(9) * t68 - t34;
t132 = -t108 * t26 + t112 * t24;
t9 = -qJD(5) * t180 + t132;
t4 = -t29 * pkin(10) + t155 + t9;
t31 = -t161 * t173 + (t170 * t185 - t142) * t112 + t66 * t108;
t8 = -t108 * t24 - t112 * t26 - t160 * t46 + t161 * t54;
t5 = -pkin(10) * t31 - t8;
t147 = -t107 * t5 + t111 * t4;
t146 = qJD(4) * t183;
t141 = t109 * t102;
t140 = t109 * t164;
t138 = t113 * t102;
t137 = t110 * t102;
t131 = -t108 * t54 + t112 * t46;
t61 = t78 * t114;
t17 = pkin(10) * t61 + t131 + t182;
t135 = -t17 - t182;
t16 = t20 * t159;
t134 = -t107 * t4 + t16;
t130 = -t108 * t93 - t112 * t92;
t100 = pkin(5) + t181;
t129 = qJD(6) * (-pkin(5) - t100);
t77 = pkin(4) * t173 - t169;
t126 = t109 * t136;
t125 = t109 * t139;
t122 = t107 * t17 + t177;
t38 = -pkin(10) * t79 + t130;
t39 = -pkin(10) * t78 + t179;
t121 = t107 * t39 - t111 * t38;
t120 = t107 * t38 + t111 * t39;
t58 = t79 * t110;
t119 = -t107 * t60 + t111 * t58;
t118 = -t107 * t58 - t111 * t60;
t36 = -t107 * t61 + t111 * t59;
t37 = -t107 * t59 - t111 * t61;
t50 = t107 * t79 + t111 * t78;
t51 = -t107 * t78 + t111 * t79;
t97 = t115 * t157;
t55 = -pkin(4) * t68 + t97;
t83 = t109 * t146;
t84 = t113 * t146;
t32 = t108 * t84 + t112 * t83 + t160 * t92 + t161 * t93;
t2 = -qJD(6) * t122 + t147;
t33 = -qJD(5) * t179 + t108 * t83 - t112 * t84;
t117 = (t108 * t159 + (-t111 * t112 + t175) * qJD(5)) * pkin(4);
t116 = (-t108 * t158 + (-t107 * t112 - t174) * qJD(5)) * pkin(4);
t96 = 0.2e1 * t137;
t67 = t110 * t164 + t141;
t65 = t110 * t165 - t138;
t57 = pkin(5) * t78 + t101;
t52 = t185 * t78;
t49 = pkin(5) * t59 + t77;
t45 = -t100 * t159 + t116;
t44 = -t100 * t158 + t117;
t41 = pkin(5) * t53 + t154;
t35 = -qJD(4) * t178 - t126 + t63;
t30 = t185 * t60 - t148;
t28 = t108 * t141 - t112 * t138 + t176;
t21 = pkin(5) * t31 + t55;
t19 = t52 * pkin(10) + t33;
t18 = -pkin(10) * t53 - t32;
t15 = qJD(6) * t51 - t107 * t52 + t111 * t53;
t14 = -qJD(6) * t50 - t107 * t53 - t111 * t52;
t13 = qJD(6) * t37 + t107 * t29 + t111 * t31;
t12 = -qJD(6) * t118 + t107 * t28 + t111 * t30;
t11 = -qJD(6) * t36 - t107 * t31 + t111 * t29;
t10 = qJD(6) * t119 - t107 * t30 + t111 * t28;
t7 = -qJD(6) * t120 - t107 * t18 + t111 * t19;
t6 = qJD(6) * t121 - t107 * t19 - t111 * t18;
t1 = (-qJD(6) * t17 - t5) * t111 + t134;
t3 = [0, 0, 0, 0, t184, qJ(2) * t184, -0.2e1 * t137, 0.2e1 * t127, 0, 0, 0, 0.2e1 * qJ(2) * t102 + 0.2e1 * qJD(2) * t110, -0.2e1 * qJ(2) * t157 + 0.2e1 * qJD(2) * t114, -0.2e1 * t105 * t137 - 0.2e1 * t106 * t140, 0.2e1 * t106 * t128 + 0.4e1 * t114 * t125, -0.2e1 * t110 * t145 - 0.2e1 * t113 * t127, 0.2e1 * t109 * t127 - 0.2e1 * t110 * t143, t96, -0.2e1 * t106 * t113 * t162 + 0.2e1 * t75 * t102 + 0.2e1 * (t35 + t126) * t110, 0.2e1 * t106 * t144 + 0.2e1 * t34 * t110 + 0.2e1 * (-t74 + t94) * t102, -0.2e1 * t61 * t29, -0.2e1 * t29 * t59 + 0.2e1 * t31 * t61, -0.2e1 * t102 * t61 + 0.2e1 * t110 * t29, -0.2e1 * t102 * t59 - 0.2e1 * t110 * t31, t96, 0.2e1 * t102 * t131 + 0.2e1 * t110 * t9 + 0.2e1 * t31 * t77 + 0.2e1 * t55 * t59, -0.2e1 * t102 * t180 + 0.2e1 * t110 * t8 + 0.2e1 * t29 * t77 - 0.2e1 * t55 * t61, 0.2e1 * t37 * t11, -0.2e1 * t11 * t36 - 0.2e1 * t13 * t37, 0.2e1 * t102 * t37 + 0.2e1 * t11 * t110, -0.2e1 * t102 * t36 - 0.2e1 * t110 * t13, t96, 0.2e1 * t2 * t110 + 0.2e1 * (-t107 * t20 + t111 * t17) * t102 + 0.2e1 * t21 * t36 + 0.2e1 * t49 * t13, 0.2e1 * t1 * t110 - 0.2e1 * t102 * t122 + 0.2e1 * t11 * t49 + 0.2e1 * t21 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166 * t164, t166 * t165, 0, 0, 0, 0, 0, t30 * t110 - t114 * t31 + (t110 * t59 - t114 * t58) * qJD(3), t28 * t110 - t114 * t29 + (-t110 * t61 + t114 * t60) * qJD(3), 0, 0, 0, 0, 0, t12 * t110 - t114 * t13 + (t110 * t36 - t114 * t119) * qJD(3), t10 * t110 - t114 * t11 + (t110 * t37 - t114 * t118) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t157, -t102, 0, -t97, -t136, -t114 * t128 - t125, -0.4e1 * t114 * t140 + t157 * t168, t67, -t65, 0 (-t109 * t169 - t113 * t124) * qJD(4) + (t109 * t123 - t94) * qJD(3) (t109 * t124 - t113 * t169) * qJD(4) + (-pkin(8) * t170 + (pkin(3) * t113 + t172) * t110) * qJD(3), t29 * t79 + t52 * t61, -t29 * t78 - t31 * t79 + t52 * t59 + t53 * t61, -t110 * t52 + t148, -t102 * t78 - t176, 0, t101 * t31 + t102 * t130 + t110 * t33 + t154 * t59 + t53 * t77 + t55 * t78, t101 * t29 - t102 * t179 + t110 * t32 - t154 * t61 - t52 * t77 + t55 * t79, t11 * t51 + t14 * t37, -t11 * t50 - t13 * t51 - t14 * t36 - t15 * t37, t102 * t51 + t110 * t14, -t102 * t50 - t110 * t15, 0, -t102 * t121 + t110 * t7 + t13 * t57 + t15 * t49 + t21 * t50 + t36 * t41, -t102 * t120 + t11 * t57 + t110 * t6 + t14 * t49 + t21 * t51 + t37 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, -t102, 0, 0, 0, 0, 0, t66, t68, 0, 0, 0, 0, 0, -t114 * t53 + t149, t114 * t52 + t157 * t79, 0, 0, 0, 0, 0, -t114 * t15 + t157 * t50, -t114 * t14 + t157 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t140, -0.2e1 * t128, 0, 0, 0, t109 * t156, t113 * t156, -0.2e1 * t79 * t52, 0.2e1 * t52 * t78 - 0.2e1 * t53 * t79, 0, 0, 0, 0.2e1 * t101 * t53 + 0.2e1 * t154 * t78, -0.2e1 * t101 * t52 + 0.2e1 * t154 * t79, 0.2e1 * t51 * t14, -0.2e1 * t14 * t50 - 0.2e1 * t15 * t51, 0, 0, 0, 0.2e1 * t15 * t57 + 0.2e1 * t41 * t50, 0.2e1 * t14 * t57 + 0.2e1 * t41 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t68, t102, t35, t34, 0, 0, t29, -t31, t102, t102 * t181 + (-t48 + (-pkin(4) * t110 - t46) * t108) * qJD(5) + t132 (-t102 * t108 - t110 * t160) * pkin(4) + t8, 0, 0, t11, -t13, t102, t45 * t110 + (-pkin(4) * t175 + t100 * t111) * t102 + t2, t44 * t110 - (pkin(4) * t174 + t100 * t107) * t102 - t111 * t5 - t17 * t158 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t65, 0, 0, 0, 0, 0, t30, t28, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t165, 0, -pkin(8) * t164, pkin(8) * t165, 0, 0, -t52, -t53, 0, t33, t32, 0, 0, t14, -t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t153, -0.2e1 * t152, 0, 0, 0, 0, 0, 0.2e1 * t45, 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t31, t102, t9, t8, 0, 0, t11, -t13, t102, t111 * t155 + (t107 * t135 - t177) * qJD(6) + t147, t16 + (-t4 - t155) * t107 + (qJD(6) * t135 - t5) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t28, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, t33, t32, 0, 0, t14, -t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t152, 0, 0, 0, 0, 0, t107 * t129 + t116, t111 * t129 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t151, -0.2e1 * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t13, t102, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, -t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
