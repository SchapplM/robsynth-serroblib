% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:40
% EndTime: 2019-03-08 22:52:47
% DurationCPUTime: 2.14s
% Computational Cost: add. (1350->246), mult. (3624->402), div. (0->0), fcn. (3048->8), ass. (0->140)
t80 = cos(qJ(4));
t140 = t80 * qJD(5);
t77 = sin(qJ(4));
t151 = t77 * qJ(5);
t169 = pkin(4) + pkin(5);
t176 = t169 * t80 + t151;
t180 = t176 * qJD(4) - t140;
t79 = sin(qJ(2));
t149 = qJD(2) * t79;
t76 = sin(pkin(6));
t124 = t76 * t149;
t82 = cos(qJ(2));
t160 = t76 * t82;
t133 = t77 * t160;
t123 = qJD(2) * t160;
t150 = cos(pkin(6));
t161 = t76 * t79;
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t39 = -t150 * t81 + t78 * t161;
t25 = -t39 * qJD(3) + t81 * t123;
t40 = t150 * t78 + t81 * t161;
t69 = qJD(4) * t80;
t10 = -qJD(4) * t133 - t80 * t124 + t25 * t77 + t40 * t69;
t138 = t81 * qJD(3);
t121 = t77 * t138;
t179 = t10 * t81 + t39 * t121;
t144 = qJD(4) * t81;
t127 = t80 * t144;
t141 = t78 * qJD(3);
t178 = t77 * t141 - t127;
t26 = t40 * qJD(3) + t78 * t123;
t14 = t26 * t77 + t39 * t69;
t177 = t169 * qJD(3);
t175 = -0.4e1 * t78;
t27 = t80 * t160 + t40 * t77;
t174 = t14 * t78 - t27 * t141 + t179;
t74 = t80 ^ 2;
t155 = t77 ^ 2 - t74;
t113 = t155 * qJD(4);
t115 = qJ(6) * t138;
t128 = t78 * t69;
t142 = t77 * qJD(6);
t173 = qJ(6) * t128 + t77 * t115 + t78 * t142;
t167 = pkin(9) * t81;
t107 = pkin(3) * t78 - t167;
t49 = t107 * qJD(3);
t166 = t78 * pkin(9);
t108 = -t81 * pkin(3) - t166;
t53 = -pkin(2) + t108;
t129 = t77 * t144;
t90 = t80 * t141 + t129;
t18 = t90 * pkin(8) - t77 * t49 - t53 * t69;
t11 = -t27 * qJD(4) + t77 * t124 + t25 * t80;
t159 = t78 * t80;
t28 = t40 * t80 - t133;
t118 = t80 * t138;
t145 = qJD(4) * t78;
t91 = -t77 * t145 + t118;
t3 = t11 * t81 - t28 * t141 + t26 * t159 + t91 * t39;
t153 = qJ(5) * t80;
t104 = pkin(4) * t77 - t153;
t100 = pkin(8) + t104;
t35 = t100 * t78;
t143 = t77 * qJD(5);
t37 = t104 * qJD(4) - t143;
t105 = t80 * pkin(4) + t151;
t50 = -pkin(3) - t105;
t171 = qJD(3) * (-t50 * t81 + t166) - qJD(4) * t35 - t37 * t78;
t170 = t105 * qJD(4) - t140;
t84 = 0.2e1 * qJD(5);
t168 = pkin(8) * t81;
t165 = t11 * qJ(5) + t28 * qJD(5);
t97 = -t169 * t77 + t153;
t92 = -pkin(8) + t97;
t12 = t92 * t138 - t180 * t78;
t164 = t12 * t77;
t163 = t12 * t80;
t158 = pkin(9) - qJ(6);
t64 = t80 * t168;
t156 = t77 * t53 + t64;
t73 = t78 ^ 2;
t154 = -t81 ^ 2 + t73;
t152 = qJ(6) * t78;
t148 = qJD(3) * t77;
t147 = qJD(3) * t80;
t146 = qJD(4) * t77;
t139 = t80 * qJD(6);
t137 = t81 * qJD(5);
t136 = qJ(5) * qJD(3);
t63 = t77 * t168;
t135 = -0.2e1 * pkin(2) * qJD(3);
t134 = -0.2e1 * pkin(3) * qJD(4);
t132 = pkin(4) * t141;
t131 = pkin(9) * t146;
t130 = pkin(9) * t69;
t47 = pkin(3) + t176;
t125 = t47 * t69;
t120 = t77 * t69;
t119 = t78 * t138;
t117 = t81 * t136;
t116 = qJ(5) * t145;
t55 = t158 * t80;
t114 = t80 * t53 - t63;
t112 = t154 * qJD(3);
t111 = 0.2e1 * t119;
t110 = -t178 * pkin(8) + t53 * t146 - t80 * t49;
t109 = t77 * t118;
t30 = -t81 * qJ(5) + t156;
t103 = t27 * t80 - t28 * t77;
t70 = t81 * pkin(4);
t31 = -t114 + t70;
t102 = -t30 * t77 + t31 * t80;
t15 = t39 * t146 - t26 * t80;
t34 = t97 * qJD(4) + t143;
t99 = -t47 * t146 + t34 * t80;
t98 = 0.2e1 * t27 * t10 + 0.2e1 * t28 * t11 + 0.2e1 * t39 * t26;
t29 = t92 * t78;
t96 = qJD(4) * t29 + t47 * t138;
t94 = -qJ(6) * t146 + t139;
t89 = -t80 * t115 + t110;
t2 = t103 * qJD(4) + t10 * t77 + t11 * t80;
t67 = t78 * t136;
t13 = t67 - t18 - t137;
t16 = t110 - t132;
t86 = t102 * qJD(4) + t13 * t80 + t16 * t77;
t85 = 0.2e1 * t67 - 0.2e1 * t137 - t18;
t71 = qJ(5) * t84;
t60 = pkin(9) * t127;
t54 = t158 * t77;
t42 = -t121 - t128;
t38 = qJD(4) * t55 - t142;
t36 = -t158 * t146 - t139;
t23 = t77 * t152 + t30;
t20 = t81 * pkin(5) + t63 + t70 + (-t53 - t152) * t80;
t17 = t100 * t138 + t170 * t78;
t7 = t13 + t173;
t6 = (-t94 - t177) * t78 + t89;
t1 = -t103 * t138 + (-t10 * t80 + t11 * t77 + (t27 * t77 + t28 * t80) * qJD(4)) * t78;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, t98; 0, 0, -t124, -t123, 0, 0, 0, 0, 0 (-t82 * t141 - t81 * t149) * t76 (-t82 * t138 + t78 * t149) * t76, 0, 0, 0, 0, 0, t174, t3, t174, -t1, -t3, t10 * t31 + t11 * t30 + t28 * t13 + t27 * t16 + t39 * t17 + t26 * t35 (-qJD(3) * t27 + t14) * t78 + t179, -t3, t1, t10 * t20 + t11 * t23 - t39 * t12 - t26 * t29 + t27 * t6 + t28 * t7; 0, 0, 0, 0, t111, -0.2e1 * t112, 0, 0, 0, t78 * t135, t81 * t135, 0.2e1 * t74 * t119 - 0.2e1 * t73 * t120, t109 * t175 + 0.2e1 * t113 * t73, 0.2e1 * t78 * t129 + 0.2e1 * t154 * t147, -0.2e1 * t112 * t77 + 0.2e1 * t127 * t78, -0.2e1 * t119, 0.2e1 * t110 * t81 + 0.2e1 * t114 * t141 + 0.2e1 * (t111 * t77 + t69 * t73) * pkin(8), -0.2e1 * t18 * t81 - 0.2e1 * t156 * t141 + 0.2e1 * (t111 * t80 - t146 * t73) * pkin(8), 0.2e1 * (t148 * t35 + t16) * t81 + 0.2e1 * (-qJD(3) * t31 + t17 * t77 + t35 * t69) * t78, 0.2e1 * t102 * t138 + 0.2e1 * (-t13 * t77 + t16 * t80 + (-t30 * t80 - t31 * t77) * qJD(4)) * t78, 0.2e1 * (-t147 * t35 - t13) * t81 + 0.2e1 * (qJD(3) * t30 + t146 * t35 - t17 * t80) * t78, 0.2e1 * t30 * t13 + 0.2e1 * t31 * t16 + 0.2e1 * t35 * t17, 0.2e1 * (-t148 * t29 + t6) * t81 + 0.2e1 * (-qJD(3) * t20 - t29 * t69 - t164) * t78, 0.2e1 * (t147 * t29 - t7) * t81 + 0.2e1 * (qJD(3) * t23 - t146 * t29 + t163) * t78, 0.2e1 * (-t20 * t80 + t23 * t77) * t138 + 0.2e1 * (-t6 * t80 + t7 * t77 + (t20 * t77 + t23 * t80) * qJD(4)) * t78, 0.2e1 * t29 * t12 + 0.2e1 * t20 * t6 + 0.2e1 * t23 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, 0, 0, 0, 0, t15, t14, t15, t2, -t14, pkin(9) * t2 + t26 * t50 + t39 * t37, t15, -t14, -t2, t10 * t54 + t11 * t55 - t26 * t47 + t27 * t38 + t28 * t36 - t39 * t34; 0, 0, 0, 0, 0, 0, t138, -t141, 0, -pkin(8) * t138, pkin(8) * t141, -t78 * t113 + t109, t120 * t175 - t155 * t138, t178, t90, 0, t60 + (-pkin(3) * t80 + pkin(8) * t77) * t145 + (t108 * t77 - t64) * qJD(3) (pkin(8) * t159 + t107 * t77) * qJD(4) + (t108 * t80 + t63) * qJD(3), t60 + (t145 * t50 - t17) * t80 - t171 * t77, t86 (-t17 + (t50 * t78 + t167) * qJD(4)) * t77 + t171 * t80, pkin(9) * t86 + t17 * t50 + t35 * t37, t163 + t38 * t81 + (-qJD(3) * t54 - t125) * t78 + (-t34 * t78 - t96) * t77, t164 - t36 * t81 + t96 * t80 + (qJD(3) * t55 + t99) * t78 (-t54 * t138 - t38 * t78 - t7 + (t55 * t78 - t20) * qJD(4)) * t80 + (t55 * t138 + t36 * t78 - t6 + (t54 * t78 + t23) * qJD(4)) * t77, t12 * t47 + t20 * t38 + t23 * t36 + t29 * t34 + t6 * t54 + t7 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t120, -0.2e1 * t113, 0, 0, 0, t77 * t134, t80 * t134, 0.2e1 * t146 * t50 - 0.2e1 * t37 * t80, 0, -0.2e1 * t37 * t77 - 0.2e1 * t50 * t69, 0.2e1 * t50 * t37, 0.2e1 * t99, 0.2e1 * t34 * t77 + 0.2e1 * t125, -0.2e1 * t36 * t80 - 0.2e1 * t38 * t77 + 0.2e1 * (-t54 * t80 + t55 * t77) * qJD(4), 0.2e1 * t47 * t34 + 0.2e1 * t55 * t36 + 0.2e1 * t54 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -t10, 0, t11, -t10 * pkin(4) + t165, -t10, t11, 0, -t10 * t169 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t42, t141, -t110, t18, -t110 + 0.2e1 * t132 (-pkin(4) * t138 - t116) * t80 + (-t117 + (pkin(4) * qJD(4) - qJD(5)) * t78) * t77, t85, -t16 * pkin(4) + t13 * qJ(5) + t30 * qJD(5) (t94 + 0.2e1 * t177) * t78 - t89, t85 + t173 (t138 * t169 + t116) * t80 + (t117 + (-qJD(4) * t169 + qJD(5)) * t78) * t77, t7 * qJ(5) + t23 * qJD(5) - t169 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t146, 0, -t130, t131, -t130, -t170, -t131, -t170 * pkin(9), -t38, t36, t180, t36 * qJ(5) + t55 * qJD(5) - t169 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t71, 0, t84, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t91, 0, t16, -t141, 0, -t91, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t130, 0, 0, -t69, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t91, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t69, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
