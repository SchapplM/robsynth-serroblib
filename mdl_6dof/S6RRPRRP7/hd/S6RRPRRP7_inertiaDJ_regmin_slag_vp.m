% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:06
% EndTime: 2019-03-09 12:20:12
% DurationCPUTime: 1.90s
% Computational Cost: add. (2507->216), mult. (5274->371), div. (0->0), fcn. (4589->6), ass. (0->131)
t78 = sin(qJ(5));
t76 = t78 ^ 2;
t81 = cos(qJ(5));
t77 = t81 ^ 2;
t139 = t76 - t77;
t163 = qJD(5) * t139;
t169 = 0.2e1 * t163;
t153 = cos(qJ(4));
t80 = sin(qJ(2));
t119 = t80 * t153;
t79 = sin(qJ(4));
t82 = cos(qJ(2));
t53 = -t82 * t79 + t119;
t115 = qJD(4) * t153;
t137 = qJD(4) * t79;
t70 = t82 * qJD(2);
t164 = qJD(2) * t119 - t79 * t70;
t28 = t80 * t115 - t82 * t137 - t164;
t97 = -t82 * t153 - t80 * t79;
t29 = (-qJD(2) + qJD(4)) * t97;
t96 = -t153 * t29 - t28 * t79;
t168 = (qJD(4) * (-t153 * t97 - t53 * t79) - t96) * t81;
t166 = -t82 * pkin(2) - t80 * qJ(3);
t58 = -pkin(1) + t166;
t48 = t82 * pkin(3) - t58;
t24 = -pkin(4) * t97 - pkin(9) * t53 + t48;
t157 = pkin(7) - pkin(8);
t60 = t157 * t80;
t61 = t157 * t82;
t32 = t153 * t61 + t79 * t60;
t167 = t78 * t24 + t81 * t32;
t165 = t76 + t77;
t158 = -pkin(2) - pkin(3);
t93 = t153 * qJ(3) + t79 * t158;
t56 = -pkin(9) + t93;
t149 = t97 * t56;
t17 = t97 * qJD(2) * t157 + t32 * qJD(4);
t111 = t153 * t158;
t55 = t79 * qJ(3) + pkin(4) - t111;
t161 = qJD(5) * (-t53 * t55 - t149) - t17;
t132 = t80 * qJD(2);
t140 = qJ(3) * t70 + t80 * qJD(3);
t33 = t158 * t132 + t140;
t12 = pkin(4) * t28 - pkin(9) * t29 + t33;
t16 = -t60 * t115 + t61 * t137 + t164 * t157;
t5 = -qJD(5) * t167 + t81 * t12 + t16 * t78;
t104 = t78 * pkin(5) - t81 * qJ(6);
t39 = t104 * qJD(5) - t78 * qJD(6);
t105 = t81 * pkin(5) + t78 * qJ(6);
t38 = t105 * qJD(5) - t81 * qJD(6);
t160 = 0.2e1 * qJD(3);
t159 = 0.2e1 * qJD(6);
t156 = pkin(5) * t28;
t155 = pkin(9) * t28;
t154 = pkin(9) * t97;
t37 = t79 * qJD(3) + t93 * qJD(4);
t152 = t37 * t78;
t151 = t37 * t81;
t150 = t39 * t53;
t148 = t53 * t29;
t147 = t53 * t78;
t146 = t53 * t81;
t144 = t81 * t29;
t27 = t37 - t39;
t143 = -t27 + t39;
t34 = t105 + t55;
t57 = -pkin(4) - t105;
t141 = t34 - t57;
t138 = qJ(6) * t28;
t31 = -t153 * t60 + t79 * t61;
t15 = t104 * t53 + t31;
t136 = qJD(5) * t15;
t135 = qJD(5) * t78;
t71 = qJD(5) * t81;
t134 = qJD(6) * t97;
t130 = -0.2e1 * pkin(1) * qJD(2);
t129 = -0.2e1 * pkin(4) * qJD(5);
t128 = t78 * t144;
t127 = pkin(9) * t135;
t126 = pkin(9) * t71;
t125 = pkin(7) * t132;
t124 = pkin(7) * t70;
t122 = t78 * t137;
t121 = t78 * t71;
t120 = t78 * t153;
t118 = t81 * t153;
t36 = qJ(3) * t137 - t153 * qJD(3) - qJD(4) * t111;
t23 = t165 * t36;
t114 = qJD(5) * (pkin(4) + t55);
t113 = qJD(5) * t153;
t109 = -pkin(4) * t29 - t155;
t108 = pkin(4) * t53 - t154;
t7 = -qJ(6) * t97 + t167;
t103 = t24 * t81 - t32 * t78;
t8 = pkin(5) * t97 - t103;
t107 = -t7 * t81 - t78 * t8;
t106 = -t7 * t78 + t8 * t81;
t101 = t28 * t56 + t36 * t97;
t99 = t78 * t29 + t53 * t71;
t98 = t53 * t135 - t144;
t94 = t153 * t53 - t79 * t97;
t4 = -t78 * t12 + t32 * t135 + t81 * t16 - t24 * t71;
t92 = t165 * t153;
t6 = t104 * t29 + t38 * t53 + t17;
t91 = -t6 + (t53 * t57 + t154) * qJD(5);
t90 = qJD(5) * t94;
t89 = t6 + (t34 * t53 + t149) * qJD(5);
t88 = -t29 * t57 - t136 - t150 + t155;
t87 = t166 * qJD(2) + qJD(3) * t82;
t2 = -t134 - t4 + t138;
t3 = -t156 - t5;
t1 = t106 * qJD(5) + t2 * t81 + t3 * t78;
t85 = -t27 * t53 - t29 * t34 + t101 + t136;
t84 = -qJD(5) * t31 + t29 * t55 + t37 * t53 - t101;
t83 = t53 * t122 - t81 * t90 + (t115 * t97 + t96) * t78;
t62 = 0.2e1 * t121;
t51 = -0.2e1 * t163;
t50 = t53 ^ 2;
t45 = t78 * t113 + t81 * t137;
t44 = -t81 * t113 + t122;
t43 = t78 * t115 + t79 * t71;
t42 = -t81 * t115 + t79 * t135;
t41 = t92 * qJD(4);
t40 = pkin(2) * t132 - t140;
t26 = -t36 * t78 + t56 * t71;
t25 = t56 * t135 + t36 * t81;
t20 = t28 * t78 - t71 * t97;
t18 = -t135 * t97 - t28 * t81;
t13 = t163 * t53 - t128;
t11 = 0.4e1 * t53 * t121 + t139 * t29;
t9 = [0, 0, 0, 0.2e1 * t80 * t70, 0.2e1 * (-t80 ^ 2 + t82 ^ 2) * qJD(2), 0, 0, 0, t80 * t130, t82 * t130, 0.2e1 * t58 * t132 - 0.2e1 * t40 * t82, 0, -0.2e1 * t40 * t80 - 0.2e1 * t58 * t70, 0.2e1 * t58 * t40, 0.2e1 * t148, -0.2e1 * t28 * t53 + 0.2e1 * t29 * t97, 0, 0, 0, 0.2e1 * t28 * t48 - 0.2e1 * t33 * t97, 0.2e1 * t29 * t48 + 0.2e1 * t33 * t53, -0.2e1 * t121 * t50 + 0.2e1 * t77 * t148, -0.4e1 * t128 * t53 + t50 * t169, 0.2e1 * t28 * t146 + 0.2e1 * t97 * t98, -0.2e1 * t28 * t147 + 0.2e1 * t97 * t99, -0.2e1 * t97 * t28, 0.2e1 * t103 * t28 + 0.2e1 * t17 * t147 + 0.2e1 * t31 * t99 - 0.2e1 * t5 * t97, 0.2e1 * t17 * t146 - 0.2e1 * t167 * t28 - 0.2e1 * t98 * t31 - 0.2e1 * t4 * t97, 0.2e1 * t6 * t147 + 0.2e1 * t15 * t99 - 0.2e1 * t28 * t8 + 0.2e1 * t3 * t97, 0.2e1 * t106 * t29 + 0.2e1 * (qJD(5) * t107 - t2 * t78 + t3 * t81) * t53, -0.2e1 * t6 * t146 + 0.2e1 * t15 * t98 - 0.2e1 * t2 * t97 + 0.2e1 * t28 * t7, 0.2e1 * t15 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t3 * t8; 0, 0, 0, 0, 0, t70, -t132, 0, -t124, t125, -t124, t87, -t125, t87 * pkin(7), 0, 0, -t29, t28, 0, t17, -t16, t13, t11, -t20, t18, 0, -t161 * t81 + t84 * t78, t161 * t78 + t84 * t81, -t78 * t85 + t81 * t89, -t1, t78 * t89 + t81 * t85, t1 * t56 + t107 * t36 + t15 * t27 + t34 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, qJ(3) * t160, 0, 0, 0, 0, 0, 0.2e1 * t37, -0.2e1 * t36, t62, t51, 0, 0, 0, -0.2e1 * t55 * t135 + 0.2e1 * t151, -0.2e1 * t55 * t71 - 0.2e1 * t152, -0.2e1 * t34 * t135 + 0.2e1 * t27 * t81, 0.2e1 * t23, 0.2e1 * t27 * t78 + 0.2e1 * t34 * t71, -0.2e1 * t23 * t56 + 0.2e1 * t27 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, t124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t94 * t135 - t168, t83, 0, -t78 * t90 + t168, -t6 * t153 + (t118 * t7 + t120 * t8) * qJD(4) + (qJD(4) * t15 + t1) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t115, 0, 0, 0, 0, 0, t45, -t44, t45, -t41, t44, -t27 * t153 - t79 * t23 + (t34 * t79 + t56 * t92) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t153 + t92) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, -t17, t16, -t13, -t11, t20, -t18, 0, -t17 * t81 + t109 * t78 + (-t108 * t81 + t31 * t78) * qJD(5), t17 * t78 + t109 * t81 + (t108 * t78 + t31 * t81) * qJD(5), -t78 * t88 + t81 * t91, t1, t78 * t91 + t81 * t88, pkin(9) * t1 + t15 * t39 + t57 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, -0.2e1 * t121, t169, 0, 0, 0, t114 * t78 - t151, t114 * t81 + t152, t141 * t135 + t143 * t81, -t23, -t141 * t71 + t143 * t78, -pkin(9) * t23 + t27 * t57 + t34 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t115, 0, 0, 0, 0, 0, -t45, t44, -t45, t41, -t44, -t153 * t39 + (pkin(9) * t92 + t57 * t79) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t51, 0, 0, 0, t78 * t129, t81 * t129, 0.2e1 * t57 * t135 - 0.2e1 * t39 * t81, 0, -0.2e1 * t39 * t78 - 0.2e1 * t57 * t71, 0.2e1 * t57 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t99, t28, t5, t4, t5 + 0.2e1 * t156, -t105 * t29 + t150, -0.2e1 * t134 - t4 + 0.2e1 * t138, -pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t135, 0, -t26, t25, -t26, t38, -t25, t104 * t36 - t38 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t42, -t43, 0, -t42 (-pkin(5) * t120 + qJ(6) * t118) * qJD(4) - t38 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t135, 0, -t126, t127, -t126, -t38, -t127, -t38 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, qJ(6) * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t98, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
