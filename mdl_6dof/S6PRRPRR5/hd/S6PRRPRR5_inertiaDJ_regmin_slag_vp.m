% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:47
% EndTime: 2019-03-08 22:20:52
% DurationCPUTime: 1.61s
% Computational Cost: add. (2053->231), mult. (5578->432), div. (0->0), fcn. (5527->12), ass. (0->131)
t101 = sin(pkin(12));
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t139 = qJD(5) * t110;
t103 = cos(pkin(12));
t144 = t110 * t103;
t158 = pkin(9) + qJ(4);
t86 = t158 * t101;
t87 = t158 * t103;
t29 = (qJD(4) * t101 + qJD(5) * t87) * t106 - qJD(4) * t144 + t86 * t139;
t95 = -pkin(4) * t103 - pkin(3);
t162 = 0.2e1 * t95;
t82 = t101 * t110 + t103 * t106;
t71 = t82 * qJD(5);
t161 = pkin(5) * t71;
t111 = cos(qJ(3));
t160 = pkin(5) * t111;
t109 = cos(qJ(6));
t141 = qJD(5) * t106;
t107 = sin(qJ(3));
t145 = t103 * t111;
t97 = t107 * qJD(3);
t132 = pkin(8) * t97;
t68 = -t107 * qJD(4) + (pkin(3) * t107 - qJ(4) * t111) * qJD(3);
t52 = t101 * t132 + t103 * t68;
t34 = (pkin(4) * t107 - pkin(9) * t145) * qJD(3) + t52;
t146 = t103 * t107;
t149 = t101 * t111;
t63 = t101 * t68;
t42 = t63 + (-pkin(8) * t146 - pkin(9) * t149) * qJD(3);
t122 = -pkin(3) * t111 - qJ(4) * t107;
t84 = -pkin(2) + t122;
t77 = t103 * t84;
t49 = -pkin(9) * t146 + t77 + (-pkin(8) * t101 - pkin(4)) * t111;
t150 = t101 * t107;
t93 = pkin(8) * t145;
t60 = t101 * t84 + t93;
t58 = -pkin(9) * t150 + t60;
t12 = -t106 * t34 - t110 * t42 - t49 * t139 + t58 * t141;
t138 = t111 * qJD(3);
t140 = qJD(5) * t107;
t151 = t101 * t106;
t40 = t82 * t138 + t139 * t146 - t140 * t151;
t9 = -pkin(10) * t40 - t12;
t159 = t109 * t9;
t157 = t106 * t49 + t110 * t58;
t155 = -t106 * t86 + t110 * t87;
t130 = t101 * t138;
t96 = pkin(8) * t138;
t74 = pkin(4) * t130 + t96;
t83 = pkin(4) * t150 + t107 * pkin(8);
t154 = pkin(5) * qJD(6);
t105 = sin(qJ(6));
t66 = t82 * t107;
t17 = -pkin(10) * t66 + t157;
t153 = t105 * t17;
t152 = t109 * t17;
t102 = sin(pkin(6));
t108 = sin(qJ(2));
t148 = t102 * t108;
t112 = cos(qJ(2));
t147 = t102 * t112;
t143 = qJD(2) * t108;
t142 = qJD(4) * t111;
t137 = -0.2e1 * pkin(2) * qJD(3);
t136 = pkin(8) * t149;
t135 = pkin(5) * t97;
t134 = t105 * t154;
t133 = t109 * t154;
t13 = -t157 * qJD(5) - t106 * t42 + t110 * t34;
t81 = -t144 + t151;
t39 = -t81 * t138 - t82 * t140;
t8 = -pkin(10) * t39 + t13 + t135;
t131 = -t105 * t9 + t109 * t8;
t129 = t102 * t143;
t128 = qJD(2) * t147;
t127 = t107 * t138;
t125 = -t106 * t58 + t110 * t49;
t67 = t81 * t107;
t16 = pkin(10) * t67 + t125 - t160;
t126 = -t16 + t160;
t124 = -t106 * t87 - t110 * t86;
t123 = 0.2e1 * (t101 ^ 2 + t103 ^ 2) * qJD(4);
t104 = cos(pkin(6));
t72 = -t104 * t111 + t107 * t148;
t56 = -qJD(3) * t72 + t111 * t128;
t35 = -t101 * t56 + t103 * t129;
t36 = t101 * t129 + t103 * t56;
t121 = -t101 * t35 + t103 * t36;
t53 = -t103 * t132 + t63;
t120 = -t101 * t52 + t103 * t53;
t119 = -t109 * t16 + t153;
t118 = t105 * t16 + t152;
t73 = t104 * t107 + t111 * t148;
t54 = -t73 * t101 - t103 * t147;
t55 = -t101 * t147 + t73 * t103;
t20 = -t106 * t55 + t110 * t54;
t21 = t106 * t54 + t110 * t55;
t117 = t105 * t21 - t109 * t20;
t116 = t105 * t20 + t109 * t21;
t37 = -pkin(10) * t82 + t124;
t38 = -pkin(10) * t81 + t155;
t115 = t105 * t38 - t109 * t37;
t114 = t105 * t37 + t109 * t38;
t32 = -t105 * t67 + t109 * t66;
t33 = -t105 * t66 - t109 * t67;
t47 = t105 * t82 + t109 * t81;
t48 = -t105 * t81 + t109 * t82;
t30 = -t82 * qJD(4) - t155 * qJD(5);
t90 = -0.2e1 * t127;
t70 = t81 * qJD(5);
t62 = pkin(5) * t81 + t95;
t59 = t77 - t136;
t57 = qJD(3) * t73 + t107 * t128;
t51 = pkin(5) * t66 + t83;
t24 = pkin(5) * t40 + t74;
t23 = pkin(10) * t70 + t30;
t22 = -pkin(10) * t71 - t29;
t19 = qJD(6) * t48 - t105 * t70 + t109 * t71;
t18 = -qJD(6) * t47 - t105 * t71 - t109 * t70;
t15 = qJD(6) * t33 + t105 * t39 + t109 * t40;
t14 = -qJD(6) * t32 - t105 * t40 + t109 * t39;
t11 = -qJD(5) * t21 - t106 * t36 + t110 * t35;
t10 = -t106 * t35 - t110 * t36 - t54 * t139 + t55 * t141;
t6 = -qJD(6) * t114 - t105 * t22 + t109 * t23;
t5 = t115 * qJD(6) - t105 * t23 - t109 * t22;
t4 = -t116 * qJD(6) + t105 * t10 + t109 * t11;
t3 = t117 * qJD(6) + t109 * t10 - t105 * t11;
t2 = -t118 * qJD(6) + t131;
t1 = t119 * qJD(6) - t105 * t8 - t159;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t35 * t54 + 0.2e1 * t36 * t55 + 0.2e1 * t57 * t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t129, -t128, 0, 0, 0, 0, 0 (-t111 * t143 - t112 * t97) * t102 (t107 * t143 - t112 * t138) * t102, t57 * t150 - t35 * t111 + (t107 * t54 + t72 * t149) * qJD(3), t57 * t146 + t36 * t111 + (-t107 * t55 + t72 * t145) * qJD(3) (-t101 * t36 - t103 * t35) * t107 + (-t101 * t55 - t103 * t54) * t138, t35 * t59 + t36 * t60 + t52 * t54 + t53 * t55 + (t107 * t57 + t72 * t138) * pkin(8), 0, 0, 0, 0, 0, -t11 * t111 + t20 * t97 + t40 * t72 + t57 * t66, -t10 * t111 - t21 * t97 + t39 * t72 - t57 * t67, 0, 0, 0, 0, 0, -t4 * t111 - t117 * t97 + t72 * t15 + t57 * t32, -t3 * t111 - t116 * t97 + t72 * t14 + t57 * t33; 0, 0, 0, 0, 0.2e1 * t127, 0.2e1 * (-t107 ^ 2 + t111 ^ 2) * qJD(3), 0, 0, 0, t107 * t137, t111 * t137, -0.2e1 * t52 * t111 + 0.2e1 * (t59 + 0.2e1 * t136) * t97, 0.2e1 * t53 * t111 + 0.2e1 * (-t60 + 0.2e1 * t93) * t97, 0.2e1 * (-t101 * t53 - t103 * t52) * t107 + 0.2e1 * (-t101 * t60 - t103 * t59) * t138, 0.2e1 * pkin(8) ^ 2 * t127 + 0.2e1 * t52 * t59 + 0.2e1 * t53 * t60, -0.2e1 * t67 * t39, -0.2e1 * t39 * t66 + 0.2e1 * t40 * t67, -0.2e1 * t111 * t39 - 0.2e1 * t67 * t97, 0.2e1 * t111 * t40 - 0.2e1 * t66 * t97, t90, -0.2e1 * t13 * t111 + 0.2e1 * t125 * t97 + 0.2e1 * t83 * t40 + 0.2e1 * t74 * t66, -0.2e1 * t12 * t111 - 0.2e1 * t157 * t97 + 0.2e1 * t83 * t39 - 0.2e1 * t74 * t67, 0.2e1 * t33 * t14, -0.2e1 * t14 * t32 - 0.2e1 * t15 * t33, -0.2e1 * t111 * t14 + 0.2e1 * t33 * t97, 0.2e1 * t111 * t15 - 0.2e1 * t32 * t97, t90, -0.2e1 * t2 * t111 - 0.2e1 * t119 * t97 + 0.2e1 * t51 * t15 + 0.2e1 * t24 * t32, -0.2e1 * t1 * t111 - 0.2e1 * t118 * t97 + 0.2e1 * t51 * t14 + 0.2e1 * t24 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, -t57 * t103, t57 * t101, t121, -pkin(3) * t57 + (-t101 * t54 + t103 * t55) * qJD(4) + t121 * qJ(4), 0, 0, 0, 0, 0, t57 * t81 + t71 * t72, t57 * t82 - t70 * t72, 0, 0, 0, 0, 0, t19 * t72 + t47 * t57, t18 * t72 + t48 * t57; 0, 0, 0, 0, 0, 0, t138, -t97, 0, -t96, t132, t101 * t142 + (t122 * t101 - t93) * qJD(3), t103 * t142 + (t122 * t103 + t136) * qJD(3), t120, -pkin(3) * t96 + (-t101 * t59 + t103 * t60) * qJD(4) + t120 * qJ(4), t39 * t82 + t67 * t70, -t39 * t81 - t40 * t82 + t66 * t70 + t67 * t71, t111 * t70 + t82 * t97, t111 * t71 - t81 * t97, 0, -t30 * t111 + t124 * t97 + t95 * t40 + t83 * t71 + t74 * t81, -t29 * t111 - t155 * t97 + t95 * t39 - t83 * t70 + t74 * t82, t14 * t48 + t18 * t33, -t14 * t47 - t15 * t48 - t18 * t32 - t19 * t33, -t111 * t18 + t48 * t97, t111 * t19 - t47 * t97, 0, -t6 * t111 - t115 * t97 + t62 * t15 + t32 * t161 + t51 * t19 + t24 * t47, -t5 * t111 - t114 * t97 + t62 * t14 + t33 * t161 + t51 * t18 + t24 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, qJ(4) * t123, -0.2e1 * t82 * t70, 0.2e1 * t70 * t81 - 0.2e1 * t71 * t82, 0, 0, 0, t71 * t162, -t70 * t162, 0.2e1 * t48 * t18, -0.2e1 * t18 * t47 - 0.2e1 * t19 * t48, 0, 0, 0, 0.2e1 * t47 * t161 + 0.2e1 * t19 * t62, 0.2e1 * t48 * t161 + 0.2e1 * t18 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t103 * t138, 0, t96, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t70, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, t97, t13, t12, 0, 0, t14, -t15, t97, t109 * t135 + (t105 * t126 - t152) * qJD(6) + t131, -t159 + (-t8 - t135) * t105 + (t109 * t126 + t153) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t71, 0, t30, t29, 0, 0, t18, -t19, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t134, -0.2e1 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t97, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
