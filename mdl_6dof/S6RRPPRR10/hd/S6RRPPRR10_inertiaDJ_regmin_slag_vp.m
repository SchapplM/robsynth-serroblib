% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:01
% EndTime: 2019-03-09 09:37:06
% DurationCPUTime: 1.60s
% Computational Cost: add. (2007->192), mult. (4606->358), div. (0->0), fcn. (4234->8), ass. (0->117)
t106 = sin(qJ(2));
t104 = sin(qJ(6));
t107 = cos(qJ(6));
t101 = sin(pkin(10));
t102 = cos(pkin(10));
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t114 = t101 * t105 - t102 * t108;
t74 = t101 * t108 + t102 * t105;
t47 = -t104 * t74 - t107 * t114;
t66 = t74 * qJD(5);
t140 = qJD(5) * t108;
t141 = qJD(5) * t105;
t67 = -t101 * t141 + t102 * t140;
t15 = -qJD(6) * t47 + t104 * t66 - t107 * t67;
t46 = -t104 * t114 + t107 * t74;
t109 = cos(qJ(2));
t95 = qJD(2) * t109;
t164 = t106 * t15 - t46 * t95;
t160 = pkin(3) + pkin(7);
t82 = (t101 ^ 2 + t102 ^ 2) * qJD(4);
t103 = -pkin(2) - qJ(4);
t157 = -pkin(8) + t103;
t80 = t157 * t101;
t81 = t157 * t102;
t27 = qJD(4) * t74 - t140 * t81 + t141 * t80;
t145 = t106 * qJ(3);
t115 = -t103 * t109 + t145;
t110 = -qJD(6) * t46 - t104 * t67 - t107 * t66;
t163 = t106 * t110 + t47 * t95;
t138 = t109 * qJD(3);
t162 = qJD(2) * t115 + qJD(4) * t106 - t138;
t161 = 0.2e1 * qJD(3);
t159 = pkin(5) * t106;
t142 = qJD(2) * t106;
t128 = t102 * t142;
t129 = t101 * t142;
t139 = qJD(5) * t109;
t42 = -t105 * t129 + t108 * t128 + t139 * t74;
t148 = t101 * t106;
t124 = pkin(2) * t142 - t106 * qJD(3);
t51 = -t109 * qJD(4) + (-qJ(3) * t109 + qJ(4) * t106) * qJD(2) + t124;
t93 = pkin(7) * t95;
t79 = pkin(3) * t95 + t93;
t30 = -t101 * t51 + t102 * t79;
t24 = (pkin(4) * t109 - pkin(8) * t148) * qJD(2) + t30;
t31 = t101 * t79 + t102 * t51;
t26 = pkin(8) * t128 + t31;
t73 = -pkin(1) - t115;
t85 = t160 * t106;
t77 = t102 * t85;
t39 = t106 * pkin(4) + t77 + (pkin(8) * t109 - t73) * t101;
t147 = t102 * t109;
t49 = t101 * t85 + t102 * t73;
t43 = -pkin(8) * t147 + t49;
t8 = -t105 * t24 - t108 * t26 - t140 * t39 + t141 * t43;
t5 = pkin(9) * t42 - t8;
t158 = t107 * t5;
t156 = t105 * t39 + t108 * t43;
t155 = t105 * t81 + t108 * t80;
t86 = t160 * t109;
t154 = pkin(5) * qJD(6);
t58 = t114 * t109;
t13 = pkin(9) * t58 + t156;
t153 = t104 * t13;
t151 = t107 * t13;
t78 = t160 * t142;
t150 = t78 * t101;
t92 = pkin(4) * t101 + qJ(3);
t143 = qJ(3) * qJD(3);
t137 = -0.2e1 * pkin(1) * qJD(2);
t65 = pkin(4) * t147 + t86;
t136 = pkin(5) * t95;
t135 = t104 * t154;
t134 = t107 * t154;
t133 = pkin(7) * t142;
t41 = t114 * t139 + t142 * t74;
t9 = -qJD(5) * t156 - t105 * t26 + t108 * t24;
t4 = -t41 * pkin(9) + t136 + t9;
t130 = -t104 * t5 + t107 * t4;
t126 = -t105 * t43 + t108 * t39;
t59 = t74 * t109;
t12 = pkin(9) * t59 + t126 + t159;
t127 = -t12 - t159;
t125 = -t105 * t80 + t108 * t81;
t122 = -pkin(2) * t109 - t145;
t14 = t101 * t31 + t102 * t30;
t121 = -t107 * t12 + t153;
t120 = t104 * t12 + t151;
t34 = pkin(9) * t114 + t125;
t35 = -pkin(9) * t74 + t155;
t119 = t104 * t35 - t107 * t34;
t118 = t104 * t34 + t107 * t35;
t37 = -t104 * t59 - t107 * t58;
t38 = t104 * t58 - t107 * t59;
t112 = -t106 * t67 - t74 * t95;
t54 = (-pkin(4) * t102 - t160) * t142;
t111 = qJD(2) * t122 + t138;
t28 = qJD(4) * t114 - t155 * qJD(5);
t88 = 0.2e1 * t106 * t95;
t83 = -pkin(1) + t122;
t61 = -qJ(3) * t95 + t124;
t55 = pkin(5) * t67 + qJD(3);
t53 = pkin(5) * t74 + t92;
t48 = -t101 * t73 + t77;
t45 = -pkin(5) * t58 + t65;
t44 = -t106 * t66 - t114 * t95;
t23 = -t42 * pkin(5) + t54;
t20 = t66 * pkin(9) + t28;
t19 = -t67 * pkin(9) - t27;
t11 = qJD(6) * t38 + t104 * t41 - t107 * t42;
t10 = -qJD(6) * t37 + t104 * t42 + t107 * t41;
t7 = -qJD(6) * t118 - t104 * t19 + t107 * t20;
t6 = qJD(6) * t119 - t104 * t20 - t107 * t19;
t2 = -qJD(6) * t120 + t130;
t1 = qJD(6) * t121 - t104 * t4 - t158;
t3 = [0, 0, 0, t88, 0.2e1 * (-t106 ^ 2 + t109 ^ 2) * qJD(2), 0, 0, 0, t106 * t137, t109 * t137, 0, 0.2e1 * t109 * t61 - 0.2e1 * t142 * t83, -0.2e1 * t106 * t61 - 0.2e1 * t83 * t95, 0.2e1 * t83 * t61, -0.2e1 * t78 * t147 + 0.2e1 * t30 * t106 + 0.2e1 * (-t102 * t106 * t86 + t109 * t48) * qJD(2), 0.2e1 * t109 * t150 - 0.2e1 * t31 * t106 + 0.2e1 * (-t109 * t49 + t148 * t86) * qJD(2), 0.2e1 * (t101 * t30 - t102 * t31) * t109 + 0.2e1 * (-t101 * t48 + t102 * t49) * t142, 0.2e1 * t30 * t48 + 0.2e1 * t31 * t49 - 0.2e1 * t78 * t86, -0.2e1 * t59 * t41, 0.2e1 * t41 * t58 - 0.2e1 * t42 * t59, 0.2e1 * t106 * t41 - 0.2e1 * t59 * t95, 0.2e1 * t106 * t42 + 0.2e1 * t58 * t95, t88, 0.2e1 * t106 * t9 + 0.2e1 * t126 * t95 - 0.2e1 * t42 * t65 - 0.2e1 * t54 * t58, 0.2e1 * t106 * t8 - 0.2e1 * t156 * t95 + 0.2e1 * t41 * t65 - 0.2e1 * t54 * t59, 0.2e1 * t38 * t10, -0.2e1 * t10 * t37 - 0.2e1 * t11 * t38, 0.2e1 * t10 * t106 + 0.2e1 * t38 * t95, -0.2e1 * t106 * t11 - 0.2e1 * t37 * t95, t88, 0.2e1 * t106 * t2 + 0.2e1 * t11 * t45 - 0.2e1 * t121 * t95 + 0.2e1 * t23 * t37, 0.2e1 * t1 * t106 + 0.2e1 * t10 * t45 - 0.2e1 * t120 * t95 + 0.2e1 * t23 * t38; 0, 0, 0, 0, 0, t95, -t142, 0, -t93, t133, t111, t93, -t133, t111 * pkin(7), -t102 * t162 - t150, t101 * t162 - t78 * t102, -t14, -t78 * qJ(3) + t86 * qJD(3) + t14 * t103 + (-t101 * t49 - t102 * t48) * qJD(4), -t114 * t41 + t59 * t66, -t114 * t42 - t41 * t74 - t58 * t66 + t59 * t67, t44, t112, 0, -qJD(3) * t58 + t106 * t28 + t125 * t95 - t42 * t92 + t54 * t74 + t65 * t67, -qJD(3) * t59 + t106 * t27 - t114 * t54 - t155 * t95 + t41 * t92 - t65 * t66, t10 * t47 + t110 * t38, -t10 * t46 - t11 * t47 - t110 * t37 + t15 * t38, t163, t164, 0, t106 * t7 + t11 * t53 - t119 * t95 - t15 * t45 + t23 * t46 + t37 * t55, t10 * t53 + t106 * t6 + t110 * t45 - t118 * t95 + t23 * t47 + t38 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, 0.2e1 * t143, t101 * t161, t102 * t161, 0.2e1 * t82, -0.2e1 * t103 * t82 + 0.2e1 * t143, 0.2e1 * t114 * t66, 0.2e1 * t114 * t67 + 0.2e1 * t66 * t74, 0, 0, 0, 0.2e1 * qJD(3) * t74 + 0.2e1 * t67 * t92, -0.2e1 * qJD(3) * t114 - 0.2e1 * t66 * t92, 0.2e1 * t47 * t110, -0.2e1 * t110 * t46 + 0.2e1 * t15 * t47, 0, 0, 0, -0.2e1 * t15 * t53 + 0.2e1 * t46 * t55, 0.2e1 * t110 * t53 + 0.2e1 * t47 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, t93, t102 * t95, -t101 * t95, 0, t14, 0, 0, 0, 0, 0, t44, t112, 0, 0, 0, 0, 0, t163, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t129, 0, -t78, 0, 0, 0, 0, 0, -t42, t41, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, t67, -t66, 0, 0, 0, 0, 0, -t15, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, t95, t9, t8, 0, 0, t10, -t11, t95, t107 * t136 + (t104 * t127 - t151) * qJD(6) + t130, -t158 + (-t4 - t136) * t104 + (t107 * t127 + t153) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t67, 0, t28, t27, 0, 0, t110, t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t67, 0, 0, 0, 0, 0, t110, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t135, -0.2e1 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t95, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
