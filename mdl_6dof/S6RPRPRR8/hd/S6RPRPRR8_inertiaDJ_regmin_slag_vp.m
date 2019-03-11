% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:31
% EndTime: 2019-03-09 03:59:35
% DurationCPUTime: 1.44s
% Computational Cost: add. (1715->177), mult. (3708->307), div. (0->0), fcn. (3615->8), ass. (0->117)
t126 = sin(pkin(10));
t100 = qJD(3) * t126;
t127 = cos(pkin(10));
t101 = qJD(3) * t127;
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t46 = t77 * t100 - t80 * t101;
t75 = sin(qJ(6));
t76 = sin(qJ(5));
t135 = t75 * t76;
t78 = cos(qJ(6));
t79 = cos(qJ(5));
t55 = -t78 * t79 + t135;
t137 = t55 * t46;
t118 = qJD(5) + qJD(6);
t56 = t75 * t79 + t78 * t76;
t150 = t118 * t56;
t103 = t126 * t80;
t53 = -t127 * t77 - t103;
t17 = t150 * t53 + t137;
t47 = -t80 * t100 - t77 * t101;
t104 = t127 * t80;
t52 = -t126 * t77 + t104;
t8 = -t150 * t52 - t47 * t55;
t140 = t52 * t47;
t142 = t46 * t53;
t86 = 0.2e1 * t140 + 0.2e1 * t142;
t149 = (-t126 * t46 + t127 * t47) * pkin(3);
t129 = t77 * pkin(3) + qJ(2);
t32 = -t53 * pkin(4) - t52 * pkin(8) + t129;
t81 = -pkin(1) - pkin(7);
t128 = qJ(4) - t81;
t57 = t128 * t77;
t38 = -t128 * t103 - t127 * t57;
t33 = t79 * t38;
t131 = t76 * t32 + t33;
t125 = qJD(5) * t76;
t132 = t79 * t47;
t88 = -t52 * t125 + t132;
t51 = t52 ^ 2;
t148 = 2 * qJD(2);
t147 = 0.2e1 * qJD(5);
t146 = t46 * pkin(5);
t145 = t53 * pkin(5);
t124 = qJD(5) * t79;
t120 = t80 * qJD(3);
t43 = -t77 * qJD(4) - t128 * t120;
t121 = t77 * qJD(3);
t85 = -t80 * qJD(4) + t128 * t121;
t23 = t126 * t85 + t127 * t43;
t59 = pkin(3) * t120 + qJD(2);
t24 = -t46 * pkin(4) - t47 * pkin(8) + t59;
t6 = -t32 * t124 + t38 * t125 - t79 * t23 - t76 * t24;
t134 = t76 * t47;
t89 = t52 * t124 + t134;
t5 = -t89 * pkin(9) - t6;
t144 = t78 * t5;
t68 = t126 * pkin(3) + pkin(8);
t143 = pkin(9) + t68;
t139 = t52 * t76;
t138 = t52 * t79;
t87 = t56 * t46;
t13 = -pkin(9) * t139 + t131;
t136 = t75 * t13;
t133 = t78 * t13;
t74 = t79 ^ 2;
t130 = t76 ^ 2 - t74;
t123 = qJD(6) * t75;
t122 = qJD(6) * t78;
t119 = qJ(2) * qJD(3);
t69 = -t127 * pkin(3) - pkin(4);
t117 = t69 * t147;
t116 = pkin(5) * t125;
t115 = pkin(5) * t123;
t114 = pkin(5) * t122;
t112 = t76 * t124;
t111 = t53 ^ 2 + t51;
t107 = -t76 * t23 + t79 * t24;
t4 = -pkin(9) * t132 - t146 + (-t33 + (pkin(9) * t52 - t32) * t76) * qJD(5) + t107;
t110 = t78 * t4 - t75 * t5;
t106 = t79 * t32 - t76 * t38;
t12 = -pkin(9) * t138 + t106 - t145;
t109 = -t12 + t145;
t108 = -0.4e1 * t76 * t138;
t105 = qJD(5) * t143;
t102 = t130 * qJD(5);
t99 = t78 * t12 - t136;
t98 = t75 * t12 + t133;
t35 = t118 * t135 - t79 * t122 - t78 * t124;
t97 = t35 * t53 - t87;
t95 = t46 * t68 + t47 * t69;
t49 = t143 * t76;
t50 = t143 * t79;
t94 = -t78 * t49 - t75 * t50;
t93 = -t75 * t49 + t78 * t50;
t92 = t52 * t69 + t53 * t68;
t22 = t126 * t43 - t127 * t85;
t37 = t128 * t104 - t126 * t57;
t91 = t53 * t124 + t76 * t46;
t90 = t53 * t125 - t79 * t46;
t83 = t22 * t52 + t23 * t53 + t37 * t47 + t38 * t46;
t58 = -t79 * pkin(5) + t69;
t45 = t79 * t105;
t44 = t76 * t105;
t34 = 0.2e1 * t142;
t27 = t55 * t52;
t26 = t56 * t52;
t19 = pkin(5) * t139 + t37;
t16 = -t93 * qJD(6) + t75 * t44 - t78 * t45;
t15 = -t94 * qJD(6) + t78 * t44 + t75 * t45;
t14 = t89 * pkin(5) + t22;
t11 = -t118 * t53 * t55 + t87;
t10 = -t123 * t139 + (t118 * t138 + t134) * t78 + t88 * t75;
t7 = -t131 * qJD(5) + t107;
t2 = -t98 * qJD(6) + t110;
t1 = -t99 * qJD(6) - t75 * t4 - t144;
t3 = [0, 0, 0, 0, t148, qJ(2) * t148, -0.2e1 * t77 * t120, 0.2e1 * (t77 ^ 2 - t80 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t77 + 0.2e1 * t80 * t119, 0.2e1 * qJD(2) * t80 - 0.2e1 * t77 * t119, 0.2e1 * t83, 0.2e1 * t129 * t59 + 0.2e1 * t37 * t22 + 0.2e1 * t38 * t23, -0.2e1 * t51 * t112 + 0.2e1 * t74 * t140, t130 * t51 * t147 + t47 * t108, -0.2e1 * t53 * t132 + 0.2e1 * t90 * t52, 0.2e1 * t53 * t134 + 0.2e1 * t91 * t52, t34, -0.2e1 * t106 * t46 + 0.2e1 * t22 * t139 + 0.2e1 * t89 * t37 - 0.2e1 * t7 * t53, 0.2e1 * t131 * t46 + 0.2e1 * t22 * t138 + 0.2e1 * t88 * t37 - 0.2e1 * t6 * t53, -0.2e1 * t27 * t8, 0.2e1 * t27 * t10 - 0.2e1 * t8 * t26, 0.2e1 * t27 * t46 - 0.2e1 * t8 * t53, 0.2e1 * t10 * t53 + 0.2e1 * t26 * t46, t34, 0.2e1 * t19 * t10 + 0.2e1 * t14 * t26 - 0.2e1 * t2 * t53 - 0.2e1 * t99 * t46, -0.2e1 * t1 * t53 - 0.2e1 * t14 * t27 + 0.2e1 * t19 * t8 + 0.2e1 * t46 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t83, 0, 0, 0, 0, 0, -t111 * t124 - t76 * t86, t111 * t125 - t79 * t86, 0, 0, 0, 0, 0, -t52 * t10 - t47 * t26 + (-t11 - t87) * t53, t47 * t27 - t52 * t8 + (t17 + t137) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t120, 0, -t81 * t121, -t81 * t120, -t149 (t126 * t23 - t127 * t22) * pkin(3), -t52 * t102 + t76 * t132, qJD(5) * t108 - t130 * t47, -t91, t90, 0, -t22 * t79 + t95 * t76 + (t37 * t76 + t92 * t79) * qJD(5), t22 * t76 + t95 * t79 + (t37 * t79 - t92 * t76) * qJD(5), t27 * t35 + t8 * t56, -t56 * t10 + t150 * t27 + t35 * t26 - t8 * t55, t97, t17, 0, t58 * t10 + t26 * t116 + t14 * t55 + t150 * t19 - t16 * t53 - t94 * t46, -t116 * t27 + t14 * t56 - t15 * t53 - t19 * t35 + t46 * t93 + t58 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t120, 0, t149, 0, 0, 0, 0, 0, t88, -t89, 0, 0, 0, 0, 0, t8, t52 * t35 - t47 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t112, -0.2e1 * t102, 0, 0, 0, t76 * t117, t79 * t117, -0.2e1 * t56 * t35, -0.2e1 * t150 * t56 + 0.2e1 * t35 * t55, 0, 0, 0, 0.2e1 * t55 * t116 + 0.2e1 * t150 * t58, 0.2e1 * t116 * t56 - 0.2e1 * t58 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, t90, t91, 0, 0, 0, 0, 0, t17, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t89, -t46, t7, t6, 0, 0, t8, -t10, -t46, -t78 * t146 + (t109 * t75 - t133) * qJD(6) + t110, -t144 + (-t4 + t146) * t75 + (t109 * t78 + t136) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t90, 0, 0, 0, 0, 0, t11, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t125, 0, -t68 * t124, t68 * t125, 0, 0, -t35, -t150, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124, 0, 0, 0, 0, 0, -t150, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t115, -0.2e1 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t10, -t46, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t150, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
