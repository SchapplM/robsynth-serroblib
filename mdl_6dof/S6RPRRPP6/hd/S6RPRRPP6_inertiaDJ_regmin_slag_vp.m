% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:33
% EndTime: 2019-03-09 04:48:36
% DurationCPUTime: 1.38s
% Computational Cost: add. (2020->216), mult. (4522->377), div. (0->0), fcn. (3735->6), ass. (0->120)
t94 = sin(qJ(3));
t87 = t94 * qJD(3);
t93 = sin(qJ(4));
t125 = t93 * t87;
t96 = cos(qJ(3));
t138 = qJD(4) * t96;
t95 = cos(qJ(4));
t127 = t95 * t138;
t58 = t125 - t127;
t142 = cos(pkin(9));
t116 = t142 * t95;
t141 = qJD(4) * t93;
t92 = sin(pkin(9));
t54 = qJD(4) * t116 - t92 * t141;
t158 = t96 * t54;
t154 = t96 * pkin(8);
t108 = t94 * pkin(3) - t154;
t67 = qJ(2) + t108;
t60 = t93 * t67;
t150 = t94 * t95;
t97 = -pkin(1) - pkin(7);
t76 = t97 * t150;
t157 = -t76 - t60;
t89 = t94 ^ 2;
t91 = t96 ^ 2;
t113 = (t89 - t91) * qJD(3);
t90 = t95 ^ 2;
t146 = t93 ^ 2 - t90;
t114 = t146 * qJD(4);
t156 = 2 * qJD(2);
t155 = 2 * qJD(6);
t134 = t96 * qJD(3);
t120 = t97 * t134;
t139 = qJD(4) * t95;
t109 = pkin(3) * t96 + pkin(8) * t94;
t62 = t109 * qJD(3) + qJD(2);
t130 = t95 * t120 + t67 * t139 + t93 * t62;
t137 = qJD(4) * t97;
t11 = -qJ(5) * t127 + (-qJD(5) * t96 + (qJ(5) * qJD(3) - t137) * t94) * t93 + t130;
t151 = t93 * t97;
t119 = pkin(4) - t151;
t123 = t95 * t87;
t135 = t95 * qJD(5);
t99 = t157 * qJD(4) + t95 * t62;
t8 = qJ(5) * t123 + (qJ(5) * t141 + t119 * qJD(3) - t135) * t96 + t99;
t4 = t142 * t11 + t92 * t8;
t153 = t92 * t93;
t152 = t93 * t96;
t117 = t142 * t93;
t64 = t92 * t95 + t117;
t53 = t64 * qJD(4);
t149 = t96 * t53;
t148 = t96 * t97;
t147 = -qJ(5) - pkin(8);
t143 = qJ(5) * t96;
t61 = t95 * t67;
t39 = t119 * t94 - t95 * t143 + t61;
t41 = -t93 * t143 - t157;
t14 = t142 * t41 + t92 * t39;
t144 = t89 + t91;
t140 = qJD(4) * t94;
t136 = t94 * qJD(6);
t133 = qJ(2) * qJD(3);
t132 = qJ(6) * t134 + t4;
t131 = -0.2e1 * pkin(3) * qJD(4);
t86 = pkin(4) * t141;
t129 = t93 * t138;
t128 = t93 * t137;
t124 = t93 * t139;
t122 = t95 * t134;
t121 = t94 * t134;
t85 = -t95 * pkin(4) - pkin(3);
t3 = -t92 * t11 + t142 * t8;
t115 = qJD(4) * t147;
t101 = -t93 * qJD(5) + t95 * t115;
t50 = t93 * t115 + t135;
t34 = -t142 * t101 + t92 * t50;
t35 = t92 * t101 + t142 * t50;
t74 = t147 * t95;
t42 = -t147 * t117 - t92 * t74;
t43 = -t142 * t74 + t147 * t153;
t118 = t42 * t34 + t43 * t35;
t65 = pkin(4) * t152 - t148;
t112 = qJD(3) * t142;
t111 = t93 * t120;
t110 = t93 * t123;
t79 = t97 * t87;
t44 = -t58 * pkin(4) + t79;
t13 = t142 * t39 - t92 * t41;
t106 = t116 - t153;
t30 = -t112 * t152 - t92 * t122 - t54 * t94;
t32 = t106 * t134 - t64 * t140;
t46 = t64 * t94;
t48 = t106 * t94;
t105 = -t30 * t42 + t32 * t43 + t46 * t34 + t48 * t35;
t31 = qJD(3) * t46 - t158;
t33 = t112 * t150 - t92 * t125 + t149;
t47 = t64 * t96;
t49 = t96 * t116 - t92 * t152;
t104 = t43 * t31 - t42 * t33 + t34 * t49 - t35 * t47;
t103 = -t30 * t49 + t48 * t31 - t32 * t47 - t46 * t33;
t102 = t106 * t32 - t30 * t64 + t46 * t54 - t48 * t53;
t100 = 0.2e1 * t106 * t35 + 0.2e1 * t34 * t64 + 0.2e1 * t42 * t54 - 0.2e1 * t43 * t53;
t98 = -0.2e1 * t46 * t30 + 0.2e1 * t48 * t32 - 0.2e1 * t121;
t83 = -t142 * pkin(4) - pkin(5);
t81 = t92 * pkin(4) + qJ(6);
t57 = t93 * t134 + t94 * t139;
t56 = -t123 - t129;
t55 = t93 * t140 - t122;
t37 = -pkin(5) * t106 - t64 * qJ(6) + t85;
t28 = t99 - t111;
t27 = t94 * t128 - t130;
t26 = t47 * pkin(5) - t49 * qJ(6) + t65;
t23 = t53 * pkin(5) - t54 * qJ(6) - t64 * qJD(6) + t86;
t12 = -t94 * pkin(5) - t13;
t10 = t94 * qJ(6) + t14;
t5 = -t31 * pkin(5) + t33 * qJ(6) - t49 * qJD(6) + t44;
t2 = -pkin(5) * t134 - t3;
t1 = t132 + t136;
t6 = [0, 0, 0, 0, t156, qJ(2) * t156, -0.2e1 * t121, 0.2e1 * t113, 0, 0, 0, 0.2e1 * qJD(2) * t94 + 0.2e1 * t96 * t133, 0.2e1 * qJD(2) * t96 - 0.2e1 * t94 * t133, -0.2e1 * t90 * t121 - 0.2e1 * t91 * t124, 0.4e1 * t110 * t96 + 0.2e1 * t91 * t114, -0.2e1 * t95 * t113 - 0.2e1 * t129 * t94, 0.2e1 * t113 * t93 - 0.2e1 * t127 * t94, 0.2e1 * t121, -0.2e1 * t91 * t95 * t137 + 0.2e1 * t61 * t134 + 0.2e1 * (t28 + t111) * t94, 0.2e1 * t91 * t128 + 0.2e1 * t27 * t94 + 0.2e1 * (-t60 + t76) * t134, 0.2e1 * t13 * t33 + 0.2e1 * t14 * t31 - 0.2e1 * t3 * t49 - 0.2e1 * t4 * t47, 0.2e1 * t13 * t3 + 0.2e1 * t14 * t4 + 0.2e1 * t65 * t44, -0.2e1 * t12 * t134 - 0.2e1 * t2 * t94 - 0.2e1 * t26 * t31 + 0.2e1 * t5 * t47, -0.2e1 * t1 * t47 + 0.2e1 * t10 * t31 - 0.2e1 * t12 * t33 + 0.2e1 * t2 * t49, 0.2e1 * t1 * t94 + 0.2e1 * t10 * t134 + 0.2e1 * t26 * t33 - 0.2e1 * t5 * t49, 0.2e1 * t10 * t1 + 0.2e1 * t12 * t2 + 0.2e1 * t26 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144 * t139, t144 * t141, t103, t13 * t30 + t14 * t32 - t3 * t46 + t4 * t48 - t44 * t96 + t65 * t87, t30 * t94 + t96 * t31 + (-t46 * t96 + t47 * t94) * qJD(3), t103, t32 * t94 - t96 * t33 + (t48 * t96 - t49 * t94) * qJD(3), t1 * t48 + t10 * t32 - t12 * t30 + t2 * t46 + t26 * t87 - t5 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t134, 0, -t79, -t120, -t96 * t114 - t110, -0.4e1 * t124 * t96 + t146 * t87, t57, -t55, 0 (-t109 * t95 - t148 * t93) * qJD(4) + (t108 * t93 - t76) * qJD(3) (t109 * t93 - t148 * t95) * qJD(4) + (-t95 * t154 + (pkin(3) * t95 + t151) * t94) * qJD(3), t106 * t4 - t13 * t54 - t14 * t53 - t3 * t64 + t104, -t13 * t34 + t14 * t35 - t3 * t42 + t4 * t43 + t44 * t85 + t65 * t86, -t106 * t5 - t134 * t42 + t23 * t47 + t26 * t53 - t37 * t31 - t34 * t94, t1 * t106 - t10 * t53 + t12 * t54 + t2 * t64 + t104, t134 * t43 - t23 * t49 - t26 * t54 + t37 * t33 + t35 * t94 - t5 * t64, t1 * t43 + t10 * t35 + t12 * t34 + t2 * t42 + t26 * t23 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t134, 0, 0, 0, 0, 0, t56, t58, t102, -pkin(4) * t129 + t85 * t87 + t105, -t106 * t87 - t149, t102, -t64 * t87 + t158, -t96 * t23 + t37 * t87 + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t124, -0.2e1 * t114, 0, 0, 0, t93 * t131, t95 * t131, t100, 0.2e1 * t85 * t86 + 0.2e1 * t118, -0.2e1 * t106 * t23 + 0.2e1 * t37 * t53, t100, -0.2e1 * t23 * t64 - 0.2e1 * t37 * t54, 0.2e1 * t37 * t23 + 0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t58, t134, t28, t27 (t142 * t33 + t31 * t92) * pkin(4) (t142 * t3 + t4 * t92) * pkin(4) (pkin(5) - t83) * t134 + t3, -qJD(6) * t47 + t81 * t31 - t83 * t33, t134 * t81 + t132 + 0.2e1 * t136, t10 * qJD(6) + t1 * t81 + t2 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t55, 0 (t142 * t30 + t32 * t92) * pkin(4), t30, 0, t32, t48 * qJD(6) - t30 * t83 + t32 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t141, 0, -pkin(8) * t139, pkin(8) * t141 (-t142 * t54 - t53 * t92) * pkin(4) (-t142 * t34 + t35 * t92) * pkin(4), -t34, qJD(6) * t106 - t81 * t53 + t83 * t54, t35, t43 * qJD(6) + t34 * t83 + t35 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t81 * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t31, 0, t33, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t53, 0, -t54, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t33, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
