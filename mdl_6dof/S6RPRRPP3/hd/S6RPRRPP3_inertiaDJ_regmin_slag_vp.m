% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:37:03
% EndTime: 2019-03-09 04:37:08
% DurationCPUTime: 1.70s
% Computational Cost: add. (1227->230), mult. (2700->358), div. (0->0), fcn. (1974->6), ass. (0->127)
t71 = cos(qJ(3));
t124 = t71 * qJD(3);
t70 = cos(qJ(4));
t60 = qJD(4) * t70;
t68 = sin(qJ(4));
t69 = sin(qJ(3));
t24 = t68 * t124 + t69 * t60;
t152 = -0.4e1 * t69;
t67 = -pkin(4) - qJ(6);
t151 = t67 * t69;
t150 = t24 * pkin(4);
t59 = qJD(5) * t70;
t131 = t68 * qJ(5);
t87 = -t70 * pkin(4) - t131;
t149 = t87 * qJD(4) + t59;
t84 = t67 * t70 - t131;
t62 = t68 ^ 2;
t64 = t70 ^ 2;
t97 = (t62 - t64) * qJD(4);
t38 = -pkin(3) + t87;
t142 = t38 * t69;
t146 = pkin(8) * t71;
t53 = sin(pkin(9)) * pkin(1) + pkin(7);
t37 = t53 * t124;
t121 = qJ(5) * qJD(4);
t100 = t69 * t121;
t101 = qJ(5) * t124;
t116 = -t68 * t100 + t70 * t101 + t69 * t59;
t93 = t116 - t150;
t8 = t37 - t93;
t148 = (t142 + t146) * qJD(4) - t8;
t54 = -cos(pkin(9)) * pkin(1) - pkin(2);
t143 = t69 * pkin(8);
t92 = -t71 * pkin(3) - t143;
t29 = t54 + t92;
t91 = pkin(3) * t69 - t146;
t36 = t91 * qJD(3);
t136 = t29 * t60 + t68 * t36;
t127 = qJD(4) * t71;
t111 = t68 * t127;
t58 = t69 * qJD(3);
t23 = t70 * t58 + t111;
t6 = t23 * t53 - t136;
t147 = pkin(5) + pkin(8);
t126 = qJD(6) * t68;
t79 = -t69 * t126 + t116;
t3 = t24 * qJ(6) + t150 + t37 - t79;
t145 = t3 * t68;
t144 = t3 * t70;
t140 = t68 * t69;
t139 = t68 * t71;
t138 = t69 * t70;
t137 = t70 * t71;
t135 = pkin(4) * t140 + t69 * t53;
t63 = t69 ^ 2;
t133 = -t71 ^ 2 + t63;
t132 = qJ(5) * t70;
t130 = qJD(3) * t68;
t129 = qJD(3) * t70;
t128 = qJD(4) * t68;
t125 = t68 * qJD(5);
t123 = t71 * qJD(5);
t122 = t71 * qJD(6);
t120 = qJ(5) * qJD(5);
t119 = pkin(5) * t137;
t118 = -0.2e1 * pkin(3) * qJD(4);
t32 = t53 * t139;
t33 = t53 * t137;
t117 = 0.2e1 * t58;
t115 = pkin(4) * t58;
t114 = pkin(8) * t128;
t113 = t53 * t128;
t112 = t69 * t128;
t109 = t70 * t127;
t28 = -pkin(3) + t84;
t108 = t28 * t60;
t107 = t62 * t124;
t105 = t68 * t60;
t104 = t69 * t124;
t103 = t53 * t58;
t102 = t70 * t124;
t99 = -t70 * t29 + t32;
t98 = t53 * t70 - qJ(5);
t96 = t133 * qJD(3);
t95 = -t68 * t103 + t53 * t109 + t29 * t128 - t70 * t36;
t94 = t68 * t102;
t19 = t68 * t29;
t11 = t71 * qJ(5) - t19 - t33;
t10 = -pkin(5) * t140 - t11;
t61 = t71 * pkin(4);
t9 = t71 * qJ(6) + t32 + t61 + (pkin(5) * t69 - t29) * t70;
t89 = t10 * t70 + t68 * t9;
t88 = -t10 * t68 + t70 * t9;
t12 = t61 + t99;
t86 = t11 * t70 - t12 * t68;
t85 = t11 * t68 + t12 * t70;
t83 = qJ(6) * t68 - t132;
t82 = -t70 * qJD(6) - t125;
t55 = pkin(4) * t128;
t15 = t83 * qJD(4) + t55 + t82;
t81 = t28 * t128 - t15 * t70;
t13 = t83 * t69 + t135;
t80 = -qJD(4) * t13 - t28 * t124;
t77 = -pkin(5) * t112 + t95;
t4 = (qJD(5) + t113) * t71 + t98 * t58 - t136;
t5 = t95 - t115;
t76 = t85 * qJD(4) - t4 * t70 + t5 * t68;
t34 = t147 * t128;
t56 = pkin(8) * t60;
t35 = pkin(5) * t60 + t56;
t41 = t147 * t68;
t42 = t147 * t70;
t75 = -t34 * t70 + t35 * t68 + (t41 * t70 - t42 * t68) * qJD(4);
t16 = -t69 * t132 + t135;
t20 = -t70 * t121 - t125 + t55;
t74 = -qJD(4) * t16 - t20 * t69 + (-t38 * t71 + t143) * qJD(3);
t73 = qJ(5) * t117 - 0.2e1 * t123 - t6;
t72 = 0.2e1 * qJD(5);
t48 = t64 * t124;
t39 = t64 * t104;
t25 = -t68 * t58 + t109;
t22 = -t102 + t112;
t21 = t48 + t107;
t14 = 0.2e1 * t39 + 0.2e1 * (t62 - 0.1e1) * t104;
t2 = -t123 + (-pkin(5) * t138 - t32) * qJD(4) + (-pkin(5) * t139 - t98 * t69) * qJD(3) + t136;
t1 = t122 + (t119 + t151) * qJD(3) + t77;
t7 = [0, 0, 0, 0, 0.2e1 * t104, -0.2e1 * t96, 0, 0, 0, t54 * t117, 0.2e1 * t54 * t124, -0.2e1 * t63 * t105 + 0.2e1 * t39, t152 * t94 + 0.2e1 * t63 * t97, 0.2e1 * t69 * t111 + 0.2e1 * t133 * t129, 0.2e1 * t109 * t69 - 0.2e1 * t68 * t96, -0.2e1 * t104, 0.2e1 * t63 * t53 * t60 + 0.2e1 * t95 * t71 + 0.2e1 * (-t99 + 0.2e1 * t32) * t58, -0.2e1 * t63 * t113 - 0.2e1 * t6 * t71 + 0.2e1 * (-t19 + t33) * t58, 0.2e1 * t85 * t124 + 0.2e1 * (t86 * qJD(4) + t4 * t68 + t5 * t70) * t69, 0.2e1 * (-t16 * t130 - t5) * t71 + 0.2e1 * (qJD(3) * t12 - t16 * t60 - t8 * t68) * t69, 0.2e1 * (-t16 * t129 + t4) * t71 + 0.2e1 * (-qJD(3) * t11 + t16 * t128 - t8 * t70) * t69, 0.2e1 * t11 * t4 + 0.2e1 * t12 * t5 + 0.2e1 * t16 * t8, 0.2e1 * t88 * t124 + 0.2e1 * (-t89 * qJD(4) + t1 * t70 - t2 * t68) * t69, 0.2e1 * (-t13 * t129 - t2) * t71 + 0.2e1 * (qJD(3) * t10 + t13 * t128 - t144) * t69, 0.2e1 * (t13 * t130 + t1) * t71 + 0.2e1 * (-qJD(3) * t9 + t13 * t60 + t145) * t69, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t13 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-qJD(3) * t86 - t8) * t71 + (qJD(3) * t16 + t76) * t69, 0, 0, 0 (qJD(3) * t89 - t3) * t71 + (qJD(3) * t13 + qJD(4) * t88 + t1 * t68 + t2 * t70) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, t124, -t58, 0, -t37, t103, -t69 * t97 + t94, t105 * t152 - t107 + t48, -t25, t23, 0 (pkin(8) * t137 + (-pkin(3) * t70 + t53 * t68) * t69) * qJD(4) + (t68 * t92 - t33) * qJD(3) (t53 * t138 + t68 * t91) * qJD(4) + (t70 * t92 + t32) * qJD(3), t76, -t148 * t70 + t74 * t68, t148 * t68 + t74 * t70, pkin(8) * t76 + t16 * t20 + t8 * t38 (t41 * t124 + t35 * t69 + t2 + (-t42 * t69 + t9) * qJD(4)) * t70 + (-t42 * t124 + t34 * t69 + t1 + (-t41 * t69 - t10) * qJD(4)) * t68, -t145 + t34 * t71 + t80 * t70 + (qJD(3) * t42 + t81) * t69, -t144 + t35 * t71 + (-qJD(3) * t41 + t108) * t69 + (t15 * t69 - t80) * t68, t1 * t41 - t10 * t34 + t13 * t15 + t2 * t42 + t3 * t28 + t9 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t124, 0, 0, 0, 0, 0, -t23, -t25, t21, t23, t25, -t71 * t20 + (t142 + (t62 + t64) * t146) * qJD(3), t21, t25, -t23 (-t15 + (t41 * t68 + t42 * t70) * qJD(3)) * t71 + (qJD(3) * t28 + t75) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t105, -0.2e1 * t97, 0, 0, 0, t68 * t118, t70 * t118, 0, -0.2e1 * t38 * t128 + 0.2e1 * t20 * t70, -0.2e1 * t20 * t68 - 0.2e1 * t38 * t60, 0.2e1 * t38 * t20, 0.2e1 * t75, -0.2e1 * t15 * t68 - 0.2e1 * t108, 0.2e1 * t81, 0.2e1 * t28 * t15 - 0.2e1 * t42 * t34 + 0.2e1 * t41 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t24, t58, -t95, t6 (-pkin(4) * t124 - t100) * t70 + (-t101 + (pkin(4) * qJD(4) - qJD(5)) * t69) * t68, t95 - 0.2e1 * t115, t73, -t5 * pkin(4) - t4 * qJ(5) - t11 * qJD(5), t84 * t124 + ((-t67 * t68 - t132) * qJD(4) + t82) * t69, -pkin(5) * t24 + t73, -0.2e1 * t122 + (-t119 - 0.2e1 * t151) * qJD(3) - t77, t2 * qJ(5) + t10 * qJD(5) - t9 * qJD(6) + t1 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t22, 0, t24, -t22, t93, 0, -t22, -t24, t24 * t67 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t128, 0, -t56, t114, t149, t56, -t114, t149 * pkin(8), qJD(4) * t84 - t126 + t59, -t34, -t35, -t34 * qJ(5) + t42 * qJD(5) - t41 * qJD(6) + t35 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0.2e1 * t120, 0, t72, 0.2e1 * qJD(6), -0.2e1 * t67 * qJD(6) + 0.2e1 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t58, 0, t5, -t22, 0, -t58, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, t56, t60, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t58, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, 0, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
