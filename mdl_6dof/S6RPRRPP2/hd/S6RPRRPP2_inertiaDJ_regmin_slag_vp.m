% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP2
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:21
% EndTime: 2019-03-09 04:33:26
% DurationCPUTime: 1.74s
% Computational Cost: add. (1222->221), mult. (2697->343), div. (0->0), fcn. (1977->6), ass. (0->127)
t63 = sin(qJ(3));
t126 = qJD(4) * t63;
t64 = cos(qJ(4));
t121 = t64 * qJD(5);
t117 = qJ(5) * qJD(3);
t65 = cos(qJ(3));
t99 = t65 * t117;
t136 = t63 * t121 + t64 * t99;
t62 = sin(qJ(4));
t131 = t62 * qJ(5);
t145 = pkin(4) + pkin(5);
t79 = t145 * t64 + t131;
t153 = -t126 * t79 + t136;
t152 = t145 * qJD(3);
t151 = -0.4e1 * t63;
t54 = qJD(4) * t64;
t108 = t63 * t54;
t123 = t62 * qJD(6);
t119 = t65 * qJD(3);
t97 = qJ(6) * t119;
t150 = qJ(6) * t108 + t63 * t123 + t62 * t97;
t105 = t62 * t119;
t26 = t105 + t108;
t98 = qJ(5) * t126;
t149 = -t26 * pkin(4) - t62 * t98 + t136;
t125 = qJD(4) * t65;
t109 = t62 * t125;
t122 = t63 * qJD(3);
t25 = t64 * t122 + t109;
t51 = -cos(pkin(9)) * pkin(1) - pkin(2);
t141 = t63 * pkin(8);
t90 = -t65 * pkin(3) - t141;
t30 = t51 + t90;
t144 = pkin(8) * t65;
t89 = pkin(3) * t63 - t144;
t36 = t89 * qJD(3);
t50 = sin(pkin(9)) * pkin(1) + pkin(7);
t6 = t25 * t50 - t30 * t54 - t62 * t36;
t133 = qJ(5) * t64;
t84 = pkin(4) * t62 - t133;
t15 = (t50 + t84) * t63;
t124 = t62 * qJD(5);
t21 = t84 * qJD(4) - t124;
t85 = -t64 * pkin(4) - t131;
t37 = -pkin(3) + t85;
t148 = qJD(3) * (-t37 * t65 + t141) - qJD(4) * t15 - t21 * t63;
t120 = t64 * qJD(6);
t127 = qJD(4) * t62;
t138 = pkin(8) - qJ(6);
t20 = -t138 * t127 - t120;
t40 = t138 * t64;
t22 = qJD(4) * t40 - t123;
t39 = t138 * t62;
t147 = -(t39 * t64 - t40 * t62) * qJD(4) - t20 * t64 - t22 * t62;
t67 = 0.2e1 * qJD(5);
t111 = t145 * t62;
t4 = (-t50 - t111) * t119 + t153;
t143 = t4 * t62;
t142 = t4 * t64;
t140 = t37 * t63;
t139 = t50 * t65;
t57 = t62 ^ 2;
t59 = t64 ^ 2;
t135 = t57 - t59;
t58 = t63 ^ 2;
t134 = -t65 ^ 2 + t58;
t132 = qJ(6) * t63;
t130 = qJD(3) * t62;
t129 = qJD(3) * t64;
t128 = qJD(4) * t58;
t118 = t65 * qJD(5);
t116 = -0.2e1 * pkin(3) * qJD(4);
t34 = t62 * t139;
t35 = t64 * t139;
t115 = 0.2e1 * qJD(3) * t51;
t114 = pkin(4) * t122;
t113 = pkin(8) * t127;
t112 = pkin(8) * t54;
t110 = t50 * t128;
t107 = t64 * t125;
t31 = pkin(3) + t79;
t106 = t31 * t54;
t104 = t62 * t54;
t103 = t63 * t119;
t102 = t50 * t122;
t101 = t64 * t119;
t100 = t50 * t119;
t96 = (t57 + t59) * t65;
t95 = -t64 * t30 + t34;
t94 = t135 * qJD(4);
t93 = t134 * qJD(3);
t92 = -t62 * t102 + t50 * t107 + t30 * t127 - t64 * t36;
t91 = t62 * t101;
t19 = t62 * t30;
t11 = -t65 * qJ(5) + t19 + t35;
t10 = t62 * t132 + t11;
t55 = t65 * pkin(4);
t9 = t65 * pkin(5) + t34 + t55 + (-t30 - t132) * t64;
t87 = t10 * t64 + t62 * t9;
t86 = t10 * t62 - t64 * t9;
t12 = t55 + t95;
t83 = t11 * t64 + t12 * t62;
t82 = -t11 * t62 + t12 * t64;
t76 = -t111 + t133;
t16 = t76 * qJD(4) + t124;
t77 = -t31 * t127 + t16 * t64;
t13 = (-t50 + t76) * t63;
t75 = qJD(4) * t13 + t31 * t119;
t74 = -qJ(6) * t127 + t120;
t72 = -t64 * t97 + t92;
t71 = t85 * qJD(4) + t121;
t52 = t63 * t117;
t3 = t52 - t6 - t118;
t5 = t92 - t114;
t69 = t82 * qJD(4) + t3 * t64 + t5 * t62;
t68 = 0.2e1 * t52 - 0.2e1 * t118 - t6;
t56 = qJ(5) * t67;
t44 = pkin(8) * t107;
t38 = t59 * t103;
t27 = -t62 * t122 + t107;
t24 = t62 * t126 - t101;
t23 = qJD(3) * t96;
t14 = 0.2e1 * t38 + 0.2e1 * (t57 - 0.1e1) * t103;
t8 = t100 - t149;
t2 = t3 + t150;
t1 = (-t74 - t152) * t63 + t72;
t7 = [0, 0, 0, 0, 0.2e1 * t103, -0.2e1 * t93, 0, 0, 0, t63 * t115, t65 * t115, -0.2e1 * t58 * t104 + 0.2e1 * t38, 0.2e1 * t135 * t128 + t91 * t151, 0.2e1 * t63 * t109 + 0.2e1 * t134 * t129, 0.2e1 * t107 * t63 - 0.2e1 * t62 * t93, -0.2e1 * t103, 0.2e1 * t64 * t110 + 0.2e1 * t92 * t65 + 0.2e1 * (-t95 + 0.2e1 * t34) * t122, -0.2e1 * t62 * t110 - 0.2e1 * t6 * t65 + 0.2e1 * (-t19 + t35) * t122, 0.2e1 * (t130 * t15 + t5) * t65 + 0.2e1 * (-qJD(3) * t12 + t15 * t54 + t8 * t62) * t63, 0.2e1 * t82 * t119 + 0.2e1 * (-qJD(4) * t83 - t3 * t62 + t5 * t64) * t63, 0.2e1 * (-t129 * t15 - t3) * t65 + 0.2e1 * (qJD(3) * t11 + t127 * t15 - t8 * t64) * t63, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t5 + 0.2e1 * t15 * t8, 0.2e1 * (-t13 * t130 + t1) * t65 + 0.2e1 * (-qJD(3) * t9 - t13 * t54 - t143) * t63, 0.2e1 * (t129 * t13 - t2) * t65 + 0.2e1 * (qJD(3) * t10 - t127 * t13 + t142) * t63, 0.2e1 * t86 * t119 + 0.2e1 * (qJD(4) * t87 - t1 * t64 + t2 * t62) * t63, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t13 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(3) * t83 - t8) * t65 + (qJD(3) * t15 + t69) * t63, 0, 0, 0 (qJD(3) * t87 + t4) * t65 + (-qJD(3) * t13 - qJD(4) * t86 + t1 * t62 + t2 * t64) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, t119, -t122, 0, -t100, t102, -t63 * t94 + t91, t104 * t151 - t135 * t119, -t27, t25, 0, t44 + (-pkin(3) * t64 + t50 * t62) * t126 + (t62 * t90 - t35) * qJD(3) (t50 * t63 * t64 + t62 * t89) * qJD(4) + (t64 * t90 + t34) * qJD(3), t44 + (t126 * t37 - t8) * t64 - t148 * t62, t69 (-t8 + (t140 + t144) * qJD(4)) * t62 + t148 * t64, pkin(8) * t69 + t15 * t21 + t8 * t37, t22 * t65 + t142 + (-qJD(3) * t39 - t106) * t63 + (-t16 * t63 - t75) * t62, -t20 * t65 + t143 + t75 * t64 + (qJD(3) * t40 + t77) * t63 (-t39 * t119 - t22 * t63 - t2 + (t40 * t63 - t9) * qJD(4)) * t64 + (t40 * t119 + t20 * t63 - t1 + (t39 * t63 + t10) * qJD(4)) * t62, t1 * t39 + t10 * t20 + t13 * t16 + t2 * t40 + t9 * t22 + t4 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t119, 0, 0, 0, 0, 0, -t25, -t27, -t25, t23, t27, -t65 * t21 + (pkin(8) * t96 + t140) * qJD(3), -t25, t27, -t23 (t16 + (t39 * t62 + t40 * t64) * qJD(3)) * t65 + (-qJD(3) * t31 - t147) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t104, -0.2e1 * t94, 0, 0, 0, t62 * t116, t64 * t116, 0.2e1 * t127 * t37 - 0.2e1 * t21 * t64, 0, -0.2e1 * t21 * t62 - 0.2e1 * t37 * t54, 0.2e1 * t37 * t21, 0.2e1 * t77, 0.2e1 * t16 * t62 + 0.2e1 * t106, 0.2e1 * t147, 0.2e1 * t31 * t16 + 0.2e1 * t40 * t20 + 0.2e1 * t39 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, t122, -t92, t6, -t92 + 0.2e1 * t114 (-pkin(4) * t119 - t98) * t64 + (-t99 + (pkin(4) * qJD(4) - qJD(5)) * t63) * t62, t68, -t5 * pkin(4) + t3 * qJ(5) + t11 * qJD(5) (t74 + 0.2e1 * t152) * t63 - t72, t68 + t150 (t119 * t145 + t98) * t64 + (t99 + (-qJD(4) * t145 + qJD(5)) * t63) * t62, t2 * qJ(5) + t10 * qJD(5) - t1 * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t24, -t26, 0, -t24, t149, -t26, -t24, 0, -t105 * t145 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t127, 0, -t112, t113, -t112, t71, -t113, t71 * pkin(8), -t22, t20, qJD(4) * t79 - t121, t20 * qJ(5) + t40 * qJD(5) - t145 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t56, 0, t67, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t24, 0, t5, -t122, 0, t24, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t112, 0, 0, -t54, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t24, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t54, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
