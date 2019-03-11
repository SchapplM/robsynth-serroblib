% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:24
% EndTime: 2019-03-08 19:05:28
% DurationCPUTime: 1.40s
% Computational Cost: add. (1331->201), mult. (4069->378), div. (0->0), fcn. (4215->14), ass. (0->140)
t81 = sin(qJ(4));
t159 = -0.4e1 * t81;
t79 = sin(qJ(6));
t80 = sin(qJ(5));
t83 = cos(qJ(6));
t84 = cos(qJ(5));
t51 = t79 * t84 + t83 * t80;
t39 = t51 * t81;
t85 = cos(qJ(4));
t103 = -t85 * pkin(4) - t81 * pkin(10);
t58 = -pkin(3) + t103;
t145 = t84 * t85;
t65 = pkin(9) * t145;
t142 = t80 * t58 + t65;
t74 = sin(pkin(7));
t86 = cos(qJ(3));
t152 = t74 * t86;
t117 = qJD(3) * t152;
t82 = sin(qJ(3));
t153 = t74 * t82;
t77 = cos(pkin(7));
t42 = t81 * t153 - t85 * t77;
t158 = -qJD(4) * t42 + t85 * t117;
t130 = t85 * qJD(4);
t112 = t84 * t130;
t134 = qJD(5) * t80;
t157 = -t81 * t134 + t112;
t71 = t84 ^ 2;
t141 = t80 ^ 2 - t71;
t108 = t141 * qJD(5);
t156 = qJD(5) + qJD(6);
t133 = qJD(5) * t84;
t102 = pkin(4) * t81 - pkin(10) * t85;
t54 = t102 * qJD(4);
t132 = qJD(5) * t85;
t120 = t80 * t132;
t68 = t81 * qJD(4);
t88 = t84 * t68 + t120;
t24 = t88 * pkin(9) - t58 * t133 - t80 * t54;
t155 = pkin(10) + pkin(11);
t154 = pkin(9) * t80;
t76 = cos(pkin(13));
t151 = t76 * t77;
t149 = t80 * t81;
t37 = -pkin(11) * t149 + t142;
t150 = t79 * t37;
t148 = t81 * t84;
t115 = t80 * t130;
t87 = t81 * t133 + t115;
t17 = -t87 * pkin(11) - t24;
t147 = t83 * t17;
t146 = t83 * t37;
t116 = t80 * t68;
t143 = pkin(9) * t116 + t84 * t54;
t70 = t81 ^ 2;
t140 = -t85 ^ 2 + t70;
t139 = qJD(3) * t82;
t138 = qJD(3) * t85;
t136 = qJD(4) * t80;
t135 = qJD(4) * t84;
t131 = qJD(6) * t79;
t129 = -0.2e1 * pkin(3) * qJD(4);
t128 = -0.2e1 * pkin(4) * qJD(5);
t127 = t86 * t151;
t126 = pkin(5) * t134;
t125 = pkin(5) * t68;
t124 = pkin(5) * t131;
t123 = qJD(6) * t83 * pkin(5);
t122 = pkin(9) * t130;
t119 = t84 * t132;
t118 = t74 * t139;
t114 = t80 * t133;
t113 = t81 * t130;
t49 = t84 * t58;
t29 = -pkin(11) * t148 + t49 + (-pkin(5) - t154) * t85;
t111 = t85 * pkin(5) - t29;
t110 = qJD(5) * t155;
t16 = (pkin(5) * t81 - pkin(11) * t145) * qJD(4) + (-t65 + (pkin(11) * t81 - t58) * t80) * qJD(5) + t143;
t109 = t83 * t16 - t79 * t17;
t107 = t140 * qJD(4);
t106 = 0.2e1 * t113;
t104 = t80 * t112;
t73 = sin(pkin(13));
t75 = sin(pkin(6));
t78 = cos(pkin(6));
t31 = t78 * t153 + (t82 * t151 + t73 * t86) * t75;
t41 = -t75 * t76 * t74 + t78 * t77;
t23 = t31 * t85 + t41 * t81;
t30 = -t78 * t152 + (t73 * t82 - t127) * t75;
t11 = -t23 * t80 + t30 * t84;
t12 = t23 * t84 + t30 * t80;
t101 = t83 * t11 - t79 * t12;
t100 = t79 * t11 + t83 * t12;
t99 = t83 * t29 - t150;
t98 = t79 * t29 + t146;
t22 = t31 * t81 - t41 * t85;
t43 = t85 * t153 + t81 * t77;
t36 = -t80 * t152 + t84 * t43;
t93 = t84 * t152 + t80 * t43;
t97 = -t79 * t36 - t83 * t93;
t96 = t83 * t36 - t79 * t93;
t61 = t155 * t80;
t62 = t155 * t84;
t95 = -t83 * t61 - t79 * t62;
t94 = -t79 * t61 + t83 * t62;
t50 = t79 * t80 - t83 * t84;
t26 = -t78 * t117 + (-qJD(3) * t127 + t139 * t73) * t75;
t9 = t23 * qJD(4) - t26 * t81;
t92 = t22 * t133 + t9 * t80;
t91 = t22 * t134 - t9 * t84;
t34 = t43 * qJD(4) + t81 * t117;
t90 = t42 * t133 + t34 * t80;
t89 = t42 * t134 - t34 * t84;
t67 = -t84 * pkin(5) - pkin(4);
t64 = -0.2e1 * t113;
t55 = (pkin(5) * t80 + pkin(9)) * t81;
t53 = t84 * t110;
t52 = t80 * t110;
t40 = t50 * t81;
t38 = t87 * pkin(5) + t122;
t33 = t156 * t51;
t32 = t156 * t50;
t27 = t31 * qJD(3);
t25 = -t142 * qJD(5) + t143;
t21 = -t94 * qJD(6) + t79 * t52 - t83 * t53;
t20 = -t95 * qJD(6) + t83 * t52 + t79 * t53;
t19 = -t131 * t149 + (t156 * t148 + t115) * t83 + t157 * t79;
t18 = -t50 * t130 - t156 * t39;
t14 = -t77 * t115 - t43 * t133 + (t84 * t139 + (t82 * t68 + (qJD(5) - t138) * t86) * t80) * t74;
t13 = t93 * qJD(5) - t80 * t118 - t84 * t158;
t10 = -qJD(4) * t22 - t26 * t85;
t8 = -t98 * qJD(6) + t109;
t7 = -t99 * qJD(6) - t79 * t16 - t147;
t6 = -t96 * qJD(6) + t79 * t13 + t83 * t14;
t5 = -t97 * qJD(6) + t83 * t13 - t79 * t14;
t4 = t11 * qJD(5) + t10 * t84 + t27 * t80;
t3 = -t12 * qJD(5) - t10 * t80 + t27 * t84;
t2 = -t100 * qJD(6) + t83 * t3 - t79 * t4;
t1 = -t101 * qJD(6) - t79 * t3 - t83 * t4;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t27, t26, 0, 0, 0, 0, 0, -t27 * t85 + t30 * t68, t30 * t130 + t27 * t81, 0, 0, 0, 0, 0 (t22 * t136 - t3) * t85 + (qJD(4) * t11 + t92) * t81 (t22 * t135 + t4) * t85 + (-qJD(4) * t12 - t91) * t81, 0, 0, 0, 0, 0, t101 * t68 + t22 * t19 - t2 * t85 + t9 * t39, -t1 * t85 - t100 * t68 + t22 * t18 - t9 * t40; 0, 0, 0, -t118, -t117, 0, 0, 0, 0, 0 (-t82 * t138 - t86 * t68) * t74 (-t86 * t130 + t81 * t139) * t74, 0, 0, 0, 0, 0 (t42 * t136 - t14) * t85 + (-qJD(4) * t93 + t90) * t81 (t42 * t135 - t13) * t85 + (-qJD(4) * t36 - t89) * t81, 0, 0, 0, 0, 0, t42 * t19 + t34 * t39 - t6 * t85 + t68 * t97, t42 * t18 - t34 * t40 - t5 * t85 - t68 * t96; 0, 0, 0, 0, 0, t106, -0.2e1 * t107, 0, 0, 0, t81 * t129, t85 * t129, 0.2e1 * t71 * t113 - 0.2e1 * t70 * t114, t104 * t159 + 0.2e1 * t70 * t108, 0.2e1 * t81 * t120 + 0.2e1 * t140 * t135, -0.2e1 * t80 * t107 + 0.2e1 * t81 * t119, t64, 0.2e1 * t49 * t68 - 0.2e1 * t25 * t85 + 0.2e1 * (t80 * t113 + t70 * t133) * pkin(9), -0.2e1 * t24 * t85 - 0.2e1 * t142 * t68 + 0.2e1 * (t84 * t106 - t70 * t134) * pkin(9), -0.2e1 * t40 * t18, -0.2e1 * t18 * t39 + 0.2e1 * t40 * t19, -0.2e1 * t18 * t85 - 0.2e1 * t40 * t68, 0.2e1 * t19 * t85 - 0.2e1 * t39 * t68, t64, 0.2e1 * t55 * t19 + 0.2e1 * t38 * t39 + 0.2e1 * t68 * t99 - 0.2e1 * t8 * t85, 0.2e1 * t55 * t18 - 0.2e1 * t38 * t40 - 0.2e1 * t68 * t98 - 0.2e1 * t7 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, t91, t92, 0, 0, 0, 0, 0, t22 * t33 + t9 * t50, -t22 * t32 + t9 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t158, 0, 0, 0, 0, 0, t89, t90, 0, 0, 0, 0, 0, t42 * t33 + t34 * t50, -t42 * t32 + t34 * t51; 0, 0, 0, 0, 0, 0, 0, t130, -t68, 0, -t122, pkin(9) * t68, -t81 * t108 + t104, t114 * t159 - t141 * t130, t116 - t119, t88, 0 (pkin(10) * t145 + (-pkin(4) * t84 + t154) * t81) * qJD(5) + (t103 * t80 - t65) * qJD(4) (pkin(9) * t148 + t102 * t80) * qJD(5) + (t103 * t84 + t85 * t154) * qJD(4), t18 * t51 + t40 * t32, -t18 * t50 - t51 * t19 + t32 * t39 + t40 * t33, t32 * t85 + t51 * t68, t33 * t85 - t50 * t68, 0, t126 * t39 + t67 * t19 - t21 * t85 + t55 * t33 + t38 * t50 + t68 * t95, -t126 * t40 + t67 * t18 - t20 * t85 - t55 * t32 + t38 * t51 - t68 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t114, -0.2e1 * t108, 0, 0, 0, t80 * t128, t84 * t128, -0.2e1 * t51 * t32, 0.2e1 * t32 * t50 - 0.2e1 * t51 * t33, 0, 0, 0, 0.2e1 * t126 * t50 + 0.2e1 * t67 * t33, 0.2e1 * t126 * t51 - 0.2e1 * t67 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, -t87, t68, t25, t24, 0, 0, t18, -t19, t68, t83 * t125 + (t111 * t79 - t146) * qJD(6) + t109, -t147 + (-t16 - t125) * t79 + (t111 * t83 + t150) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, -t134, 0, -pkin(10) * t133, pkin(10) * t134, 0, 0, -t32, -t33, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t124, -0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, t68, t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
