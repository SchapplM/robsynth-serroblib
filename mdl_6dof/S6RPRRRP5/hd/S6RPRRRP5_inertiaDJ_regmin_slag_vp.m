% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:25
% EndTime: 2019-03-09 06:12:29
% DurationCPUTime: 1.70s
% Computational Cost: add. (4676->192), mult. (10077->321), div. (0->0), fcn. (10543->8), ass. (0->126)
t148 = sin(qJ(3));
t149 = cos(qJ(3));
t76 = sin(pkin(10));
t77 = cos(pkin(10));
t160 = t148 * t77 + t149 * t76;
t138 = pkin(7) + qJ(2);
t61 = t138 * t76;
t62 = t138 * t77;
t92 = t148 * t62 + t149 * t61;
t38 = -pkin(8) * t160 - t92;
t91 = t148 * t76 - t149 * t77;
t93 = t148 * t61 - t149 * t62;
t39 = -t91 * pkin(8) - t93;
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t29 = t79 * t38 + t81 * t39;
t41 = t160 * t79 + t81 * t91;
t42 = t160 * t81 - t79 * t91;
t67 = -t77 * pkin(2) - pkin(1);
t44 = t91 * pkin(3) + t67;
t30 = t41 * pkin(4) - t42 * pkin(9) + t44;
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t159 = t80 * t29 + t78 * t30;
t134 = t160 * qJD(3);
t53 = t91 * qJD(3);
t31 = -t41 * qJD(4) - t79 * t134 - t81 * t53;
t71 = qJD(5) * t80;
t97 = t78 * t31 + t42 * t71;
t74 = t78 ^ 2;
t75 = t80 ^ 2;
t133 = t74 - t75;
t108 = t133 * qJD(5);
t113 = t134 * pkin(3);
t32 = t42 * qJD(4) + t81 * t134 - t79 * t53;
t17 = t32 * pkin(4) - t31 * pkin(9) + t113;
t128 = qJD(4) * t81;
t129 = qJD(4) * t79;
t157 = t91 * qJD(2) + t92 * qJD(3);
t83 = -t134 * pkin(8) - t157;
t85 = -qJD(2) * t160 + t93 * qJD(3);
t84 = t53 * pkin(8) + t85;
t9 = -t38 * t128 + t39 * t129 - t79 * t84 - t81 * t83;
t5 = -qJD(5) * t159 + t80 * t17 + t78 * t9;
t11 = t41 * qJ(6) + t159;
t100 = -t78 * t29 + t80 * t30;
t12 = -t41 * pkin(5) - t100;
t103 = t11 * t80 + t12 * t78;
t125 = t41 * qJD(6);
t131 = t32 * qJ(6);
t126 = qJD(5) * t78;
t4 = t29 * t126 - t78 * t17 - t30 * t71 + t80 * t9;
t2 = t125 - t4 + t131;
t152 = t32 * pkin(5);
t3 = -t152 - t5;
t158 = t103 * qJD(5) + t2 * t78 - t3 * t80;
t156 = 0.2e1 * t67;
t155 = 0.2e1 * qJD(6);
t154 = pkin(9) * t32;
t153 = pkin(9) * t41;
t151 = t81 * pkin(3);
t82 = t79 * t83 - t81 * t84;
t10 = t29 * qJD(4) + t82;
t28 = -t81 * t38 + t79 * t39;
t150 = t10 * t78 + t28 * t71;
t68 = t79 * pkin(3) + pkin(9);
t147 = t32 * t68;
t146 = t41 * t68;
t145 = t42 * t31;
t144 = t42 * t78;
t143 = t42 * t80;
t141 = t78 * t32;
t140 = t80 * t31;
t139 = t80 * t32;
t120 = pkin(3) * t129;
t124 = t78 * qJD(6);
t51 = -pkin(5) * t126 + qJ(6) * t71 + t124;
t43 = -t51 + t120;
t136 = -t43 + t51;
t69 = -pkin(4) - t151;
t135 = t78 * t120 + t69 * t71;
t132 = t74 + t75;
t130 = t80 * qJ(6);
t104 = t78 * pkin(5) - t130;
t18 = t104 * t42 + t28;
t127 = qJD(5) * t18;
t123 = t80 * qJD(6);
t122 = pkin(4) * t126;
t121 = pkin(4) * t71;
t119 = pkin(3) * t128;
t118 = pkin(9) * t126;
t117 = pkin(9) * t71;
t116 = t42 * t126;
t114 = t78 * t71;
t109 = -0.4e1 * t78 * t143;
t107 = 0.2e1 * (t76 ^ 2 + t77 ^ 2) * qJD(2);
t105 = -t80 * pkin(5) - t78 * qJ(6);
t102 = -t11 * t78 + t12 * t80;
t99 = -t42 * t69 + t146;
t49 = t132 * t119;
t98 = -t80 * t120 + t69 * t126;
t60 = -pkin(4) + t105;
t96 = t116 - t140;
t22 = t41 * t71 + t141;
t95 = t41 * t126 - t139;
t94 = t31 * t60 - t42 * t51 - t154;
t6 = t97 * pkin(5) + qJ(6) * t116 - t42 * t123 + t39 * t128 + t38 * t129 - t31 * t130 + t82;
t89 = -t6 + (t42 * t60 - t153) * qJD(5);
t55 = t60 - t151;
t88 = -t6 + (t42 * t55 - t146) * qJD(5);
t50 = t105 * qJD(5) + t123;
t87 = -t41 * t119 + t31 * t55 + t42 * t43 - t147;
t1 = t102 * qJD(5) + t2 * t80 + t3 * t78;
t86 = t31 * t69 - t147 + (-t41 * t81 + t42 * t79) * qJD(4) * pkin(3);
t64 = 0.2e1 * t114;
t57 = -0.2e1 * t108;
t54 = t60 * t126;
t48 = t55 * t126;
t46 = t78 * t119 + t68 * t71;
t45 = -t80 * t119 + t68 * t126;
t40 = t42 ^ 2;
t23 = t28 * t126;
t19 = -t42 * t108 + t78 * t140;
t16 = t18 * t126;
t13 = qJD(5) * t109 - t133 * t31;
t7 = [0, 0, 0, 0, 0, t107, qJ(2) * t107, -0.2e1 * t160 * t53, -0.2e1 * t134 * t160 + 0.2e1 * t53 * t91, 0, 0, 0, t134 * t156, -t53 * t156, 0.2e1 * t145, -0.2e1 * t31 * t41 - 0.2e1 * t42 * t32, 0, 0, 0, 0.2e1 * t41 * t113 + 0.2e1 * t44 * t32, 0.2e1 * t42 * t113 + 0.2e1 * t44 * t31, -0.2e1 * t40 * t114 + 0.2e1 * t75 * t145, 0.2e1 * t40 * t108 + t31 * t109, 0.2e1 * t42 * t139 - 0.2e1 * t96 * t41, -0.2e1 * t42 * t141 - 0.2e1 * t97 * t41, 0.2e1 * t41 * t32, 0.2e1 * t10 * t144 + 0.2e1 * t100 * t32 + 0.2e1 * t97 * t28 + 0.2e1 * t5 * t41, 0.2e1 * t10 * t143 - 0.2e1 * t159 * t32 - 0.2e1 * t96 * t28 + 0.2e1 * t4 * t41, -0.2e1 * t12 * t32 + 0.2e1 * t6 * t144 + 0.2e1 * t97 * t18 - 0.2e1 * t3 * t41, 0.2e1 * t102 * t31 - 0.2e1 * t158 * t42, 0.2e1 * t11 * t32 - 0.2e1 * t6 * t143 + 0.2e1 * t96 * t18 + 0.2e1 * t2 * t41, 0.2e1 * t11 * t2 + 0.2e1 * t12 * t3 + 0.2e1 * t18 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t53, 0, 0, 0, 0, 0, t32, t31, 0, 0, 0, 0, 0, -t95, -t22, -t95, -t132 * t31, t22, t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t134, 0, t85, t157, 0, 0, t31, -t32, 0, -t10, t9, t19, t13, t22, -t95, 0, t23 + (-t99 * qJD(5) - t10) * t80 + t86 * t78, t99 * t126 + t86 * t80 + t150, t87 * t78 + t88 * t80 + t16, t1, t88 * t78 + (-t87 - t127) * t80, t1 * t68 + t103 * t119 + t18 * t43 + t6 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t120, -0.2e1 * t119, t64, t57, 0, 0, 0, 0.2e1 * t98, 0.2e1 * t135, -0.2e1 * t43 * t80 + 0.2e1 * t48, 0.2e1 * t49, -0.2e1 * t43 * t78 - 0.2e1 * t55 * t71, 0.2e1 * t55 * t43 + 0.2e1 * t49 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t32, 0, -t10, t9, t19, t13, t22, -t95, 0, t23 + (-pkin(4) * t31 - t154) * t78 + (-t10 + (-pkin(4) * t42 - t153) * qJD(5)) * t80, t96 * pkin(4) + t95 * pkin(9) + t150, t94 * t78 + t89 * t80 + t16, t1, t89 * t78 + (-t94 - t127) * t80, t1 * pkin(9) - t18 * t51 + t6 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, t64, t57, 0, 0, 0, t98 - t122, -t121 + t135, t136 * t80 + t48 + t54, t49, t136 * t78 + (-t55 - t60) * t71, pkin(9) * t49 + t43 * t60 - t55 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t57, 0, 0, 0, -0.2e1 * t122, -0.2e1 * t121, 0.2e1 * t51 * t80 + 0.2e1 * t54, 0, 0.2e1 * t51 * t78 - 0.2e1 * t60 * t71, -0.2e1 * t60 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t97, t32, t5, t4, t5 + 0.2e1 * t152, t105 * t31 + (t104 * qJD(5) - t124) * t42, 0.2e1 * t125 - t4 + 0.2e1 * t131, -t3 * pkin(5) + t2 * qJ(6) + t11 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t71, -t126, 0, t71, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t126, 0, -t46, t45, -t46, t50, -t45, -t104 * t119 + t50 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t126, 0, -t117, t118, -t117, t50, -t118, t50 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, qJ(6) * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t96, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
