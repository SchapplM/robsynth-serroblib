% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:43
% EndTime: 2019-03-09 10:47:48
% DurationCPUTime: 1.37s
% Computational Cost: add. (2054->175), mult. (4437->319), div. (0->0), fcn. (4147->8), ass. (0->118)
t77 = cos(qJ(6));
t72 = t77 ^ 2;
t74 = sin(qJ(6));
t129 = t74 ^ 2 - t72;
t145 = qJD(6) * t129;
t149 = 0.2e1 * t145;
t76 = sin(qJ(2));
t79 = cos(qJ(2));
t148 = -t79 * pkin(2) - t76 * qJ(3);
t75 = sin(qJ(4));
t134 = t76 * t75;
t78 = cos(qJ(4));
t95 = t78 * t79 + t134;
t35 = (qJD(2) - qJD(4)) * t95;
t52 = -t75 * t79 + t76 * t78;
t140 = pkin(7) - pkin(8);
t112 = t140 * t79;
t105 = t78 * t112;
t83 = (t140 * t134 + t105) * qJD(2);
t147 = -t35 * qJ(5) - t52 * qJD(5) + t83;
t113 = t140 * t76;
t146 = (t75 * t113 + t105) * qJD(4);
t128 = cos(pkin(10));
t73 = sin(pkin(10));
t29 = t128 * t95 + t52 * t73;
t122 = t76 * qJD(2);
t111 = t78 * t122;
t109 = qJD(4) * t140;
t104 = t76 * t109;
t106 = t75 * t112;
t117 = -qJD(2) * t106 - t78 * t104 + t140 * t111;
t120 = -qJ(5) + t140;
t127 = qJD(4) * t75;
t126 = qJD(4) * t78;
t67 = t79 * qJD(2);
t131 = -t76 * t126 - t75 * t67;
t11 = (t111 + t131) * qJ(5) - qJD(5) * t134 + (-t78 * qJD(5) - t120 * t127) * t79 - t117;
t3 = t11 * t73 - t128 * (-t146 + t147);
t30 = t128 * t52 - t73 * t95;
t141 = pkin(2) + pkin(3);
t114 = t78 * t141;
t53 = -t75 * qJ(3) - pkin(4) - t114;
t54 = t78 * qJ(3) - t75 * t141;
t33 = t128 * t53 - t73 * t54;
t31 = pkin(5) - t33;
t34 = t128 * t54 + t73 * t53;
t32 = -pkin(9) + t34;
t144 = -t3 + (t29 * t32 - t30 * t31) * qJD(6);
t143 = t75 * qJD(3) + t54 * qJD(4);
t142 = 0.2e1 * qJD(3);
t139 = t3 * t77;
t93 = t79 * t127 + t131;
t87 = t93 + t111;
t17 = t128 * t35 + t73 * t87;
t138 = t17 * t77;
t38 = qJ(3) * t127 - t78 * qJD(3) + qJD(4) * t114;
t21 = t128 * t143 - t38 * t73;
t137 = t21 * t74;
t136 = t21 * t77;
t135 = t30 * t74;
t16 = -t128 * t87 + t35 * t73;
t133 = t77 * t16;
t130 = qJ(3) * t67 + t76 * qJD(3);
t125 = qJD(6) * t74;
t124 = qJD(6) * t77;
t121 = -0.2e1 * pkin(1) * qJD(2);
t55 = -pkin(1) + t148;
t119 = t74 * t138;
t63 = -t128 * pkin(4) - pkin(5);
t118 = 0.2e1 * qJD(6) * t63;
t116 = pkin(7) * t122;
t115 = pkin(7) * t67;
t110 = t74 * t124;
t47 = t79 * pkin(3) - t55;
t108 = qJD(6) * (t31 - t63);
t103 = t79 * t109;
t84 = t95 * pkin(4) + t47;
t12 = t29 * pkin(5) - t30 * pkin(9) + t84;
t24 = t95 * t120;
t81 = -t52 * qJ(5) + t78 * t113 - t106;
t14 = t128 * t24 + t73 * t81;
t102 = t12 * t77 - t14 * t74;
t101 = t12 * t74 + t14 * t77;
t62 = pkin(4) * t73 + pkin(9);
t100 = -t16 * t62 + t17 * t63;
t22 = -t128 * t38 - t143 * t73;
t99 = t21 * t30 - t22 * t29;
t50 = t128 * t75 + t73 * t78;
t90 = t128 * t78 - t73 * t75;
t97 = t29 * t50 + t30 * t90;
t96 = t29 * t62 - t30 * t63;
t9 = t29 * t124 + t16 * t74;
t92 = t30 * t124 + t17 * t74;
t91 = -t30 * t125 + t138;
t41 = t50 * qJD(4);
t42 = t90 * qJD(4);
t86 = -t16 * t50 - t17 * t90 - t29 * t42 + t30 * t41;
t85 = t148 * qJD(2) + t79 * qJD(3);
t13 = -t128 * t81 + t24 * t73;
t80 = -qJD(6) * t13 - t16 * t32 + t17 * t31 + t99;
t20 = -t93 * pkin(4) + (-pkin(4) * t78 - t141) * t122 + t130;
t57 = 0.2e1 * t110;
t51 = -0.2e1 * t145;
t39 = pkin(2) * t122 - t130;
t36 = -t141 * t122 + t130;
t27 = t30 ^ 2;
t26 = t125 * t90 + t41 * t77;
t25 = -t124 * t90 + t41 * t74;
t19 = -t83 + t146;
t18 = t75 * t103 + t117;
t8 = t29 * t125 - t133;
t7 = t145 * t30 - t119;
t6 = 0.4e1 * t30 * t110 + t129 * t17;
t5 = t16 * pkin(5) - t17 * pkin(9) + t20;
t4 = t128 * t11 + (-t78 * t103 - t75 * t104 + t147) * t73;
t2 = -t101 * qJD(6) - t74 * t4 + t77 * t5;
t1 = -t102 * qJD(6) - t77 * t4 - t74 * t5;
t10 = [0, 0, 0, 0.2e1 * t76 * t67, 0.2e1 * (-t76 ^ 2 + t79 ^ 2) * qJD(2), 0, 0, 0, t76 * t121, t79 * t121, 0.2e1 * t55 * t122 - 0.2e1 * t39 * t79, 0, -0.2e1 * t39 * t76 - 0.2e1 * t55 * t67, 0.2e1 * t55 * t39, 0.2e1 * t52 * t35, -0.2e1 * t35 * t95 + 0.2e1 * t52 * t87, 0, 0, 0, 0.2e1 * t36 * t95 - 0.2e1 * t47 * t87, 0.2e1 * t35 * t47 + 0.2e1 * t36 * t52, 0.2e1 * t13 * t17 - 0.2e1 * t14 * t16 - 0.2e1 * t29 * t4 + 0.2e1 * t3 * t30, 0.2e1 * t13 * t3 + 0.2e1 * t14 * t4 + 0.2e1 * t84 * t20, 0.2e1 * t17 * t30 * t72 - 0.2e1 * t27 * t110, -0.4e1 * t30 * t119 + t27 * t149, 0.2e1 * t30 * t133 + 0.2e1 * t91 * t29, -0.2e1 * t16 * t135 - 0.2e1 * t92 * t29, 0.2e1 * t29 * t16, 0.2e1 * t102 * t16 + 0.2e1 * t92 * t13 + 0.2e1 * t3 * t135 + 0.2e1 * t2 * t29, 0.2e1 * t1 * t29 - 0.2e1 * t101 * t16 + 0.2e1 * t91 * t13 + 0.2e1 * t30 * t139; 0, 0, 0, 0, 0, t67, -t122, 0, -t115, t116, -t115, t85, -t116, t85 * pkin(7), 0, 0, -t35, -t87, 0, t19, -t18, -t16 * t34 - t17 * t33 + t99, t13 * t21 + t14 * t22 - t3 * t33 + t34 * t4, t7, t6, -t9, t8, 0, -t144 * t77 + t80 * t74, t144 * t74 + t80 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, qJ(3) * t142, 0, 0, 0, 0, 0, 0.2e1 * t143, -0.2e1 * t38, 0, -0.2e1 * t21 * t33 + 0.2e1 * t22 * t34, t57, t51, 0, 0, 0, -0.2e1 * t31 * t125 + 0.2e1 * t136, -0.2e1 * t31 * t124 - 0.2e1 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, t115, 0, 0, 0, 0, 0, 0, 0, t86, t13 * t41 + t14 * t42 - t3 * t90 + t4 * t50, 0, 0, 0, 0, 0, -t97 * t124 + t86 * t74, t97 * t125 + t86 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t126, 0, -t21 * t90 + t22 * t50 - t33 * t41 + t34 * t42, 0, 0, 0, 0, 0, t26, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t41 * t90 + 0.2e1 * t42 * t50, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t87, 0, -t19, t18 (-t128 * t17 - t16 * t73) * pkin(4) (-t128 * t3 + t4 * t73) * pkin(4), -t7, -t6, t9, -t8, 0, -t139 + t100 * t74 + (t13 * t74 - t96 * t77) * qJD(6), t3 * t74 + t100 * t77 + (t13 * t77 + t96 * t74) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t38, 0 (-t128 * t21 + t22 * t73) * pkin(4), -0.2e1 * t110, t149, 0, 0, 0, t74 * t108 - t136, t77 * t108 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t126, 0 (-t128 * t41 + t42 * t73) * pkin(4), 0, 0, 0, 0, 0, -t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t51, 0, 0, 0, t74 * t118, t77 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t92, t16, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t125, 0, -t32 * t124 - t22 * t74, t32 * t125 - t22 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t124 - t42 * t74, t50 * t125 - t42 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t125, 0, -t62 * t124, t62 * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
