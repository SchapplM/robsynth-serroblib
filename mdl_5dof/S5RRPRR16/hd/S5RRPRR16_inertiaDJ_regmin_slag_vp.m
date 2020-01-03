% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR16_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:26
% EndTime: 2019-12-31 20:47:31
% DurationCPUTime: 1.47s
% Computational Cost: add. (1162->207), mult. (3309->403), div. (0->0), fcn. (2969->8), ass. (0->129)
t78 = cos(qJ(5));
t76 = sin(qJ(4));
t79 = cos(qJ(4));
t95 = pkin(4) * t79 + pkin(9) * t76;
t151 = t78 * t95;
t80 = cos(qJ(2));
t100 = -pkin(1) * t80 - pkin(2);
t73 = sin(pkin(5));
t77 = sin(qJ(2));
t142 = t73 * t77;
t62 = pkin(7) * t142;
t74 = cos(pkin(5));
t23 = pkin(3) * t142 + t62 + (-pkin(8) + t100) * t74;
t133 = qJ(3) * t77;
t81 = -pkin(2) - pkin(8);
t31 = (t81 * t80 - pkin(1) - t133) * t73;
t150 = t76 * t23 + t79 * t31;
t141 = t73 * t80;
t145 = pkin(1) * t74;
t65 = t77 * t145;
t149 = pkin(7) * t141 + t65;
t71 = t78 ^ 2;
t75 = sin(qJ(5));
t136 = t75 ^ 2 - t71;
t99 = t136 * qJD(5);
t130 = qJD(3) * t77;
t132 = qJD(2) * t77;
t105 = t73 * t132;
t59 = pkin(2) * t105;
t21 = t59 + (-t130 + (pkin(8) * t77 - qJ(3) * t80) * qJD(2)) * t73;
t146 = pkin(3) + pkin(7);
t32 = (t146 * t141 + t65) * qJD(2);
t8 = -qJD(4) * t150 - t76 * t21 + t79 * t32;
t148 = 0.2e1 * t73;
t147 = 0.2e1 * qJD(3);
t45 = t79 * t141 + t74 * t76;
t127 = qJD(4) * t45;
t26 = t76 * t105 - t127;
t131 = qJD(2) * t80;
t61 = t73 * t131;
t115 = t76 * t141;
t46 = t74 * t79 - t115;
t88 = t78 * t142 - t46 * t75;
t11 = t88 * qJD(5) + t26 * t78 + t75 * t61;
t144 = t11 * t75;
t43 = t149 * qJD(2);
t143 = t43 * t74;
t140 = t76 * t81;
t139 = t78 * t79;
t138 = t79 * t81;
t70 = t76 ^ 2;
t72 = t79 ^ 2;
t135 = t70 - t72;
t134 = t70 + t72;
t129 = qJD(4) * t88;
t29 = t75 * t142 + t46 * t78;
t128 = qJD(4) * t29;
t126 = qJD(4) * t76;
t125 = qJD(4) * t78;
t124 = qJD(4) * t79;
t123 = qJD(4) * t81;
t122 = qJD(5) * t75;
t121 = qJD(5) * t78;
t120 = qJD(5) * t79;
t119 = qJD(5) * t81;
t118 = qJ(3) * qJD(4);
t117 = -0.2e1 * pkin(4) * qJD(5);
t116 = t81 * t142;
t114 = t75 * t140;
t113 = t75 * t138;
t112 = t78 * t140;
t39 = -t74 * qJ(3) - t149;
t68 = t73 ^ 2;
t111 = t68 * t131;
t110 = t75 * t127;
t109 = t45 * t125;
t108 = t75 * t120;
t107 = t75 * t119;
t106 = t78 * t120;
t104 = t75 * t121;
t103 = t76 * t125;
t102 = t78 * t124;
t101 = t76 * t124;
t98 = t135 * qJD(4);
t30 = pkin(3) * t141 - t39;
t97 = t81 * t61;
t96 = t75 * t103;
t94 = t76 * pkin(4) - t79 * pkin(9);
t93 = -pkin(2) * t80 - t133;
t13 = pkin(9) * t142 + t150;
t15 = t45 * pkin(4) - t46 * pkin(9) + t30;
t6 = t78 * t13 + t75 * t15;
t91 = t79 * t23 - t76 * t31;
t89 = t29 * t75 - t78 * t88;
t60 = t131 * t145;
t42 = pkin(7) * t105 - t60;
t53 = qJ(3) + t94;
t36 = t75 * t53 + t112;
t12 = -pkin(4) * t142 - t91;
t4 = -pkin(4) * t61 - t8;
t87 = t12 * t121 + t4 * t75;
t86 = t12 * t122 - t4 * t78;
t27 = -qJD(4) * t115 - t79 * t105 + t74 * t124;
t85 = t45 * t121 + t75 * t27;
t84 = t45 * t122 - t78 * t27;
t7 = -t23 * t124 + t31 * t126 - t79 * t21 - t76 * t32;
t48 = -t103 - t108;
t66 = t74 * qJD(3);
t22 = -t146 * t105 + t60 + t66;
t83 = pkin(9) * t61 - t7;
t82 = t27 * pkin(4) - t26 * pkin(9) + t22;
t52 = 0.2e1 * t77 * t111;
t50 = t75 * t126 - t106;
t49 = t76 * t121 + t75 * t124;
t47 = t76 * t122 - t102;
t41 = t100 * t74 + t62;
t40 = (-pkin(1) + t93) * t73;
t38 = (-t77 * t126 + t79 * t131) * t73;
t37 = (-t77 * t124 - t76 * t131) * t73;
t35 = t78 * t53 - t114;
t34 = t42 - t66;
t33 = t59 + (-qJ(3) * t131 - t130) * t73;
t17 = t78 * qJD(3) - t36 * qJD(5) + (-t113 + t151) * qJD(4);
t16 = t76 * t107 - t75 * (t95 * qJD(4) + qJD(3)) - t53 * t121 - t81 * t102;
t10 = t29 * qJD(5) + t26 * t75 - t78 * t61;
t5 = -t75 * t13 + t78 * t15;
t2 = -t6 * qJD(5) - t75 * t83 + t78 * t82;
t1 = -t15 * t121 + t13 * t122 - t75 * t82 - t78 * t83;
t3 = [0, 0, 0, t52, 0.2e1 * (-t77 ^ 2 + t80 ^ 2) * t68 * qJD(2), 0.2e1 * t74 * t61, -0.2e1 * t74 * t105, 0, -0.2e1 * t68 * pkin(1) * t132 - 0.2e1 * t143, -0.2e1 * pkin(1) * t111 + 0.2e1 * t42 * t74, (-t34 * t80 + t43 * t77 + (t39 * t77 + t41 * t80) * qJD(2)) * t148, 0.2e1 * t143 + 0.2e1 * (-t40 * t132 + t33 * t80) * t73, -0.2e1 * t34 * t74 + 0.2e1 * (-t40 * t131 - t33 * t77) * t73, 0.2e1 * t40 * t33 + 0.2e1 * t39 * t34 + 0.2e1 * t41 * t43, 0.2e1 * t46 * t26, -0.2e1 * t26 * t45 - 0.2e1 * t46 * t27, (t46 * t131 + t26 * t77) * t148, (-t45 * t131 - t27 * t77) * t148, t52, 0.2e1 * t22 * t45 + 0.2e1 * t30 * t27 + 0.2e1 * (t91 * t131 + t8 * t77) * t73, 0.2e1 * t22 * t46 + 0.2e1 * t30 * t26 + 0.2e1 * (-t131 * t150 + t7 * t77) * t73, 0.2e1 * t29 * t11, -0.2e1 * t29 * t10 + 0.2e1 * t11 * t88, 0.2e1 * t11 * t45 + 0.2e1 * t29 * t27, -0.2e1 * t10 * t45 + 0.2e1 * t27 * t88, 0.2e1 * t45 * t27, 0.2e1 * t12 * t10 + 0.2e1 * t2 * t45 + 0.2e1 * t5 * t27 - 0.2e1 * t4 * t88, 0.2e1 * t1 * t45 + 0.2e1 * t12 * t11 - 0.2e1 * t6 * t27 + 0.2e1 * t4 * t29; 0, 0, 0, 0, 0, t61, -t105, 0, -t43, t42, (t93 * qJD(2) + qJD(3) * t80) * t73, t43, -t42 + 0.2e1 * t66, -t43 * pkin(2) - t34 * qJ(3) - t39 * qJD(3), -t46 * t126 + t26 * t79, -t26 * t76 - t79 * t27 + (t45 * t76 - t46 * t79) * qJD(4), t38, t37, 0, t79 * t97 + qJ(3) * t27 + qJD(3) * t45 + t22 * t76 + (-t76 * t116 + t30 * t79) * qJD(4), -t76 * t97 + qJ(3) * t26 + qJD(3) * t46 + t22 * t79 + (-t79 * t116 - t30 * t76) * qJD(4), t11 * t139 + t29 * t48, t89 * t126 + (-t10 * t78 - t144 + (-t29 * t78 - t75 * t88) * qJD(5)) * t79, (t11 - t109) * t76 + (-t84 + t128) * t79, (-t10 + t110) * t76 + (-t85 + t129) * t79, t45 * t124 + t27 * t76, t17 * t45 + t35 * t27 + (t2 + (-t12 * t75 - t81 * t88) * qJD(4)) * t76 + (qJD(4) * t5 - t10 * t81 + t87) * t79, t16 * t45 - t36 * t27 + (t1 + (-t12 * t78 + t29 * t81) * qJD(4)) * t76 + (-qJD(4) * t6 - t11 * t81 - t86) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, qJ(3) * t147, -0.2e1 * t101, 0.2e1 * t98, 0, 0, 0, 0.2e1 * qJD(3) * t76 + 0.2e1 * t79 * t118, 0.2e1 * qJD(3) * t79 - 0.2e1 * t76 * t118, -0.2e1 * t71 * t101 - 0.2e1 * t72 * t104, 0.2e1 * t72 * t99 + 0.4e1 * t79 * t96, -0.2e1 * t76 * t108 - 0.2e1 * t135 * t125, -0.2e1 * t76 * t106 + 0.2e1 * t75 * t98, 0.2e1 * t101, -0.2e1 * t72 * t78 * t119 + 0.2e1 * t17 * t76 + 0.2e1 * (t35 + 0.2e1 * t114) * t124, 0.2e1 * t72 * t107 + 0.2e1 * t16 * t76 + 0.2e1 * (-t36 + 0.2e1 * t112) * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, t43, 0, 0, 0, 0, 0, t38, t37, 0, 0, 0, 0, 0, (-t10 - t110) * t79 + (-t85 - t129) * t76, (-t11 - t109) * t79 + (t84 + t128) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t121, t134 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, t61, t8, t7, t29 * t121 + t144, -qJD(5) * t89 - t75 * t10 + t11 * t78, t85, -t84, 0, -pkin(4) * t10 - pkin(9) * t85 + t86, -pkin(4) * t11 + pkin(9) * t84 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t124, 0, -t76 * t123, -t79 * t123, -t79 * t99 - t96, -0.4e1 * t79 * t104 + t136 * t126, t49, -t47, 0, (-t113 - t151) * qJD(5) + (t75 * t94 - t112) * qJD(4), (-t78 * t138 + t75 * t95) * qJD(5) + (-pkin(9) * t139 + (t78 * pkin(4) + t75 * t81) * t76) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t124, 0, 0, 0, 0, 0, t48, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t104, -0.2e1 * t99, 0, 0, 0, t75 * t117, t78 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t27, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t50, t124, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t122, 0, -pkin(9) * t121, pkin(9) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
