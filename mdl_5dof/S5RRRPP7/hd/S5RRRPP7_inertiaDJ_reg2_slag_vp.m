% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:06:00
% DurationCPUTime: 1.82s
% Computational Cost: add. (959->209), mult. (2529->343), div. (0->0), fcn. (1748->4), ass. (0->127)
t75 = cos(qJ(3));
t126 = t75 * qJD(4);
t73 = sin(qJ(3));
t134 = t73 * qJ(4);
t146 = pkin(3) + pkin(4);
t152 = t146 * t75 + t134;
t154 = t152 * qJD(3) - t126;
t76 = cos(qJ(2));
t129 = qJD(3) * t76;
t113 = t75 * t129;
t74 = sin(qJ(2));
t65 = t74 * qJD(2);
t34 = -t73 * t65 + t113;
t153 = t146 * qJD(2);
t69 = t73 ^ 2;
t71 = t75 ^ 2;
t138 = t69 - t71;
t151 = t138 * qJD(3);
t124 = t76 * qJD(2);
t105 = qJ(5) * t124;
t66 = qJD(3) * t75;
t114 = t74 * t66;
t127 = t73 * qJD(5);
t150 = qJ(5) * t114 + t73 * t105 + t74 * t127;
t115 = t73 * t129;
t32 = t75 * t65 + t115;
t141 = t74 * pkin(7);
t99 = -t76 * pkin(2) - t141;
t44 = -pkin(1) + t99;
t144 = pkin(7) * t76;
t98 = pkin(2) * t74 - t144;
t88 = t98 * qJD(2);
t7 = t32 * pkin(6) - t44 * t66 - t73 * t88;
t136 = qJ(4) * t75;
t95 = pkin(3) * t73 - t136;
t91 = pkin(6) + t95;
t27 = t91 * t74;
t128 = t73 * qJD(4);
t29 = t95 * qJD(3) - t128;
t96 = t75 * pkin(3) + t134;
t42 = -pkin(2) - t96;
t149 = qJD(2) * (-t42 * t76 + t141) - qJD(3) * t27 - t29 * t74;
t147 = t96 * qJD(3) - t126;
t108 = t75 * t124;
t101 = t73 * t108;
t70 = t74 ^ 2;
t11 = -0.4e1 * t74 * t101 + 0.2e1 * t70 * t151;
t41 = -0.2e1 * t151;
t78 = 0.2e1 * qJD(4);
t145 = pkin(6) * t76;
t89 = -t146 * t73 + t136;
t85 = -pkin(6) + t89;
t3 = t85 * t124 - t154 * t74;
t143 = t3 * t73;
t142 = t3 * t75;
t140 = pkin(7) - qJ(5);
t61 = t75 * t145;
t26 = t73 * t44 + t61;
t137 = -t76 ^ 2 + t70;
t135 = qJ(5) * t74;
t133 = qJD(2) * t73;
t132 = qJD(2) * t75;
t131 = qJD(3) * t73;
t130 = qJD(3) * t74;
t125 = t75 * qJD(5);
t123 = t76 * qJD(4);
t122 = qJ(4) * qJD(2);
t60 = t73 * t145;
t121 = -0.2e1 * pkin(1) * qJD(2);
t120 = -0.2e1 * pkin(2) * qJD(3);
t119 = pkin(3) * t65;
t118 = pkin(7) * t131;
t117 = pkin(7) * t66;
t116 = pkin(6) * t124;
t39 = pkin(2) + t152;
t112 = t39 * t66;
t110 = t73 * t66;
t109 = t74 * t124;
t107 = t76 * t122;
t106 = qJ(4) * t130;
t47 = t140 * t75;
t25 = t75 * t44 - t60;
t104 = t137 * qJD(2);
t103 = 0.2e1 * t109;
t102 = t34 * pkin(6) + t44 * t131 - t75 * t88;
t100 = t70 * t110;
t19 = -t76 * qJ(4) + t26;
t67 = t76 * pkin(3);
t21 = -t25 + t67;
t94 = -t19 * t73 + t21 * t75;
t93 = -t25 * t75 - t26 * t73;
t24 = t89 * qJD(3) + t128;
t90 = -t39 * t131 + t24 * t75;
t15 = t85 * t74;
t87 = qJD(3) * t15 + t39 * t124;
t86 = -qJ(5) * t131 + t125;
t84 = -t75 * t105 + t102;
t63 = t74 * t122;
t4 = t63 - t7 - t123;
t5 = t102 - t119;
t81 = t94 * qJD(3) + t4 * t75 + t5 * t73;
t80 = t93 * qJD(3) + t102 * t73 - t7 * t75;
t79 = 0.2e1 * t63 - 0.2e1 * t123 - t7;
t68 = qJ(4) * t78;
t55 = -0.2e1 * t109;
t54 = -0.2e1 * t110;
t53 = 0.2e1 * t110;
t52 = pkin(7) * t113;
t46 = t140 * t73;
t33 = t73 * t124 + t114;
t31 = t73 * t130 - t108;
t30 = qJD(3) * t47 - t127;
t28 = -t140 * t131 - t125;
t23 = 0.2e1 * t71 * t109 - 0.2e1 * t100;
t22 = 0.2e1 * t69 * t109 + 0.2e1 * t100;
t20 = t138 * t130 - t101;
t18 = -t73 * t104 + t74 * t113;
t17 = t74 * t115 + t137 * t132;
t16 = 0.4e1 * t74 * t110 + t138 * t124;
t14 = 0.2e1 * t18;
t13 = 0.2e1 * t17;
t10 = t73 * t135 + t19;
t9 = t76 * pkin(4) + t60 + t67 + (-t44 - t135) * t75;
t6 = t91 * t124 + t147 * t74;
t2 = t4 + t150;
t1 = (-t86 - t153) * t74 + t84;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, -0.2e1 * t104, 0, t55, 0, 0, t74 * t121, t76 * t121, 0, 0, t23, t11, t13, t22, t14, t55, 0.2e1 * t25 * t65 + 0.2e1 * t102 * t76 + 0.2e1 * (t103 * t73 + t70 * t66) * pkin(6), -0.2e1 * t26 * t65 - 0.2e1 * t7 * t76 + 0.2e1 * (t103 * t75 - t70 * t131) * pkin(6), 0.2e1 * t93 * t124 + 0.2e1 * (t7 * t73 + t75 * t102 + (t25 * t73 - t26 * t75) * qJD(3)) * t74, 0.2e1 * pkin(6) ^ 2 * t109 - 0.2e1 * t102 * t25 - 0.2e1 * t26 * t7, t23, t13, -t11, t55, -0.2e1 * t18, t22, 0.2e1 * (t133 * t27 + t5) * t76 + 0.2e1 * (-qJD(2) * t21 + t27 * t66 + t6 * t73) * t74, 0.2e1 * t94 * t124 + 0.2e1 * (-t4 * t73 + t5 * t75 + (-t19 * t75 - t21 * t73) * qJD(3)) * t74, 0.2e1 * (-t132 * t27 - t4) * t76 + 0.2e1 * (qJD(2) * t19 + t131 * t27 - t6 * t75) * t74, 0.2e1 * t19 * t4 + 0.2e1 * t21 * t5 + 0.2e1 * t27 * t6, t23, -t11, -0.2e1 * t17, t22, t14, t55, 0.2e1 * (-t133 * t15 + t1) * t76 + 0.2e1 * (-qJD(2) * t9 - t15 * t66 - t143) * t74, 0.2e1 * (t132 * t15 - t2) * t76 + 0.2e1 * (qJD(2) * t10 - t131 * t15 + t142) * t74, 0.2e1 * (t10 * t73 - t75 * t9) * t124 + 0.2e1 * (-t1 * t75 + t2 * t73 + (t10 * t75 + t73 * t9) * qJD(3)) * t74, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t15 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, -t65, 0, -t116, pkin(6) * t65, 0, 0, -t20, -t16, -t34, t20, t32, 0, t52 + (-pkin(2) * t75 + pkin(6) * t73) * t130 + (t73 * t99 - t61) * qJD(2), (pkin(6) * t74 * t75 + t73 * t98) * qJD(3) + (t75 * t99 + t60) * qJD(2), t80, -pkin(2) * t116 + t80 * pkin(7), -t20, -t34, t16, 0, -t32, t20, t52 + (t130 * t42 - t6) * t75 - t149 * t73, t81, (-t6 + (t42 * t74 + t144) * qJD(3)) * t73 + t149 * t75, t81 * pkin(7) + t27 * t29 + t6 * t42, -t20, t16, t34, t20, t32, 0, t142 + t30 * t76 + (-qJD(2) * t46 - t112) * t74 + (-t24 * t74 - t87) * t73, -t28 * t76 + t143 + t87 * t75 + (qJD(2) * t47 + t90) * t74, (-t46 * t124 - t30 * t74 - t2 + (t47 * t74 - t9) * qJD(3)) * t75 + (t47 * t124 + t28 * t74 - t1 + (t46 * t74 + t10) * qJD(3)) * t73, t1 * t46 + t10 * t28 + t15 * t24 + t2 * t47 + t3 * t39 + t9 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t41, 0, t54, 0, 0, t73 * t120, t75 * t120, 0, 0, t53, 0, -t41, 0, 0, t54, 0.2e1 * t131 * t42 - 0.2e1 * t29 * t75, 0, -0.2e1 * t29 * t73 - 0.2e1 * t42 * t66, 0.2e1 * t42 * t29, t53, -t41, 0, t54, 0, 0, 0.2e1 * t90, 0.2e1 * t24 * t73 + 0.2e1 * t112, -0.2e1 * t28 * t75 - 0.2e1 * t30 * t73 + 0.2e1 * (-t46 * t75 + t47 * t73) * qJD(3), 0.2e1 * t39 * t24 + 0.2e1 * t47 * t28 + 0.2e1 * t46 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t33, t65, -t102, t7, 0, 0, 0, -t31, 0, t65, t33, 0, -t102 + 0.2e1 * t119, (-pkin(3) * t124 - t106) * t75 + (-t107 + (pkin(3) * qJD(3) - qJD(4)) * t74) * t73, t79, -t5 * pkin(3) + t4 * qJ(4) + t19 * qJD(4), 0, 0, t31, 0, -t33, t65, (t86 + 0.2e1 * t153) * t74 - t84, t79 + t150, (t124 * t146 + t106) * t75 + (t107 + (-qJD(3) * t146 + qJD(4)) * t74) * t73, t2 * qJ(4) + t10 * qJD(4) - t1 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t131, 0, -t117, t118, 0, 0, 0, t66, 0, 0, t131, 0, -t117, -t147, -t118, -t147 * pkin(7), 0, 0, -t66, 0, -t131, 0, -t30, t28, t154, t28 * qJ(4) + t47 * qJD(4) - t146 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t68, 0, 0, 0, 0, 0, 0, 0, t78, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t31, 0, t5, 0, 0, 0, 0, 0, 0, -t65, 0, t31, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t31, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t66, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
