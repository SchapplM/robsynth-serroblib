% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP6
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
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:49
% EndTime: 2019-03-09 06:16:54
% DurationCPUTime: 1.72s
% Computational Cost: add. (3373->214), mult. (7662->376), div. (0->0), fcn. (7669->8), ass. (0->122)
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t140 = qJD(4) * t108;
t103 = sin(pkin(10));
t109 = cos(qJ(3));
t104 = cos(pkin(10));
t107 = sin(qJ(3));
t145 = t107 * t104;
t77 = t109 * t103 + t145;
t132 = t77 * t140;
t76 = t107 * t103 - t109 * t104;
t65 = t76 * qJD(3);
t168 = -t106 * t65 + t132;
t160 = cos(qJ(5));
t125 = t160 * qJD(5);
t167 = t160 * qJD(4) + t125;
t129 = t160 * t108;
t105 = sin(qJ(5));
t147 = t105 * t106;
t78 = -t129 + t147;
t102 = t108 ^ 2;
t144 = t106 ^ 2 - t102;
t120 = t144 * qJD(4);
t166 = qJD(4) + qJD(5);
t95 = -t104 * pkin(2) - pkin(1);
t165 = 0.2e1 * t95;
t164 = pkin(8) + pkin(9);
t163 = pkin(9) * t77;
t130 = t160 * t106;
t79 = t105 * t108 + t130;
t54 = t166 * t79;
t162 = t54 * pkin(5);
t161 = t76 * pkin(4);
t159 = t77 * t65;
t53 = -t167 * t108 + t166 * t147;
t158 = t79 * t53;
t157 = pkin(7) + qJ(2);
t50 = t76 * pkin(3) - t77 * pkin(8) + t95;
t46 = t108 * t50;
t85 = t157 * t103;
t149 = t107 * t85;
t86 = t157 * t104;
t56 = t109 * t86 - t149;
t23 = -t106 * t56 - t108 * t163 + t161 + t46;
t150 = t106 * t77;
t51 = t108 * t56;
t153 = t106 * t50 + t51;
t25 = -pkin(9) * t150 + t153;
t24 = t160 * t25;
t156 = t105 * t23 + t24;
t141 = qJD(4) * t106;
t133 = t77 * t141;
t134 = t77 * t147;
t17 = -t65 * t130 - t105 * t133 - qJD(5) * t134 + (-t105 * t65 + t167 * t77) * t108;
t43 = t79 * t77;
t155 = -t79 * t17 + t53 * t43;
t87 = t164 * t106;
t88 = t164 * t108;
t152 = -t105 * t87 + t160 * t88;
t148 = t108 * t65;
t146 = t106 * t108;
t143 = qJD(2) * t109;
t142 = qJD(3) * t109;
t139 = qJD(5) * t105;
t138 = -0.2e1 * pkin(3) * qJD(4);
t137 = t160 * pkin(4);
t136 = pkin(4) * t141;
t135 = pkin(4) * t139;
t98 = -t108 * pkin(4) - pkin(3);
t131 = qJD(4) * t164;
t128 = t106 * t140;
t33 = t85 * t142 - t104 * t143 + (qJD(2) * t103 + qJD(3) * t86) * t107;
t66 = t77 * qJD(3);
t49 = t66 * pkin(3) + t65 * pkin(8);
t13 = -t106 * t49 + t108 * t33 - t50 * t140 + t56 * t141;
t12 = -pkin(9) * t168 - t13;
t122 = t106 * t33 + t108 * t49;
t9 = pkin(9) * t148 + t66 * pkin(4) + (-t51 + (-t50 + t163) * t106) * qJD(4) + t122;
t127 = -t105 * t12 + t160 * t9;
t124 = -t105 * t25 + t160 * t23;
t123 = -t105 * t88 - t160 * t87;
t55 = t107 * t86 + t109 * t85;
t121 = -0.4e1 * t77 * t146;
t34 = qJD(2) * t145 - qJD(3) * t149 + t103 * t143 + t86 * t142;
t119 = pkin(4) * t125;
t118 = 0.2e1 * (t103 ^ 2 + t104 ^ 2) * qJD(2);
t35 = pkin(4) * t150 + t55;
t117 = pkin(3) * t65 - pkin(8) * t66;
t116 = pkin(3) * t77 + pkin(8) * t76;
t16 = t54 * t77 - t78 * t65;
t44 = t77 * t129 - t134;
t115 = -t78 * t16 + t44 * t54;
t114 = t34 * t77 - t55 * t65;
t113 = t53 * t76 - t79 * t66;
t112 = -t65 * t76 + t77 * t66;
t27 = t168 * pkin(4) + t34;
t110 = t106 * t66 + t76 * t140;
t3 = -t105 * t9 - t160 * t12 - t23 * t125 + t25 * t139;
t83 = t106 * t131;
t84 = t108 * t131;
t30 = t105 * t84 + t87 * t125 + t88 * t139 + t160 * t83;
t4 = -t156 * qJD(5) + t127;
t31 = -t152 * qJD(5) + t105 * t83 - t160 * t84;
t97 = t137 + pkin(5);
t73 = t77 ^ 2;
t61 = t78 * pkin(5) + t98;
t52 = 0.2e1 * t76 * t66;
t48 = t136 + t162;
t47 = t108 * t66 - t76 * t141;
t40 = -t78 * qJ(6) + t152;
t39 = -t79 * qJ(6) + t123;
t28 = -t54 * t76 - t78 * t66;
t26 = t43 * pkin(5) + t35;
t22 = t53 * qJ(6) - t79 * qJD(6) + t31;
t21 = -t54 * qJ(6) - t78 * qJD(6) - t30;
t14 = -t153 * qJD(4) + t122;
t10 = t17 * pkin(5) + t27;
t6 = -t43 * qJ(6) + t156;
t5 = t76 * pkin(5) - t44 * qJ(6) + t124;
t2 = -t17 * qJ(6) - t43 * qJD(6) - t3;
t1 = t66 * pkin(5) + t16 * qJ(6) - t44 * qJD(6) + t4;
t7 = [0, 0, 0, 0, 0, t118, qJ(2) * t118, -0.2e1 * t159, -0.2e1 * t112, 0, 0, 0, t66 * t165, -t65 * t165, -0.2e1 * t102 * t159 - 0.2e1 * t73 * t128, 0.2e1 * t73 * t120 - t65 * t121, 0.2e1 * t112 * t108 - 0.2e1 * t76 * t133, -0.2e1 * t112 * t106 - 0.2e1 * t76 * t132, t52, 0.2e1 * t55 * t132 + 0.2e1 * t14 * t76 + 0.2e1 * t46 * t66 + 0.2e1 * (-t56 * t66 + t114) * t106, 0.2e1 * t114 * t108 + 0.2e1 * t13 * t76 - 0.2e1 * t55 * t133 - 0.2e1 * t153 * t66, -0.2e1 * t44 * t16, 0.2e1 * t16 * t43 - 0.2e1 * t44 * t17, -0.2e1 * t16 * t76 + 0.2e1 * t44 * t66, -0.2e1 * t17 * t76 - 0.2e1 * t43 * t66, t52, 0.2e1 * t124 * t66 + 0.2e1 * t35 * t17 + 0.2e1 * t27 * t43 + 0.2e1 * t4 * t76, -0.2e1 * t156 * t66 - 0.2e1 * t35 * t16 + 0.2e1 * t27 * t44 + 0.2e1 * t3 * t76, -0.2e1 * t1 * t44 + 0.2e1 * t5 * t16 - 0.2e1 * t6 * t17 - 0.2e1 * t2 * t43, 0.2e1 * t5 * t1 + 0.2e1 * t26 * t10 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t65, 0, 0, 0, 0, 0, t47, -t110, 0, 0, 0, 0, 0, t28, t113, t115 + t155, -t1 * t78 + t2 * t79 - t5 * t54 - t6 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78 * t54 - 0.2e1 * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t66, 0, -t34, t33, -t77 * t120 - t65 * t146, qJD(4) * t121 + t144 * t65, t110, t47, 0, -t34 * t108 + t117 * t106 + (t106 * t55 - t116 * t108) * qJD(4), t34 * t106 + t117 * t108 + (t116 * t106 + t108 * t55) * qJD(4), -t16 * t79 - t44 * t53, -t115 + t155, -t113, t28, 0, t123 * t66 + t43 * t136 + t98 * t17 + t27 * t78 + t31 * t76 + t35 * t54, t44 * t136 - t152 * t66 - t98 * t16 + t27 * t79 + t30 * t76 - t35 * t53, -t1 * t79 + t39 * t16 - t40 * t17 - t2 * t78 - t21 * t43 - t22 * t44 + t5 * t53 - t6 * t54, t1 * t39 + t10 * t61 + t2 * t40 + t6 * t21 + t5 * t22 + t26 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t21 - t78 * t22 - t54 * t39 - t53 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t128, -0.2e1 * t120, 0, 0, 0, t106 * t138, t108 * t138, -0.2e1 * t158, 0.2e1 * t78 * t53 - 0.2e1 * t79 * t54, 0, 0, 0, 0.2e1 * t78 * t136 + 0.2e1 * t98 * t54, 0.2e1 * t79 * t136 - 0.2e1 * t98 * t53, -0.2e1 * t21 * t78 - 0.2e1 * t22 * t79 + 0.2e1 * t39 * t53 - 0.2e1 * t40 * t54, 0.2e1 * t40 * t21 + 0.2e1 * t39 * t22 + 0.2e1 * t61 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 - t148, -t168, t66, t14, t13, 0, 0, -t16, -t17, t66, t66 * t137 + (-t24 + (-t23 - t161) * t105) * qJD(5) + t127 (-t105 * t66 - t76 * t125) * pkin(4) + t3, t97 * t16 + (-t105 * t17 + (t105 * t44 - t160 * t43) * qJD(5)) * pkin(4), t1 * t97 + (t105 * t2 + (-t105 * t5 + t160 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t140, 0, 0, 0, 0, 0, -t54, t53, 0, -t54 * t97 + (-t105 * t53 + (t105 * t78 + t160 * t79) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t141, 0, -pkin(8) * t140, pkin(8) * t141, 0, 0, -t53, -t54, 0, t31, t30, t97 * t53 + (-t105 * t54 + (t105 * t79 - t160 * t78) * qJD(5)) * pkin(4), t22 * t97 + (t105 * t21 + (-t105 * t39 + t160 * t40) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t135, -0.2e1 * t119, 0, 0.2e1 * (t137 - t97) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, t66, t4, t3, pkin(5) * t16, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t53, 0, -t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t54, 0, t31, t30, pkin(5) * t53, t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t119, 0, -pkin(5) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
