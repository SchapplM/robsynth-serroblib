% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:11
% EndTime: 2019-03-08 22:34:16
% DurationCPUTime: 1.45s
% Computational Cost: add. (1007->191), mult. (2632->350), div. (0->0), fcn. (2353->10), ass. (0->137)
t152 = pkin(4) + pkin(8);
t77 = sin(qJ(3));
t137 = t77 * qJ(4);
t81 = cos(qJ(3));
t83 = -pkin(3) - pkin(9);
t93 = -t81 * t83 + t137;
t41 = -pkin(2) - t93;
t56 = t152 * t77;
t76 = sin(qJ(5));
t47 = t76 * t56;
t80 = cos(qJ(5));
t141 = t80 * t41 + t47;
t72 = t81 ^ 2;
t105 = qJD(3) * (t77 ^ 2 - t72);
t69 = t76 ^ 2;
t140 = -t80 ^ 2 + t69;
t104 = t140 * qJD(5);
t156 = qJD(5) + qJD(6);
t57 = t152 * t81;
t126 = t57 * qJD(5);
t132 = qJD(4) * t81;
t155 = qJD(3) * t93 - t126 - t132;
t154 = t156 * t81;
t153 = 0.2e1 * qJD(4);
t151 = t77 * pkin(5);
t130 = qJD(5) * t80;
t131 = qJD(5) * t76;
t125 = t77 * qJD(3);
t103 = pkin(3) * t125 - t77 * qJD(4);
t138 = qJ(4) * t81;
t27 = (pkin(9) * t77 - t138) * qJD(3) + t103;
t66 = t81 * qJD(3);
t64 = pkin(8) * t66;
t50 = pkin(4) * t66 + t64;
t8 = -t56 * t130 + t41 * t131 - t80 * t27 - t76 * t50;
t111 = t80 * t125;
t129 = qJD(5) * t81;
t117 = t76 * t129;
t87 = t111 + t117;
t7 = t87 * pkin(10) - t8;
t79 = cos(qJ(6));
t150 = t79 * t7;
t149 = pkin(10) - t83;
t73 = sin(pkin(6));
t78 = sin(qJ(2));
t148 = t73 * t78;
t82 = cos(qJ(2));
t147 = t73 * t82;
t142 = t80 * t81;
t18 = -pkin(10) * t142 + t141;
t75 = sin(qJ(6));
t146 = t75 * t18;
t145 = t75 * t76;
t144 = t79 * t18;
t143 = t79 * t80;
t136 = qJD(2) * t78;
t135 = qJD(2) * t82;
t134 = qJD(3) * t76;
t133 = qJD(3) * t80;
t128 = qJD(5) * t83;
t127 = qJD(6) * t75;
t124 = qJ(4) * qJD(5);
t123 = -0.2e1 * pkin(2) * qJD(3);
t122 = t77 * t148;
t121 = pkin(5) * t66;
t120 = pkin(5) * t127;
t119 = qJD(6) * t79 * pkin(5);
t118 = pkin(8) * t125;
t116 = t80 * t129;
t115 = t73 * t136;
t114 = t73 * t135;
t113 = t76 * t125;
t112 = t77 * t66;
t110 = t76 * t130;
t106 = -t76 * t27 + t80 * t50;
t107 = pkin(10) * t81 - t41;
t6 = (-pkin(10) * t76 * t77 + pkin(5) * t81) * qJD(3) + (t107 * t80 - t47) * qJD(5) + t106;
t109 = t79 * t6 - t75 * t7;
t52 = t149 * t80;
t48 = t80 * t56;
t17 = t107 * t76 + t151 + t48;
t108 = -t17 - t151;
t102 = t76 * t111;
t101 = -t81 * pkin(3) - t137;
t100 = t79 * t17 - t146;
t99 = t75 * t17 + t144;
t74 = cos(pkin(6));
t35 = -t74 * t81 + t122;
t23 = t76 * t147 + t35 * t80;
t91 = t80 * t147 - t35 * t76;
t98 = t79 * t23 + t75 * t91;
t97 = t75 * t23 - t79 * t91;
t51 = t149 * t76;
t96 = -t79 * t51 - t75 * t52;
t95 = -t75 * t51 + t79 * t52;
t45 = t75 * t80 + t79 * t76;
t94 = -t143 + t145;
t36 = t81 * t148 + t74 * t77;
t20 = -t76 * t127 - t75 * t131 + t156 * t143;
t90 = -t20 * t77 - t45 * t66;
t21 = -qJD(3) * t122 + (qJD(3) * t74 + t114) * t81;
t89 = t36 * t130 + t21 * t76;
t88 = -t36 * t131 + t21 * t80;
t49 = t152 * t125;
t86 = -t49 + (-t77 * t83 - t138) * qJD(5);
t85 = t101 * qJD(3) + t132;
t22 = qJD(3) * t36 + t77 * t114;
t84 = t21 * t81 + t22 * t77 + (t35 * t81 - t36 * t77) * qJD(3);
t63 = t76 * pkin(5) + qJ(4);
t60 = pkin(5) * t130 + qJD(4);
t59 = 0.2e1 * t112;
t53 = -pkin(2) + t101;
t44 = qJD(5) * t52;
t43 = t149 * t131;
t38 = -t77 * t130 - t76 * t66;
t37 = -t77 * t131 + t80 * t66;
t34 = pkin(5) * t142 + t57;
t33 = -qJ(4) * t66 + t103;
t31 = t45 * t81;
t30 = -t79 * t142 + t81 * t145;
t29 = (t82 * t125 + t81 * t136) * t73;
t28 = (t77 * t136 - t82 * t66) * t73;
t25 = -pkin(5) * t117 + (-pkin(5) * t80 - t152) * t125;
t19 = t156 * t45;
t16 = -t19 * t77 - t66 * t94;
t15 = t79 * t111 - t75 * t113 + t45 * t154;
t14 = t45 * t125 + t94 * t154;
t13 = -t96 * qJD(6) + t79 * t43 + t75 * t44;
t12 = t95 * qJD(6) - t75 * t43 + t79 * t44;
t11 = t23 * qJD(5) + t80 * t115 + t22 * t76;
t10 = t91 * qJD(5) - t76 * t115 + t22 * t80;
t9 = -t141 * qJD(5) + t106;
t4 = -t97 * qJD(6) + t79 * t10 - t75 * t11;
t3 = -t98 * qJD(6) - t75 * t10 - t79 * t11;
t2 = -t99 * qJD(6) + t109;
t1 = -t100 * qJD(6) - t75 * t6 - t150;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t73 ^ 2 * t78 * t135 + 0.2e1 * t36 * t21 + 0.2e1 * t35 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t115, -t114, 0, 0, 0, 0, 0, -t29, t28, t84, t29, -t28 (t53 * t136 - t33 * t82) * t73 + t84 * pkin(8), 0, 0, 0, 0, 0 (-t36 * t133 + t10) * t77 + (qJD(3) * t23 + t88) * t81 (t36 * t134 - t11) * t77 + (qJD(3) * t91 - t89) * t81, 0, 0, 0, 0, 0, -t36 * t15 - t21 * t30 + t4 * t77 + t98 * t66, t36 * t14 - t21 * t31 + t3 * t77 - t97 * t66; 0, 0, 0, 0, t59, -0.2e1 * t105, 0, 0, 0, t77 * t123, t81 * t123, 0, -0.2e1 * t53 * t125 + 0.2e1 * t33 * t81, -0.2e1 * t33 * t77 - 0.2e1 * t53 * t66, 0.2e1 * t53 * t33, 0.2e1 * t110 * t72 - 0.2e1 * t112 * t69, -0.4e1 * t81 * t102 - 0.2e1 * t72 * t104, 0.2e1 * t105 * t76 - 0.2e1 * t116 * t77, 0.2e1 * t105 * t80 + 0.2e1 * t117 * t77, t59, 0.2e1 * (-t57 * t133 + t9) * t77 + 0.2e1 * ((-t76 * t41 + t48) * qJD(3) - t49 * t80 - t76 * t126) * t81, 0.2e1 * (t57 * t134 + t8) * t77 + 0.2e1 * (-t141 * qJD(3) - t80 * t126 + t49 * t76) * t81, -0.2e1 * t31 * t14, 0.2e1 * t14 * t30 - 0.2e1 * t31 * t15, 0.2e1 * t14 * t77 - 0.2e1 * t31 * t66, 0.2e1 * t15 * t77 + 0.2e1 * t30 * t66, t59, 0.2e1 * t100 * t66 - 0.2e1 * t34 * t15 + 0.2e1 * t2 * t77 - 0.2e1 * t25 * t30, 0.2e1 * t1 * t77 + 0.2e1 * t34 * t14 - 0.2e1 * t25 * t31 - 0.2e1 * t99 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21, 0, t22, t21, -t22 * pkin(3) + t21 * qJ(4) + t36 * qJD(4), 0, 0, 0, 0, 0, t89, t88, 0, 0, 0, 0, 0, t36 * t20 + t21 * t45, -t36 * t19 - t21 * t94; 0, 0, 0, 0, 0, 0, t66, -t125, 0, -t64, t118, t85, t64, -t118, t85 * pkin(8), t104 * t81 + t102, 0.4e1 * t81 * t110 - t140 * t125, t37, t38, 0, -t155 * t80 + t86 * t76, t155 * t76 + t86 * t80, -t14 * t94 + t31 * t19, -t14 * t45 - t15 * t94 - t19 * t30 + t31 * t20, t16, t90, 0, t13 * t77 - t63 * t15 + t34 * t20 + t25 * t45 - t60 * t30 - t95 * t66, t12 * t77 + t63 * t14 - t34 * t19 - t25 * t94 - t60 * t31 - t96 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, qJ(4) * t153, -0.2e1 * t110, 0.2e1 * t104, 0, 0, 0, 0.2e1 * qJD(4) * t76 + 0.2e1 * t80 * t124, 0.2e1 * qJD(4) * t80 - 0.2e1 * t76 * t124, 0.2e1 * t94 * t19, 0.2e1 * t19 * t45 + 0.2e1 * t20 * t94, 0, 0, 0, 0.2e1 * t63 * t20 + 0.2e1 * t60 * t45, -0.2e1 * t63 * t19 - 0.2e1 * t60 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, t64, 0, 0, 0, 0, 0, t37, t38, 0, 0, 0, 0, 0, t16, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 - t116, t87, t66, t9, t8, 0, 0, t14, t15, t66, t79 * t121 + (t108 * t75 - t144) * qJD(6) + t109, -t150 + (-t6 - t121) * t75 + (t108 * t79 + t146) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t130, 0, -t76 * t128, -t80 * t128, 0, 0, -t19, -t20, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t130, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t120, -0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, t66, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
