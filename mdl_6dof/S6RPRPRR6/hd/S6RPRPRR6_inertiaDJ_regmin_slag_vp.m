% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:53
% EndTime: 2019-03-09 03:52:57
% DurationCPUTime: 1.34s
% Computational Cost: add. (3232->193), mult. (7591->351), div. (0->0), fcn. (7948->10), ass. (0->118)
t109 = cos(pkin(11));
t112 = sin(qJ(5));
t107 = sin(pkin(11));
t114 = cos(qJ(5));
t148 = t114 * t107;
t87 = t112 * t109 + t148;
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t113 = sin(qJ(3));
t149 = t113 * t110;
t161 = cos(qJ(3));
t88 = t161 * t108 + t149;
t55 = t87 * t88;
t144 = qJD(5) * t114;
t147 = t114 * t109;
t158 = pkin(8) + qJ(4);
t89 = t158 * t107;
t91 = t158 * t109;
t47 = (t107 * qJD(4) + qJD(5) * t91) * t112 - qJD(4) * t147 + t89 * t144;
t145 = qJD(5) * t112;
t136 = t107 * t145;
t75 = -t109 * t144 + t136;
t165 = -0.2e1 * t75;
t102 = -t110 * pkin(2) - pkin(1);
t164 = 0.2e1 * t102;
t77 = t88 * qJD(3);
t163 = t77 * pkin(5);
t117 = -t113 * t108 + t161 * t110;
t162 = t117 * pkin(5);
t160 = cos(qJ(6));
t159 = pkin(7) + qJ(2);
t59 = -pkin(3) * t117 - t88 * qJ(4) + t102;
t90 = t159 * t108;
t150 = t113 * t90;
t92 = t159 * t110;
t65 = t161 * t92 - t150;
t37 = -t107 * t65 + t109 * t59;
t26 = -t109 * t88 * pkin(8) - pkin(4) * t117 + t37;
t152 = t107 * t88;
t38 = t107 * t59 + t109 * t65;
t29 = -pkin(8) * t152 + t38;
t157 = t112 * t26 + t114 * t29;
t76 = t117 * qJD(3);
t45 = t77 * pkin(3) - t76 * qJ(4) - t88 * qJD(4);
t134 = qJD(3) * t161;
t135 = qJD(2) * t161;
t49 = t90 * t134 - t110 * t135 + (qJD(2) * t108 + qJD(3) * t92) * t113;
t25 = t107 * t45 - t109 * t49;
t155 = -t112 * t89 + t114 * t91;
t154 = t107 * t144 + t109 * t145;
t153 = t107 * t76;
t151 = t109 * t76;
t146 = t107 ^ 2 + t109 ^ 2;
t111 = sin(qJ(6));
t143 = qJD(6) * t111;
t142 = t160 * pkin(5);
t31 = t76 * t148 - t88 * t136 + (t112 * t76 + t88 * t144) * t109;
t24 = t107 * t49 + t109 * t45;
t17 = t77 * pkin(4) - pkin(8) * t151 + t24;
t20 = -pkin(8) * t153 + t25;
t6 = -t112 * t17 - t114 * t20 - t26 * t144 + t29 * t145;
t5 = -t31 * pkin(9) - t6;
t141 = t160 * t5;
t140 = pkin(5) * t143;
t139 = pkin(5) * t154;
t11 = -t55 * pkin(9) + t157;
t138 = t160 * t11;
t101 = -t109 * pkin(4) - pkin(3);
t124 = t112 * t107 - t147;
t30 = -qJD(5) * t55 - t124 * t76;
t7 = -t157 * qJD(5) - t112 * t20 + t114 * t17;
t4 = -t30 * pkin(9) + t163 + t7;
t137 = -t111 * t5 + t160 * t4;
t133 = -t112 * t29 + t114 * t26;
t132 = -t112 * t91 - t114 * t89;
t50 = qJD(2) * t149 - qJD(3) * t150 + t108 * t135 + t92 * t134;
t131 = qJD(6) * t142;
t130 = 0.2e1 * t146 * qJD(4);
t129 = 0.2e1 * (t108 ^ 2 + t110 ^ 2) * qJD(2);
t39 = pkin(4) * t153 + t50;
t64 = t113 * t92 + t161 * t90;
t115 = t160 * t124;
t35 = qJD(6) * t115 + t111 * t154 + t87 * t143 + t160 * t75;
t61 = -t111 * t124 + t160 * t87;
t128 = -t117 * t35 - t61 * t77;
t127 = t50 * t88 + t64 * t76;
t126 = -t117 * t75 - t87 * t77;
t125 = -t24 * t107 + t25 * t109;
t51 = pkin(4) * t152 + t64;
t56 = t124 * t88;
t10 = t56 * pkin(9) + t133 - t162;
t122 = -t160 * t10 + t111 * t11;
t121 = t111 * t10 + t138;
t53 = -t87 * pkin(9) + t132;
t54 = -t124 * pkin(9) + t155;
t120 = t111 * t54 - t160 * t53;
t119 = t111 * t53 + t160 * t54;
t118 = t111 * t56 - t160 * t55;
t33 = -t111 * t55 - t160 * t56;
t116 = -pkin(3) * t76 - qJ(4) * t77 + qJD(4) * t117;
t48 = -t87 * qJD(4) - t155 * qJD(5);
t66 = t124 * pkin(5) + t101;
t63 = -0.2e1 * t117 * t77;
t60 = t111 * t87 + t115;
t42 = t117 * t154 - t124 * t77;
t41 = t75 * pkin(9) + t48;
t40 = -t154 * pkin(9) - t47;
t36 = t61 * qJD(6) - t111 * t75 + t160 * t154;
t34 = t55 * pkin(5) + t51;
t18 = t117 * t36 - t60 * t77;
t14 = t31 * pkin(5) + t39;
t13 = -t119 * qJD(6) - t111 * t40 + t160 * t41;
t12 = t120 * qJD(6) - t111 * t41 - t160 * t40;
t9 = t33 * qJD(6) + t111 * t30 + t160 * t31;
t8 = t118 * qJD(6) - t111 * t31 + t160 * t30;
t2 = -t121 * qJD(6) + t137;
t1 = t122 * qJD(6) - t111 * t4 - t141;
t3 = [0, 0, 0, 0, 0, t129, qJ(2) * t129, 0.2e1 * t88 * t76, 0.2e1 * t117 * t76 - 0.2e1 * t88 * t77, 0, 0, 0, t77 * t164, t76 * t164, 0.2e1 * t127 * t107 - 0.2e1 * t117 * t24 + 0.2e1 * t37 * t77, 0.2e1 * t127 * t109 + 0.2e1 * t117 * t25 - 0.2e1 * t38 * t77, 0.2e1 * (-t24 * t88 - t37 * t76) * t109 + 0.2e1 * (-t25 * t88 - t38 * t76) * t107, 0.2e1 * t37 * t24 + 0.2e1 * t38 * t25 + 0.2e1 * t64 * t50, -0.2e1 * t56 * t30, -0.2e1 * t30 * t55 + 0.2e1 * t56 * t31, -0.2e1 * t117 * t30 - 0.2e1 * t56 * t77, 0.2e1 * t117 * t31 - 0.2e1 * t55 * t77, t63, -0.2e1 * t117 * t7 + 0.2e1 * t133 * t77 + 0.2e1 * t51 * t31 + 0.2e1 * t39 * t55, -0.2e1 * t117 * t6 - 0.2e1 * t157 * t77 + 0.2e1 * t51 * t30 - 0.2e1 * t39 * t56, 0.2e1 * t33 * t8, 0.2e1 * t118 * t8 - 0.2e1 * t33 * t9, -0.2e1 * t117 * t8 + 0.2e1 * t33 * t77, 0.2e1 * t117 * t9 + 0.2e1 * t118 * t77, t63, -0.2e1 * t117 * t2 - 0.2e1 * t118 * t14 - 0.2e1 * t122 * t77 + 0.2e1 * t34 * t9, -0.2e1 * t1 * t117 - 0.2e1 * t121 * t77 + 0.2e1 * t14 * t33 + 0.2e1 * t34 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, t109 * t77, -t107 * t77, -t146 * t76, t25 * t107 + t24 * t109, 0, 0, 0, 0, 0, t42, t126, 0, 0, 0, 0, 0, t18, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t77, 0, -t50, t49, t116 * t107 - t50 * t109, t50 * t107 + t116 * t109, t125, -t50 * pkin(3) + (-t107 * t37 + t109 * t38) * qJD(4) + t125 * qJ(4), t30 * t87 + t56 * t75, -t30 * t124 + t56 * t154 - t87 * t31 + t75 * t55, -t126, t42, 0, t101 * t31 - t117 * t48 + t39 * t124 + t132 * t77 + t51 * t154, t101 * t30 - t117 * t47 - t155 * t77 + t39 * t87 - t51 * t75, -t33 * t35 + t8 * t61, -t118 * t35 - t33 * t36 - t8 * t60 - t61 * t9, -t128, t18, 0, -t117 * t13 - t118 * t139 - t120 * t77 + t14 * t60 + t34 * t36 + t66 * t9, -t117 * t12 - t119 * t77 + t33 * t139 + t14 * t61 - t34 * t35 + t66 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, qJ(4) * t130, t87 * t165, 0.2e1 * t75 * t124 - 0.2e1 * t87 * t154, 0, 0, 0, 0.2e1 * t101 * t154, t101 * t165, -0.2e1 * t61 * t35, 0.2e1 * t35 * t60 - 0.2e1 * t61 * t36, 0, 0, 0, 0.2e1 * t60 * t139 + 0.2e1 * t66 * t36, 0.2e1 * t61 * t139 - 0.2e1 * t66 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t151, 0, t50, 0, 0, 0, 0, 0, t31, t30, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t75, 0, 0, 0, 0, 0, t36, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, t77, t7, t6, 0, 0, t8, -t9, t77, t77 * t142 + (-t138 + (-t10 + t162) * t111) * qJD(6) + t137, -t141 + (-t4 - t163) * t111 + (t117 * t142 + t122) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t75, 0, 0, 0, 0, 0, -t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t154, 0, t48, t47, 0, 0, -t35, -t36, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t140, -0.2e1 * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t77, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
