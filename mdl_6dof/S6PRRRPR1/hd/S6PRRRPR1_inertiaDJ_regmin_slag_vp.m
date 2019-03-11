% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:22
% EndTime: 2019-03-08 23:03:26
% DurationCPUTime: 1.44s
% Computational Cost: add. (2332->165), mult. (5786->326), div. (0->0), fcn. (5801->12), ass. (0->122)
t156 = qJD(3) + qJD(4);
t150 = pkin(8) + pkin(9);
t79 = sin(qJ(3));
t127 = t150 * t79;
t83 = cos(qJ(3));
t153 = t150 * t83;
t78 = sin(qJ(4));
t82 = cos(qJ(4));
t91 = t82 * t127 + t153 * t78;
t75 = sin(pkin(6));
t80 = sin(qJ(2));
t145 = t75 * t80;
t76 = cos(pkin(6));
t98 = t79 * t145 - t76 * t83;
t99 = t83 * t145 + t76 * t79;
t154 = t78 * t99 + t82 * t98;
t113 = t78 * t127;
t31 = t156 * (-t153 * t82 + t113);
t142 = t82 * t83;
t67 = t78 * t79;
t117 = -t67 + t142;
t43 = t156 * t117;
t57 = t78 * t83 + t82 * t79;
t152 = -t43 * qJ(5) - t57 * qJD(5) + t31;
t81 = cos(qJ(6));
t73 = t81 ^ 2;
t77 = sin(qJ(6));
t140 = t77 ^ 2 - t73;
t115 = t140 * qJD(6);
t35 = -t78 * t98 + t82 * t99;
t84 = cos(qJ(2));
t136 = qJD(2) * t84;
t121 = t75 * t136;
t151 = qJD(3) * t98 - t83 * t121;
t30 = t156 * t91;
t138 = cos(pkin(12));
t131 = t83 * qJD(3);
t132 = t79 * qJD(3);
t134 = qJD(4) * t82;
t135 = qJD(4) * t78;
t114 = t78 * t131 + t82 * t132 + t79 * t134 + t83 * t135;
t20 = -t114 * qJ(5) + t117 * qJD(5) - t30;
t74 = sin(pkin(12));
t10 = -t138 * t152 + t20 * t74;
t34 = -t67 * qJ(5) + (qJ(5) + t150) * t142 - t113;
t88 = -t57 * qJ(5) - t91;
t24 = -t138 * t88 + t34 * t74;
t71 = qJD(6) * t81;
t149 = t10 * t77 + t24 * t71;
t40 = t74 * t117 + t138 * t57;
t148 = t40 * t77;
t147 = t40 * t81;
t146 = t74 * t78;
t144 = t75 * t84;
t29 = -t74 * t114 + t138 * t43;
t143 = t81 * t29;
t68 = pkin(3) * t82 + pkin(4);
t53 = -pkin(3) * t146 + t138 * t68;
t49 = -pkin(5) - t53;
t116 = t138 * t78;
t139 = pkin(3) * qJD(4);
t51 = (t74 * t82 + t116) * t139;
t141 = t49 * t71 + t51 * t77;
t54 = pkin(3) * t116 + t74 * t68;
t137 = qJD(2) * t80;
t133 = qJD(6) * t77;
t130 = -0.2e1 * pkin(2) * qJD(3);
t70 = pkin(3) * t132;
t129 = pkin(3) * t135;
t128 = pkin(3) * t134;
t66 = -t138 * pkin(4) - pkin(5);
t124 = t66 * t133;
t123 = t66 * t71;
t122 = t75 * t137;
t120 = t77 * t71;
t69 = -t83 * pkin(3) - pkin(2);
t119 = -0.4e1 * t77 * t147;
t118 = t49 * t133 - t51 * t81;
t39 = -t138 * t117 + t57 * t74;
t48 = -t117 * pkin(4) + t69;
t23 = t39 * pkin(5) - t40 * pkin(10) + t48;
t25 = t138 * t34 + t74 * t88;
t109 = t23 * t81 - t25 * t77;
t108 = t23 * t77 + t25 * t81;
t28 = t138 * t114 + t43 * t74;
t65 = pkin(4) * t74 + pkin(10);
t107 = -t28 * t65 + t29 * t66;
t50 = pkin(10) + t54;
t106 = t39 * t50 - t40 * t49;
t52 = (t138 * t82 - t146) * t139;
t105 = -t39 * t52 + t40 * t51;
t104 = t39 * t65 - t40 * t66;
t27 = t138 * t35 - t154 * t74;
t101 = t81 * t144 + t77 * t27;
t100 = t77 * t144 - t81 * t27;
t17 = t28 * t77 + t39 * t71;
t97 = t77 * t29 + t40 * t71;
t96 = -t40 * t133 + t143;
t36 = t114 * pkin(4) + t70;
t90 = -t28 * t50 + t29 * t49 + t105;
t87 = -qJD(3) * t99 - t79 * t121;
t85 = -t57 * t121 - t156 * t35;
t59 = 0.2e1 * t120;
t56 = -0.2e1 * t115;
t37 = t40 ^ 2;
t26 = t138 * t154 + t35 * t74;
t21 = t24 * t133;
t19 = t154 * qJD(4) + t82 * t151 - t78 * t87;
t16 = -t39 * t133 + t28 * t81;
t14 = -t40 * t115 + t77 * t143;
t13 = t28 * pkin(5) - t29 * pkin(10) + t36;
t12 = qJD(6) * t119 - t140 * t29;
t11 = t138 * t20 + t152 * t74;
t9 = -t138 * t19 + t74 * t85;
t8 = -t138 * t85 - t19 * t74;
t6 = t26 * t133 - t8 * t81;
t5 = t26 * t71 + t8 * t77;
t4 = qJD(6) * t100 + t81 * t122 - t77 * t9;
t3 = qJD(6) * t101 - t77 * t122 - t81 * t9;
t2 = -t108 * qJD(6) - t11 * t77 + t13 * t81;
t1 = -t109 * qJD(6) - t11 * t81 - t13 * t77;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t75 ^ 2 * t80 * t136 + 0.2e1 * t26 * t8 + 0.2e1 * t27 * t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t122, -t121, 0, 0, 0, 0, 0 (-t84 * t132 - t83 * t137) * t75 (-t84 * t131 + t79 * t137) * t75, 0, 0, 0, 0, 0 (-t84 * t114 - t117 * t137) * t75 (t57 * t137 - t43 * t84) * t75, t26 * t29 - t27 * t28 - t39 * t9 + t8 * t40, t26 * t10 + t27 * t11 + t8 * t24 + t9 * t25 + (t48 * t137 - t36 * t84) * t75, 0, 0, 0, 0, 0, -t101 * t28 + t8 * t148 + t26 * t97 + t4 * t39, t100 * t28 + t8 * t147 + t26 * t96 + t3 * t39; 0, 0, 0, 0, 0.2e1 * t79 * t131, 0.2e1 * (-t79 ^ 2 + t83 ^ 2) * qJD(3), 0, 0, 0, t79 * t130, t83 * t130, 0.2e1 * t57 * t43, -0.2e1 * t57 * t114 + 0.2e1 * t43 * t117, 0, 0, 0, 0.2e1 * t69 * t114 - 0.2e1 * t117 * t70, 0.2e1 * t43 * t69 + 0.2e1 * t57 * t70, 0.2e1 * t10 * t40 - 0.2e1 * t11 * t39 + 0.2e1 * t24 * t29 - 0.2e1 * t25 * t28, 0.2e1 * t10 * t24 + 0.2e1 * t11 * t25 + 0.2e1 * t36 * t48, 0.2e1 * t29 * t40 * t73 - 0.2e1 * t37 * t120, 0.2e1 * t37 * t115 + t29 * t119, 0.2e1 * t28 * t147 + 0.2e1 * t39 * t96, -0.2e1 * t28 * t148 - 0.2e1 * t39 * t97, 0.2e1 * t39 * t28, 0.2e1 * t10 * t148 + 0.2e1 * t109 * t28 + 0.2e1 * t2 * t39 + 0.2e1 * t97 * t24, 0.2e1 * t1 * t39 + 0.2e1 * t10 * t147 - 0.2e1 * t108 * t28 + 0.2e1 * t96 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t151, 0, 0, 0, 0, 0, t85, t19, 0, t26 * t51 + t27 * t52 - t53 * t8 + t54 * t9, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, t131, -t132, 0, -pkin(8) * t131, pkin(8) * t132, 0, 0, t43, -t114, 0, t31, t30, -t28 * t54 - t29 * t53 + t105, -t10 * t53 + t11 * t54 + t24 * t51 + t25 * t52, t14, t12, t17, t16, 0, t21 + (-t106 * qJD(6) - t10) * t81 + t90 * t77, t106 * t133 + t90 * t81 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t129, -0.2e1 * t128, 0, -0.2e1 * t51 * t53 + 0.2e1 * t52 * t54, t59, t56, 0, 0, 0, 0.2e1 * t118, 0.2e1 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t19, 0 (-t138 * t8 + t74 * t9) * pkin(4), 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t114, 0, t31, t30 (-t138 * t29 - t28 * t74) * pkin(4) (-t138 * t10 + t11 * t74) * pkin(4), t14, t12, t17, t16, 0, t21 + t107 * t77 + (-t104 * qJD(6) - t10) * t81, t104 * t133 + t107 * t81 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t128, 0 (-t138 * t51 + t52 * t74) * pkin(4), t59, t56, 0, 0, 0, t118 + t124, t123 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t56, 0, 0, 0, 0.2e1 * t124, 0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t97, t28, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t133, 0, -t50 * t71 - t52 * t77, t50 * t133 - t52 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t133, 0, -t65 * t71, t65 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
