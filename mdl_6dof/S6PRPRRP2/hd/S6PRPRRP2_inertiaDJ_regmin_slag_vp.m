% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:06
% EndTime: 2019-03-08 20:03:09
% DurationCPUTime: 1.11s
% Computational Cost: add. (1091->177), mult. (2967->294), div. (0->0), fcn. (2805->10), ass. (0->107)
t61 = cos(qJ(5));
t58 = sin(qJ(5));
t77 = t61 * pkin(5) + t58 * qJ(6);
t122 = t77 * qJD(5) - t61 * qJD(6);
t59 = sin(qJ(4));
t125 = t122 * t59;
t54 = sin(pkin(11));
t55 = sin(pkin(6));
t56 = cos(pkin(11));
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t27 = (t54 * t63 + t56 * t60) * t55;
t57 = cos(pkin(6));
t62 = cos(qJ(4));
t21 = t27 * t62 + t57 * t59;
t109 = qJD(2) * t55;
t71 = t54 * t60 - t56 * t63;
t25 = t71 * t109;
t11 = t21 * qJD(4) - t25 * t59;
t20 = t27 * t59 - t57 * t62;
t49 = qJD(5) * t61;
t4 = t11 * t58 + t20 * t49;
t124 = -0.4e1 * t59;
t48 = -t56 * pkin(2) - pkin(3);
t119 = t59 * pkin(9);
t80 = -t62 * pkin(4) - t119;
t37 = t48 + t80;
t29 = t58 * t37;
t47 = t54 * pkin(2) + pkin(8);
t114 = t47 * t62;
t39 = t61 * t114;
t69 = t39 + t29;
t50 = t58 ^ 2;
t52 = t61 ^ 2;
t83 = (t50 - t52) * qJD(5);
t76 = pkin(5) * t58 - qJ(6) * t61;
t68 = t47 + t76;
t23 = t68 * t59;
t30 = t76 * qJD(5) - t58 * qJD(6);
t41 = -pkin(4) - t77;
t123 = (-t41 * t62 + t119) * qJD(4) - qJD(5) * t23 - t30 * t59;
t120 = pkin(9) * t62;
t79 = pkin(4) * t59 - t120;
t40 = t79 * qJD(4);
t102 = t59 * qJD(4);
t86 = t47 * t102;
t10 = -t69 * qJD(5) + t61 * t40 + t58 * t86;
t121 = 0.2e1 * qJD(6);
t117 = t41 * t59;
t116 = t47 * t58;
t115 = t47 * t61;
t113 = t37 * t49 + t58 * t40;
t51 = t59 ^ 2;
t110 = -t62 ^ 2 + t51;
t108 = qJD(4) * t20;
t107 = qJD(4) * t23;
t106 = qJD(4) * t61;
t105 = qJD(5) * t58;
t104 = qJD(5) * t59;
t103 = qJD(5) * t62;
t100 = t62 * qJD(4);
t99 = -0.2e1 * pkin(4) * qJD(5);
t98 = 0.2e1 * qJD(4) * t48;
t97 = pkin(5) * t102;
t96 = pkin(9) * t105;
t95 = pkin(9) * t49;
t94 = t47 * t105;
t93 = t58 * t103;
t92 = t61 * t103;
t90 = t50 * t100;
t89 = t58 * t100;
t88 = t58 * t49;
t87 = t59 * t100;
t85 = t61 * t102;
t84 = t61 * t100;
t82 = t110 * qJD(4);
t81 = t58 * t84;
t26 = t71 * t55;
t13 = t21 * t58 - t26 * t61;
t14 = t21 * t61 + t26 * t58;
t75 = t13 * t61 - t14 * t58;
t74 = t13 * t58 + t14 * t61;
t18 = -t62 * qJ(6) + t69;
t19 = -t61 * t37 + (pkin(5) + t116) * t62;
t73 = t18 * t61 + t19 * t58;
t72 = -t18 * t58 + t19 * t61;
t5 = t20 * t105 - t11 * t61;
t32 = t85 + t93;
t12 = -t25 * t62 - t108;
t24 = qJD(2) * t27;
t2 = t14 * qJD(5) + t12 * t58 - t24 * t61;
t67 = -t13 * t102 + t2 * t62 + t20 * t89 + t4 * t59;
t3 = -t13 * qJD(5) + t12 * t61 + t24 * t58;
t65 = t75 * qJD(5) + t2 * t58 + t3 * t61;
t6 = (-qJD(6) - t94) * t62 + (qJ(6) - t115) * t102 + t113;
t7 = -t10 - t97;
t64 = t72 * qJD(5) + t7 * t58 + t6 * t61;
t45 = t52 * t100;
t43 = pkin(9) * t92;
t42 = t52 * t87;
t34 = -t58 * t102 + t92;
t33 = t59 * t49 + t89;
t31 = t58 * t104 - t84;
t15 = t68 * t100 + t125;
t9 = t32 * t47 - t113;
t1 = (t20 * t106 + t3) * t62 + (-qJD(4) * t14 - t5) * t59;
t8 = [0, 0, 0, 0, 0.2e1 * t26 * t24 - 0.2e1 * t27 * t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t20 * t11 + 0.2e1 * t13 * t2 + 0.2e1 * t14 * t3; 0, 0, -t60 * t109, -t63 * t109 (-t24 * t56 - t25 * t54) * pkin(2), 0, 0, 0, 0, 0, t26 * t102 - t24 * t62, t26 * t100 + t24 * t59, 0, 0, 0, 0, 0, t67, t1, t67, t75 * t100 + (-t74 * qJD(5) + t2 * t61 - t3 * t58) * t59, -t1, t11 * t23 + t13 * t7 + t14 * t6 + t20 * t15 + t3 * t18 + t2 * t19; 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t82, 0, 0, 0, t59 * t98, t62 * t98, -0.2e1 * t51 * t88 + 0.2e1 * t42, t81 * t124 + 0.2e1 * t51 * t83, 0.2e1 * t110 * t106 + 0.2e1 * t59 * t93, -0.2e1 * t58 * t82 + 0.2e1 * t59 * t92, -0.2e1 * t87, 0.2e1 * t37 * t85 - 0.2e1 * t10 * t62 + 0.2e1 * (t51 * t49 + t58 * t87) * t47, -0.2e1 * t51 * t94 - 0.2e1 * t9 * t62 + 0.2e1 * (-t29 + t39) * t102, 0.2e1 * (t58 * t107 + t7) * t62 + 0.2e1 * (-qJD(4) * t19 + t15 * t58 + t23 * t49) * t59, 0.2e1 * t72 * t100 + 0.2e1 * (-t73 * qJD(5) - t58 * t6 + t61 * t7) * t59, 0.2e1 * (-t23 * t106 - t6) * t62 + 0.2e1 * (qJD(4) * t18 + t23 * t105 - t15 * t61) * t59, 0.2e1 * t23 * t15 + 0.2e1 * t18 * t6 + 0.2e1 * t19 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t74 * qJD(4) - t11) * t62 + (t65 + t108) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t73 * qJD(4) - t15) * t62 + (t64 + t107) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42 + 0.2e1 * (t50 - 0.1e1) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, t5, t4, t5, t65, -t4, pkin(9) * t65 + t11 * t41 + t20 * t30; 0, 0, 0, 0, 0, 0, 0, t100, -t102, 0, -t47 * t100, t86, -t59 * t83 + t81, t88 * t124 + t45 - t90, -t34, t32, 0, t43 + (-pkin(4) * t61 + t116) * t104 + (t80 * t58 - t39) * qJD(4) (t59 * t115 + t79 * t58) * qJD(5) + (t58 * t114 + t80 * t61) * qJD(4), t43 + (t41 * t104 - t15) * t61 - t123 * t58, t64 (-t15 + (t117 + t120) * qJD(5)) * t58 + t123 * t61, pkin(9) * t64 + t15 * t41 + t23 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t100, 0, 0, 0, 0, 0, -t32, -t34, -t32, t45 + t90, t34, -t62 * t30 + (t117 + (t50 + t52) * t120) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t88, -0.2e1 * t83, 0, 0, 0, t58 * t99, t61 * t99, 0.2e1 * t41 * t105 - 0.2e1 * t30 * t61, 0, -0.2e1 * t30 * t58 - 0.2e1 * t41 * t49, 0.2e1 * t41 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(5) + t3 * qJ(6) + t14 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, t102, t10, t9, t10 + 0.2e1 * t97 (-pkin(5) * t100 - qJ(6) * t104) * t61 + (-qJ(6) * t100 + (pkin(5) * qJD(5) - qJD(6)) * t59) * t58 (-0.2e1 * qJD(6) - t94) * t62 + (0.2e1 * qJ(6) - t115) * t102 + t113, -t7 * pkin(5) + t6 * qJ(6) + t18 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t31, -t33, 0, -t31, -t100 * t76 - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t105, 0, -t95, t96, -t95, -t122, -t96, -t122 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, qJ(6) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t31, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
