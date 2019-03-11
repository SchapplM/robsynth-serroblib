% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:53
% EndTime: 2019-03-08 19:32:57
% DurationCPUTime: 1.06s
% Computational Cost: add. (867->152), mult. (2538->301), div. (0->0), fcn. (2497->12), ass. (0->97)
t61 = sin(pkin(12));
t64 = cos(pkin(12));
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t111 = -t67 * t61 + t70 * t64;
t40 = t111 * qJD(6);
t55 = -t64 * pkin(5) - pkin(4);
t110 = 0.2e1 * t55;
t68 = sin(qJ(4));
t109 = pkin(4) * t68;
t108 = t61 * t68;
t71 = cos(qJ(4));
t107 = t61 * t71;
t106 = t64 * t68;
t105 = t64 * t71;
t102 = pkin(9) + qJ(5);
t39 = -t68 * qJD(5) + (-qJ(5) * t71 + t109) * qJD(4);
t62 = sin(pkin(11));
t53 = t62 * pkin(2) + pkin(8);
t56 = t68 * qJD(4);
t90 = t53 * t56;
t20 = t64 * t39 + t61 * t90;
t65 = cos(pkin(11));
t54 = -t65 * pkin(2) - pkin(3);
t81 = -t71 * pkin(4) - t68 * qJ(5);
t42 = t54 + t81;
t46 = t53 * t105;
t23 = t61 * t42 + t46;
t101 = t61 ^ 2 + t64 ^ 2;
t100 = t68 ^ 2 - t71 ^ 2;
t63 = sin(pkin(6));
t99 = qJD(2) * t63;
t98 = qJD(5) * t71;
t96 = t71 * qJD(4);
t95 = t53 * t107;
t94 = 0.2e1 * qJD(4) * t54;
t93 = t61 * t96;
t92 = t64 * t96;
t91 = t68 * t96;
t47 = t53 * t96;
t89 = t101 * t71;
t88 = t101 * qJD(5);
t87 = 0.2e1 * t91;
t86 = 0.2e1 * t88;
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t33 = (t62 * t72 + t65 * t69) * t63;
t27 = qJD(2) * t33;
t66 = cos(pkin(6));
t24 = t33 * t68 - t66 * t71;
t74 = t62 * t69 - t65 * t72;
t28 = t74 * t99;
t8 = -t24 * qJD(4) - t28 * t71;
t5 = t27 * t64 - t8 * t61;
t6 = t27 * t61 + t8 * t64;
t85 = -t5 * t61 + t6 * t64;
t25 = t33 * t71 + t66 * t68;
t32 = t74 * t63;
t10 = t25 * t64 + t32 * t61;
t9 = -t25 * t61 + t32 * t64;
t84 = t10 * t64 - t61 * t9;
t83 = t70 * t10 + t67 * t9;
t82 = t67 * t10 - t70 * t9;
t37 = t64 * t42;
t15 = -pkin(9) * t106 + t37 + (-t53 * t61 - pkin(5)) * t71;
t18 = -pkin(9) * t108 + t23;
t80 = t70 * t15 - t67 * t18;
t79 = t67 * t15 + t70 * t18;
t29 = t61 * t39;
t21 = -t64 * t90 + t29;
t78 = -t20 * t61 + t21 * t64;
t22 = t37 - t95;
t77 = -t22 * t61 + t23 * t64;
t50 = t102 * t61;
t51 = t102 * t64;
t76 = -t70 * t50 - t67 * t51;
t75 = -t67 * t50 + t70 * t51;
t45 = t70 * t61 + t67 * t64;
t41 = t45 * qJD(6);
t73 = -t111 * t56 - t71 * t41;
t38 = (pkin(5) * t61 + t53) * t68;
t35 = t111 * t68;
t34 = t45 * t68;
t31 = pkin(5) * t93 + t47;
t19 = -t40 * t71 + t45 * t56;
t17 = t68 * t40 + t45 * t96;
t16 = t68 * t41 + t67 * t93 - t70 * t92;
t14 = t29 + (-pkin(9) * t107 - t53 * t106) * qJD(4);
t13 = -t45 * qJD(5) - t75 * qJD(6);
t12 = -qJD(5) * t111 - t76 * qJD(6);
t11 = (pkin(5) * t68 - pkin(9) * t105) * qJD(4) + t20;
t7 = t25 * qJD(4) - t28 * t68;
t4 = -t79 * qJD(6) + t70 * t11 - t67 * t14;
t3 = -t80 * qJD(6) - t67 * t11 - t70 * t14;
t2 = -t83 * qJD(6) + t70 * t5 - t67 * t6;
t1 = t82 * qJD(6) - t67 * t5 - t70 * t6;
t26 = [0, 0, 0, 0, 0.2e1 * t32 * t27 - 0.2e1 * t33 * t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t6 + 0.2e1 * t24 * t7 + 0.2e1 * t9 * t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t69 * t99, -t72 * t99 (-t27 * t65 - t28 * t62) * pkin(2), 0, 0, 0, 0, 0, -t27 * t71 + t32 * t56, t27 * t68 + t32 * t96, t7 * t108 - t5 * t71 + (t24 * t107 + t68 * t9) * qJD(4), t7 * t106 + t6 * t71 + (-t10 * t68 + t24 * t105) * qJD(4) (-t5 * t64 - t6 * t61) * t68 + (-t10 * t61 - t64 * t9) * t96, t10 * t21 + t9 * t20 + t5 * t22 + t6 * t23 + (t24 * t96 + t68 * t7) * t53, 0, 0, 0, 0, 0, t24 * t17 - t2 * t71 + t7 * t34 - t82 * t56, -t1 * t71 - t24 * t16 + t7 * t35 - t83 * t56; 0, 0, 0, 0, 0, t87, -0.2e1 * t100 * qJD(4), 0, 0, 0, t68 * t94, t71 * t94, -0.2e1 * t20 * t71 + 0.2e1 * (t22 + 0.2e1 * t95) * t56, 0.2e1 * t21 * t71 + 0.2e1 * (-t23 + 0.2e1 * t46) * t56, 0.2e1 * (-t20 * t64 - t21 * t61) * t68 + 0.2e1 * (-t22 * t64 - t23 * t61) * t96, 0.2e1 * t53 ^ 2 * t91 + 0.2e1 * t22 * t20 + 0.2e1 * t23 * t21, -0.2e1 * t35 * t16, 0.2e1 * t16 * t34 - 0.2e1 * t35 * t17, 0.2e1 * t16 * t71 + 0.2e1 * t35 * t56, 0.2e1 * t71 * t17 - 0.2e1 * t34 * t56, -0.2e1 * t91, 0.2e1 * t38 * t17 + 0.2e1 * t31 * t34 - 0.2e1 * t4 * t71 + 0.2e1 * t80 * t56, -0.2e1 * t38 * t16 - 0.2e1 * t3 * t71 + 0.2e1 * t31 * t35 - 0.2e1 * t79 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t71 + t85 * t68 + (t24 * t68 + t84 * t71) * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t68 + (t100 * t53 + t77 * t71) * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t101) * t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7 * t64, t7 * t61, t85, -t7 * pkin(4) + t85 * qJ(5) + t84 * qJD(5), 0, 0, 0, 0, 0, -t111 * t7 + t24 * t41, t24 * t40 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, t96, -t56, 0, -t47, t90, t61 * t98 + (t81 * t61 - t46) * qJD(4), t64 * t98 + (t81 * t64 + t95) * qJD(4), t78, -pkin(4) * t47 + t78 * qJ(5) + t77 * qJD(5), -t16 * t45 + t35 * t40, -t111 * t16 - t45 * t17 - t40 * t34 - t35 * t41, t19, -t73, 0, -t111 * t31 - t13 * t71 + t55 * t17 + t38 * t41 + t76 * t56, -t12 * t71 - t55 * t16 + t31 * t45 + t38 * t40 - t75 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t96, -t64 * t56, t61 * t56, qJD(4) * t89, t68 * t88 + (qJ(5) * t89 - t109) * qJD(4), 0, 0, 0, 0, 0, t73, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, qJ(5) * t86, 0.2e1 * t45 * t40, 0.2e1 * t111 * t40 - 0.2e1 * t45 * t41, 0, 0, 0, t41 * t110, t40 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t92, 0, t47, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, t56, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t26;
