% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:43:25
% EndTime: 2019-05-06 15:43:29
% DurationCPUTime: 1.11s
% Computational Cost: add. (1450->148), mult. (3341->285), div. (0->0), fcn. (3887->10), ass. (0->96)
t102 = cos(qJ(4));
t63 = sin(pkin(11));
t65 = cos(pkin(11));
t66 = cos(pkin(6));
t64 = sin(pkin(6));
t69 = sin(qJ(2));
t98 = t64 * t69;
t42 = t63 * t98 - t66 * t65;
t43 = t66 * t63 + t65 * t98;
t68 = sin(qJ(4));
t29 = t102 * t43 - t68 * t42;
t110 = -0.2e1 * t29;
t109 = 0.2e1 * t29;
t50 = -t102 * t65 + t68 * t63;
t51 = t102 * t63 + t68 * t65;
t57 = -t65 * pkin(3) - pkin(2);
t75 = -t51 * qJ(5) + t57;
t31 = t50 * pkin(4) + t75;
t108 = -0.2e1 * t31;
t107 = 0.2e1 * t57;
t106 = 0.2e1 * qJ(5);
t105 = pkin(4) + pkin(10);
t104 = pkin(1) * t69;
t71 = cos(qJ(2));
t103 = pkin(1) * t71;
t28 = t102 * t42 + t68 * t43;
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t97 = t64 * t71;
t20 = t67 * t28 - t70 * t97;
t101 = t20 * t70;
t100 = t51 * t50;
t59 = t64 ^ 2;
t99 = t59 * t71;
t96 = t66 * t69;
t95 = t67 * t29;
t94 = t67 * t50;
t93 = t67 * t51;
t92 = t67 * t105;
t24 = t70 * t29;
t91 = t70 * t50;
t44 = t70 * t51;
t90 = t70 * t67;
t89 = t70 * t105;
t88 = pkin(9) + qJ(3);
t83 = pkin(8) * t97;
t38 = t83 + (qJ(3) + t104) * t66;
t39 = (-pkin(2) * t71 - qJ(3) * t69 - pkin(1)) * t64;
t22 = -t63 * t38 + t65 * t39;
t14 = -pkin(3) * t97 - t43 * pkin(9) + t22;
t23 = t65 * t38 + t63 * t39;
t17 = -t42 * pkin(9) + t23;
t9 = t102 * t17 + t68 * t14;
t87 = t63 ^ 2 + t65 ^ 2;
t86 = qJ(5) * t50;
t85 = 0.2e1 * t100;
t84 = 0.2e1 * t97;
t52 = t88 * t63;
t53 = t88 * t65;
t32 = t102 * t52 + t68 * t53;
t82 = t32 * t97;
t33 = t102 * t53 - t68 * t52;
t81 = t33 * t97;
t80 = qJ(3) * t97;
t79 = qJ(5) * t97;
t78 = -t102 * t14 + t68 * t17;
t77 = -t22 * t63 + t23 * t65;
t76 = t105 * t51 + t86;
t6 = t79 - t9;
t55 = pkin(4) * t97;
t7 = t55 + t78;
t54 = pkin(8) * t98;
t41 = t54 + (-pkin(2) - t103) * t66;
t30 = t42 * pkin(3) + t41;
t74 = -t29 * qJ(5) + t30;
t62 = t70 ^ 2;
t61 = t67 ^ 2;
t48 = t51 ^ 2;
t47 = t50 ^ 2;
t46 = pkin(1) * t96 + t83;
t45 = t66 * t103 - t54;
t27 = t29 ^ 2;
t26 = -t50 * pkin(5) + t33;
t25 = t51 * pkin(5) + t32;
t21 = t105 * t50 + t75;
t19 = t70 * t28 + t67 * t97;
t18 = t29 * t51;
t12 = t70 * t21 + t67 * t25;
t11 = -t67 * t21 + t70 * t25;
t10 = t28 * pkin(4) + t74;
t5 = t105 * t28 + t74;
t4 = -t28 * pkin(5) - t6;
t3 = t29 * pkin(5) + pkin(10) * t97 + t7;
t2 = t67 * t3 + t70 * t5;
t1 = t70 * t3 - t67 * t5;
t8 = [1, 0, 0, t59 * t69 ^ 2, 0.2e1 * t69 * t99, 0.2e1 * t64 * t96, t66 * t84, t66 ^ 2, 0.2e1 * pkin(1) * t99 + 0.2e1 * t45 * t66, -0.2e1 * t59 * t104 - 0.2e1 * t46 * t66, -0.2e1 * t22 * t97 + 0.2e1 * t41 * t42, 0.2e1 * t23 * t97 + 0.2e1 * t41 * t43, -0.2e1 * t22 * t43 - 0.2e1 * t23 * t42, t22 ^ 2 + t23 ^ 2 + t41 ^ 2, t27, t28 * t110, t97 * t110, t28 * t84, t59 * t71 ^ 2, 0.2e1 * t30 * t28 + 0.2e1 * t78 * t97, 0.2e1 * t30 * t29 + 0.2e1 * t9 * t97, 0.2e1 * t6 * t28 + 0.2e1 * t7 * t29, -0.2e1 * t10 * t28 - 0.2e1 * t7 * t97, -0.2e1 * t10 * t29 + 0.2e1 * t6 * t97, t10 ^ 2 + t6 ^ 2 + t7 ^ 2, t20 ^ 2, 0.2e1 * t20 * t19, t20 * t109, t19 * t109, t27, 0.2e1 * t1 * t29 - 0.2e1 * t4 * t19, -0.2e1 * t2 * t29 + 0.2e1 * t4 * t20; 0, 0, 0, 0, 0, t98, t97, t66, t45, -t46, -pkin(2) * t42 - t41 * t65 + t63 * t80, -pkin(2) * t43 + t41 * t63 + t65 * t80 (-t42 * t65 + t43 * t63) * qJ(3) + t77, -t41 * pkin(2) + t77 * qJ(3), t18, -t51 * t28 - t29 * t50, -t51 * t97, t50 * t97, 0, t57 * t28 + t30 * t50 + t82, t57 * t29 + t30 * t51 + t81, -t33 * t28 + t32 * t29 + t6 * t50 + t7 * t51, -t10 * t50 - t31 * t28 - t82, -t10 * t51 - t31 * t29 - t81, t10 * t31 + t7 * t32 - t6 * t33, t20 * t94 (t19 * t67 + t101) * t50, t20 * t51 + t29 * t94, t19 * t51 + t29 * t91, t18, t1 * t51 + t11 * t29 - t26 * t19 - t4 * t91, -t12 * t29 - t2 * t51 + t26 * t20 + t4 * t94; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t65, -0.2e1 * pkin(2) * t63, 0.2e1 * t87 * qJ(3), t87 * qJ(3) ^ 2 + pkin(2) ^ 2, t48, -0.2e1 * t100, 0, 0, 0, t50 * t107, t51 * t107, 0.2e1 * t32 * t51 - 0.2e1 * t33 * t50, t50 * t108, t51 * t108, t31 ^ 2 + t32 ^ 2 + t33 ^ 2, t61 * t47, 0.2e1 * t47 * t90, t67 * t85, t70 * t85, t48, 0.2e1 * t11 * t51 - 0.2e1 * t26 * t91, -0.2e1 * t12 * t51 + 0.2e1 * t26 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t43, 0, t41, 0, 0, 0, 0, 0, t28, t29, 0, -t28, -t29, t10, 0, 0, 0, 0, 0, -t95, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t63, 0, -pkin(2), 0, 0, 0, 0, 0, t50, t51, 0, -t50, -t51, t31, 0, 0, 0, 0, 0, -t93, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, -t97, -t78, -t9, -pkin(4) * t29 - qJ(5) * t28, 0.2e1 * t55 + t78, -0.2e1 * t79 + t9, -t7 * pkin(4) - t6 * qJ(5), t101, t70 * t19 - t20 * t67, t24, -t95, 0, -qJ(5) * t19 - t29 * t89 + t4 * t67, qJ(5) * t20 + t29 * t92 + t4 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, 0, -t32, -t33, -pkin(4) * t51 - t86, t32, t33, -t32 * pkin(4) + t33 * qJ(5), t50 * t90 (-t61 + t62) * t50, t44, -t93, 0, t26 * t67 - t76 * t70, t26 * t70 + t76 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t106, pkin(4) ^ 2 + qJ(5) ^ 2, t62, -0.2e1 * t90, 0, 0, 0, t67 * t106, t70 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t97, 0, t7, 0, 0, 0, 0, 0, t24, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, t32, 0, 0, 0, 0, 0, t44, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t19, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t91, t51, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t67, 0, -t89, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
