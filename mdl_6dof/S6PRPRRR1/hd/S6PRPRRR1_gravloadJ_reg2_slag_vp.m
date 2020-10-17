% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:07:46
% EndTime: 2019-05-05 00:07:49
% DurationCPUTime: 0.64s
% Computational Cost: add. (629->120), mult. (1321->184), div. (0->0), fcn. (1671->14), ass. (0->74)
t48 = sin(pkin(12));
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t80 = cos(pkin(12));
t36 = -t58 * t48 - t55 * t80;
t52 = cos(pkin(6));
t34 = t36 * t52;
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t66 = -t55 * t48 + t58 * t80;
t20 = -t51 * t34 + t49 * t66;
t21 = -t49 * t34 - t51 * t66;
t50 = sin(pkin(6));
t33 = t36 * t50;
t7 = g(1) * t21 - g(2) * t20 + g(3) * t33;
t54 = sin(qJ(4));
t94 = t21 * t54;
t93 = t33 * t54;
t47 = qJ(4) + qJ(5);
t46 = cos(t47);
t53 = sin(qJ(6));
t92 = t46 * t53;
t56 = cos(qJ(6));
t91 = t46 * t56;
t90 = t49 * t50;
t89 = t49 * t55;
t88 = t50 * t51;
t87 = t50 * t54;
t57 = cos(qJ(4));
t86 = t50 * t57;
t85 = t50 * t58;
t84 = t52 * t55;
t83 = t52 * t57;
t82 = t52 * t58;
t79 = t49 * t86;
t78 = t51 * t86;
t77 = t51 * t82;
t32 = t66 * t50;
t41 = pkin(2) * t85;
t44 = t57 * pkin(4) + pkin(3);
t59 = -pkin(9) - pkin(8);
t76 = t32 * t44 + t33 * t59 + t41;
t45 = sin(t47);
t10 = -t20 * t45 - t46 * t88;
t11 = t20 * t46 - t45 * t88;
t75 = t10 * pkin(5) + t11 * pkin(10);
t12 = t21 * t45 + t46 * t90;
t13 = -t21 * t46 + t45 * t90;
t74 = t12 * pkin(5) + t13 * pkin(10);
t25 = t33 * t45 + t52 * t46;
t26 = -t33 * t46 + t52 * t45;
t73 = t25 * pkin(5) + t26 * pkin(10);
t39 = pkin(2) * t77;
t71 = -pkin(2) * t89 + t39;
t70 = pkin(5) * t46 + pkin(10) * t45;
t69 = -t20 * t54 - t78;
t68 = -t49 * t82 - t51 * t55;
t62 = t66 * t52;
t19 = t49 * t36 + t51 * t62;
t67 = t19 * t44 - t20 * t59 + t71;
t65 = g(1) * t12 + g(2) * t10 + g(3) * t25;
t5 = g(1) * t13 + g(2) * t11 + g(3) * t26;
t22 = t51 * t36 - t49 * t62;
t64 = g(1) * t22 + g(2) * t19 + g(3) * t32;
t63 = t68 * pkin(2);
t61 = t21 * t59 + t22 * t44 + t63;
t60 = -g(1) * t68 - g(3) * t85;
t42 = pkin(4) * t83;
t38 = pkin(4) * t79;
t31 = -g(3) * t52 + (-g(1) * t49 + g(2) * t51) * t50;
t6 = t64 * t45;
t2 = t65 * t56;
t1 = t65 * t53;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t77 - t89) + t60, -g(1) * (t49 * t84 - t51 * t58) - g(2) * (-t49 * t58 - t51 * t84) + g(3) * t50 * t55, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t7, 0, -g(2) * t39 + (g(2) * t89 + t60) * pkin(2), 0, 0, 0, 0, 0, 0, -t64 * t57, t64 * t54, t7, -g(1) * (t22 * pkin(3) - t21 * pkin(8) + t63) - g(2) * (t19 * pkin(3) + pkin(8) * t20 + t71) - g(3) * (t32 * pkin(3) - t33 * pkin(8) + t41) 0, 0, 0, 0, 0, 0, -t64 * t46, t6, t7, -g(1) * t61 - g(2) * t67 - g(3) * t76, 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t53 + t22 * t91) - g(2) * (t19 * t91 + t20 * t53) - g(3) * (t32 * t91 - t33 * t53) -g(1) * (-t21 * t56 - t22 * t92) - g(2) * (-t19 * t92 + t20 * t56) - g(3) * (-t32 * t92 - t33 * t56) -t6, -g(1) * (t70 * t22 + t61) - g(2) * (t70 * t19 + t67) - g(3) * (t70 * t32 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t79 + t94) - g(2) * t69 - g(3) * (t83 + t93) -g(1) * (t21 * t57 - t49 * t87) - g(2) * (-t20 * t57 + t51 * t87) - g(3) * (t33 * t57 - t52 * t54) 0, 0, 0, 0, 0, 0, 0, 0, -t65, t5, 0, -g(1) * t38 - g(3) * t42 + (g(2) * t78 - t7 * t54) * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (pkin(4) * t94 + t38 + t74) - g(2) * (t69 * pkin(4) + t75) - g(3) * (pkin(4) * t93 + t42 + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t74 - g(2) * t75 - g(3) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t53 - t22 * t56) - g(2) * (-t11 * t53 - t19 * t56) - g(3) * (-t26 * t53 - t32 * t56) -g(1) * (-t13 * t56 + t22 * t53) - g(2) * (-t11 * t56 + t19 * t53) - g(3) * (-t26 * t56 + t32 * t53) 0, 0;];
taug_reg  = t3;
