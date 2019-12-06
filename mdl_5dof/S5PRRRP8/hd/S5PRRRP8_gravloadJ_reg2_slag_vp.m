% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t57 = cos(qJ(3));
t54 = sin(qJ(3));
t87 = pkin(8) * t54;
t89 = -pkin(3) * t57 - t87;
t52 = sin(pkin(5));
t55 = sin(qJ(2));
t86 = t52 * t55;
t85 = t52 * t57;
t58 = cos(qJ(2));
t84 = t52 * t58;
t53 = sin(qJ(4));
t83 = t53 * t57;
t56 = cos(qJ(4));
t82 = t56 * t57;
t81 = t57 * t58;
t80 = pkin(2) * t84 + pkin(7) * t86;
t79 = cos(pkin(5));
t78 = cos(pkin(9));
t77 = t53 * t84;
t51 = sin(pkin(9));
t65 = t79 * t78;
t36 = t51 * t55 - t58 * t65;
t37 = t51 * t58 + t55 * t65;
t76 = -t36 * pkin(2) + pkin(7) * t37;
t71 = t51 * t79;
t38 = t55 * t78 + t58 * t71;
t39 = -t55 * t71 + t58 * t78;
t75 = -t38 * pkin(2) + pkin(7) * t39;
t70 = t52 * t78;
t17 = -t37 * t54 - t57 * t70;
t18 = t37 * t57 - t54 * t70;
t74 = t17 * pkin(3) + t18 * pkin(8);
t19 = -t39 * t54 + t51 * t85;
t20 = t51 * t52 * t54 + t39 * t57;
t73 = t19 * pkin(3) + t20 * pkin(8);
t40 = -t54 * t86 + t57 * t79;
t41 = t54 * t79 + t55 * t85;
t72 = t40 * pkin(3) + t41 * pkin(8);
t69 = t52 * pkin(3) * t81 + t84 * t87 + t80;
t68 = pkin(4) * t56 + qJ(5) * t53;
t67 = t89 * t36 + t76;
t66 = t38 * t89 + t75;
t21 = t41 * t53 + t56 * t84;
t6 = t18 * t53 - t36 * t56;
t8 = t20 * t53 - t38 * t56;
t1 = g(1) * t8 + g(2) * t6 + g(3) * t21;
t22 = t41 * t56 - t77;
t7 = t18 * t56 + t36 * t53;
t9 = t20 * t56 + t38 * t53;
t64 = g(1) * t9 + g(2) * t7 + g(3) * t22;
t11 = -t36 * t83 - t37 * t56;
t13 = -t38 * t83 - t39 * t56;
t23 = -t56 * t86 + t57 * t77;
t63 = g(1) * t13 + g(2) * t11 + g(3) * t23;
t62 = g(1) * t19 + g(2) * t17 + g(3) * t40;
t61 = g(1) * t20 + g(2) * t18 + g(3) * t41;
t60 = -g(1) * t38 - g(2) * t36 + g(3) * t84;
t59 = g(1) * t39 + g(2) * t37 + g(3) * t86;
t24 = (t53 * t55 + t56 * t81) * t52;
t14 = -t38 * t82 + t39 * t53;
t12 = -t36 * t82 + t37 * t53;
t10 = t60 * t54;
t4 = t62 * t56;
t3 = t62 * t53;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t24;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * t57, t10, -t59, -g(1) * t75 - g(2) * t76 - g(3) * t80, 0, 0, 0, 0, 0, 0, t2, t63, -t10, -g(1) * t66 - g(2) * t67 - g(3) * t69, 0, 0, 0, 0, 0, 0, t2, -t10, -t63, -g(1) * (pkin(4) * t14 + qJ(5) * t13 + t66) - g(2) * (pkin(4) * t12 + qJ(5) * t11 + t67) - g(3) * (pkin(4) * t24 + qJ(5) * t23 + t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t61, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t61, -g(1) * t73 - g(2) * t74 - g(3) * t72, 0, 0, 0, 0, 0, 0, -t4, -t61, -t3, -g(1) * (t19 * t68 + t73) - g(2) * (t17 * t68 + t74) - g(3) * (t40 * t68 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t64, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t64, -g(1) * (-pkin(4) * t8 + qJ(5) * t9) - g(2) * (-pkin(4) * t6 + qJ(5) * t7) - g(3) * (-pkin(4) * t21 + qJ(5) * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
