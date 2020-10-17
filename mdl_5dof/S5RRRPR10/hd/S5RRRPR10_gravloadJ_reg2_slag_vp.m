% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:51
% EndTime: 2019-12-31 21:29:53
% DurationCPUTime: 0.60s
% Computational Cost: add. (435->122), mult. (838->199), div. (0->0), fcn. (1000->12), ass. (0->72)
t51 = sin(qJ(2));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t76 = cos(pkin(5));
t67 = t56 * t76;
t25 = t52 * t51 - t55 * t67;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t26 = t51 * t67 + t52 * t55;
t46 = qJ(3) + pkin(10);
t43 = sin(t46);
t44 = cos(t46);
t47 = sin(pkin(5));
t83 = t47 * t56;
t8 = t26 * t44 - t43 * t83;
t95 = -t25 * t53 + t8 * t49;
t94 = t25 * t49 + t8 * t53;
t93 = g(3) * t47;
t68 = t52 * t76;
t28 = -t51 * t68 + t56 * t55;
t50 = sin(qJ(3));
t90 = t28 * t50;
t89 = t44 * t49;
t88 = t44 * t53;
t87 = t47 * t51;
t86 = t47 * t52;
t54 = cos(qJ(3));
t85 = t47 * t54;
t84 = t47 * t55;
t48 = -qJ(4) - pkin(8);
t82 = t48 * t51;
t81 = t49 * t55;
t80 = t53 * t55;
t42 = t54 * pkin(3) + pkin(2);
t79 = -t25 * t42 - t26 * t48;
t27 = t56 * t51 + t55 * t68;
t78 = -t27 * t42 - t28 * t48;
t77 = t56 * pkin(1) + pkin(7) * t86;
t75 = t50 * t87;
t74 = t50 * t86;
t73 = t52 * t85;
t37 = t50 * t83;
t72 = t54 * t83;
t71 = -t52 * pkin(1) + pkin(7) * t83;
t70 = -t26 * t43 - t44 * t83;
t69 = t26 * t54 - t37;
t66 = t76 * t54;
t11 = t28 * t43 - t44 * t86;
t65 = g(1) * t70 + g(2) * t11;
t64 = pkin(3) * t74 - t27 * t48 + t28 * t42 + t77;
t63 = pkin(4) * t44 + pkin(9) * t43;
t6 = g(1) * t25 - g(2) * t27;
t62 = g(1) * t56 + g(2) * t52;
t61 = t26 * t50 + t72;
t60 = pkin(3) * t37 + t25 * t48 - t26 * t42 + t71;
t19 = -t43 * t87 + t76 * t44;
t59 = g(1) * t11 - g(2) * t70 - g(3) * t19;
t12 = t28 * t44 + t43 * t86;
t20 = t76 * t43 + t44 * t87;
t58 = g(1) * t12 + g(2) * t8 + g(3) * t20;
t4 = -g(1) * t27 - g(2) * t25 + g(3) * t84;
t57 = g(1) * t28 + g(2) * t26 + g(3) * t87;
t41 = pkin(3) * t66;
t34 = pkin(3) * t73;
t29 = t42 * t84;
t14 = t28 * t54 + t74;
t13 = t73 - t90;
t3 = t12 * t53 + t27 * t49;
t2 = -t12 * t49 + t27 * t53;
t1 = t4 * t43;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t56, t62, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t28, -t6, -t62 * t47, -g(1) * t71 - g(2) * t77, 0, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t14, -g(1) * t61 - g(2) * t13, t6, -g(1) * (-t26 * pkin(2) - t25 * pkin(8) + t71) - g(2) * (t28 * pkin(2) + t27 * pkin(8) + t77), 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t12, t65, t6, -g(1) * t60 - g(2) * t64, 0, 0, 0, 0, 0, 0, g(1) * t94 - g(2) * t3, -g(1) * t95 - g(2) * t2, -t65, -g(1) * (-pkin(4) * t8 + pkin(9) * t70 + t60) - g(2) * (t12 * pkin(4) + t11 * pkin(9) + t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t57, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t54, t4 * t50, -t57, -g(1) * (-t27 * pkin(2) + t28 * pkin(8)) - g(2) * (-t25 * pkin(2) + t26 * pkin(8)) - (pkin(2) * t55 + pkin(8) * t51) * t93, 0, 0, 0, 0, 0, 0, -t4 * t44, t1, -t57, -g(1) * t78 - g(2) * t79 - g(3) * (-t47 * t82 + t29), 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t88 + t28 * t49) - g(2) * (-t25 * t88 + t26 * t49) - (t44 * t80 + t49 * t51) * t93, -g(1) * (t27 * t89 + t28 * t53) - g(2) * (t25 * t89 + t26 * t53) - (-t44 * t81 + t51 * t53) * t93, -t1, -g(1) * (-t63 * t27 + t78) - g(2) * (-t63 * t25 + t79) - g(3) * t29 - (t63 * t55 - t82) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t61 - g(3) * (t66 - t75), g(1) * t14 + g(2) * t69 - g(3) * (-t76 * t50 - t51 * t85), 0, 0, 0, 0, 0, 0, 0, 0, t59, t58, 0, -g(1) * t34 - g(3) * t41 + (g(2) * t72 + t57 * t50) * pkin(3), 0, 0, 0, 0, 0, 0, t59 * t53, -t59 * t49, -t58, -g(1) * (-pkin(3) * t90 - t11 * pkin(4) + t12 * pkin(9) + t34) - g(2) * (-t61 * pkin(3) + pkin(4) * t70 + t8 * pkin(9)) - g(3) * (-pkin(3) * t75 + t19 * pkin(4) + t20 * pkin(9) + t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t95 - g(3) * (-t20 * t49 - t47 * t80), g(1) * t3 + g(2) * t94 - g(3) * (-t20 * t53 + t47 * t81), 0, 0;];
taug_reg = t5;
