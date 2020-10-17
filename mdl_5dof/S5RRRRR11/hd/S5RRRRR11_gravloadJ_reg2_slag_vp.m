% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:42
% EndTime: 2019-12-31 22:43:45
% DurationCPUTime: 0.76s
% Computational Cost: add. (484->135), mult. (1115->219), div. (0->0), fcn. (1364->12), ass. (0->74)
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t47 = cos(qJ(2));
t72 = cos(pkin(5));
t89 = cos(qJ(1));
t58 = t72 * t89;
t23 = t43 * t58 + t44 * t47;
t42 = sin(qJ(3));
t46 = cos(qJ(3));
t40 = sin(pkin(5));
t68 = t40 * t89;
t11 = t23 * t46 - t42 * t68;
t22 = t44 * t43 - t47 * t58;
t39 = qJ(4) + qJ(5);
t36 = sin(t39);
t37 = cos(t39);
t96 = t11 * t36 - t22 * t37;
t95 = t11 * t37 + t22 * t36;
t41 = sin(qJ(4));
t45 = cos(qJ(4));
t94 = t11 * t41 - t22 * t45;
t93 = t11 * t45 + t22 * t41;
t80 = t40 * t46;
t21 = t72 * t42 + t43 * t80;
t63 = t44 * t72;
t25 = -t43 * t63 + t89 * t47;
t81 = t40 * t44;
t15 = t25 * t46 + t42 * t81;
t24 = t89 * t43 + t47 * t63;
t7 = -t15 * t41 + t24 * t45;
t79 = t40 * t47;
t92 = g(2) * t94 - g(3) * (-t21 * t41 - t45 * t79) - g(1) * t7;
t90 = g(3) * t40;
t84 = t36 * t46;
t83 = t37 * t46;
t82 = t40 * t43;
t78 = t41 * t43;
t77 = t41 * t46;
t76 = t45 * t46;
t75 = t46 * t47;
t74 = pkin(2) * t79 + pkin(8) * t82;
t73 = t89 * pkin(1) + pkin(7) * t81;
t71 = t25 * pkin(2) + t73;
t70 = pkin(4) * t41 + pkin(8);
t69 = g(3) * t74;
t67 = -t44 * pkin(1) + pkin(7) * t68;
t16 = t22 * pkin(2);
t66 = t23 * pkin(8) - t16;
t18 = t24 * pkin(2);
t65 = t25 * pkin(8) - t18;
t64 = -t23 * t42 - t46 * t68;
t62 = -t23 * pkin(2) + t67;
t61 = pkin(3) * t46 + pkin(9) * t42;
t14 = t25 * t42 - t44 * t80;
t60 = g(1) * t64 + g(2) * t14;
t59 = g(1) * t22 - g(2) * t24;
t35 = t45 * pkin(4) + pkin(3);
t48 = -pkin(10) - pkin(9);
t57 = t35 * t46 - t42 * t48;
t56 = t24 * pkin(8) + t71;
t55 = g(1) * t89 + g(2) * t44;
t54 = -t22 * pkin(8) + t62;
t20 = -t42 * t82 + t72 * t46;
t53 = g(1) * t14 - g(2) * t64 - g(3) * t20;
t52 = g(1) * t15 + g(2) * t11 + g(3) * t21;
t51 = -g(1) * t24 - g(2) * t22 + g(3) * t79;
t50 = g(1) * t25 + g(2) * t23 + g(3) * t82;
t9 = t51 * t42;
t8 = t15 * t45 + t24 * t41;
t6 = t15 * t37 + t24 * t36;
t5 = -t15 * t36 + t24 * t37;
t2 = g(1) * t6 + g(2) * t95 - g(3) * (-t21 * t37 + t36 * t79);
t1 = -g(1) * t5 + g(2) * t96 - g(3) * (-t21 * t36 - t37 * t79);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t89, t55, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t23 - g(2) * t25, -t59, -t55 * t40, -g(1) * t67 - g(2) * t73, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t15, t60, t59, -g(1) * t54 - g(2) * t56, 0, 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t8, -g(1) * t94 - g(2) * t7, -t60, -g(1) * (-pkin(3) * t11 + pkin(9) * t64 + t54) - g(2) * (t15 * pkin(3) + t14 * pkin(9) + t56), 0, 0, 0, 0, 0, 0, g(1) * t95 - g(2) * t6, -g(1) * t96 - g(2) * t5, -t60, -g(1) * (-t11 * t35 - t70 * t22 - t48 * t64 + t62) - g(2) * (-t14 * t48 + t15 * t35 + t70 * t24 + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t46, t9, -t50, -g(1) * t65 - g(2) * t66 - t69, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t76 + t25 * t41) - g(2) * (-t22 * t76 + t23 * t41) - (t45 * t75 + t78) * t90, -g(1) * (t24 * t77 + t25 * t45) - g(2) * (t22 * t77 + t23 * t45) - (-t41 * t75 + t43 * t45) * t90, -t9, -g(1) * (-t61 * t24 + t65) - g(2) * (-t61 * t22 + t66) - g(3) * (t61 * t79 + t74), 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t83 + t25 * t36) - g(2) * (-t22 * t83 + t23 * t36) - (t36 * t43 + t37 * t75) * t90, -g(1) * (t24 * t84 + t25 * t37) - g(2) * (t22 * t84 + t23 * t37) - (-t36 * t75 + t37 * t43) * t90, -t9, -g(1) * (-t57 * t24 + t70 * t25 - t18) - g(2) * (-t57 * t22 + t70 * t23 - t16) - t69 - (pkin(4) * t78 + t57 * t47) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t52, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t45, -t53 * t41, -t52, -g(1) * (-t14 * pkin(3) + t15 * pkin(9)) - g(2) * (pkin(3) * t64 + t11 * pkin(9)) - g(3) * (t20 * pkin(3) + t21 * pkin(9)), 0, 0, 0, 0, 0, 0, t53 * t37, -t53 * t36, -t52, -g(1) * (-t14 * t35 - t15 * t48) - g(2) * (-t11 * t48 + t35 * t64) - g(3) * (t20 * t35 - t21 * t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, g(1) * t8 + g(2) * t93 - g(3) * (-t21 * t45 + t41 * t79), 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t92 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t3;
