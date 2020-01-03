% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t102 = cos(qJ(3));
t53 = cos(qJ(1));
t101 = sin(qJ(1));
t95 = cos(pkin(11));
t97 = cos(pkin(5));
t81 = t97 * t95;
t92 = sin(pkin(11));
t67 = t101 * t92 - t53 * t81;
t93 = sin(pkin(6));
t94 = sin(pkin(5));
t78 = t94 * t93;
t96 = cos(pkin(6));
t108 = t53 * t78 + t67 * t96;
t80 = t97 * t92;
t36 = t101 * t95 + t53 * t80;
t50 = sin(qJ(3));
t16 = t108 * t102 + t36 * t50;
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t19 = -t36 * t102 + t108 * t50;
t79 = t94 * t96;
t28 = -t53 * t79 + t67 * t93;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t7 = t19 * t52 - t28 * t49;
t112 = t16 * t51 + t7 * t48;
t111 = -t16 * t48 + t7 * t51;
t107 = t19 * t49 + t28 * t52;
t61 = t101 * t81 + t53 * t92;
t106 = t101 * t79 + t61 * t93;
t104 = -t101 * t78 + t61 * t96;
t103 = t95 * t79 + t93 * t97;
t100 = t48 * t52;
t99 = t51 * t52;
t82 = t101 * t94;
t98 = t53 * pkin(1) + qJ(2) * t82;
t91 = -t16 * pkin(3) - pkin(9) * t19;
t37 = -t101 * t80 + t53 * t95;
t20 = t104 * t102 + t37 * t50;
t21 = t37 * t102 - t104 * t50;
t90 = -t20 * pkin(3) + t21 * pkin(9);
t77 = t94 * t92;
t25 = -t103 * t102 + t50 * t77;
t26 = t102 * t77 + t103 * t50;
t89 = -t25 * pkin(3) + t26 * pkin(9);
t88 = t53 * t94;
t8 = -t106 * t52 + t21 * t49;
t87 = g(1) * t107 + g(2) * t8;
t86 = -t101 * pkin(1) + qJ(2) * t88;
t85 = -pkin(4) * t52 - pkin(10) * t49;
t84 = -g(1) * t16 + g(2) * t20;
t35 = -t95 * t78 + t97 * t96;
t14 = -t26 * t49 + t35 * t52;
t73 = g(1) * t8 - g(2) * t107 - g(3) * t14;
t15 = t26 * t52 + t35 * t49;
t9 = t106 * t49 + t21 * t52;
t72 = g(1) * t9 - g(2) * t7 + g(3) * t15;
t71 = g(1) * t20 + g(2) * t16 + g(3) * t25;
t70 = g(1) * t21 - g(2) * t19 + g(3) * t26;
t58 = -t36 * pkin(2) - t28 * pkin(8) + t86;
t57 = t37 * pkin(2) + t106 * pkin(8) + t98;
t56 = t19 * pkin(3) - pkin(9) * t16 + t58;
t55 = t21 * pkin(3) + t20 * pkin(9) + t57;
t32 = -g(1) * t82 + g(2) * t88 - g(3) * t97;
t3 = t20 * t48 + t9 * t51;
t2 = t20 * t51 - t9 * t48;
t1 = t71 * t49;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t101 - g(2) * t53, g(1) * t53 + g(2) * t101, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t37, -g(1) * t67 + g(2) * t61, -g(1) * t88 - g(2) * t82, -g(1) * t86 - g(2) * t98, 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, t84, g(1) * t28 - g(2) * t106, -g(1) * t58 - g(2) * t57, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t87, -t84, -g(1) * t56 - g(2) * t55, 0, 0, 0, 0, 0, 0, -g(1) * t111 - g(2) * t3, g(1) * t112 - g(2) * t2, -t87, -g(1) * (t7 * pkin(4) + pkin(10) * t107 + t56) - g(2) * (t9 * pkin(4) + t8 * pkin(10) + t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t52, -t1, -t70, -g(1) * t90 - g(2) * t91 - g(3) * t89, 0, 0, 0, 0, 0, 0, -g(1) * (-t20 * t99 + t21 * t48) - g(2) * (-t16 * t99 - t19 * t48) - g(3) * (-t25 * t99 + t26 * t48), -g(1) * (t20 * t100 + t21 * t51) - g(2) * (t16 * t100 - t19 * t51) - g(3) * (t25 * t100 + t26 * t51), t1, -g(1) * (t85 * t20 + t90) - g(2) * (t85 * t16 + t91) - g(3) * (t85 * t25 + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t51, -t73 * t48, -t72, -g(1) * (-t8 * pkin(4) + t9 * pkin(10)) - g(2) * (pkin(4) * t107 - pkin(10) * t7) - g(3) * (t14 * pkin(4) + t15 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t112 - g(3) * (-t15 * t48 + t25 * t51), g(1) * t3 - g(2) * t111 - g(3) * (-t15 * t51 - t25 * t48), 0, 0;];
taug_reg = t4;
