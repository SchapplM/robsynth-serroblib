% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t27 = g(1) * t54 + g(2) * t51;
t50 = sin(qJ(2));
t107 = t27 * t50;
t53 = cos(qJ(2));
t55 = -pkin(9) - pkin(8);
t49 = sin(qJ(4));
t91 = t51 * t49;
t77 = pkin(4) * t91;
t94 = t50 * t51;
t106 = t53 * t77 + t55 * t94;
t85 = t54 * t49;
t37 = pkin(4) * t85;
t83 = t54 * t55;
t105 = t53 * t37 + t50 * t83;
t14 = g(3) * t50 + t27 * t53;
t103 = pkin(2) * t50;
t102 = pkin(4) * t49;
t101 = g(1) * t51;
t97 = g(3) * t53;
t52 = cos(qJ(4));
t96 = t52 * pkin(4);
t44 = t53 * pkin(2);
t95 = t53 * pkin(8);
t48 = qJ(4) + qJ(5);
t40 = sin(t48);
t93 = t51 * t40;
t41 = cos(t48);
t92 = t51 * t41;
t90 = t51 * t52;
t89 = t53 * t54;
t88 = t53 * t55;
t87 = t54 * t40;
t86 = t54 * t41;
t84 = t54 * t52;
t76 = t50 * t90;
t82 = pkin(4) * t76 + t37;
t42 = t50 * qJ(3);
t81 = t44 + t42;
t80 = t54 * pkin(1) + t51 * pkin(7);
t79 = qJ(3) * t53;
t78 = t52 * t97;
t75 = t50 * t85;
t74 = t50 * t84;
t39 = pkin(3) + t96;
t45 = t54 * pkin(7);
t73 = t54 * t39 + t51 * t88 + t45;
t72 = -pkin(1) - t44;
t10 = t50 * t87 + t92;
t9 = -t50 * t86 + t93;
t71 = -t9 * pkin(5) + t10 * qJ(6);
t11 = t50 * t92 + t87;
t12 = -t50 * t93 + t86;
t70 = t11 * pkin(5) - t12 * qJ(6);
t69 = pkin(2) * t89 + t54 * t42 + t80;
t33 = t51 * t79;
t68 = -pkin(2) * t94 + t33;
t35 = t54 * t79;
t67 = -t54 * t103 + t35;
t66 = pkin(4) * t74 - t77;
t65 = g(1) * t11 + g(2) * t9;
t64 = -g(2) * t54 + t101;
t63 = -pkin(5) * t41 - qJ(6) * t40;
t62 = pkin(5) * t40 - qJ(6) * t41;
t61 = t50 * t102 + t81 - t88;
t60 = t72 - t42;
t1 = g(1) * t9 - g(2) * t11 + t41 * t97;
t2 = -g(1) * t10 + g(2) * t12 + t40 * t97;
t59 = pkin(4) * t75 + t51 * t39 - t53 * t83 + t69;
t56 = ((-qJ(3) - t102) * t50 + t72) * t101;
t20 = t64 * t53;
t19 = t64 * t50;
t18 = -t50 * t91 + t84;
t17 = t76 + t85;
t16 = t75 + t90;
t15 = t74 - t91;
t13 = -t97 + t107;
t6 = t14 * t41;
t5 = t14 * t40;
t4 = -g(1) * t12 - g(2) * t10;
t3 = [0, 0, 0, 0, 0, 0, t64, t27, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, -t27, -g(1) * (-t51 * pkin(1) + t45) - g(2) * t80, 0, 0, 0, 0, 0, 0, -t27, -t20, t19, -g(1) * t45 - g(2) * t69 - t60 * t101, 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16, g(1) * t17 - g(2) * t15, t20, -g(1) * (t54 * pkin(3) + t45) - g(2) * (pkin(8) * t89 + t69) + (-g(1) * (t60 - t95) - g(2) * pkin(3)) * t51, 0, 0, 0, 0, 0, 0, t4, t65, t20, -g(1) * t73 - g(2) * t59 - t56, 0, 0, 0, 0, 0, 0, t4, t20, -t65, -g(1) * (t12 * pkin(5) + t11 * qJ(6) + t73) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t59) - t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -g(1) * t67 - g(2) * t68 - g(3) * t81, 0, 0, 0, 0, 0, 0, -t14 * t49, -t14 * t52, t13, -g(1) * t35 - g(2) * t33 - g(3) * (t81 + t95) + (pkin(2) + pkin(8)) * t107, 0, 0, 0, 0, 0, 0, -t5, -t6, t13, -g(1) * (t67 + t105) - g(2) * (t68 + t106) - g(3) * t61, 0, 0, 0, 0, 0, 0, -t5, t13, t6, -g(1) * (t35 + t105) - g(2) * (t33 + t106) - g(3) * (t50 * t62 + t61) + t27 * (-t62 * t53 + t103); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17 + t78, g(1) * t16 - g(2) * t18 - t49 * t97, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, pkin(4) * t78 - g(1) * t66 - g(2) * t82, 0, 0, 0, 0, 0, 0, t1, 0, t2, -g(1) * (t66 + t71) - g(2) * (t70 + t82) - (t63 - t96) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, t2, -g(1) * t71 - g(2) * t70 - t63 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t3;
