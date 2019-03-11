% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t49 = sin(qJ(1));
t50 = cos(qJ(3));
t51 = cos(qJ(2));
t47 = sin(qJ(3));
t52 = cos(qJ(1));
t86 = t52 * t47;
t21 = -t49 * t50 + t51 * t86;
t85 = t52 * t50;
t22 = t49 * t47 + t51 * t85;
t43 = pkin(10) + qJ(6);
t36 = sin(t43);
t37 = cos(t43);
t1 = t21 * t37 - t22 * t36;
t48 = sin(qJ(2));
t100 = g(3) * t48;
t89 = t49 * t51;
t19 = t47 * t89 + t85;
t88 = t50 * t51;
t20 = t49 * t88 - t86;
t113 = t19 * t37 - t20 * t36;
t63 = t36 * t50 - t37 * t47;
t114 = g(1) * t1 + g(2) * t113 - t63 * t100;
t24 = g(1) * t52 + g(2) * t49;
t111 = t24 * t48;
t109 = g(3) * t51 - t111;
t2 = t21 * t36 + t22 * t37;
t62 = t36 * t47 + t37 * t50;
t66 = t19 * t36 + t20 * t37;
t110 = g(1) * t2 + g(2) * t66 + t62 * t100;
t104 = -pkin(3) - pkin(4);
t103 = g(1) * t49;
t38 = t48 * pkin(8);
t40 = t51 * pkin(2);
t45 = cos(pkin(10));
t33 = t45 * pkin(5) + pkin(4);
t98 = -pkin(3) - t33;
t44 = sin(pkin(10));
t97 = t19 * t44;
t95 = t44 * t47;
t94 = t44 * t50;
t93 = t47 * t48;
t92 = t48 * t49;
t91 = t48 * t50;
t90 = t48 * t52;
t87 = t51 * t52;
t84 = t40 + t38;
t83 = t52 * pkin(1) + t49 * pkin(7);
t82 = qJ(4) * t47;
t81 = qJ(5) * t52;
t80 = -pkin(1) - t40;
t79 = -pkin(2) - t82;
t78 = pkin(5) * t44 + qJ(4);
t76 = -t19 * t45 + t20 * t44;
t75 = t21 * t44 + t22 * t45;
t15 = t19 * pkin(3);
t74 = t20 * qJ(4) - t15;
t17 = t21 * pkin(3);
t73 = t22 * qJ(4) - t17;
t72 = pkin(3) * t88 + t51 * t82 + t84;
t71 = pkin(2) * t87 + pkin(8) * t90 + t83;
t41 = t52 * pkin(7);
t70 = -t20 * pkin(3) - t19 * qJ(4) + t41;
t69 = t22 * pkin(3) + t71;
t68 = g(1) * t19 - g(2) * t21;
t67 = -g(2) * t52 + t103;
t65 = t20 * t45 + t97;
t64 = t21 * t45 - t22 * t44;
t61 = -t45 * t47 + t94;
t60 = t45 * t50 + t95;
t58 = g(3) * t60;
t56 = t21 * qJ(4) + t69;
t55 = (t80 - t38) * t103;
t4 = g(1) * t21 + g(2) * t19 + g(3) * t93;
t54 = g(1) * t22 + g(2) * t20 + g(3) * t91;
t46 = -pkin(9) - qJ(5);
t30 = pkin(8) * t87;
t27 = pkin(8) * t89;
t25 = qJ(4) * t91;
t23 = g(1) * t92 - g(2) * t90;
t14 = t24 * t51 + t100;
t7 = t109 * t50;
t6 = t109 * t47;
t5 = g(1) * t20 - g(2) * t22;
t3 = [0, 0, 0, 0, 0, 0, t67, t24, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t51, -t23, -t24, -g(1) * (-t49 * pkin(1) + t41) - g(2) * t83, 0, 0, 0, 0, 0, 0, t5, -t68, t23, -g(1) * t41 - g(2) * t71 - t55, 0, 0, 0, 0, 0, 0, t5, t23, t68, -g(1) * t70 - g(2) * t56 - t55, 0, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t75, -g(1) * t76 - g(2) * t64, -t23, -g(1) * (-t20 * pkin(4) + t70) - g(2) * (t22 * pkin(4) - t48 * t81 + t56) - ((-pkin(8) + qJ(5)) * t48 + t80) * t103, 0, 0, 0, 0, 0, 0, g(1) * t66 - g(2) * t2, g(1) * t113 - g(2) * t1, -t23, -g(1) * (-pkin(5) * t97 - t20 * t33 + t70) - g(2) * (t78 * t21 + t22 * t33 + t46 * t90 + t69) - ((-pkin(8) - t46) * t48 + t80) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t14, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t14, -g(1) * (-pkin(2) * t90 + t30) - g(2) * (-pkin(2) * t92 + t27) - g(3) * t84, 0, 0, 0, 0, 0, 0, -t7, -t14, -t6, -g(1) * t30 - g(2) * t27 - g(3) * t72 + (pkin(3) * t50 - t79) * t111, 0, 0, 0, 0, 0, 0, t60 * t111 - t51 * t58, t109 * t61, t14, -g(1) * (-t51 * t81 + t30) - g(2) * (-qJ(5) * t89 + t27) - g(3) * (pkin(4) * t88 + t72) + (g(3) * qJ(5) + t24 * (-t104 * t50 - t79)) * t48, 0, 0, 0, 0, 0, 0, -t109 * t62, t109 * t63, t14, -g(1) * (t46 * t87 + t30) - g(2) * (t46 * t89 + t27) - g(3) * (t51 * pkin(5) * t95 + t33 * t88 + t72) + (-g(3) * t46 + t24 * (t78 * t47 - t98 * t50 + pkin(2))) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t54, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, -t54, -g(1) * t73 - g(2) * t74 - g(3) * (-pkin(3) * t93 + t25) 0, 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t76 - t61 * t100, -g(1) * t75 - g(2) * t65 - t48 * t58, 0, -g(1) * (-t21 * pkin(4) + t73) - g(2) * (-t19 * pkin(4) + t74) - g(3) * (t104 * t93 + t25) 0, 0, 0, 0, 0, 0, t114, -t110, 0, -g(1) * (-t21 * t33 + t78 * t22 - t17) - g(2) * (-t19 * t33 + t78 * t20 - t15) - g(3) * t25 - (pkin(5) * t94 + t98 * t47) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t110, 0, 0;];
taug_reg  = t3;
