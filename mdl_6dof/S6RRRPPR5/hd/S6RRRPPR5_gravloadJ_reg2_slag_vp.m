% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t86 = cos(pkin(6));
t55 = sin(pkin(6));
t60 = sin(qJ(2));
t94 = t55 * t60;
t112 = -t59 * t94 + t86 * t62;
t107 = cos(qJ(1));
t63 = cos(qJ(2));
t61 = sin(qJ(1));
t79 = t61 * t86;
t31 = t107 * t63 - t60 * t79;
t92 = t55 * t62;
t16 = -t31 * t59 + t61 * t92;
t74 = t86 * t107;
t29 = t60 * t74 + t61 * t63;
t53 = qJ(3) + pkin(11);
t48 = sin(t53);
t50 = cos(t53);
t82 = t55 * t107;
t11 = t29 * t50 - t48 * t82;
t28 = t61 * t60 - t63 * t74;
t52 = pkin(12) + qJ(6);
t47 = sin(t52);
t49 = cos(t52);
t111 = t11 * t47 - t28 * t49;
t110 = t11 * t49 + t28 * t47;
t46 = t62 * pkin(3) + pkin(2);
t91 = t55 * t63;
t32 = t46 * t91;
t109 = g(3) * t32;
t108 = g(3) * t55;
t54 = sin(pkin(12));
t104 = t28 * t54;
t103 = t29 * t54;
t30 = t107 * t60 + t63 * t79;
t102 = t30 * t54;
t101 = t31 * t54;
t99 = t47 * t50;
t98 = t49 * t50;
t97 = t50 * t54;
t56 = cos(pkin(12));
t96 = t50 * t56;
t95 = t50 * t63;
t93 = t55 * t61;
t57 = -qJ(4) - pkin(9);
t90 = t57 * t60;
t89 = -t28 * t46 - t29 * t57;
t88 = -t30 * t46 - t31 * t57;
t87 = t107 * pkin(1) + pkin(8) * t93;
t84 = t59 * t93;
t81 = -t61 * pkin(1) + pkin(8) * t82;
t39 = t59 * t82;
t80 = t29 * t62 - t39;
t77 = t16 * pkin(3);
t76 = pkin(3) * t84 - t30 * t57 + t31 * t46 + t87;
t10 = t29 * t48 + t50 * t82;
t14 = t31 * t48 - t50 * t93;
t75 = -g(1) * t10 + g(2) * t14;
t9 = g(1) * t28 - g(2) * t30;
t73 = pkin(4) * t50 + qJ(5) * t48;
t45 = t56 * pkin(5) + pkin(4);
t58 = -pkin(10) - qJ(5);
t72 = t45 * t50 - t48 * t58;
t71 = t112 * pkin(3);
t70 = g(1) * t107 + g(2) * t61;
t69 = pkin(3) * t39 + t28 * t57 - t29 * t46 + t81;
t22 = t48 * t94 - t50 * t86;
t68 = g(1) * t14 + g(2) * t10 + g(3) * t22;
t15 = t31 * t50 + t48 * t93;
t23 = t48 * t86 + t50 * t94;
t67 = g(1) * t15 + g(2) * t11 + g(3) * t23;
t66 = t29 * t59 + t62 * t82;
t7 = -g(1) * t30 - g(2) * t28 + g(3) * t91;
t65 = g(1) * t31 + g(2) * t29 + g(3) * t94;
t64 = t66 * pkin(3);
t17 = t31 * t62 + t84;
t6 = t7 * t48;
t5 = t15 * t49 + t30 * t47;
t4 = -t15 * t47 + t30 * t49;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t61 - g(2) * t107, t70, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t31, -t9, -t70 * t55, -g(1) * t81 - g(2) * t87, 0, 0, 0, 0, 0, 0, g(1) * t80 - g(2) * t17, -g(1) * t66 - g(2) * t16, t9, -g(1) * (-t29 * pkin(2) - t28 * pkin(9) + t81) - g(2) * (t31 * pkin(2) + t30 * pkin(9) + t87) 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t15, t75, t9, -g(1) * t69 - g(2) * t76, 0, 0, 0, 0, 0, 0, -g(1) * (-t11 * t56 - t104) - g(2) * (t15 * t56 + t102) -g(1) * (t11 * t54 - t28 * t56) - g(2) * (-t15 * t54 + t30 * t56) -t75, -g(1) * (-pkin(4) * t11 - qJ(5) * t10 + t69) - g(2) * (t15 * pkin(4) + t14 * qJ(5) + t76) 0, 0, 0, 0, 0, 0, g(1) * t110 - g(2) * t5, -g(1) * t111 - g(2) * t4, -t75, -g(1) * (-pkin(5) * t104 + t10 * t58 - t11 * t45 + t69) - g(2) * (pkin(5) * t102 - t14 * t58 + t15 * t45 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t65, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t62, t7 * t59, -t65, -g(1) * (-t30 * pkin(2) + t31 * pkin(9)) - g(2) * (-t28 * pkin(2) + t29 * pkin(9)) - (pkin(2) * t63 + pkin(9) * t60) * t108, 0, 0, 0, 0, 0, 0, -t7 * t50, t6, -t65, -g(1) * t88 - g(2) * t89 - g(3) * (-t55 * t90 + t32) 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t96 + t101) - g(2) * (-t28 * t96 + t103) - (t54 * t60 + t56 * t95) * t108, -g(1) * (t30 * t97 + t31 * t56) - g(2) * (t28 * t97 + t29 * t56) - (-t54 * t95 + t56 * t60) * t108, -t6, -g(1) * (-t30 * t73 + t88) - g(2) * (-t28 * t73 + t89) - t109 - (t63 * t73 - t90) * t108, 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t98 + t31 * t47) - g(2) * (-t28 * t98 + t29 * t47) - (t47 * t60 + t49 * t95) * t108, -g(1) * (t30 * t99 + t31 * t49) - g(2) * (t28 * t99 + t29 * t49) - (-t47 * t95 + t49 * t60) * t108, -t6, -g(1) * (pkin(5) * t101 - t30 * t72 + t88) - g(2) * (pkin(5) * t103 - t28 * t72 + t89) - t109 - (t72 * t63 + (pkin(5) * t54 - t57) * t60) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t66 - g(3) * t112, g(1) * t17 + g(2) * t80 - g(3) * (-t59 * t86 - t60 * t92) 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, 0, -g(1) * t77 + g(2) * t64 - g(3) * t71, 0, 0, 0, 0, 0, 0, t68 * t56, -t68 * t54, -t67, -g(1) * (-t14 * pkin(4) + t15 * qJ(5) + t77) - g(2) * (-t10 * pkin(4) + t11 * qJ(5) - t64) - g(3) * (-t22 * pkin(4) + t23 * qJ(5) + t71) 0, 0, 0, 0, 0, 0, t68 * t49, -t68 * t47, -t67, -g(1) * (-t14 * t45 - t15 * t58 + t77) - g(2) * (-t10 * t45 - t11 * t58 - t64) - g(3) * (-t22 * t45 - t23 * t58 + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t111 - g(3) * (-t23 * t47 - t49 * t91) g(1) * t5 + g(2) * t110 - g(3) * (-t23 * t49 + t47 * t91) 0, 0;];
taug_reg  = t1;
