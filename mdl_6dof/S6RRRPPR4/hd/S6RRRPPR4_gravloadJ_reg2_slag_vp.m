% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:46:09
% EndTime: 2019-05-07 04:46:12
% DurationCPUTime: 0.69s
% Computational Cost: add. (499->123), mult. (742->166), div. (0->0), fcn. (805->10), ass. (0->79)
t46 = qJ(3) + pkin(10);
t40 = sin(t46);
t41 = cos(t46);
t55 = cos(qJ(1));
t86 = t55 * t41;
t51 = sin(qJ(1));
t54 = cos(qJ(2));
t88 = t51 * t54;
t15 = t40 * t88 + t86;
t87 = t55 * t40;
t16 = t41 * t88 - t87;
t48 = sin(qJ(6));
t52 = cos(qJ(6));
t113 = t15 * t52 - t16 * t48;
t17 = -t51 * t41 + t54 * t87;
t18 = t51 * t40 + t54 * t86;
t3 = t17 * t52 - t18 * t48;
t65 = t40 * t52 - t41 * t48;
t50 = sin(qJ(2));
t97 = g(3) * t50;
t115 = g(1) * t3 + g(2) * t113 + t65 * t97;
t99 = g(2) * t51;
t29 = g(1) * t55 + t99;
t114 = -g(3) * t54 + t29 * t50;
t111 = -t16 * pkin(4) - t15 * qJ(5);
t4 = t17 * t48 + t18 * t52;
t64 = t40 * t48 + t41 * t52;
t67 = t15 * t48 + t16 * t52;
t109 = g(1) * t4 + g(2) * t67 + t64 * t97;
t104 = -pkin(4) - pkin(5);
t49 = sin(qJ(3));
t103 = pkin(3) * t49;
t102 = pkin(4) * t41;
t101 = g(1) * t51;
t47 = -qJ(4) - pkin(8);
t95 = pkin(9) + t47;
t93 = t41 * t50;
t92 = t50 * t51;
t91 = t50 * t55;
t90 = t51 * t49;
t53 = cos(qJ(3));
t89 = t51 * t53;
t39 = t53 * pkin(3) + pkin(2);
t31 = t54 * t39;
t85 = t55 * t49;
t84 = t55 * t53;
t83 = t55 * pkin(1) + t51 * pkin(7);
t82 = qJ(5) * t40;
t80 = t47 * t91;
t79 = t54 * t85;
t43 = t55 * pkin(7);
t78 = pkin(3) * t85 + t47 * t92 + t43;
t77 = t95 * t55;
t76 = -pkin(1) - t31;
t74 = -t39 - t82;
t73 = pkin(3) * t90 + t55 * t31 + t83;
t72 = g(3) * (t31 + (t102 + t82) * t54);
t71 = t54 * pkin(2) + t50 * pkin(8);
t69 = g(1) * t15 - g(2) * t17;
t68 = -g(2) * t55 + t101;
t21 = t49 * t88 + t84;
t62 = t29 * t54;
t60 = t18 * pkin(4) + t17 * qJ(5) + t73;
t2 = g(1) * t17 + g(2) * t15 + t40 * t97;
t59 = g(1) * t18 + g(2) * t16 + g(3) * t93;
t58 = t76 * t51 + t78;
t36 = pkin(3) * t89;
t57 = -pkin(3) * t79 - t17 * pkin(4) + t18 * qJ(5) + t36;
t20 = t62 + t97;
t56 = -t21 * pkin(3) - t15 * pkin(4) + t16 * qJ(5);
t27 = qJ(5) * t93;
t25 = g(1) * t92 - g(2) * t91;
t24 = t54 * t84 + t90;
t23 = -t79 + t89;
t22 = -t53 * t88 + t85;
t7 = t114 * t41;
t6 = t114 * t40;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t68, t29, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t54, -t25, -t29, -g(1) * (-t51 * pkin(1) + t43) - g(2) * t83, 0, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, -g(1) * t21 - g(2) * t23, t25, -g(1) * t43 - g(2) * (t71 * t55 + t83) - (-pkin(1) - t71) * t101, 0, 0, 0, 0, 0, 0, t5, -t69, t25, -g(1) * t58 - g(2) * (t73 - t80) 0, 0, 0, 0, 0, 0, t5, t25, t69, -g(1) * (t58 + t111) - g(2) * (t60 - t80) 0, 0, 0, 0, 0, 0, g(1) * t67 - g(2) * t4, g(1) * t113 - g(2) * t3, -t25, -g(1) * (-t16 * pkin(5) + t111 + t78) - g(2) * (t18 * pkin(5) - t50 * t77 + t60) - (t50 * pkin(9) + t76) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t20, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t53, -t114 * t49, -t20, -g(3) * t71 + t29 * (pkin(2) * t50 - pkin(8) * t54) 0, 0, 0, 0, 0, 0, t7, -t6, -t20, -g(3) * (-t50 * t47 + t31) + t29 * (t39 * t50 + t47 * t54) 0, 0, 0, 0, 0, 0, t7, -t20, t6, -t72 + t47 * t62 + (g(3) * t47 + t29 * (-t74 + t102)) * t50, 0, 0, 0, 0, 0, 0, t114 * t64, t114 * t65, t20, -t72 + (-g(3) * pkin(5) * t41 + g(1) * t77 + t95 * t99) * t54 + (g(3) * t95 + t29 * (-t104 * t41 - t74)) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 + g(2) * t21 + t49 * t97, g(1) * t24 - g(2) * t22 + t53 * t97, 0, 0, 0, 0, 0, 0, 0, 0, t2, t59, 0, -g(1) * t36 + (g(2) * t84 + t20 * t49) * pkin(3), 0, 0, 0, 0, 0, 0, t2, 0, -t59, -g(1) * t57 - g(2) * t56 - g(3) * (t27 + (-pkin(4) * t40 - t103) * t50) 0, 0, 0, 0, 0, 0, t115, -t109, 0, -g(1) * (-t17 * pkin(5) + t57) - g(2) * (-t15 * pkin(5) + t56) - g(3) * t27 - (t104 * t40 - t103) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t109, 0, 0;];
taug_reg  = t1;
