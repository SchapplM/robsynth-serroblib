% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:06:46
% EndTime: 2019-05-07 15:06:50
% DurationCPUTime: 1.18s
% Computational Cost: add. (781->187), mult. (1540->287), div. (0->0), fcn. (1881->14), ass. (0->85)
t53 = sin(qJ(2));
t100 = cos(qJ(1));
t83 = cos(pkin(6));
t69 = t83 * t100;
t98 = sin(qJ(1));
t99 = cos(qJ(2));
t24 = t53 * t69 + t98 * t99;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t49 = sin(pkin(6));
t79 = t49 * t100;
t12 = t24 * t54 - t52 * t79;
t23 = t53 * t98 - t69 * t99;
t47 = pkin(12) + qJ(5);
t43 = qJ(6) + t47;
t38 = sin(t43);
t39 = cos(t43);
t109 = t12 * t38 - t23 * t39;
t108 = t12 * t39 + t23 * t38;
t41 = sin(t47);
t42 = cos(t47);
t107 = t12 * t41 - t23 * t42;
t106 = t12 * t42 + t23 * t41;
t87 = t49 * t53;
t22 = t52 * t83 + t54 * t87;
t78 = t49 * t99;
t68 = t83 * t98;
t26 = t100 * t99 - t53 * t68;
t77 = t49 * t98;
t16 = t26 * t54 + t52 * t77;
t25 = t100 * t53 + t68 * t99;
t8 = -t16 * t41 + t25 * t42;
t105 = -g(3) * (-t22 * t41 - t42 * t78) + g(2) * t107 - g(1) * t8;
t103 = g(3) * t49;
t48 = sin(pkin(12));
t102 = t48 * pkin(4);
t28 = pkin(5) * t41 + t102;
t101 = pkin(9) + t28;
t50 = cos(pkin(12));
t40 = t50 * pkin(4) + pkin(3);
t93 = t38 * t54;
t92 = t39 * t54;
t91 = t41 * t54;
t90 = t42 * t54;
t89 = t48 * t53;
t88 = t48 * t54;
t86 = t50 * t54;
t51 = -pkin(10) - qJ(4);
t85 = pkin(2) * t78 + pkin(9) * t87;
t84 = t100 * pkin(1) + pkin(8) * t77;
t82 = t26 * pkin(2) + t84;
t81 = pkin(9) + t102;
t80 = g(3) * t85;
t76 = t52 * t99;
t75 = t54 * t99;
t17 = t23 * pkin(2);
t74 = t24 * pkin(9) - t17;
t19 = t25 * pkin(2);
t73 = t26 * pkin(9) - t19;
t72 = -pkin(1) * t98 + pkin(8) * t79;
t11 = t24 * t52 + t54 * t79;
t15 = t26 * t52 - t54 * t77;
t71 = -g(1) * t11 + g(2) * t15;
t70 = g(1) * t23 - g(2) * t25;
t67 = -pkin(3) * t54 - qJ(4) * t52;
t27 = pkin(5) * t42 + t40;
t46 = -pkin(11) + t51;
t66 = -t27 * t54 + t46 * t52;
t65 = -t40 * t54 + t51 * t52;
t64 = t25 * pkin(9) + t82;
t63 = -t24 * pkin(2) + t72;
t21 = t52 * t87 - t54 * t83;
t61 = g(1) * t15 + g(2) * t11 + g(3) * t21;
t60 = g(1) * t16 + g(2) * t12 + g(3) * t22;
t59 = g(1) * t100 + g(2) * t98;
t58 = -t23 * pkin(9) + t63;
t57 = g(1) * t26 + g(2) * t24 + g(3) * t87;
t55 = -g(1) * t25 - g(2) * t23 + g(3) * t78;
t10 = t55 * t52;
t9 = t16 * t42 + t25 * t41;
t7 = t16 * t39 + t25 * t38;
t6 = -t16 * t38 + t25 * t39;
t2 = g(1) * t7 + g(2) * t108 - g(3) * (-t22 * t39 + t38 * t78);
t1 = -g(1) * t6 + g(2) * t109 - g(3) * (-t22 * t38 - t39 * t78);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t98 - g(2) * t100, t59, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t24 - g(2) * t26, -t70, -t59 * t49, -g(1) * t72 - g(2) * t84, 0, 0, 0, 0, 0, 0, g(1) * t12 - g(2) * t16, t71, t70, -g(1) * t58 - g(2) * t64, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t50 - t23 * t48) - g(2) * (t16 * t50 + t25 * t48) -g(1) * (t12 * t48 - t23 * t50) - g(2) * (-t16 * t48 + t25 * t50) -t71, -g(1) * (-pkin(3) * t12 - qJ(4) * t11 + t58) - g(2) * (t16 * pkin(3) + t15 * qJ(4) + t64) 0, 0, 0, 0, 0, 0, g(1) * t106 - g(2) * t9, -g(1) * t107 - g(2) * t8, -t71, -g(1) * (t11 * t51 - t12 * t40 - t23 * t81 + t63) - g(2) * (-t15 * t51 + t16 * t40 + t25 * t81 + t82) 0, 0, 0, 0, 0, 0, g(1) * t108 - g(2) * t7, -g(1) * t109 - g(2) * t6, -t71, -g(1) * (-t101 * t23 + t11 * t46 - t12 * t27 + t63) - g(2) * (t101 * t25 - t15 * t46 + t16 * t27 + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t57, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t54, t10, -t57, -g(1) * t73 - g(2) * t74 - t80, 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t86 + t26 * t48) - g(2) * (-t23 * t86 + t24 * t48) - (t50 * t75 + t89) * t103, -g(1) * (t25 * t88 + t26 * t50) - g(2) * (t23 * t88 + t24 * t50) - (-t48 * t75 + t50 * t53) * t103, -t10, -g(1) * (t25 * t67 + t73) - g(2) * (t23 * t67 + t74) - g(3) * ((pkin(3) * t75 + qJ(4) * t76) * t49 + t85) 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t90 + t26 * t41) - g(2) * (-t23 * t90 + t24 * t41) - (t41 * t53 + t42 * t75) * t103, -g(1) * (t25 * t91 + t26 * t42) - g(2) * (t23 * t91 + t24 * t42) - (-t41 * t75 + t42 * t53) * t103, -t10, -g(1) * (t25 * t65 + t26 * t81 - t19) - g(2) * (t23 * t65 + t24 * t81 - t17) - t80 - (pkin(4) * t89 + t40 * t75 - t51 * t76) * t103, 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t92 + t26 * t38) - g(2) * (-t23 * t92 + t24 * t38) - (t38 * t53 + t39 * t75) * t103, -g(1) * (t25 * t93 + t26 * t39) - g(2) * (t23 * t93 + t24 * t39) - (-t38 * t75 + t39 * t53) * t103, -t10, -g(1) * (t101 * t26 + t25 * t66 - t19) - g(2) * (t101 * t24 + t23 * t66 - t17) - t80 - (t27 * t75 + t28 * t53 - t46 * t76) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t50, -t61 * t48, -t60, -g(1) * (-t15 * pkin(3) + t16 * qJ(4)) - g(2) * (-t11 * pkin(3) + t12 * qJ(4)) - g(3) * (-t21 * pkin(3) + t22 * qJ(4)) 0, 0, 0, 0, 0, 0, t61 * t42, -t61 * t41, -t60, -g(1) * (-t15 * t40 - t16 * t51) - g(2) * (-t11 * t40 - t12 * t51) - g(3) * (-t21 * t40 - t22 * t51) 0, 0, 0, 0, 0, 0, t61 * t39, -t61 * t38, -t60, -g(1) * (-t15 * t27 - t16 * t46) - g(2) * (-t11 * t27 - t12 * t46) - g(3) * (-t21 * t27 - t22 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, g(1) * t9 + g(2) * t106 - g(3) * (-t22 * t42 + t41 * t78) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t105 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
