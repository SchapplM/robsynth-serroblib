% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t108 = sin(pkin(7));
t120 = cos(qJ(1));
t107 = sin(pkin(13));
t118 = sin(qJ(1));
t110 = cos(pkin(13));
t112 = cos(pkin(6));
t94 = t112 * t110;
t79 = t118 * t107 - t120 * t94;
t109 = sin(pkin(6));
t111 = cos(pkin(7));
t91 = t111 * t109;
t125 = t79 * t108 - t120 * t91;
t119 = cos(qJ(3));
t90 = t109 * t108;
t131 = t79 * t111 + t120 * t90;
t92 = t112 * t107;
t41 = t118 * t110 + t120 * t92;
t59 = sin(qJ(3));
t25 = -t41 * t119 + t131 * t59;
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t13 = -t125 * t58 + t25 * t61;
t22 = t131 * t119 + t41 * t59;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t139 = t13 * t57 + t22 * t60;
t138 = t13 * t60 - t22 * t57;
t56 = qJ(5) + qJ(6);
t53 = sin(t56);
t54 = cos(t56);
t137 = t13 * t53 + t22 * t54;
t136 = t13 * t54 - t22 * t53;
t130 = t125 * t61 + t25 * t58;
t74 = t120 * t107 + t118 * t94;
t129 = t74 * t108 + t118 * t91;
t126 = t74 * t111 - t118 * t90;
t124 = t112 * t108 + t110 * t91;
t89 = t109 * t107;
t32 = t119 * t89 + t124 * t59;
t40 = -t110 * t90 + t112 * t111;
t21 = t32 * t61 + t40 * t58;
t31 = -t124 * t119 + t59 * t89;
t42 = t120 * t110 - t118 * t92;
t27 = t42 * t119 - t126 * t59;
t15 = t129 * t58 + t27 * t61;
t26 = t126 * t119 + t42 * t59;
t8 = -t15 * t57 + t26 * t60;
t122 = -g(1) * t8 - g(2) * t139 - g(3) * (-t21 * t57 + t31 * t60);
t117 = t53 * t61;
t116 = t54 * t61;
t115 = t57 * t61;
t114 = t60 * t61;
t96 = t109 * t118;
t113 = t120 * pkin(1) + qJ(2) * t96;
t106 = t57 * pkin(5) + pkin(10);
t16 = t22 * pkin(3);
t105 = -pkin(10) * t25 - t16;
t18 = t26 * pkin(3);
t104 = t27 * pkin(10) - t18;
t30 = t31 * pkin(3);
t103 = t32 * pkin(10) - t30;
t97 = t120 * t109;
t102 = -t118 * pkin(1) + qJ(2) * t97;
t101 = -pkin(4) * t61 - pkin(11) * t58;
t14 = -t129 * t61 + t27 * t58;
t100 = g(1) * t130 + g(2) * t14;
t99 = -g(1) * t22 + g(2) * t26;
t52 = t60 * pkin(5) + pkin(4);
t62 = -pkin(12) - pkin(11);
t95 = -t52 * t61 + t58 * t62;
t20 = -t32 * t58 + t40 * t61;
t88 = g(1) * t14 - g(2) * t130 - g(3) * t20;
t87 = g(1) * t15 - g(2) * t13 + g(3) * t21;
t86 = g(1) * t26 + g(2) * t22 + g(3) * t31;
t85 = g(1) * t27 - g(2) * t25 + g(3) * t32;
t70 = -t41 * pkin(2) - t125 * pkin(9) + t102;
t68 = t25 * pkin(3) + t70;
t67 = t42 * pkin(2) + t129 * pkin(9) + t113;
t66 = t27 * pkin(3) + t67;
t65 = -pkin(10) * t22 + t68;
t64 = t26 * pkin(10) + t66;
t37 = -g(1) * t96 + g(2) * t97 - g(3) * t112;
t9 = t15 * t60 + t26 * t57;
t7 = t86 * t58;
t6 = t15 * t54 + t26 * t53;
t5 = -t15 * t53 + t26 * t54;
t2 = g(1) * t6 - g(2) * t136 - g(3) * (-t21 * t54 - t31 * t53);
t1 = -g(1) * t5 - g(2) * t137 - g(3) * (-t21 * t53 + t31 * t54);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t118 - g(2) * t120, g(1) * t120 + g(2) * t118, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t42, -g(1) * t79 + g(2) * t74, -g(1) * t97 - g(2) * t96, -g(1) * t102 - g(2) * t113, 0, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t27, t99, g(1) * t125 - g(2) * t129, -g(1) * t70 - g(2) * t67, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, t100, -t99, -g(1) * t65 - g(2) * t64, 0, 0, 0, 0, 0, 0, -g(1) * t138 - g(2) * t9, g(1) * t139 - g(2) * t8, -t100, -g(1) * (t13 * pkin(4) + pkin(11) * t130 + t65) - g(2) * (t15 * pkin(4) + t14 * pkin(11) + t64) 0, 0, 0, 0, 0, 0, -g(1) * t136 - g(2) * t6, g(1) * t137 - g(2) * t5, -t100, -g(1) * (-t106 * t22 + t13 * t52 - t130 * t62 + t68) - g(2) * (t106 * t26 - t14 * t62 + t15 * t52 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t85, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t61, -t7, -t85, -g(1) * t104 - g(2) * t105 - g(3) * t103, 0, 0, 0, 0, 0, 0, -g(1) * (-t26 * t114 + t27 * t57) - g(2) * (-t22 * t114 - t25 * t57) - g(3) * (-t31 * t114 + t32 * t57) -g(1) * (t26 * t115 + t27 * t60) - g(2) * (t22 * t115 - t25 * t60) - g(3) * (t31 * t115 + t32 * t60) t7, -g(1) * (t101 * t26 + t104) - g(2) * (t101 * t22 + t105) - g(3) * (t101 * t31 + t103) 0, 0, 0, 0, 0, 0, -g(1) * (-t26 * t116 + t27 * t53) - g(2) * (-t22 * t116 - t25 * t53) - g(3) * (-t31 * t116 + t32 * t53) -g(1) * (t26 * t117 + t27 * t54) - g(2) * (t22 * t117 - t25 * t54) - g(3) * (t31 * t117 + t32 * t54) t7, -g(1) * (t106 * t27 + t95 * t26 - t18) - g(2) * (-t106 * t25 + t95 * t22 - t16) - g(3) * (t106 * t32 + t95 * t31 - t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t87, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t60, -t88 * t57, -t87, -g(1) * (-pkin(4) * t14 + pkin(11) * t15) - g(2) * (pkin(4) * t130 - pkin(11) * t13) - g(3) * (pkin(4) * t20 + pkin(11) * t21) 0, 0, 0, 0, 0, 0, t88 * t54, -t88 * t53, -t87, -g(1) * (-t14 * t52 - t15 * t62) - g(2) * (t13 * t62 + t130 * t52) - g(3) * (t20 * t52 - t21 * t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, g(1) * t9 - g(2) * t138 - g(3) * (-t21 * t60 - t31 * t57) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t122 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
