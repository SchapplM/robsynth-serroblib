% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t74 = sin(qJ(4));
t128 = pkin(4) * t74;
t78 = cos(qJ(4));
t129 = -pkin(10) * t78 + qJ(3) + t128;
t113 = cos(pkin(6));
t80 = cos(qJ(1));
t104 = t80 * t113;
t75 = sin(qJ(2));
t76 = sin(qJ(1));
t79 = cos(qJ(2));
t54 = t75 * t104 + t76 * t79;
t73 = sin(qJ(5));
t77 = cos(qJ(5));
t72 = sin(pkin(6));
t119 = t72 * t80;
t53 = -t79 * t104 + t76 * t75;
t95 = t78 * t119 - t53 * t74;
t10 = t54 * t77 + t73 * t95;
t11 = -t54 * t73 + t77 * t95;
t127 = t53 * pkin(9);
t105 = t76 * t113;
t55 = t79 * t105 + t80 * t75;
t126 = t55 * pkin(9);
t123 = t72 * t75;
t122 = t72 * t76;
t121 = t72 * t78;
t120 = t72 * t79;
t118 = t73 * t74;
t117 = t74 * t77;
t116 = t75 * t77;
t115 = pkin(2) * t120 + qJ(3) * t123;
t114 = t80 * pkin(1) + pkin(8) * t122;
t112 = t73 * t123;
t111 = pkin(9) * t120 + t115;
t110 = -t76 * pkin(1) + pkin(8) * t119;
t28 = t74 * t122 - t55 * t78;
t29 = t76 * t121 + t55 * t74;
t109 = -t28 * pkin(4) + t29 * pkin(10);
t94 = t74 * t119 + t53 * t78;
t108 = pkin(4) * t94 - pkin(10) * t95;
t51 = -t113 * t74 - t78 * t120;
t52 = t113 * t78 - t74 * t120;
t107 = t51 * pkin(4) + t52 * pkin(10);
t47 = t53 * pkin(2);
t103 = t54 * qJ(3) - t47;
t49 = t55 * pkin(2);
t56 = -t75 * t105 + t80 * t79;
t102 = t56 * qJ(3) - t49;
t8 = t29 * t73 - t56 * t77;
t101 = g(1) * t10 + g(2) * t8;
t100 = g(1) * t94 + g(2) * t28;
t99 = g(1) * t53 - g(2) * t55;
t21 = g(1) * t54 - g(2) * t56;
t98 = g(1) * t80 + g(2) * t76;
t97 = pkin(5) * t77 + qJ(6) * t73;
t96 = t56 * pkin(2) + t55 * qJ(3) + t114;
t26 = -t72 * t116 + t52 * t73;
t1 = g(1) * t8 - g(2) * t10 + g(3) * t26;
t27 = t52 * t77 + t112;
t9 = t29 * t77 + t56 * t73;
t93 = g(1) * t9 - g(2) * t11 + g(3) * t27;
t92 = -t54 * pkin(2) - t53 * qJ(3) + t110;
t17 = t54 * t118 + t53 * t77;
t19 = t56 * t118 + t55 * t77;
t34 = t74 * t112 - t77 * t120;
t91 = g(1) * t19 + g(2) * t17 + g(3) * t34;
t90 = g(1) * t28 - g(2) * t94 - g(3) * t51;
t89 = g(1) * t29 - g(2) * t95 + g(3) * t52;
t15 = -g(1) * t55 - g(2) * t53 + g(3) * t120;
t88 = g(1) * t56 + g(2) * t54 + g(3) * t123;
t87 = -t75 * pkin(10) * t121 + t123 * t128 + t111;
t86 = pkin(3) * t122 + t56 * pkin(9) + t96;
t85 = t129 * t54 - t127 - t47;
t84 = t129 * t56 - t126 - t49;
t83 = pkin(3) * t119 - t54 * pkin(9) + t92;
t82 = t29 * pkin(4) + t28 * pkin(10) + t86;
t81 = pkin(4) * t95 + pkin(10) * t94 + t83;
t57 = t98 * t72;
t35 = (t74 * t116 + t73 * t79) * t72;
t20 = t56 * t117 - t55 * t73;
t18 = t54 * t117 - t53 * t73;
t14 = t88 * t78;
t5 = t90 * t77;
t4 = t90 * t73;
t3 = -g(1) * t11 - g(2) * t9;
t2 = -g(1) * t20 - g(2) * t18 - g(3) * t35;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t76 - g(2) * t80, t98, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t99, -t57, -g(1) * t110 - g(2) * t114, 0, 0, 0, 0, 0, 0, -t57, -t21, t99, -g(1) * t92 - g(2) * t96, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t29, t100, t21, -g(1) * t83 - g(2) * t86, 0, 0, 0, 0, 0, 0, t3, t101, -t100, -g(1) * t81 - g(2) * t82, 0, 0, 0, 0, 0, 0, t3, -t100, -t101, -g(1) * (t11 * pkin(5) + t10 * qJ(6) + t81) - g(2) * (t9 * pkin(5) + t8 * qJ(6) + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t88, -g(1) * t102 - g(2) * t103 - g(3) * t115, 0, 0, 0, 0, 0, 0, -t88 * t74, -t14, -t15, -g(1) * (t102 - t126) - g(2) * (t103 - t127) - g(3) * t111, 0, 0, 0, 0, 0, 0, t2, t91, t14, -g(1) * t84 - g(2) * t85 - g(3) * t87, 0, 0, 0, 0, 0, 0, t2, t14, -t91, -g(1) * (t20 * pkin(5) + t19 * qJ(6) + t84) - g(2) * (t18 * pkin(5) + t17 * qJ(6) + t85) - g(3) * (t35 * pkin(5) + t34 * qJ(6) + t87); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t89, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t89, -g(1) * t109 - g(2) * t108 - g(3) * t107, 0, 0, 0, 0, 0, 0, t5, -t89, t4, -g(1) * (-t97 * t28 + t109) - g(2) * (t94 * t97 + t108) - g(3) * (t97 * t51 + t107); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t93, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t93, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (pkin(5) * t10 - qJ(6) * t11) - g(3) * (-t26 * pkin(5) + t27 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
