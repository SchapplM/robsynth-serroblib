% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t160 = pkin(5) * t96 + qJ(6) * t92;
t136 = cos(pkin(6));
t91 = sin(pkin(6));
t94 = sin(qJ(2));
t145 = t91 * t94;
t93 = sin(qJ(3));
t97 = cos(qJ(3));
t159 = t136 * t97 - t93 * t145;
t143 = t91 * t97;
t95 = sin(qJ(1));
t124 = t95 * t136;
t151 = cos(qJ(1));
t98 = cos(qJ(2));
t68 = -t94 * t124 + t151 * t98;
t39 = t95 * t143 - t68 * t93;
t90 = qJ(3) + qJ(4);
t87 = sin(t90);
t88 = cos(t90);
t158 = pkin(4) * t88 + pkin(11) * t87;
t131 = t91 * t151;
t115 = t136 * t151;
t66 = t94 * t115 + t95 * t98;
t34 = -t87 * t131 + t66 * t88;
t65 = -t98 * t115 + t95 * t94;
t11 = t34 * t92 - t65 * t96;
t12 = t34 * t96 + t65 * t92;
t126 = -t88 * t131 - t66 * t87;
t157 = t160 * t126;
t144 = t91 * t95;
t37 = -t88 * t144 + t68 * t87;
t156 = t160 * t37;
t54 = t136 * t88 - t87 * t145;
t155 = t160 * t54;
t147 = t88 * t92;
t146 = t88 * t96;
t142 = t91 * t98;
t141 = t96 * t98;
t86 = t97 * pkin(3) + pkin(2);
t99 = -pkin(10) - pkin(9);
t140 = -t65 * t86 - t66 * t99;
t67 = t98 * t124 + t151 * t94;
t139 = -t67 * t86 - t68 * t99;
t138 = t151 * pkin(1) + pkin(8) * t144;
t134 = t92 * t142;
t133 = t93 * t144;
t130 = -t95 * pkin(1) + pkin(8) * t131;
t129 = pkin(4) * t126 + t34 * pkin(11);
t38 = t87 * t144 + t68 * t88;
t128 = -t37 * pkin(4) + t38 * pkin(11);
t55 = t136 * t87 + t88 * t145;
t127 = t54 * pkin(4) + t55 * pkin(11);
t81 = t93 * t131;
t125 = t66 * t97 - t81;
t122 = -t158 * t65 + t140;
t121 = -t158 * t67 + t139;
t120 = t39 * pkin(3);
t119 = t86 * t142 - t99 * t145;
t118 = pkin(3) * t133 - t67 * t99 + t68 * t86 + t138;
t15 = t38 * t92 - t67 * t96;
t117 = -g(1) * t11 + g(2) * t15;
t116 = g(1) * t126 + g(2) * t37;
t26 = g(1) * t65 - g(2) * t67;
t114 = t159 * pkin(3);
t113 = g(1) * t151 + g(2) * t95;
t112 = pkin(3) * t81 + t65 * t99 - t66 * t86 + t130;
t111 = t158 * t142 + t119;
t31 = t91 * t141 + t55 * t92;
t1 = g(1) * t15 + g(2) * t11 + g(3) * t31;
t16 = t38 * t96 + t67 * t92;
t32 = t55 * t96 - t134;
t110 = g(1) * t16 + g(2) * t12 + g(3) * t32;
t17 = -t65 * t147 - t66 * t96;
t19 = -t67 * t147 - t68 * t96;
t43 = t88 * t134 - t96 * t145;
t109 = g(1) * t19 + g(2) * t17 + g(3) * t43;
t6 = g(1) * t37 - g(2) * t126 - g(3) * t54;
t8 = g(1) * t38 + g(2) * t34 + g(3) * t55;
t108 = t97 * t131 + t66 * t93;
t107 = t120 + t128;
t106 = -g(1) * t67 - g(2) * t65 + g(3) * t142;
t105 = g(1) * t68 + g(2) * t66 + g(3) * t145;
t104 = t38 * pkin(4) + t37 * pkin(11) + t118;
t103 = t114 + t127;
t102 = t108 * pkin(3);
t101 = -pkin(4) * t34 + pkin(11) * t126 + t112;
t100 = -t102 + t129;
t44 = (t88 * t141 + t92 * t94) * t91;
t40 = t68 * t97 + t133;
t20 = -t67 * t146 + t68 * t92;
t18 = -t65 * t146 + t66 * t92;
t10 = t106 * t87;
t5 = t6 * t96;
t4 = t6 * t92;
t3 = g(1) * t12 - g(2) * t16;
t2 = -g(1) * t20 - g(2) * t18 - g(3) * t44;
t7 = [0, 0, 0, 0, 0, 0, g(1) * t95 - g(2) * t151, t113, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t66 - g(2) * t68, -t26, -t113 * t91, -g(1) * t130 - g(2) * t138, 0, 0, 0, 0, 0, 0, g(1) * t125 - g(2) * t40, -g(1) * t108 - g(2) * t39, t26, -g(1) * (-t66 * pkin(2) - t65 * pkin(9) + t130) - g(2) * (t68 * pkin(2) + t67 * pkin(9) + t138) 0, 0, 0, 0, 0, 0, g(1) * t34 - g(2) * t38, t116, t26, -g(1) * t112 - g(2) * t118, 0, 0, 0, 0, 0, 0, t3, t117, -t116, -g(1) * t101 - g(2) * t104, 0, 0, 0, 0, 0, 0, t3, -t116, -t117, -g(1) * (-pkin(5) * t12 - qJ(6) * t11 + t101) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, t105, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t97, t106 * t93, -t105, -g(1) * (-t67 * pkin(2) + t68 * pkin(9)) - g(2) * (-t65 * pkin(2) + t66 * pkin(9)) - g(3) * (pkin(2) * t98 + pkin(9) * t94) * t91, 0, 0, 0, 0, 0, 0, -t106 * t88, t10, -t105, -g(1) * t139 - g(2) * t140 - g(3) * t119, 0, 0, 0, 0, 0, 0, t2, t109, -t10, -g(1) * t121 - g(2) * t122 - g(3) * t111, 0, 0, 0, 0, 0, 0, t2, -t10, -t109, -g(1) * (t20 * pkin(5) + t19 * qJ(6) + t121) - g(2) * (t18 * pkin(5) + t17 * qJ(6) + t122) - g(3) * (t44 * pkin(5) + t43 * qJ(6) + t111); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t39 + g(2) * t108 - g(3) * t159, g(1) * t40 + g(2) * t125 - g(3) * (-t136 * t93 - t94 * t143) 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, -g(1) * t120 + g(2) * t102 - g(3) * t114, 0, 0, 0, 0, 0, 0, t5, -t4, -t8, -g(1) * t107 - g(2) * t100 - g(3) * t103, 0, 0, 0, 0, 0, 0, t5, -t8, t4, -g(1) * (t107 - t156) - g(2) * (t100 + t157) - g(3) * (t103 + t155); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t8, -g(1) * t128 - g(2) * t129 - g(3) * t127, 0, 0, 0, 0, 0, 0, t5, -t8, t4, -g(1) * (t128 - t156) - g(2) * (t129 + t157) - g(3) * (t127 + t155); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t110, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t110, -g(1) * (-t15 * pkin(5) + t16 * qJ(6)) - g(2) * (-t11 * pkin(5) + t12 * qJ(6)) - g(3) * (-t31 * pkin(5) + t32 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t7;
