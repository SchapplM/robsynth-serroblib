% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t145 = cos(qJ(1));
t82 = sin(pkin(6));
t123 = t82 * t145;
t144 = cos(qJ(3));
t128 = cos(pkin(6));
t108 = t128 * t145;
t143 = sin(qJ(1));
t86 = sin(qJ(2));
t88 = cos(qJ(2));
t65 = t86 * t108 + t143 * t88;
t85 = sin(qJ(3));
t39 = -t85 * t123 + t65 * t144;
t64 = -t88 * t108 + t143 * t86;
t81 = sin(pkin(11));
t83 = cos(pkin(11));
t15 = t39 * t81 - t64 * t83;
t16 = t39 * t83 + t64 * t81;
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t154 = t15 * t87 - t16 * t84;
t153 = t15 * t84 + t16 * t87;
t131 = qJ(5) * t81;
t122 = t82 * t144;
t38 = t145 * t122 + t65 * t85;
t142 = t38 * t83;
t152 = -pkin(4) * t142 - t38 * t131;
t121 = t82 * t143;
t107 = t128 * t143;
t67 = -t86 * t107 + t145 * t88;
t42 = -t144 * t121 + t67 * t85;
t141 = t42 * t83;
t151 = -pkin(4) * t141 - t42 * t131;
t137 = t82 * t86;
t62 = -t128 * t144 + t85 * t137;
t140 = t62 * t83;
t150 = -pkin(4) * t140 - t62 * t131;
t99 = g(1) * t42 + g(2) * t38 + g(3) * t62;
t149 = pkin(10) * t85;
t136 = t82 * t88;
t135 = -pkin(10) + qJ(4);
t134 = pkin(2) * t136 + pkin(9) * t137;
t133 = t145 * pkin(1) + pkin(8) * t121;
t132 = qJ(4) * t85;
t130 = t38 * qJ(4);
t129 = t42 * qJ(4);
t127 = t85 * t136;
t126 = t64 * t144;
t66 = t88 * t107 + t145 * t86;
t125 = t66 * t144;
t124 = t81 * t144;
t120 = t83 * t144;
t119 = t88 * t144;
t118 = -t64 * pkin(2) + t65 * pkin(9);
t117 = -t66 * pkin(2) + t67 * pkin(9);
t32 = t38 * pkin(3);
t116 = t39 * qJ(4) - t32;
t34 = t42 * pkin(3);
t43 = t85 * t121 + t67 * t144;
t115 = t43 * qJ(4) - t34;
t57 = t62 * pkin(3);
t63 = t86 * t122 + t128 * t85;
t114 = t63 * qJ(4) - t57;
t112 = t82 * t119;
t113 = pkin(3) * t112 + qJ(4) * t127 + t134;
t111 = -t143 * pkin(1) + pkin(8) * t123;
t19 = t43 * t81 - t66 * t83;
t110 = -g(1) * t15 + g(2) * t19;
t12 = -g(1) * t38 + g(2) * t42;
t109 = g(1) * t64 - g(2) * t66;
t104 = -pkin(3) * t126 - t64 * t132 + t118;
t103 = -pkin(3) * t125 - t66 * t132 + t117;
t102 = t67 * pkin(2) + t66 * pkin(9) + t133;
t101 = t43 * pkin(3) + t102;
t24 = -t64 * t124 - t65 * t83;
t26 = -t66 * t124 - t67 * t83;
t45 = t81 * t112 - t83 * t137;
t100 = g(1) * t26 + g(2) * t24 + g(3) * t45;
t10 = g(1) * t43 + g(2) * t39 + g(3) * t63;
t46 = (t83 * t119 + t81 * t86) * t82;
t98 = t46 * pkin(4) + t45 * qJ(5) + t113;
t97 = g(1) * t145 + g(2) * t143;
t96 = -t65 * pkin(2) - t64 * pkin(9) + t111;
t95 = -g(1) * t66 - g(2) * t64 + g(3) * t136;
t94 = g(1) * t67 + g(2) * t65 + g(3) * t137;
t25 = -t64 * t120 + t65 * t81;
t93 = t25 * pkin(4) + t24 * qJ(5) + t104;
t27 = -t66 * t120 + t67 * t81;
t92 = t27 * pkin(4) + t26 * qJ(5) + t103;
t91 = -pkin(3) * t39 + t96;
t20 = t43 * t83 + t66 * t81;
t90 = t20 * pkin(4) + t19 * qJ(5) + t101;
t89 = -pkin(4) * t16 - qJ(5) * t15 + t91;
t37 = -t81 * t136 + t63 * t83;
t36 = t83 * t136 + t63 * t81;
t21 = t95 * t85;
t7 = t99 * t83;
t6 = t99 * t81;
t5 = g(1) * t16 - g(2) * t20;
t4 = -g(1) * t27 - g(2) * t25 - g(3) * t46;
t3 = t19 * t84 + t20 * t87;
t2 = t19 * t87 - t20 * t84;
t1 = -g(1) * t19 - g(2) * t15 - g(3) * t36;
t8 = [0, 0, 0, 0, 0, 0, g(1) * t143 - g(2) * t145, t97, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t67, -t109, -t97 * t82, -g(1) * t111 - g(2) * t133, 0, 0, 0, 0, 0, 0, g(1) * t39 - g(2) * t43, t12, t109, -g(1) * t96 - g(2) * t102, 0, 0, 0, 0, 0, 0, t5, t110, -t12, -g(1) * (t91 - t130) - g(2) * (t101 + t129) 0, 0, 0, 0, 0, 0, t5, -t12, -t110, -g(1) * (t89 - t130) - g(2) * (t90 + t129) 0, 0, 0, 0, 0, 0, g(1) * t153 - g(2) * t3, g(1) * t154 - g(2) * t2, t12, -g(1) * (-pkin(5) * t16 - t135 * t38 + t89) - g(2) * (t20 * pkin(5) + t135 * t42 + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t94, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t125 + g(2) * t126 - g(3) * t112, t21, -t94, -g(1) * t117 - g(2) * t118 - g(3) * t134, 0, 0, 0, 0, 0, 0, t4, t100, -t21, -g(1) * t103 - g(2) * t104 - g(3) * t113, 0, 0, 0, 0, 0, 0, t4, -t21, -t100, -g(1) * t92 - g(2) * t93 - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t26 * t84 + t27 * t87) - g(2) * (t24 * t84 + t25 * t87) - g(3) * (t45 * t84 + t46 * t87) -g(1) * (t26 * t87 - t27 * t84) - g(2) * (t24 * t87 - t25 * t84) - g(3) * (t45 * t87 - t46 * t84) t21, -g(1) * (t27 * pkin(5) + t66 * t149 + t92) - g(2) * (t25 * pkin(5) + t64 * t149 + t93) - g(3) * (t46 * pkin(5) - pkin(10) * t127 + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t10, -g(1) * t115 - g(2) * t116 - g(3) * t114, 0, 0, 0, 0, 0, 0, t7, -t10, t6, -g(1) * (t115 + t151) - g(2) * (t116 + t152) - g(3) * (t114 + t150) 0, 0, 0, 0, 0, 0, t99 * (t81 * t84 + t83 * t87) t99 * (t81 * t87 - t83 * t84) t10, -g(1) * (-pkin(5) * t141 + t135 * t43 + t151 - t34) - g(2) * (-pkin(5) * t142 + t135 * t39 + t152 - t32) - g(3) * (-pkin(5) * t140 + t135 * t63 + t150 - t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t154 - g(3) * (t36 * t87 - t37 * t84) g(1) * t3 + g(2) * t153 - g(3) * (-t36 * t84 - t37 * t87) 0, 0;];
taug_reg  = t8;
