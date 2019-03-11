% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t103 = cos(pkin(11));
t65 = sin(pkin(11));
t70 = sin(qJ(2));
t74 = cos(qJ(2));
t52 = -t103 * t70 - t65 * t74;
t66 = sin(pkin(6));
t75 = cos(qJ(1));
t117 = t66 * t75;
t67 = cos(pkin(6));
t106 = t52 * t67;
t71 = sin(qJ(1));
t84 = t103 * t74 - t65 * t70;
t28 = t106 * t75 - t71 * t84;
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t14 = t117 * t73 - t28 * t69;
t79 = t84 * t67;
t29 = t71 * t52 + t75 * t79;
t68 = sin(qJ(6));
t72 = cos(qJ(6));
t134 = t14 * t68 - t29 * t72;
t133 = t14 * t72 + t29 * t68;
t104 = qJ(5) * t69;
t122 = t29 * t73;
t132 = pkin(4) * t122 + t104 * t29;
t32 = t52 * t75 - t71 * t79;
t121 = t32 * t73;
t131 = pkin(4) * t121 + t104 * t32;
t43 = t84 * t66;
t120 = t43 * t73;
t130 = pkin(4) * t120 + t104 * t43;
t118 = t66 * t74;
t108 = t75 * t70;
t111 = t71 * t74;
t48 = -t111 * t67 - t108;
t129 = -g(1) * t48 - g(3) * t118;
t33 = t106 * t71 + t75 * t84;
t128 = pkin(5) + pkin(9);
t126 = t29 * pkin(9);
t125 = t32 * pkin(9);
t119 = t66 * t71;
t116 = t68 * t69;
t115 = t69 * t72;
t112 = t71 * t70;
t107 = t75 * t74;
t105 = pkin(2) * t118 + pkin(3) * t43;
t101 = t67 * t107;
t15 = -t117 * t69 - t28 * t73;
t45 = t67 * t70 * pkin(2) + (-pkin(8) - qJ(3)) * t66;
t64 = pkin(2) * t74 + pkin(1);
t100 = -t71 * t45 + t64 * t75;
t97 = -pkin(4) * t14 + qJ(5) * t15;
t18 = -t119 * t73 + t33 * t69;
t19 = t119 * t69 + t33 * t73;
t96 = -pkin(4) * t18 + qJ(5) * t19;
t44 = t52 * t66;
t35 = -t44 * t69 - t67 * t73;
t36 = -t44 * t73 + t67 * t69;
t95 = -pkin(4) * t35 + qJ(5) * t36;
t94 = -pkin(9) * t44 + t105;
t93 = pkin(3) * t33 + t100;
t92 = -g(1) * t14 + g(2) * t18;
t91 = -g(1) * t15 + g(2) * t19;
t90 = g(1) * t29 - g(2) * t32;
t89 = g(1) * t75 + g(2) * t71;
t88 = g(1) * t71 - g(2) * t75;
t87 = -t75 * t45 - t71 * t64;
t53 = pkin(2) * t101;
t86 = -pkin(2) * t112 + pkin(3) * t29 + t53;
t85 = pkin(3) * t28 + t87;
t2 = g(1) * t18 + g(2) * t14 + g(3) * t35;
t83 = g(1) * t19 + g(2) * t15 + g(3) * t36;
t6 = -g(1) * t33 + g(2) * t28 + g(3) * t44;
t82 = g(1) * t32 + g(2) * t29 + g(3) * t43;
t81 = -pkin(9) * t28 + t86;
t80 = pkin(4) * t19 + t18 * qJ(5) + t93;
t78 = pkin(2) * t48 + pkin(3) * t32;
t77 = -pkin(4) * t15 - qJ(5) * t14 + t85;
t76 = pkin(9) * t33 + t78;
t50 = t89 * t66;
t49 = -t112 * t67 + t107;
t47 = -t108 * t67 - t111;
t46 = -t101 + t112;
t42 = -g(3) * t67 - t66 * t88;
t8 = t18 * t68 - t32 * t72;
t7 = t18 * t72 + t32 * t68;
t4 = t82 * t73;
t3 = t82 * t69;
t1 = [0, 0, 0, 0, 0, 0, t88, t89, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t47 - g(2) * t49, -g(1) * t46 - g(2) * t48, -t50, -g(1) * (-pkin(1) * t71 + pkin(8) * t117) - g(2) * (pkin(1) * t75 + pkin(8) * t119) 0, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t33, t90, -t50, -g(1) * t87 - g(2) * t100, 0, 0, 0, 0, 0, 0, -t91, t92, -t90, -g(1) * (t85 + t126) - g(2) * (t93 - t125) 0, 0, 0, 0, 0, 0, -t90, t91, -t92, -g(1) * (t77 + t126) - g(2) * (t80 - t125) 0, 0, 0, 0, 0, 0, g(1) * t134 - g(2) * t8, g(1) * t133 - g(2) * t7, -t91, -g(1) * (-pkin(10) * t15 + t128 * t29 + t77) - g(2) * (t19 * pkin(10) - t128 * t32 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t46 + t129, g(3) * t66 * t70 + g(1) * t49 - g(2) * t47, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t6, 0, -g(2) * t53 + (g(2) * t112 + t129) * pkin(2), 0, 0, 0, 0, 0, 0, -t4, t3, t6, -g(1) * t76 - g(2) * t81 - g(3) * t94, 0, 0, 0, 0, 0, 0, t6, t4, -t3, -g(1) * (t76 + t131) - g(2) * (t81 + t132) - g(3) * (t94 + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t116 * t32 + t33 * t72) - g(2) * (t116 * t29 - t28 * t72) - g(3) * (t116 * t43 - t44 * t72) -g(1) * (t115 * t32 - t33 * t68) - g(2) * (t115 * t29 + t28 * t68) - g(3) * (t115 * t43 + t44 * t68) -t4, -g(1) * (pkin(10) * t121 + t128 * t33 + t131 + t78) - g(2) * (pkin(10) * t122 - t128 * t28 + t132 + t86) - g(3) * (pkin(10) * t120 - t128 * t44 + t105 + t130); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t83, -g(1) * t96 - g(2) * t97 - g(3) * t95, 0, 0, 0, 0, 0, 0, -t83 * t68, -t83 * t72, t2, -g(1) * (-pkin(10) * t18 + t96) - g(2) * (-pkin(10) * t14 + t97) - g(3) * (-pkin(10) * t35 + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t133 - g(3) * (t35 * t72 + t43 * t68) g(1) * t8 + g(2) * t134 - g(3) * (-t35 * t68 + t43 * t72) 0, 0;];
taug_reg  = t1;
