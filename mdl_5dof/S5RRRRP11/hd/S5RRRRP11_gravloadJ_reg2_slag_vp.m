% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:32
% EndTime: 2019-12-31 22:18:34
% DurationCPUTime: 0.59s
% Computational Cost: add. (466->116), mult. (1223->178), div. (0->0), fcn. (1511->10), ass. (0->78)
t69 = sin(qJ(3));
t117 = pkin(9) * t69;
t73 = cos(qJ(3));
t119 = -pkin(3) * t73 - t117;
t116 = cos(qJ(1));
t67 = sin(pkin(5));
t102 = t67 * t116;
t70 = sin(qJ(2));
t71 = sin(qJ(1));
t74 = cos(qJ(2));
t104 = cos(pkin(5));
t89 = t104 * t116;
t50 = t70 * t89 + t71 * t74;
t26 = -t69 * t102 + t50 * t73;
t49 = t71 * t70 - t74 * t89;
t68 = sin(qJ(4));
t72 = cos(qJ(4));
t8 = t26 * t68 - t49 * t72;
t9 = t26 * t72 + t49 * t68;
t113 = t67 * t70;
t112 = t67 * t71;
t111 = t67 * t73;
t110 = t67 * t74;
t109 = t68 * t73;
t108 = t72 * t73;
t107 = t73 * t74;
t106 = pkin(2) * t110 + pkin(8) * t113;
t105 = t116 * pkin(1) + pkin(7) * t112;
t103 = t68 * t110;
t101 = -t71 * pkin(1) + pkin(7) * t102;
t100 = -t49 * pkin(2) + t50 * pkin(8);
t94 = t71 * t104;
t51 = t116 * t70 + t74 * t94;
t52 = t116 * t74 - t70 * t94;
t99 = -t51 * pkin(2) + t52 * pkin(8);
t95 = -t73 * t102 - t50 * t69;
t98 = pkin(3) * t95 + t26 * pkin(9);
t29 = -t71 * t111 + t52 * t69;
t30 = t69 * t112 + t52 * t73;
t97 = -t29 * pkin(3) + t30 * pkin(9);
t47 = t104 * t73 - t69 * t113;
t48 = t104 * t69 + t70 * t111;
t96 = t47 * pkin(3) + t48 * pkin(9);
t93 = t67 * pkin(3) * t107 + t110 * t117 + t106;
t12 = t30 * t68 - t51 * t72;
t92 = -g(1) * t8 + g(2) * t12;
t91 = g(1) * t95 + g(2) * t29;
t90 = g(1) * t49 - g(2) * t51;
t88 = pkin(4) * t72 + qJ(5) * t68;
t87 = t119 * t49 + t100;
t86 = t119 * t51 + t99;
t85 = t52 * pkin(2) + t51 * pkin(8) + t105;
t84 = g(1) * t116 + g(2) * t71;
t83 = -t50 * pkin(2) - t49 * pkin(8) + t101;
t23 = t72 * t110 + t48 * t68;
t1 = g(1) * t12 + g(2) * t8 + g(3) * t23;
t13 = t30 * t72 + t51 * t68;
t24 = t48 * t72 - t103;
t82 = g(1) * t13 + g(2) * t9 + g(3) * t24;
t15 = -t49 * t109 - t50 * t72;
t17 = -t51 * t109 - t52 * t72;
t31 = t73 * t103 - t72 * t113;
t81 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t80 = g(1) * t29 - g(2) * t95 - g(3) * t47;
t79 = g(1) * t30 + g(2) * t26 + g(3) * t48;
t78 = -g(1) * t51 - g(2) * t49 + g(3) * t110;
t77 = g(1) * t52 + g(2) * t50 + g(3) * t113;
t76 = t30 * pkin(3) + t29 * pkin(9) + t85;
t75 = -pkin(3) * t26 + pkin(9) * t95 + t83;
t32 = (t72 * t107 + t68 * t70) * t67;
t18 = -t51 * t108 + t52 * t68;
t16 = -t49 * t108 + t50 * t68;
t14 = t78 * t69;
t5 = t80 * t72;
t4 = t80 * t68;
t3 = g(1) * t9 - g(2) * t13;
t2 = -g(1) * t18 - g(2) * t16 - g(3) * t32;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t116, t84, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t50 - g(2) * t52, -t90, -t84 * t67, -g(1) * t101 - g(2) * t105, 0, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t30, t91, t90, -g(1) * t83 - g(2) * t85, 0, 0, 0, 0, 0, 0, t3, t92, -t91, -g(1) * t75 - g(2) * t76, 0, 0, 0, 0, 0, 0, t3, -t91, -t92, -g(1) * (-pkin(4) * t9 - qJ(5) * t8 + t75) - g(2) * (t13 * pkin(4) + t12 * qJ(5) + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, 0, 0, 0, 0, 0, 0, 0, -t78 * t73, t14, -t77, -g(1) * t99 - g(2) * t100 - g(3) * t106, 0, 0, 0, 0, 0, 0, t2, t81, -t14, -g(1) * t86 - g(2) * t87 - g(3) * t93, 0, 0, 0, 0, 0, 0, t2, -t14, -t81, -g(1) * (t18 * pkin(4) + t17 * qJ(5) + t86) - g(2) * (t16 * pkin(4) + t15 * qJ(5) + t87) - g(3) * (t32 * pkin(4) + t31 * qJ(5) + t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t79, -g(1) * t97 - g(2) * t98 - g(3) * t96, 0, 0, 0, 0, 0, 0, t5, -t79, t4, -g(1) * (-t88 * t29 + t97) - g(2) * (t88 * t95 + t98) - g(3) * (t88 * t47 + t96); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t82, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t82, -g(1) * (-t12 * pkin(4) + t13 * qJ(5)) - g(2) * (-t8 * pkin(4) + t9 * qJ(5)) - g(3) * (-t23 * pkin(4) + t24 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t6;
