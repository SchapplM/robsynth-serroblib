% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:37:12
% EndTime: 2019-05-06 13:37:16
% DurationCPUTime: 0.95s
% Computational Cost: add. (778->150), mult. (1693->236), div. (0->0), fcn. (2132->14), ass. (0->91)
t104 = cos(pkin(11));
t65 = sin(pkin(11));
t71 = sin(qJ(2));
t75 = cos(qJ(2));
t46 = -t71 * t104 - t75 * t65;
t72 = sin(qJ(1));
t76 = cos(qJ(1));
t67 = cos(pkin(6));
t83 = t75 * t104 - t71 * t65;
t78 = t83 * t67;
t24 = t72 * t46 + t76 * t78;
t69 = sin(qJ(6));
t73 = cos(qJ(6));
t66 = sin(pkin(6));
t111 = t66 * t76;
t38 = t46 * t67;
t25 = -t76 * t38 + t72 * t83;
t64 = qJ(4) + pkin(12);
t62 = sin(t64);
t63 = cos(t64);
t8 = -t62 * t111 + t25 * t63;
t124 = t24 * t73 + t8 * t69;
t123 = -t24 * t69 + t8 * t73;
t26 = -t72 * t38 - t76 * t83;
t37 = t46 * t66;
t4 = g(1) * t26 - g(2) * t25 + g(3) * t37;
t112 = t66 * t75;
t106 = t76 * t71;
t107 = t72 * t75;
t42 = -t67 * t107 - t106;
t122 = -g(1) * t42 - g(3) * t112;
t70 = sin(qJ(4));
t117 = t26 * t70;
t116 = t37 * t70;
t115 = t63 * t69;
t114 = t63 * t73;
t113 = t66 * t72;
t74 = cos(qJ(4));
t110 = t67 * t74;
t108 = t72 * t71;
t105 = t76 * t75;
t102 = t70 * t113;
t101 = t74 * t113;
t56 = t70 * t111;
t100 = t74 * t111;
t99 = t67 * t105;
t36 = t83 * t66;
t57 = pkin(2) * t112;
t60 = t74 * pkin(4) + pkin(3);
t68 = -qJ(5) - pkin(9);
t98 = t36 * t60 + t37 * t68 + t57;
t97 = -t63 * t111 - t25 * t62;
t96 = t25 * t74 - t56;
t39 = t67 * t71 * pkin(2) + (-pkin(8) - qJ(3)) * t66;
t61 = t75 * pkin(2) + pkin(1);
t95 = -t72 * t39 + t76 * t61;
t53 = pkin(2) * t99;
t93 = -pkin(2) * t108 + t53;
t11 = -t63 * t113 - t26 * t62;
t92 = g(1) * t97 + g(2) * t11;
t91 = pkin(5) * t63 + pkin(10) * t62;
t27 = t76 * t46 - t72 * t78;
t90 = g(1) * t24 - g(2) * t27;
t89 = g(1) * t76 + g(2) * t72;
t88 = g(1) * t72 - g(2) * t76;
t87 = -t76 * t39 - t72 * t61;
t86 = t25 * t70 + t100;
t85 = t24 * t60 - t25 * t68 + t93;
t84 = pkin(4) * t102 - t26 * t60 + t27 * t68 + t95;
t29 = t37 * t62 + t67 * t63;
t82 = g(1) * t11 - g(2) * t97 - g(3) * t29;
t12 = t62 * t113 - t26 * t63;
t30 = -t37 * t63 + t67 * t62;
t81 = g(1) * t12 + g(2) * t8 + g(3) * t30;
t5 = g(1) * t27 + g(2) * t24 + g(3) * t36;
t80 = t42 * pkin(2);
t79 = pkin(4) * t56 - t24 * t68 - t25 * t60 + t87;
t77 = t26 * t68 + t27 * t60 + t80;
t58 = pkin(4) * t110;
t51 = pkin(4) * t101;
t44 = t89 * t66;
t43 = -t67 * t108 + t105;
t41 = -t67 * t106 - t107;
t40 = -t99 + t108;
t35 = -g(3) * t67 - t88 * t66;
t14 = -t26 * t74 + t102;
t13 = t101 + t117;
t3 = t12 * t73 - t27 * t69;
t2 = -t12 * t69 - t27 * t73;
t1 = t5 * t62;
t6 = [0, 0, 0, 0, 0, 0, t88, t89, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t41 - g(2) * t43, -g(1) * t40 - g(2) * t42, -t44, -g(1) * (-t72 * pkin(1) + pkin(8) * t111) - g(2) * (t76 * pkin(1) + pkin(8) * t113) 0, 0, 0, 0, 0, 0, g(1) * t25 + g(2) * t26, t90, -t44, -g(1) * t87 - g(2) * t95, 0, 0, 0, 0, 0, 0, g(1) * t96 - g(2) * t14, -g(1) * t86 - g(2) * t13, -t90, -g(1) * (-t25 * pkin(3) + t24 * pkin(9) + t87) - g(2) * (-pkin(3) * t26 - t27 * pkin(9) + t95) 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t12, t92, -t90, -g(1) * t79 - g(2) * t84, 0, 0, 0, 0, 0, 0, g(1) * t123 - g(2) * t3, -g(1) * t124 - g(2) * t2, -t92, -g(1) * (-pkin(5) * t8 + pkin(10) * t97 + t79) - g(2) * (t12 * pkin(5) + t11 * pkin(10) + t84); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t40 + t122, g(3) * t66 * t71 + g(1) * t43 - g(2) * t41, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t4, 0, -g(2) * t53 + (g(2) * t108 + t122) * pkin(2), 0, 0, 0, 0, 0, 0, -t5 * t74, t5 * t70, t4, -g(1) * (t27 * pkin(3) - t26 * pkin(9) + t80) - g(2) * (t24 * pkin(3) + pkin(9) * t25 + t93) - g(3) * (t36 * pkin(3) - t37 * pkin(9) + t57) 0, 0, 0, 0, 0, 0, -t5 * t63, t1, t4, -g(1) * t77 - g(2) * t85 - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t27 * t114 - t26 * t69) - g(2) * (t24 * t114 + t25 * t69) - g(3) * (t36 * t114 - t37 * t69) -g(1) * (-t27 * t115 - t26 * t73) - g(2) * (-t24 * t115 + t25 * t73) - g(3) * (-t36 * t115 - t37 * t73) -t1, -g(1) * (t91 * t27 + t77) - g(2) * (t91 * t24 + t85) - g(3) * (t91 * t36 + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t86 - g(3) * (t110 + t116) g(1) * t14 + g(2) * t96 - g(3) * (t37 * t74 - t67 * t70) 0, 0, 0, 0, 0, 0, 0, 0, t82, t81, 0, -g(1) * t51 - g(3) * t58 + (g(2) * t100 - t4 * t70) * pkin(4), 0, 0, 0, 0, 0, 0, t82 * t73, -t82 * t69, -t81, -g(1) * (pkin(4) * t117 - t11 * pkin(5) + t12 * pkin(10) + t51) - g(2) * (-t86 * pkin(4) + pkin(5) * t97 + t8 * pkin(10)) - g(3) * (pkin(4) * t116 + t29 * pkin(5) + t30 * pkin(10) + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t124 - g(3) * (-t30 * t69 - t36 * t73) g(1) * t3 + g(2) * t123 - g(3) * (-t30 * t73 + t36 * t69) 0, 0;];
taug_reg  = t6;
