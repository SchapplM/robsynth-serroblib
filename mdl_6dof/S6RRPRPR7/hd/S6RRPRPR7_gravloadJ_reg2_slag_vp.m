% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t53 = sin(qJ(1));
t57 = cos(qJ(1));
t30 = g(1) * t57 + g(2) * t53;
t52 = sin(qJ(2));
t98 = t30 * t52;
t50 = sin(qJ(6));
t78 = qJ(4) + pkin(10);
t42 = sin(t78);
t56 = cos(qJ(2));
t73 = cos(t78);
t15 = t52 * t42 + t56 * t73;
t68 = t52 * t73;
t16 = -t56 * t42 + t68;
t5 = t16 * t53;
t83 = t56 * t57;
t7 = t42 * t83 - t57 * t68;
t62 = g(1) * t7 - g(2) * t5 + g(3) * t15;
t101 = t62 * t50;
t54 = cos(qJ(6));
t100 = t62 * t54;
t55 = cos(qJ(4));
t41 = t55 * pkin(4) + pkin(3);
t43 = t52 * qJ(3);
t81 = t56 * pkin(2) + t43;
t99 = t56 * t41 + t81;
t96 = pkin(4) * t53;
t95 = g(1) * t53;
t92 = g(3) * t16;
t91 = t56 * pkin(3);
t90 = pkin(2) + t41;
t51 = sin(qJ(4));
t89 = t52 * t51;
t88 = t52 * t55;
t87 = t52 * t57;
t85 = t56 * t51;
t84 = t56 * t55;
t46 = t57 * pkin(7);
t49 = -qJ(5) - pkin(8);
t82 = t57 * t49 + t46;
t80 = t57 * pkin(1) + t53 * pkin(7);
t79 = qJ(3) * t56;
t37 = pkin(4) * t89;
t77 = pkin(4) * t84;
t76 = t55 * t87;
t75 = t51 * t83;
t74 = pkin(2) * t83 + t57 * t43 + t80;
t72 = -g(1) * t5 - g(2) * t7;
t29 = -g(2) * t57 + t95;
t6 = t15 * t53;
t71 = t6 * t50 - t57 * t54;
t70 = t57 * t50 + t6 * t54;
t69 = t85 - t88;
t19 = t84 + t89;
t25 = t85 * t96;
t67 = t5 * pkin(5) + t6 * pkin(9) - t25;
t28 = pkin(4) * t75;
t8 = t15 * t57;
t66 = -t7 * pkin(5) + t8 * pkin(9) - t28;
t65 = -t15 * pkin(5) + t16 * pkin(9) - t37;
t64 = -pkin(1) - t81;
t63 = t57 * t37 + t41 * t83 + t53 * t49 + t74;
t2 = g(1) * t8 + g(2) * t6 + t92;
t11 = t69 * t53;
t13 = t75 - t76;
t61 = g(1) * t13 + g(2) * t11 + g(3) * t19;
t12 = t19 * t53;
t14 = t19 * t57;
t60 = g(1) * t14 + g(2) * t12 - g(3) * t69;
t59 = t90 * t98;
t58 = (-pkin(1) - t90 * t56 + (-pkin(4) * t51 - qJ(3)) * t52) * t95;
t36 = t57 * t79;
t34 = t53 * t79;
t27 = pkin(4) * t76;
t24 = t88 * t96;
t18 = t29 * t56;
t17 = t29 * t52;
t10 = g(3) * t52 + t30 * t56;
t9 = -g(3) * t56 + t98;
t4 = -t53 * t50 + t8 * t54;
t3 = -t8 * t50 - t53 * t54;
t1 = [0, 0, 0, 0, 0, 0, t29, t30, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t30, -g(1) * (-t53 * pkin(1) + t46) - g(2) * t80, 0, 0, 0, 0, 0, 0, t18, -t30, t17, -g(1) * t46 - g(2) * t74 - t64 * t95, 0, 0, 0, 0, 0, 0, g(1) * t12 - g(2) * t14, -g(1) * t11 + g(2) * t13, t30, -g(1) * (-t57 * pkin(8) + t46) - g(2) * (pkin(3) * t83 + t74) + (-g(1) * (t64 - t91) + g(2) * pkin(8)) * t53, 0, 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, -t72, t30, -g(1) * t82 - g(2) * t63 - t58, 0, 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t4, -g(1) * t71 - g(2) * t3, t72, -g(1) * (-t6 * pkin(5) + t5 * pkin(9) + t82) - g(2) * (t8 * pkin(5) + t7 * pkin(9) + t63) - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t10, -g(1) * (-pkin(2) * t87 + t36) - g(2) * (-t53 * t52 * pkin(2) + t34) - g(3) * t81, 0, 0, 0, 0, 0, 0, -t61, -t60, 0, -g(1) * t36 - g(2) * t34 - g(3) * (t81 + t91) + (pkin(2) + pkin(3)) * t98, 0, 0, 0, 0, 0, 0, -t62, -t2, 0, -g(1) * (t28 + t36) - g(2) * (t25 + t34) - g(3) * (t37 + t99) + t59, 0, 0, 0, 0, 0, 0, -t100, t101, t2, -g(1) * (t36 - t66) - g(2) * (t34 - t67) - g(3) * (-t65 + t99) + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, 0, 0, 0, 0, 0, 0, 0, 0, t62, t2, 0, -g(1) * (-t28 + t27) - g(2) * (-t25 + t24) - g(3) * (-t37 - t77) 0, 0, 0, 0, 0, 0, t100, -t101, -t2, -g(1) * (t27 + t66) - g(2) * (t24 + t67) - g(3) * (t65 - t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t71 + t50 * t92, g(1) * t4 + g(2) * t70 + t54 * t92, 0, 0;];
taug_reg  = t1;
