% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t103 = sin(qJ(1));
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t104 = cos(qJ(1));
t84 = cos(pkin(6));
t72 = t84 * t104;
t26 = t103 * t59 + t56 * t72;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t52 = sin(pkin(6));
t80 = t52 * t104;
t14 = t26 * t58 - t55 * t80;
t25 = t103 * t56 - t59 * t72;
t51 = qJ(4) + pkin(12);
t45 = sin(t51);
t46 = cos(t51);
t115 = t14 * t45 - t25 * t46;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t114 = t14 * t54 - t25 * t57;
t113 = t14 * t57 + t25 * t54;
t112 = t14 * t46 + t25 * t45;
t47 = qJ(6) + t51;
t42 = sin(t47);
t43 = cos(t47);
t111 = t14 * t42 - t25 * t43;
t110 = t14 * t43 + t25 * t42;
t71 = t84 * t103;
t28 = t104 * t59 - t56 * t71;
t79 = t52 * t103;
t18 = t28 * t58 + t55 * t79;
t27 = t104 * t56 + t59 * t71;
t10 = -t18 * t54 + t27 * t57;
t92 = t52 * t56;
t24 = t84 * t55 + t58 * t92;
t91 = t52 * t59;
t109 = g(2) * t114 - g(3) * (-t24 * t54 - t57 * t91) - g(1) * t10;
t107 = g(3) * t52;
t106 = t54 * pkin(4);
t30 = pkin(5) * t45 + t106;
t105 = pkin(9) + t30;
t96 = t42 * t58;
t95 = t43 * t58;
t94 = t45 * t58;
t93 = t46 * t58;
t90 = t54 * t56;
t89 = t54 * t58;
t88 = t57 * t58;
t87 = t58 * t59;
t53 = -qJ(5) - pkin(10);
t86 = pkin(2) * t91 + pkin(9) * t92;
t48 = t57 * pkin(4);
t31 = pkin(5) * t46 + t48;
t85 = t104 * pkin(1) + pkin(8) * t79;
t83 = t28 * pkin(2) + t85;
t82 = pkin(9) + t106;
t81 = g(3) * t86;
t19 = t25 * pkin(2);
t78 = t26 * pkin(9) - t19;
t21 = t27 * pkin(2);
t77 = t28 * pkin(9) - t21;
t76 = -t103 * pkin(1) + pkin(8) * t80;
t75 = pkin(3) * t58 + pkin(10) * t55;
t13 = t26 * t55 + t58 * t80;
t17 = t28 * t55 - t58 * t79;
t74 = -g(1) * t13 + g(2) * t17;
t73 = g(1) * t25 - g(2) * t27;
t29 = pkin(3) + t31;
t50 = -pkin(11) + t53;
t70 = t29 * t58 - t50 * t55;
t44 = t48 + pkin(3);
t69 = t44 * t58 - t53 * t55;
t68 = t27 * pkin(9) + t83;
t67 = -t26 * pkin(2) + t76;
t23 = t55 * t92 - t84 * t58;
t66 = g(1) * t17 + g(2) * t13 + g(3) * t23;
t65 = g(1) * t18 + g(2) * t14 + g(3) * t24;
t64 = g(1) * t104 + g(2) * t103;
t63 = -t25 * pkin(9) + t67;
t62 = -g(1) * t27 - g(2) * t25 + g(3) * t91;
t61 = g(1) * t28 + g(2) * t26 + g(3) * t92;
t12 = t62 * t55;
t11 = t18 * t57 + t27 * t54;
t9 = t18 * t46 + t27 * t45;
t8 = -t18 * t45 + t27 * t46;
t7 = t18 * t43 + t27 * t42;
t6 = -t18 * t42 + t27 * t43;
t2 = g(1) * t7 + g(2) * t110 - g(3) * (-t24 * t43 + t42 * t91);
t1 = -g(1) * t6 + g(2) * t111 - g(3) * (-t24 * t42 - t43 * t91);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t103 - g(2) * t104, t64, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t28, -t73, -t64 * t52, -g(1) * t76 - g(2) * t85, 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t74, t73, -g(1) * t63 - g(2) * t68, 0, 0, 0, 0, 0, 0, g(1) * t113 - g(2) * t11, -g(1) * t114 - g(2) * t10, -t74, -g(1) * (-pkin(3) * t14 - pkin(10) * t13 + t63) - g(2) * (t18 * pkin(3) + t17 * pkin(10) + t68) 0, 0, 0, 0, 0, 0, g(1) * t112 - g(2) * t9, -g(1) * t115 - g(2) * t8, -t74, -g(1) * (t13 * t53 - t14 * t44 - t82 * t25 + t67) - g(2) * (-t17 * t53 + t18 * t44 + t82 * t27 + t83) 0, 0, 0, 0, 0, 0, g(1) * t110 - g(2) * t7, -g(1) * t111 - g(2) * t6, -t74, -g(1) * (-t105 * t25 + t13 * t50 - t14 * t29 + t67) - g(2) * (t105 * t27 - t17 * t50 + t18 * t29 + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t61, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t58, t12, -t61, -g(1) * t77 - g(2) * t78 - t81, 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t88 + t28 * t54) - g(2) * (-t25 * t88 + t26 * t54) - (t57 * t87 + t90) * t107, -g(1) * (t27 * t89 + t28 * t57) - g(2) * (t25 * t89 + t26 * t57) - (-t54 * t87 + t56 * t57) * t107, -t12, -g(1) * (-t75 * t27 + t77) - g(2) * (-t75 * t25 + t78) - g(3) * (t75 * t91 + t86) 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t93 + t28 * t45) - g(2) * (-t25 * t93 + t26 * t45) - (t45 * t56 + t46 * t87) * t107, -g(1) * (t27 * t94 + t28 * t46) - g(2) * (t25 * t94 + t26 * t46) - (-t45 * t87 + t46 * t56) * t107, -t12, -g(1) * (-t69 * t27 + t82 * t28 - t21) - g(2) * (-t69 * t25 + t82 * t26 - t19) - t81 - (pkin(4) * t90 + t69 * t59) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t95 + t28 * t42) - g(2) * (-t25 * t95 + t26 * t42) - (t42 * t56 + t43 * t87) * t107, -g(1) * (t27 * t96 + t28 * t43) - g(2) * (t25 * t96 + t26 * t43) - (-t42 * t87 + t43 * t56) * t107, -t12, -g(1) * (t105 * t28 - t27 * t70 - t21) - g(2) * (t105 * t26 - t25 * t70 - t19) - t81 - (t30 * t56 + t59 * t70) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t57, -t66 * t54, -t65, -g(1) * (-t17 * pkin(3) + t18 * pkin(10)) - g(2) * (-t13 * pkin(3) + t14 * pkin(10)) - g(3) * (-t23 * pkin(3) + t24 * pkin(10)) 0, 0, 0, 0, 0, 0, t66 * t46, -t66 * t45, -t65, -g(1) * (-t17 * t44 - t18 * t53) - g(2) * (-t13 * t44 - t14 * t53) - g(3) * (-t23 * t44 - t24 * t53) 0, 0, 0, 0, 0, 0, t66 * t43, -t66 * t42, -t65, -g(1) * (-t17 * t29 - t18 * t50) - g(2) * (-t13 * t29 - t14 * t50) - g(3) * (-t23 * t29 - t24 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, g(1) * t11 + g(2) * t113 - g(3) * (-t24 * t57 + t54 * t91) 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t115 - g(3) * (-t24 * t45 - t46 * t91) g(1) * t9 + g(2) * t112 - g(3) * (-t24 * t46 + t45 * t91) 0, t109 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t18 * t30 + t27 * t31) - g(2) * (-t14 * t30 + t25 * t31) - g(3) * (-t24 * t30 - t31 * t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
