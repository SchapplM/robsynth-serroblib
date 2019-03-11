% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t53 = qJ(3) + qJ(4);
t46 = cos(t53);
t57 = cos(qJ(3));
t49 = t57 * pkin(3);
t34 = pkin(4) * t46 + t49;
t32 = pkin(2) + t34;
t58 = cos(qJ(2));
t27 = t58 * t32;
t60 = -pkin(9) - pkin(8);
t52 = -pkin(10) + t60;
t55 = sin(qJ(2));
t91 = t55 * t52;
t105 = t27 - t91;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t37 = g(1) * t59 + g(2) * t56;
t54 = sin(qJ(3));
t81 = t59 * t57;
t88 = t56 * t58;
t23 = t54 * t88 + t81;
t82 = t59 * t54;
t25 = t56 * t57 - t58 * t82;
t96 = g(3) * t55;
t104 = -g(1) * t25 + g(2) * t23 + t54 * t96;
t62 = -g(3) * t58 + t37 * t55;
t45 = sin(t53);
t103 = pkin(4) * t45;
t47 = qJ(5) + t53;
t42 = sin(t47);
t102 = pkin(5) * t42;
t100 = g(1) * t56;
t94 = t54 * pkin(3);
t93 = t42 * t55;
t43 = cos(t47);
t92 = t43 * t55;
t90 = t55 * t60;
t89 = t56 * t46;
t87 = t58 * t59;
t33 = t94 + t103;
t31 = t59 * t33;
t86 = t59 * t42;
t85 = t59 * t43;
t84 = t59 * t45;
t83 = t59 * t46;
t80 = -t58 * t31 + t56 * t34;
t79 = t59 * pkin(1) + t56 * pkin(7);
t77 = t58 * t84;
t11 = t42 * t88 + t85;
t12 = t43 * t88 - t86;
t76 = -t11 * pkin(5) + t12 * qJ(6);
t75 = -t33 * t88 - t59 * t34;
t13 = -t56 * t43 + t58 * t86;
t14 = t56 * t42 + t58 * t85;
t74 = -t13 * pkin(5) + t14 * qJ(6);
t73 = t58 * pkin(2) + t55 * pkin(8);
t71 = g(1) * t11 - g(2) * t13;
t70 = -g(2) * t59 + t100;
t69 = pkin(5) * t43 + qJ(6) * t42;
t44 = t49 + pkin(2);
t67 = t58 * t44 - t90;
t15 = t45 * t88 + t83;
t64 = t32 * t87 + t56 * t33 - t59 * t91 + t79;
t1 = g(1) * t13 + g(2) * t11 + g(3) * t93;
t3 = g(1) * t14 + g(2) * t12 + g(3) * t92;
t50 = t59 * pkin(7);
t63 = t31 + t50 + (-pkin(1) - t105) * t56;
t19 = t37 * t58 + t96;
t40 = pkin(4) * t89;
t35 = qJ(6) * t92;
t28 = t70 * t55;
t26 = t56 * t54 + t58 * t81;
t24 = -t57 * t88 + t82;
t18 = t56 * t45 + t58 * t83;
t17 = -t77 + t89;
t16 = -t46 * t88 + t84;
t8 = t62 * t43;
t7 = t62 * t42;
t6 = g(1) * t12 - g(2) * t14;
t5 = g(1) * t18 - g(2) * t16 + t46 * t96;
t4 = -g(1) * t17 + g(2) * t15 + t45 * t96;
t2 = [0, 0, 0, 0, 0, 0, t70, t37, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t58, -t28, -t37, -g(1) * (-t56 * pkin(1) + t50) - g(2) * t79, 0, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t26, -g(1) * t23 - g(2) * t25, t28, -g(1) * t50 - g(2) * (t73 * t59 + t79) - (-pkin(1) - t73) * t100, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t28, -g(1) * (pkin(3) * t82 + t50) - g(2) * (t44 * t87 - t59 * t90 + t79) + (-g(1) * (-pkin(1) - t67) - g(2) * t94) * t56, 0, 0, 0, 0, 0, 0, t6, -t71, t28, -g(1) * t63 - g(2) * t64, 0, 0, 0, 0, 0, 0, t6, t28, t71, -g(1) * (-t12 * pkin(5) - t11 * qJ(6) + t63) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t19, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t57, -t62 * t54, -t19, -g(3) * t73 + t37 * (pkin(2) * t55 - pkin(8) * t58) 0, 0, 0, 0, 0, 0, t62 * t46, -t62 * t45, -t19, -g(3) * t67 + t37 * (t44 * t55 + t58 * t60) 0, 0, 0, 0, 0, 0, t8, -t7, -t19, -g(3) * t105 + t37 * (t32 * t55 + t52 * t58) 0, 0, 0, 0, 0, 0, t8, -t19, t7, -g(3) * t27 + (-g(3) * t69 + t37 * t52) * t58 + (g(3) * t52 + t37 * (t32 + t69)) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, g(1) * t26 - g(2) * t24 + t57 * t96, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, t104 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t80 - g(2) * t75 + t33 * t96, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (t74 + t80) - g(2) * (t75 + t76) - g(3) * (t35 + (-t33 - t102) * t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t40 + (g(2) * t83 + t19 * t45) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t77 + t40 + t74) - g(2) * (-t15 * pkin(4) + t76) - g(3) * (t35 + (-t102 - t103) * t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t74 - g(2) * t76 - g(3) * (-pkin(5) * t93 + t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
