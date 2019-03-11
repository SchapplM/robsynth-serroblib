% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t51 = qJ(3) + qJ(4);
t47 = sin(t51);
t48 = cos(t51);
t82 = cos(pkin(6));
t53 = sin(pkin(6));
t56 = sin(qJ(2));
t91 = t53 * t56;
t99 = -t47 * t91 + t82 * t48;
t59 = cos(qJ(2));
t52 = sin(pkin(11));
t75 = t52 * t82;
t81 = cos(pkin(11));
t32 = -t56 * t75 + t81 * t59;
t92 = t52 * t53;
t98 = -t32 * t47 + t48 * t92;
t69 = t82 * t81;
t30 = t52 * t59 + t56 * t69;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t74 = t53 * t81;
t90 = t53 * t58;
t97 = -g(1) * (-t32 * t55 + t52 * t90) - g(2) * (-t30 * t55 - t58 * t74) - g(3) * (-t55 * t91 + t82 * t58);
t60 = -pkin(9) - pkin(8);
t96 = g(3) * t53;
t46 = pkin(12) + t51;
t43 = cos(t46);
t54 = sin(qJ(6));
t94 = t43 * t54;
t57 = cos(qJ(6));
t93 = t43 * t57;
t89 = t53 * t59;
t88 = t54 * t59;
t87 = t57 * t59;
t29 = t52 * t56 - t59 * t69;
t49 = t58 * pkin(3);
t38 = pkin(4) * t48 + t49;
t36 = pkin(2) + t38;
t50 = -qJ(5) + t60;
t86 = -t29 * t36 - t30 * t50;
t31 = t81 * t56 + t59 * t75;
t85 = -t31 * t36 - t32 * t50;
t37 = -t55 * pkin(3) - pkin(4) * t47;
t84 = t32 * t37 + t38 * t92;
t83 = t37 * t91 + t82 * t38;
t42 = sin(t46);
t13 = -t30 * t42 - t43 * t74;
t14 = t30 * t43 - t42 * t74;
t78 = t13 * pkin(5) + t14 * pkin(10);
t15 = -t32 * t42 + t43 * t92;
t16 = t32 * t43 + t42 * t92;
t77 = t15 * pkin(5) + t16 * pkin(10);
t22 = -t42 * t91 + t82 * t43;
t23 = t82 * t42 + t43 * t91;
t76 = t22 * pkin(5) + t23 * pkin(10);
t72 = t98 * pkin(4);
t71 = pkin(5) * t43 + pkin(10) * t42;
t70 = t99 * pkin(4);
t68 = t30 * t37 - t38 * t74;
t67 = g(1) * t15 + g(2) * t13 + g(3) * t22;
t5 = g(1) * t16 + g(2) * t14 + g(3) * t23;
t65 = -t30 * t47 - t48 * t74;
t9 = -g(1) * t31 - g(2) * t29 + g(3) * t89;
t64 = g(1) * t32 + g(2) * t30 + g(3) * t91;
t61 = t65 * pkin(4);
t45 = t49 + pkin(2);
t26 = t36 * t89;
t8 = t9 * t42;
t7 = -g(1) * (-t32 * t48 - t47 * t92) - g(2) * (-t30 * t48 + t47 * t74) - g(3) * (-t82 * t47 - t48 * t91);
t6 = -g(1) * t98 - g(2) * t65 - g(3) * t99;
t2 = t67 * t57;
t1 = t67 * t54;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t64, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t58, t9 * t55, -t64, -g(1) * (-t31 * pkin(2) + t32 * pkin(8)) - g(2) * (-t29 * pkin(2) + t30 * pkin(8)) - (pkin(2) * t59 + pkin(8) * t56) * t96, 0, 0, 0, 0, 0, 0, -t9 * t48, t9 * t47, -t64, -g(1) * (-t31 * t45 - t32 * t60) - g(2) * (-t29 * t45 - t30 * t60) - (t45 * t59 - t56 * t60) * t96, 0, 0, 0, 0, 0, 0, -t9 * t43, t8, -t64, -g(1) * t85 - g(2) * t86 - g(3) * (-t50 * t91 + t26) 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t93 + t32 * t54) - g(2) * (-t29 * t93 + t30 * t54) - (t43 * t87 + t54 * t56) * t96, -g(1) * (t31 * t94 + t32 * t57) - g(2) * (t29 * t94 + t30 * t57) - (-t43 * t88 + t56 * t57) * t96, -t8, -g(1) * (-t71 * t31 + t85) - g(2) * (-t71 * t29 + t86) - g(3) * t26 - (-t50 * t56 + t71 * t59) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -g(1) * (-t32 * t58 - t55 * t92) - g(2) * (-t30 * t58 + t55 * t74) - g(3) * (-t82 * t55 - t56 * t90) 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t97 * pkin(3), 0, 0, 0, 0, 0, 0, -t67, t5, 0, -g(1) * t84 - g(2) * t68 - g(3) * t83, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t77 + t84) - g(2) * (t68 + t78) - g(3) * (t76 + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t5, 0, -g(1) * t72 - g(2) * t61 - g(3) * t70, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t72 + t77) - g(2) * (t61 + t78) - g(3) * (t70 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t54 + t31 * t57) - g(2) * (-t14 * t54 + t29 * t57) - g(3) * (-t23 * t54 - t53 * t87) -g(1) * (-t16 * t57 - t31 * t54) - g(2) * (-t14 * t57 - t29 * t54) - g(3) * (-t23 * t57 + t53 * t88) 0, 0;];
taug_reg  = t3;
