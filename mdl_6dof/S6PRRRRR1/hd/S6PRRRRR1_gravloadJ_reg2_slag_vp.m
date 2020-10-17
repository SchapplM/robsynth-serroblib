% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:25:52
% EndTime: 2019-05-05 10:25:54
% DurationCPUTime: 0.71s
% Computational Cost: add. (791->133), mult. (1027->196), div. (0->0), fcn. (1220->14), ass. (0->72)
t50 = qJ(3) + qJ(4);
t45 = sin(t50);
t46 = cos(t50);
t82 = cos(pkin(6));
t52 = sin(pkin(6));
t55 = sin(qJ(2));
t91 = t52 * t55;
t99 = -t45 * t91 + t82 * t46;
t58 = cos(qJ(2));
t51 = sin(pkin(12));
t75 = t51 * t82;
t81 = cos(pkin(12));
t31 = -t55 * t75 + t81 * t58;
t92 = t51 * t52;
t98 = -t31 * t45 + t46 * t92;
t69 = t82 * t81;
t29 = t51 * t58 + t55 * t69;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t74 = t52 * t81;
t90 = t52 * t57;
t97 = -g(1) * (-t31 * t54 + t51 * t90) - g(2) * (-t29 * t54 - t57 * t74) - g(3) * (-t54 * t91 + t82 * t57);
t59 = -pkin(9) - pkin(8);
t96 = g(3) * t52;
t47 = qJ(5) + t50;
t43 = cos(t47);
t53 = sin(qJ(6));
t94 = t43 * t53;
t56 = cos(qJ(6));
t93 = t43 * t56;
t89 = t52 * t58;
t88 = t53 * t58;
t87 = t56 * t58;
t28 = t51 * t55 - t58 * t69;
t48 = t57 * pkin(3);
t37 = pkin(4) * t46 + t48;
t35 = pkin(2) + t37;
t49 = -pkin(10) + t59;
t86 = -t28 * t35 - t29 * t49;
t30 = t81 * t55 + t58 * t75;
t85 = -t30 * t35 - t31 * t49;
t36 = -t54 * pkin(3) - pkin(4) * t45;
t84 = t31 * t36 + t37 * t92;
t83 = t36 * t91 + t82 * t37;
t42 = sin(t47);
t12 = -t29 * t42 - t43 * t74;
t13 = t29 * t43 - t42 * t74;
t78 = t12 * pkin(5) + t13 * pkin(11);
t14 = -t31 * t42 + t43 * t92;
t15 = t31 * t43 + t42 * t92;
t77 = t14 * pkin(5) + t15 * pkin(11);
t21 = -t42 * t91 + t82 * t43;
t22 = t82 * t42 + t43 * t91;
t76 = t21 * pkin(5) + t22 * pkin(11);
t72 = t98 * pkin(4);
t71 = pkin(5) * t43 + pkin(11) * t42;
t70 = t99 * pkin(4);
t68 = t29 * t36 - t37 * t74;
t67 = g(1) * t14 + g(2) * t12 + g(3) * t21;
t5 = g(1) * t15 + g(2) * t13 + g(3) * t22;
t65 = -t29 * t45 - t46 * t74;
t64 = -g(1) * t30 - g(2) * t28 + g(3) * t89;
t63 = g(1) * t31 + g(2) * t29 + g(3) * t91;
t60 = t65 * pkin(4);
t44 = t48 + pkin(2);
t25 = t35 * t89;
t8 = t64 * t42;
t7 = -g(1) * (-t31 * t46 - t45 * t92) - g(2) * (-t29 * t46 + t45 * t74) - g(3) * (-t82 * t45 - t46 * t91);
t6 = -g(1) * t98 - g(2) * t65 - g(3) * t99;
t2 = t67 * t56;
t1 = t67 * t53;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t63, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t57, t64 * t54, -t63, -g(1) * (-t30 * pkin(2) + t31 * pkin(8)) - g(2) * (-t28 * pkin(2) + t29 * pkin(8)) - (pkin(2) * t58 + pkin(8) * t55) * t96, 0, 0, 0, 0, 0, 0, -t64 * t46, t64 * t45, -t63, -g(1) * (-t30 * t44 - t31 * t59) - g(2) * (-t28 * t44 - t29 * t59) - (t44 * t58 - t55 * t59) * t96, 0, 0, 0, 0, 0, 0, -t64 * t43, t8, -t63, -g(1) * t85 - g(2) * t86 - g(3) * (-t49 * t91 + t25) 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t93 + t31 * t53) - g(2) * (-t28 * t93 + t29 * t53) - (t43 * t87 + t53 * t55) * t96, -g(1) * (t30 * t94 + t31 * t56) - g(2) * (t28 * t94 + t29 * t56) - (-t43 * t88 + t55 * t56) * t96, -t8, -g(1) * (-t71 * t30 + t85) - g(2) * (-t71 * t28 + t86) - g(3) * t25 - (-t49 * t55 + t71 * t58) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -g(1) * (-t31 * t57 - t54 * t92) - g(2) * (-t29 * t57 + t54 * t74) - g(3) * (-t82 * t54 - t55 * t90) 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t97 * pkin(3), 0, 0, 0, 0, 0, 0, -t67, t5, 0, -g(1) * t84 - g(2) * t68 - g(3) * t83, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t77 + t84) - g(2) * (t68 + t78) - g(3) * (t76 + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t5, 0, -g(1) * t72 - g(2) * t60 - g(3) * t70, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t72 + t77) - g(2) * (t60 + t78) - g(3) * (t70 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t77 - g(2) * t78 - g(3) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t53 + t30 * t56) - g(2) * (-t13 * t53 + t28 * t56) - g(3) * (-t22 * t53 - t52 * t87) -g(1) * (-t15 * t56 - t30 * t53) - g(2) * (-t13 * t56 - t28 * t53) - g(3) * (-t22 * t56 + t52 * t88) 0, 0;];
taug_reg  = t3;
