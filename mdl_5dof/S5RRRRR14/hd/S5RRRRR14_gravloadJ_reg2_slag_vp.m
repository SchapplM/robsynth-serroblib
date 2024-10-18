% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:38
% EndTime: 2024-09-27 18:43:38
% DurationCPUTime: 0.55s
% Computational Cost: add. (682->90), mult. (471->135), div. (0->0), fcn. (450->22), ass. (0->83)
t59 = qJ(1) + qJ(2);
t52 = sin(t59);
t54 = cos(t59);
t63 = sin(qJ(3));
t61 = cos(pkin(5));
t66 = cos(qJ(3));
t90 = t61 * t66;
t20 = -t52 * t90 - t54 * t63;
t69 = -t52 * t63 + t54 * t90;
t60 = sin(pkin(5));
t93 = g(3) * t60;
t98 = -g(1) * t20 - g(2) * t69 - t66 * t93;
t96 = pkin(8) + pkin(9);
t95 = pkin(8) * t60;
t64 = sin(qJ(1));
t92 = t64 * pkin(1);
t49 = t66 * pkin(3) + pkin(2);
t91 = t61 * t63;
t62 = sin(qJ(4));
t89 = t62 * t63;
t88 = t54 * pkin(2) + t52 * t95;
t58 = qJ(3) + qJ(4);
t50 = pkin(5) + t58;
t45 = qJ(5) + t50;
t86 = sin(t45) / 0.2e1;
t85 = sin(t50) / 0.2e1;
t84 = -t52 * pkin(2) + t54 * t95;
t65 = cos(qJ(4));
t48 = t65 * pkin(4) + pkin(3);
t15 = -t60 * (pkin(10) + t96) + (t66 * t62 * pkin(4) + t63 * t48) * t61;
t53 = cos(t58);
t31 = pkin(4) * t53 + t49;
t83 = -t52 * t15 + t54 * t31;
t28 = pkin(3) * t91 - t60 * t96;
t82 = -t52 * t28 + t54 * t49;
t81 = pkin(5) - t58;
t80 = -qJ(5) + t81;
t30 = g(1) * t54 + g(2) * t52;
t29 = g(1) * t52 - g(2) * t54;
t67 = cos(qJ(1));
t79 = g(1) * t64 - g(2) * t67;
t78 = sin(t81);
t77 = -t54 * t15 - t52 * t31;
t68 = sin(t80);
t24 = t86 - t68 / 0.2e1;
t55 = qJ(5) + t58;
t47 = cos(t55);
t76 = t54 * t24 + t52 * t47;
t75 = t52 * t24 - t54 * t47;
t26 = t85 - t78 / 0.2e1;
t74 = t54 * t26 + t52 * t53;
t73 = t52 * t26 - t54 * t53;
t72 = -t54 * t28 - t52 * t49;
t70 = -pkin(4) * t89 + t48 * t66;
t57 = t67 * pkin(1);
t51 = sin(t58);
t46 = sin(t55);
t43 = cos(t81);
t42 = cos(t80);
t41 = cos(t50) / 0.2e1;
t38 = cos(t45) / 0.2e1;
t32 = -t63 * pkin(3) - pkin(4) * t51;
t27 = t43 / 0.2e1 + t41;
t25 = t42 / 0.2e1 + t38;
t22 = t30 * t60;
t21 = -t52 * t91 + t54 * t66;
t19 = -t52 * t66 - t54 * t91;
t16 = t70 * t61;
t14 = -t52 * t27 - t54 * t51;
t13 = -t54 * t27 + t52 * t51;
t12 = -t52 * t25 - t54 * t46;
t11 = -t54 * t25 + t52 * t46;
t10 = -g(1) * t19 - g(2) * t21;
t9 = g(1) * t69 - g(2) * t20;
t8 = g(1) * t74 + g(2) * t73;
t7 = -g(1) * t13 - g(2) * t14;
t6 = g(1) * t76 + g(2) * t75;
t5 = -g(1) * t11 - g(2) * t12;
t4 = -g(1) * t73 + g(2) * t74 - g(3) * (t41 - t43 / 0.2e1);
t3 = -g(1) * t14 + g(2) * t13 - g(3) * (t85 + t78 / 0.2e1);
t2 = -g(1) * t75 + g(2) * t76 - g(3) * (t38 - t42 / 0.2e1);
t1 = -g(1) * t12 + g(2) * t11 - g(3) * (t86 + t68 / 0.2e1);
t17 = [0, 0, 0, 0, 0, 0, t79, g(1) * t67 + g(2) * t64, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, t79 * pkin(1), 0, 0, 0, 0, 0, 0, t10, t9, -t22, -g(1) * (t84 - t92) - g(2) * (t57 + t88), 0, 0, 0, 0, 0, 0, t8, t7, -t22, -g(1) * (t72 - t92) - g(2) * (t57 + t82), 0, 0, 0, 0, 0, 0, t6, t5, -t22, -g(1) * (t77 - t92) - g(2) * (t57 + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, -t22, -g(1) * t84 - g(2) * t88, 0, 0, 0, 0, 0, 0, t8, t7, -t22, -g(1) * t72 - g(2) * t82, 0, 0, 0, 0, 0, 0, t6, t5, -t22, -g(1) * t77 - g(2) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, g(1) * t21 - g(2) * t19 + t63 * t93, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t98 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t52 * t16 + t54 * t32) - g(2) * (t54 * t16 + t52 * t32) - t70 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, (t30 * t51 + (t29 * t61 - t93) * (t65 * t66 - t89)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t17;
