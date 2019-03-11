% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t59 = cos(qJ(4));
t56 = sin(qJ(4));
t91 = pkin(4) * t56;
t92 = -pkin(9) * t59 + qJ(3) + t91;
t52 = sin(pkin(10));
t54 = cos(pkin(10));
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t80 = cos(pkin(6));
t72 = t60 * t80;
t37 = t52 * t57 - t54 * t72;
t90 = t37 * pkin(8);
t39 = t52 * t72 + t54 * t57;
t89 = t39 * pkin(8);
t53 = sin(pkin(6));
t88 = t53 * t56;
t87 = t53 * t57;
t86 = t53 * t59;
t85 = t53 * t60;
t55 = sin(qJ(5));
t84 = t55 * t56;
t58 = cos(qJ(5));
t83 = t56 * t58;
t82 = t57 * t58;
t81 = pkin(2) * t85 + qJ(3) * t87;
t79 = t55 * t87;
t78 = pkin(8) * t85 + t81;
t19 = t39 * t59 - t52 * t88;
t20 = t39 * t56 + t52 * t86;
t77 = t19 * pkin(4) + t20 * pkin(9);
t21 = t37 * t59 + t54 * t88;
t22 = -t37 * t56 + t54 * t86;
t76 = t21 * pkin(4) - t22 * pkin(9);
t41 = -t80 * t56 - t59 * t85;
t42 = -t56 * t85 + t80 * t59;
t75 = t41 * pkin(4) + t42 * pkin(9);
t73 = t57 * t80;
t34 = t37 * pkin(2);
t38 = t52 * t60 + t54 * t73;
t71 = t38 * qJ(3) - t34;
t35 = t39 * pkin(2);
t40 = -t52 * t73 + t54 * t60;
t70 = t40 * qJ(3) - t35;
t69 = pkin(5) * t58 + qJ(6) * t55;
t23 = t42 * t55 - t53 * t82;
t6 = t20 * t55 - t40 * t58;
t8 = -t22 * t55 - t38 * t58;
t1 = g(1) * t6 + g(2) * t8 + g(3) * t23;
t24 = t42 * t58 + t79;
t7 = t20 * t58 + t40 * t55;
t9 = -t22 * t58 + t38 * t55;
t68 = g(1) * t7 + g(2) * t9 + g(3) * t24;
t13 = t37 * t58 + t38 * t84;
t15 = t39 * t58 + t40 * t84;
t25 = t56 * t79 - t58 * t85;
t67 = g(1) * t15 + g(2) * t13 + g(3) * t25;
t66 = g(1) * t19 + g(2) * t21 + g(3) * t41;
t65 = g(1) * t20 - g(2) * t22 + g(3) * t42;
t11 = -g(1) * t39 - g(2) * t37 + g(3) * t85;
t64 = g(1) * t40 + g(2) * t38 + g(3) * t87;
t63 = -t57 * pkin(9) * t86 + t87 * t91 + t78;
t62 = t92 * t38 - t34 - t90;
t61 = t92 * t40 - t35 - t89;
t26 = (t55 * t60 + t56 * t82) * t53;
t16 = -t39 * t55 + t40 * t83;
t14 = -t37 * t55 + t38 * t83;
t10 = t64 * t59;
t4 = t66 * t58;
t3 = t66 * t55;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t64, -g(1) * t70 - g(2) * t71 - g(3) * t81, 0, 0, 0, 0, 0, 0, -t64 * t56, -t10, -t11, -g(1) * (t70 - t89) - g(2) * (t71 - t90) - g(3) * t78, 0, 0, 0, 0, 0, 0, t2, t67, t10, -g(1) * t61 - g(2) * t62 - g(3) * t63, 0, 0, 0, 0, 0, 0, t2, t10, -t67, -g(1) * (t16 * pkin(5) + t15 * qJ(6) + t61) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t62) - g(3) * (t26 * pkin(5) + t25 * qJ(6) + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t65, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t65, -g(1) * t77 - g(2) * t76 - g(3) * t75, 0, 0, 0, 0, 0, 0, -t4, -t65, -t3, -g(1) * (t69 * t19 + t77) - g(2) * (t69 * t21 + t76) - g(3) * (t69 * t41 + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t68, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t68, -g(1) * (-t6 * pkin(5) + t7 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - g(3) * (-t23 * pkin(5) + t24 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
