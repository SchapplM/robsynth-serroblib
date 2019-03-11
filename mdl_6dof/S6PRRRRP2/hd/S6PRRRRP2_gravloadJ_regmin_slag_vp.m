% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t53 = sin(pkin(11));
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t73 = cos(pkin(11));
t74 = cos(pkin(6));
t69 = t74 * t73;
t41 = t53 * t60 + t57 * t69;
t71 = t53 * t74;
t43 = -t57 * t71 + t73 * t60;
t52 = qJ(3) + qJ(4);
t50 = sin(t52);
t51 = cos(t52);
t54 = sin(pkin(6));
t70 = t54 * t73;
t79 = t54 * t57;
t80 = t53 * t54;
t65 = g(3) * (-t50 * t79 + t74 * t51) + g(2) * (-t41 * t50 - t51 * t70) + g(1) * (-t43 * t50 + t51 * t80);
t24 = t41 * t51 - t50 * t70;
t26 = t43 * t51 + t50 * t80;
t35 = t74 * t50 + t51 * t79;
t7 = g(1) * t26 + g(2) * t24 + g(3) * t35;
t40 = t53 * t57 - t60 * t69;
t42 = t73 * t57 + t60 * t71;
t86 = g(1) * t42 + g(2) * t40;
t55 = sin(qJ(5));
t82 = t51 * t55;
t58 = cos(qJ(5));
t81 = t51 * t58;
t59 = cos(qJ(3));
t78 = t54 * t59;
t77 = t54 * t60;
t76 = t58 * t60;
t72 = t55 * t77;
t68 = t59 * pkin(3) + pkin(4) * t51 + pkin(10) * t50 + pkin(2);
t11 = t26 * t55 - t42 * t58;
t27 = t35 * t55 + t54 * t76;
t9 = t24 * t55 - t40 * t58;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t27;
t10 = t24 * t58 + t40 * t55;
t12 = t26 * t58 + t42 * t55;
t28 = t35 * t58 - t72;
t67 = g(1) * t12 + g(2) * t10 + g(3) * t28;
t13 = -t40 * t82 - t41 * t58;
t15 = -t42 * t82 - t43 * t58;
t31 = t51 * t72 - t58 * t79;
t66 = g(1) * t15 + g(2) * t13 + g(3) * t31;
t64 = g(3) * t77 - t86;
t63 = -t7 * pkin(10) - t65 * (pkin(5) * t58 + qJ(6) * t55 + pkin(4));
t56 = sin(qJ(3));
t62 = -g(1) * (-t43 * t56 + t53 * t78) - g(2) * (-t41 * t56 - t59 * t70) - g(3) * (-t56 * t79 + t74 * t59);
t61 = -pkin(9) - pkin(8);
t32 = (t51 * t76 + t55 * t57) * t54;
t16 = -t42 * t81 + t43 * t55;
t14 = -t40 * t81 + t41 * t55;
t8 = t64 * t50;
t4 = t65 * t58;
t3 = t65 * t55;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t32;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t64, g(1) * t43 + g(2) * t41 + g(3) * t79, 0, 0, 0, 0, 0, -t64 * t59, t64 * t56, 0, 0, 0, 0, 0, -t64 * t51, t8, 0, 0, 0, 0, 0, t2, t66, t2, -t8, -t66, -g(1) * (t16 * pkin(5) + t15 * qJ(6) - t43 * t61) - g(2) * (t14 * pkin(5) + t13 * qJ(6) - t41 * t61) + t86 * t68 + (-t32 * pkin(5) - t31 * qJ(6) - (-t57 * t61 + t68 * t60) * t54) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -g(1) * (-t43 * t59 - t56 * t80) - g(2) * (-t41 * t59 + t56 * t70) - g(3) * (-t74 * t56 - t57 * t78) 0, 0, 0, 0, 0, -t65, t7, 0, 0, 0, 0, 0, -t4, t3, -t4, -t7, -t3, t62 * pkin(3) + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t7, 0, 0, 0, 0, 0, -t4, t3, -t4, -t7, -t3, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t67, t1, 0, -t67, -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t27 * pkin(5) + t28 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
