% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:24:51
% EndTime: 2019-05-05 00:24:54
% DurationCPUTime: 0.75s
% Computational Cost: add. (641->127), mult. (1550->205), div. (0->0), fcn. (1984->14), ass. (0->68)
t43 = sin(qJ(2));
t73 = sin(pkin(12));
t75 = cos(pkin(12));
t86 = cos(qJ(2));
t26 = -t43 * t75 - t86 * t73;
t38 = sin(pkin(11));
t39 = cos(pkin(11));
t40 = cos(pkin(6));
t71 = t40 * t86;
t87 = -t38 * t43 + t39 * t71;
t25 = t43 * t73 - t86 * t75;
t77 = t26 * t40;
t16 = -t39 * t25 + t38 * t77;
t11 = t38 * t25 + t39 * t77;
t37 = qJ(5) + qJ(6);
t35 = sin(t37);
t45 = cos(qJ(4));
t85 = t35 * t45;
t36 = cos(t37);
t84 = t36 * t45;
t80 = t40 * t43;
t41 = sin(qJ(5));
t79 = t41 * t45;
t44 = cos(qJ(5));
t78 = t44 * t45;
t74 = sin(pkin(6));
t58 = t74 * t73;
t59 = t75 * t74;
t23 = t43 * t58 - t86 * t59;
t62 = t74 * t86;
t32 = pkin(2) * t62;
t76 = -t23 * pkin(3) + t32;
t72 = pkin(5) * t41 + pkin(8);
t42 = sin(qJ(4));
t70 = t42 * t74;
t68 = t45 * t74;
t66 = t87 * pkin(2);
t24 = t43 * t59 + t86 * t58;
t65 = t24 * pkin(8) + t76;
t64 = pkin(4) * t45 + pkin(9) * t42;
t34 = t44 * pkin(5) + pkin(4);
t46 = -pkin(10) - pkin(9);
t61 = t34 * t45 - t42 * t46;
t49 = t25 * t40;
t12 = t38 * t26 - t39 * t49;
t60 = t12 * pkin(3) + t66;
t17 = -t24 * t42 + t40 * t45;
t5 = t11 * t42 - t39 * t68;
t7 = -t16 * t42 + t38 * t68;
t57 = g(1) * t7 + g(2) * t5 + g(3) * t17;
t18 = t24 * t45 + t40 * t42;
t6 = -t11 * t45 - t39 * t70;
t8 = t16 * t45 + t38 * t70;
t56 = g(1) * t8 + g(2) * t6 + g(3) * t18;
t55 = -g(1) * t16 + g(2) * t11 - g(3) * t24;
t15 = t39 * t26 + t38 * t49;
t54 = g(1) * t15 + g(2) * t12 - g(3) * t23;
t53 = -t38 * t71 - t39 * t43;
t52 = -t11 * pkin(8) + t60;
t51 = t53 * pkin(2);
t50 = t15 * pkin(3) + t51;
t48 = pkin(8) * t16 + t50;
t47 = -g(1) * (-t15 * t44 - t8 * t41) - g(2) * (-t12 * t44 - t6 * t41) - g(3) * (-t18 * t41 + t23 * t44);
t22 = -g(3) * t40 + (-g(1) * t38 + g(2) * t39) * t74;
t4 = t54 * t42;
t2 = -g(1) * (t15 * t35 - t8 * t36) - g(2) * (t12 * t35 - t6 * t36) - g(3) * (-t18 * t36 - t23 * t35);
t1 = -g(1) * (-t15 * t36 - t8 * t35) - g(2) * (-t12 * t36 - t6 * t35) - g(3) * (-t18 * t35 + t23 * t36);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t87 - g(3) * t62, -g(1) * (t38 * t80 - t39 * t86) - g(2) * (-t38 * t86 - t39 * t80) + g(3) * t74 * t43, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0, -g(1) * t51 - g(2) * t66 - g(3) * t32, 0, 0, 0, 0, 0, 0, -t54 * t45, t4, t55, -g(1) * t48 - g(2) * t52 - g(3) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t78 + t16 * t41) - g(2) * (-t11 * t41 + t12 * t78) - g(3) * (-t23 * t78 + t24 * t41) -g(1) * (-t15 * t79 + t16 * t44) - g(2) * (-t11 * t44 - t12 * t79) - g(3) * (t23 * t79 + t24 * t44) -t4, -g(1) * (t64 * t15 + t48) - g(2) * (t64 * t12 + t52) - g(3) * (-t64 * t23 + t65) 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t84 + t16 * t35) - g(2) * (-t11 * t35 + t12 * t84) - g(3) * (-t23 * t84 + t24 * t35) -g(1) * (-t15 * t85 + t16 * t36) - g(2) * (-t11 * t36 - t12 * t85) - g(3) * (t23 * t85 + t24 * t36) -t4, -g(1) * (t61 * t15 + t16 * t72 + t50) - g(2) * (-t72 * t11 + t61 * t12 + t60) - g(3) * (-t61 * t23 + t72 * t24 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t44, t57 * t41, -t56, -g(1) * (t7 * pkin(4) + t8 * pkin(9)) - g(2) * (t5 * pkin(4) + t6 * pkin(9)) - g(3) * (t17 * pkin(4) + t18 * pkin(9)) 0, 0, 0, 0, 0, 0, -t57 * t36, t57 * t35, -t56, -g(1) * (t7 * t34 - t8 * t46) - g(2) * (t5 * t34 - t6 * t46) - g(3) * (t17 * t34 - t18 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -g(1) * (t15 * t41 - t8 * t44) - g(2) * (t12 * t41 - t6 * t44) - g(3) * (-t18 * t44 - t23 * t41) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t47 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
