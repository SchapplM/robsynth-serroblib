% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:35
% EndTime: 2019-12-31 21:40:37
% DurationCPUTime: 0.69s
% Computational Cost: add. (426->133), mult. (990->212), div. (0->0), fcn. (1205->12), ass. (0->68)
t42 = sin(qJ(2));
t43 = sin(qJ(1));
t69 = cos(pkin(5));
t82 = cos(qJ(2));
t54 = t69 * t82;
t83 = cos(qJ(1));
t19 = t43 * t42 - t83 * t54;
t36 = pkin(10) + qJ(5);
t33 = sin(t36);
t34 = cos(t36);
t58 = t42 * t69;
t20 = t43 * t82 + t83 * t58;
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t38 = sin(pkin(5));
t65 = t38 * t83;
t8 = t20 * t44 - t41 * t65;
t86 = -t19 * t34 + t8 * t33;
t85 = t19 * t33 + t8 * t34;
t84 = g(3) * t38;
t79 = t33 * t44;
t78 = t34 * t44;
t37 = sin(pkin(10));
t77 = t37 * t42;
t76 = t37 * t44;
t75 = t38 * t42;
t74 = t38 * t43;
t73 = t38 * t44;
t39 = cos(pkin(10));
t72 = t39 * t44;
t64 = t38 * t82;
t71 = pkin(2) * t64 + pkin(8) * t75;
t70 = t83 * pkin(1) + pkin(7) * t74;
t22 = -t43 * t58 + t83 * t82;
t68 = t22 * pkin(2) + t70;
t67 = pkin(4) * t37 + pkin(8);
t66 = g(3) * t71;
t63 = t41 * t82;
t62 = t44 * t82;
t61 = -t43 * pkin(1) + pkin(7) * t65;
t13 = t19 * pkin(2);
t60 = t20 * pkin(8) - t13;
t21 = t83 * t42 + t43 * t54;
t15 = t21 * pkin(2);
t59 = t22 * pkin(8) - t15;
t57 = -t20 * pkin(2) + t61;
t11 = t22 * t41 - t43 * t73;
t7 = t20 * t41 + t44 * t65;
t56 = -g(1) * t7 + g(2) * t11;
t55 = g(1) * t19 - g(2) * t21;
t53 = -pkin(3) * t44 - qJ(4) * t41;
t32 = t39 * pkin(4) + pkin(3);
t40 = -pkin(9) - qJ(4);
t52 = -t32 * t44 + t40 * t41;
t51 = t21 * pkin(8) + t68;
t50 = g(1) * t83 + g(2) * t43;
t49 = -t19 * pkin(8) + t57;
t17 = t41 * t75 - t69 * t44;
t48 = g(1) * t11 + g(2) * t7 + g(3) * t17;
t12 = t22 * t44 + t41 * t74;
t18 = t69 * t41 + t42 * t73;
t47 = g(1) * t12 + g(2) * t8 + g(3) * t18;
t46 = g(1) * t22 + g(2) * t20 + g(3) * t75;
t45 = -g(1) * t21 - g(2) * t19 + g(3) * t64;
t6 = t45 * t41;
t5 = t12 * t34 + t21 * t33;
t4 = -t12 * t33 + t21 * t34;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t43 - g(2) * t83, t50, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t20 - g(2) * t22, -t55, -t50 * t38, -g(1) * t61 - g(2) * t70, 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t12, t56, t55, -g(1) * t49 - g(2) * t51, 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t37 - t39 * t8) - g(2) * (t12 * t39 + t21 * t37), -g(1) * (-t19 * t39 + t37 * t8) - g(2) * (-t12 * t37 + t21 * t39), -t56, -g(1) * (-pkin(3) * t8 - qJ(4) * t7 + t49) - g(2) * (t12 * pkin(3) + t11 * qJ(4) + t51), 0, 0, 0, 0, 0, 0, g(1) * t85 - g(2) * t5, -g(1) * t86 - g(2) * t4, -t56, -g(1) * (-t67 * t19 - t32 * t8 + t40 * t7 + t57) - g(2) * (-t11 * t40 + t12 * t32 + t67 * t21 + t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t46, 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t44, t6, -t46, -g(1) * t59 - g(2) * t60 - t66, 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t72 + t22 * t37) - g(2) * (-t19 * t72 + t20 * t37) - (t39 * t62 + t77) * t84, -g(1) * (t21 * t76 + t22 * t39) - g(2) * (t19 * t76 + t20 * t39) - (-t37 * t62 + t39 * t42) * t84, -t6, -g(1) * (t53 * t21 + t59) - g(2) * (t53 * t19 + t60) - g(3) * ((pkin(3) * t62 + qJ(4) * t63) * t38 + t71), 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t78 + t22 * t33) - g(2) * (-t19 * t78 + t20 * t33) - (t33 * t42 + t34 * t62) * t84, -g(1) * (t21 * t79 + t22 * t34) - g(2) * (t19 * t79 + t20 * t34) - (-t33 * t62 + t34 * t42) * t84, -t6, -g(1) * (t52 * t21 + t67 * t22 - t15) - g(2) * (t52 * t19 + t67 * t20 - t13) - t66 - (pkin(4) * t77 + t32 * t62 - t40 * t63) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t39, -t48 * t37, -t47, -g(1) * (-t11 * pkin(3) + t12 * qJ(4)) - g(2) * (-t7 * pkin(3) + t8 * qJ(4)) - g(3) * (-t17 * pkin(3) + t18 * qJ(4)), 0, 0, 0, 0, 0, 0, t48 * t34, -t48 * t33, -t47, -g(1) * (-t11 * t32 - t12 * t40) - g(2) * (-t7 * t32 - t8 * t40) - g(3) * (-t17 * t32 - t18 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t86 - g(3) * (-t18 * t33 - t34 * t64), g(1) * t5 + g(2) * t85 - g(3) * (-t18 * t34 + t33 * t64), 0, 0;];
taug_reg = t1;
