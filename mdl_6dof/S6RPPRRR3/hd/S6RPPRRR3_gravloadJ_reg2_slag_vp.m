% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t32 = cos(qJ(5));
t21 = t32 * pkin(5) + pkin(4);
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t35 = -pkin(9) - pkin(8);
t39 = t30 * t21 + t33 * t35;
t27 = qJ(1) + pkin(10);
t23 = cos(t27);
t22 = sin(t27);
t61 = g(1) * t22;
t65 = g(2) * t23 - t61;
t29 = sin(qJ(5));
t51 = t29 * t30;
t10 = t22 * t32 + t23 * t51;
t57 = g(3) * t33;
t8 = -t22 * t51 + t23 * t32;
t64 = -g(1) * t8 - g(2) * t10 + t29 * t57;
t37 = -g(3) * t30 - t33 * t65;
t62 = -pkin(2) - pkin(7);
t56 = t33 * pkin(8);
t55 = t23 * t29;
t54 = t23 * t30;
t28 = qJ(5) + qJ(6);
t24 = sin(t28);
t53 = t24 * t30;
t25 = cos(t28);
t52 = t25 * t30;
t49 = t30 * t32;
t34 = cos(qJ(1));
t46 = t34 * pkin(1) + t23 * pkin(2) + t22 * qJ(3);
t31 = sin(qJ(1));
t45 = -t31 * pkin(1) + t23 * qJ(3);
t44 = t23 * pkin(7) + t46;
t42 = t30 * pkin(4) - t56;
t14 = g(1) * t23 + g(2) * t22;
t41 = g(1) * t31 - g(2) * t34;
t38 = g(2) * t44;
t12 = t14 * t33;
t11 = -t22 * t29 + t23 * t49;
t9 = t22 * t49 + t55;
t7 = -g(2) * t54 + t30 * t61 + t57;
t6 = -t22 * t24 + t23 * t52;
t5 = t22 * t25 + t23 * t53;
t4 = t22 * t52 + t23 * t24;
t3 = -t22 * t53 + t23 * t25;
t2 = g(1) * t4 - g(2) * t6 + t25 * t57;
t1 = -g(1) * t3 - g(2) * t5 + t24 * t57;
t13 = [0, 0, 0, 0, 0, 0, t41, g(1) * t34 + g(2) * t31, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t14, 0, t41 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t65, -t14, -g(1) * (-t22 * pkin(2) + t45) - g(2) * t46, 0, 0, 0, 0, 0, 0, -t14 * t30, -t12, -t65, -g(1) * (t62 * t22 + t45) - t38, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9, g(1) * t10 - g(2) * t8, t12, -g(1) * (pkin(4) * t54 - t23 * t56 + t45) - t38 + (-g(1) * t62 - g(2) * t42) * t22, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t12, -g(1) * (t39 * t23 + t45) - g(2) * (pkin(5) * t55 + t44) + (-g(1) * (-pkin(5) * t29 + t62) - g(2) * t39) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t32, t37 * t29, -t7, g(3) * t42 + t65 * (pkin(4) * t33 + pkin(8) * t30) 0, 0, 0, 0, 0, 0, -t37 * t25, t37 * t24, -t7, g(3) * t39 + t65 * (t21 * t33 - t30 * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, g(1) * t9 - g(2) * t11 + t32 * t57, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t64 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t13;
