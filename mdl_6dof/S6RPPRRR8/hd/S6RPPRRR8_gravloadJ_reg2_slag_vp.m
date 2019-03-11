% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t32 = -pkin(7) - qJ(3);
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t30 = sin(pkin(10));
t64 = pkin(3) * t30;
t68 = t34 * t32 + t36 * t64;
t67 = -g(1) * t34 + g(2) * t36;
t28 = pkin(10) + qJ(4);
t20 = sin(t28);
t33 = sin(qJ(5));
t51 = t36 * t33;
t35 = cos(qJ(5));
t54 = t34 * t35;
t11 = t20 * t51 + t54;
t21 = cos(t28);
t58 = g(3) * t21;
t50 = t36 * t35;
t55 = t34 * t33;
t9 = -t20 * t55 + t50;
t66 = -g(1) * t9 - g(2) * t11 + t33 * t58;
t39 = -g(3) * t20 - t21 * t67;
t63 = pkin(5) * t33;
t29 = qJ(5) + qJ(6);
t22 = sin(t29);
t57 = t34 * t22;
t23 = cos(t29);
t56 = t34 * t23;
t53 = t36 * t22;
t52 = t36 * t23;
t49 = t36 * pkin(1) + t34 * qJ(2);
t47 = t34 * t64 + t49;
t25 = t36 * qJ(2);
t46 = -t34 * pkin(1) + t25;
t44 = t20 * pkin(4) - t21 * pkin(8);
t14 = g(1) * t36 + g(2) * t34;
t19 = t35 * pkin(5) + pkin(4);
t37 = -pkin(9) - pkin(8);
t42 = t20 * t19 + t21 * t37;
t41 = t46 + t68;
t40 = -t36 * t32 + t47;
t12 = t20 * t50 - t55;
t10 = t20 * t54 + t51;
t8 = t14 * t21;
t7 = t20 * t52 - t57;
t6 = t20 * t53 + t56;
t5 = t20 * t56 + t53;
t4 = -t20 * t57 + t52;
t3 = -t67 * t20 + t58;
t2 = g(1) * t5 - g(2) * t7 + t23 * t58;
t1 = -g(1) * t4 - g(2) * t6 + t22 * t58;
t13 = [0, 0, 0, 0, 0, 0, -t67, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t14, -g(1) * t46 - g(2) * t49, 0, 0, 0, 0, 0, 0, -t14 * t30, -t14 * cos(pkin(10)) -t67, -g(1) * (t25 + (-pkin(1) - qJ(3)) * t34) - g(2) * (t36 * qJ(3) + t49) 0, 0, 0, 0, 0, 0, -t14 * t20, -t8, -t67, -g(1) * t41 - g(2) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t8, -g(1) * (t44 * t36 + t41) - g(2) * (t44 * t34 + t40) 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t5, g(1) * t6 - g(2) * t4, t8, -g(1) * (t25 + t68) - g(2) * t47 + (-g(1) * t42 - g(2) * (-t32 + t63)) * t36 + (-g(1) * (-pkin(1) - t63) - g(2) * t42) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t35, t39 * t33, -t3, g(3) * t44 + t67 * (pkin(4) * t21 + pkin(8) * t20) 0, 0, 0, 0, 0, 0, -t39 * t23, t39 * t22, -t3, g(3) * t42 + t67 * (t19 * t21 - t20 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, g(1) * t10 - g(2) * t12 + t35 * t58, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t66 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t13;
