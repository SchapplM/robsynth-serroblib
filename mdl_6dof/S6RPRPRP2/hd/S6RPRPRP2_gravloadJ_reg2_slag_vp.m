% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = qJ(1) + pkin(9);
t25 = sin(t31);
t27 = cos(t31);
t12 = g(1) * t27 + g(2) * t25;
t30 = qJ(3) + pkin(10);
t24 = sin(t30);
t63 = t12 * t24;
t20 = t24 * pkin(8);
t26 = cos(t30);
t62 = -t26 * pkin(4) - t20;
t34 = sin(qJ(3));
t61 = pkin(3) * t34;
t58 = g(3) * t24;
t35 = sin(qJ(1));
t57 = t35 * pkin(1);
t33 = sin(qJ(5));
t56 = t25 * t33;
t36 = cos(qJ(5));
t55 = t25 * t36;
t54 = t26 * t27;
t53 = t27 * t33;
t52 = t27 * t36;
t37 = cos(qJ(3));
t28 = t37 * pkin(3);
t23 = t28 + pkin(2);
t38 = cos(qJ(1));
t29 = t38 * pkin(1);
t51 = t27 * t23 + t29;
t50 = t28 - t62;
t49 = pkin(4) * t54 + t27 * t20 + t51;
t7 = t26 * t56 + t52;
t9 = t26 * t53 - t55;
t48 = g(1) * t7 - g(2) * t9;
t47 = -pkin(4) * t24 - t61;
t11 = g(1) * t25 - g(2) * t27;
t46 = g(1) * t35 - g(2) * t38;
t32 = -qJ(4) - pkin(7);
t45 = -t27 * t32 - t57;
t44 = pkin(5) * t36 + qJ(6) * t33;
t1 = g(1) * t9 + g(2) * t7 + t33 * t58;
t10 = t26 * t52 + t56;
t8 = t26 * t55 - t53;
t42 = g(1) * t10 + g(2) * t8 + t36 * t58;
t41 = -g(3) * t26 + t63;
t40 = -g(3) * t37 + t12 * t34;
t39 = (-g(1) * (-t23 + t62) + g(2) * t32) * t25;
t15 = pkin(8) * t54;
t13 = t25 * t26 * pkin(8);
t6 = t11 * t24;
t5 = t12 * t26 + t58;
t4 = t41 * t36;
t3 = t41 * t33;
t2 = g(1) * t8 - g(2) * t10;
t14 = [0, 0, 0, 0, 0, 0, t46, g(1) * t38 + g(2) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t46 * pkin(1), 0, 0, 0, 0, 0, 0, t11 * t37, -t11 * t34, -t12, -g(1) * (-t25 * pkin(2) + t27 * pkin(7) - t57) - g(2) * (t27 * pkin(2) + t25 * pkin(7) + t29) 0, 0, 0, 0, 0, 0, t11 * t26, -t6, -t12, -g(1) * (-t25 * t23 + t45) - g(2) * (-t25 * t32 + t51) 0, 0, 0, 0, 0, 0, t2, -t48, t6, -g(1) * t45 - g(2) * t49 + t39, 0, 0, 0, 0, 0, 0, t2, t6, t48, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t45) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t49) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(3) * t34 + t12 * t37, 0, 0, 0, 0, 0, 0, 0, 0, t41, t5, 0, t40 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (t47 * t27 + t15) - g(2) * (t47 * t25 + t13) - g(3) * t50, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * (-t27 * t61 + t15) - g(2) * (-t25 * t61 + t13) - g(3) * (t44 * t26 + t50) + (pkin(4) + t44) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t42, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t42, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t33 + qJ(6) * t36) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
