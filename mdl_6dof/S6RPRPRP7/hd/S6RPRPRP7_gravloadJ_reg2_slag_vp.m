% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = -qJ(4) - pkin(7);
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t32 = sin(qJ(3));
t60 = pkin(3) * t32;
t61 = t33 * t30 + t36 * t60;
t54 = g(2) * t36;
t14 = g(1) * t33 - t54;
t28 = qJ(3) + pkin(9);
t22 = sin(t28);
t34 = cos(qJ(5));
t48 = t34 * t36;
t31 = sin(qJ(5));
t50 = t33 * t31;
t10 = -t22 * t50 + t48;
t49 = t33 * t34;
t51 = t31 * t36;
t12 = t22 * t51 + t49;
t23 = cos(t28);
t52 = g(3) * t23;
t1 = -g(1) * t10 - g(2) * t12 + t31 * t52;
t8 = -g(3) * t22 + t14 * t23;
t35 = cos(qJ(3));
t59 = pkin(3) * t35;
t58 = pkin(5) * t31;
t47 = t36 * pkin(1) + t33 * qJ(2);
t45 = t33 * t60 + t47;
t25 = t36 * qJ(2);
t44 = -t33 * pkin(1) + t25;
t43 = pkin(4) * t23 + pkin(8) * t22;
t42 = pkin(4) * t22 - pkin(8) * t23;
t15 = g(1) * t36 + g(2) * t33;
t21 = pkin(5) * t34 + pkin(4);
t29 = -qJ(6) - pkin(8);
t41 = t21 * t23 - t22 * t29;
t40 = t21 * t22 + t23 * t29;
t39 = t44 + t61;
t38 = -t30 * t36 + t45;
t37 = g(3) * t32 - t14 * t35;
t18 = t33 * t59;
t13 = t22 * t48 - t50;
t11 = t22 * t49 + t51;
t9 = t15 * t23;
t7 = t14 * t22 + t52;
t6 = t8 * t34;
t5 = t8 * t31;
t4 = -g(1) * t13 - g(2) * t11;
t3 = g(1) * t12 - g(2) * t10;
t2 = g(1) * t11 - g(2) * t13 + t34 * t52;
t16 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, -g(1) * t44 - g(2) * t47, 0, 0, 0, 0, 0, 0, -t15 * t32, -t15 * t35, t14, -g(1) * (t25 + (-pkin(1) - pkin(7)) * t33) - g(2) * (pkin(7) * t36 + t47) 0, 0, 0, 0, 0, 0, -t15 * t22, -t9, t14, -g(1) * t39 - g(2) * t38, 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * (t42 * t36 + t39) - g(2) * (t42 * t33 + t38) 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * (t25 + t61) - g(2) * t45 + (-g(1) * t40 - g(2) * (-t30 + t58)) * t36 + (-g(1) * (-pkin(1) - t58) - g(2) * t40) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(3) * t35 + t14 * t32, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, t37 * pkin(3), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(1) * (t43 * t33 + t18) - g(3) * (-t42 - t60) - (-t43 - t59) * t54, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(1) * (t41 * t33 + t18) - g(3) * (-t40 - t60) - (-t41 - t59) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t16;
