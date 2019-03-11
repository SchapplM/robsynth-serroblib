% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t35 = cos(pkin(11));
t22 = t35 * pkin(4) + pkin(3);
t36 = -pkin(8) - qJ(4);
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t43 = t39 * t22 - t37 * t36;
t32 = pkin(11) + qJ(5);
t25 = cos(t32);
t15 = pkin(5) * t25 + t22;
t31 = -pkin(9) + t36;
t45 = t39 * t15 - t37 * t31;
t33 = qJ(1) + pkin(10);
t24 = sin(t33);
t26 = cos(t33);
t14 = g(1) * t26 + g(2) * t24;
t23 = sin(t32);
t64 = g(3) * t37;
t61 = t24 * t39;
t7 = t23 * t61 + t26 * t25;
t59 = t26 * t39;
t9 = -t23 * t59 + t24 * t25;
t70 = -g(1) * t9 + g(2) * t7 + t23 * t64;
t11 = -g(3) * t39 + t14 * t37;
t67 = g(1) * t24;
t34 = sin(pkin(11));
t62 = t34 * pkin(4);
t60 = t26 * t34;
t58 = t34 * t39;
t57 = t35 * t39;
t40 = cos(qJ(1));
t51 = t40 * pkin(1) + t26 * pkin(2) + t24 * pkin(7);
t38 = sin(qJ(1));
t50 = -t38 * pkin(1) + t26 * pkin(7);
t49 = -g(2) * t26 + t67;
t48 = g(1) * t38 - g(2) * t40;
t47 = t39 * pkin(3) + t37 * qJ(4);
t27 = qJ(6) + t32;
t21 = cos(t27);
t20 = sin(t27);
t16 = pkin(5) * t23 + t62;
t13 = t49 * t37;
t12 = t14 * t39 + t64;
t10 = t24 * t23 + t25 * t59;
t8 = t26 * t23 - t25 * t61;
t6 = t24 * t20 + t21 * t59;
t5 = -t20 * t59 + t24 * t21;
t4 = t26 * t20 - t21 * t61;
t3 = t20 * t61 + t26 * t21;
t2 = g(1) * t6 - g(2) * t4 + t21 * t64;
t1 = -g(1) * t5 + g(2) * t3 + t20 * t64;
t17 = [0, 0, 0, 0, 0, 0, t48, g(1) * t40 + g(2) * t38, 0, 0, 0, 0, 0, 0, 0, 0, t49, t14, 0, t48 * pkin(1), 0, 0, 0, 0, 0, 0, t49 * t39, -t13, -t14, -g(1) * (-t24 * pkin(2) + t50) - g(2) * t51, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t57 + t60) - g(2) * (t24 * t34 + t26 * t57) -g(1) * (t24 * t58 + t26 * t35) - g(2) * (t24 * t35 - t26 * t58) t13, -g(1) * t50 - g(2) * (t47 * t26 + t51) - (-pkin(2) - t47) * t67, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t13, -g(1) * (pkin(4) * t60 + t50) - g(2) * (t43 * t26 + t51) + (-g(1) * (-pkin(2) - t43) - g(2) * t62) * t24, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t13, -g(1) * (t26 * t16 + t50) - g(2) * (t45 * t26 + t51) + (-g(1) * (-pkin(2) - t45) - g(2) * t16) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t35, -t11 * t34, -t12, -g(3) * t47 + t14 * (pkin(3) * t37 - qJ(4) * t39) 0, 0, 0, 0, 0, 0, t11 * t25, -t11 * t23, -t12, -g(3) * t43 + t14 * (t22 * t37 + t36 * t39) 0, 0, 0, 0, 0, 0, t11 * t21, -t11 * t20, -t12, -g(3) * t45 + t14 * (t15 * t37 + t31 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, g(1) * t10 - g(2) * t8 + t25 * t64, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t70 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t17;
