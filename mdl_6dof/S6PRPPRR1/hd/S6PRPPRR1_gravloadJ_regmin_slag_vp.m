% Calculate minimal parameter regressor of gravitation load for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(pkin(11));
t35 = sin(qJ(2));
t37 = cos(qJ(2));
t46 = cos(pkin(11));
t42 = -t35 * t28 + t37 * t46;
t26 = pkin(12) + qJ(5);
t25 = cos(t26);
t34 = sin(qJ(6));
t55 = t25 * t34;
t36 = cos(qJ(6));
t54 = t25 * t36;
t29 = sin(pkin(10));
t30 = sin(pkin(6));
t53 = t29 * t30;
t52 = t29 * t35;
t32 = cos(pkin(10));
t51 = t30 * t32;
t50 = t30 * t37;
t33 = cos(pkin(6));
t49 = t33 * t35;
t48 = t33 * t37;
t45 = t32 * t48;
t20 = -t37 * t28 - t35 * t46;
t18 = t20 * t33;
t7 = -t32 * t18 + t29 * t42;
t8 = -t29 * t18 - t32 * t42;
t43 = -t29 * t48 - t32 * t35;
t17 = t20 * t30;
t24 = sin(t26);
t41 = g(1) * (t24 * t8 + t25 * t53) + g(2) * (-t7 * t24 - t25 * t51) + g(3) * (t17 * t24 + t33 * t25);
t16 = t42 * t30;
t39 = t42 * t33;
t6 = t29 * t20 + t32 * t39;
t9 = t32 * t20 - t29 * t39;
t40 = g(1) * t9 + g(2) * t6 + g(3) * t16;
t38 = -g(1) * t43 - g(3) * t50;
t21 = pkin(2) * t45;
t15 = -g(3) * t33 + (-g(1) * t29 + g(2) * t32) * t30;
t12 = -t17 * t25 + t33 * t24;
t4 = t24 * t53 - t25 * t8;
t2 = -t24 * t51 + t7 * t25;
t1 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(2) * (t45 - t52) + t38, -g(1) * (t29 * t49 - t32 * t37) - g(2) * (-t29 * t37 - t32 * t49) + g(3) * t30 * t35, -g(2) * t21 + (g(2) * t52 + t38) * pkin(2), -t40 * cos(pkin(12)) t40 * sin(pkin(12)) g(1) * t8 - g(2) * t7 + g(3) * t17, -g(1) * (t43 * pkin(2) + t9 * pkin(3) - t8 * qJ(4)) - g(2) * (-pkin(2) * t52 + t6 * pkin(3) + qJ(4) * t7 + t21) - g(3) * (pkin(2) * t50 + t16 * pkin(3) - t17 * qJ(4)) 0, 0, 0, 0, 0, -t40 * t25, t40 * t24, 0, 0, 0, 0, 0, -g(1) * (-t8 * t34 + t9 * t54) - g(2) * (t34 * t7 + t6 * t54) - g(3) * (t16 * t54 - t17 * t34) -g(1) * (-t8 * t36 - t9 * t55) - g(2) * (t36 * t7 - t6 * t55) - g(3) * (-t16 * t55 - t17 * t36); 0, 0, 0, 0, t15, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, g(1) * t4 + g(2) * t2 + g(3) * t12, 0, 0, 0, 0, 0, -t41 * t36, t41 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t34 - t9 * t36) - g(2) * (-t2 * t34 - t6 * t36) - g(3) * (-t12 * t34 - t16 * t36) -g(1) * (t9 * t34 - t4 * t36) - g(2) * (-t2 * t36 + t6 * t34) - g(3) * (-t12 * t36 + t16 * t34);];
taug_reg  = t1;
