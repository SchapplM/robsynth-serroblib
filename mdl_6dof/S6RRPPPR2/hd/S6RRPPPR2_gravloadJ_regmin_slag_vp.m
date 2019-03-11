% Calculate minimal parameter regressor of gravitation load for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = qJ(2) + pkin(9);
t20 = sin(t26);
t22 = cos(t26);
t58 = t22 * pkin(3) + t20 * qJ(4);
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t12 = g(1) * t33 + g(2) * t31;
t57 = t12 * t20;
t30 = sin(qJ(2));
t55 = pkin(2) * t30;
t52 = g(3) * t22;
t29 = -qJ(3) - pkin(7);
t51 = pkin(4) - t29;
t25 = pkin(10) + qJ(6);
t19 = sin(t25);
t50 = t31 * t19;
t21 = cos(t25);
t49 = t31 * t21;
t27 = sin(pkin(10));
t48 = t31 * t27;
t28 = cos(pkin(10));
t47 = t31 * t28;
t46 = t33 * t19;
t45 = t33 * t21;
t44 = t33 * t27;
t43 = t33 * t28;
t42 = t33 * t29;
t40 = qJ(4) * t22;
t39 = t22 * qJ(5);
t32 = cos(qJ(2));
t23 = t32 * pkin(2);
t38 = t23 + t58;
t18 = t23 + pkin(1);
t14 = t33 * t18;
t37 = g(2) * (t58 * t33 + t14);
t36 = -pkin(3) * t20 - t55;
t11 = g(1) * t31 - g(2) * t33;
t35 = -t18 - t58;
t2 = -g(3) * t20 - t12 * t22;
t34 = -g(3) * t32 + t12 * t30;
t10 = t33 * t40;
t8 = t31 * t40;
t7 = t11 * t22;
t6 = -t20 * t50 + t45;
t5 = t20 * t49 + t46;
t4 = t20 * t46 + t49;
t3 = t20 * t45 - t50;
t1 = -t52 + t57;
t9 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t32, -t11 * t30, -t12, -g(1) * (-t31 * t18 - t42) - g(2) * (-t31 * t29 + t14) -t12, -t7, t11 * t20, g(1) * t42 - t37 + (-g(1) * t35 + g(2) * t29) * t31, -g(1) * (-t20 * t48 + t43) - g(2) * (t20 * t44 + t47) -g(1) * (-t20 * t47 - t44) - g(2) * (t20 * t43 - t48) t7, -t37 + (-g(1) * t51 - g(2) * t39) * t33 + (-g(1) * (t35 - t39) - g(2) * t51) * t31, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t30 + t12 * t32, 0, t34 * pkin(2), 0, -t1, t2, -g(1) * (t36 * t33 + t10) - g(2) * (t36 * t31 + t8) - g(3) * t38, t2 * t27, t2 * t28, t1, -g(1) * (-t33 * t55 + t10) - g(2) * (-t31 * t55 + t8) - g(3) * (t38 + t39) + (pkin(3) + qJ(5)) * t57, 0, 0, 0, 0, 0, t2 * t19, t2 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, -t11, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t21 * t52, g(1) * t4 - g(2) * t6 - t19 * t52;];
taug_reg  = t9;
