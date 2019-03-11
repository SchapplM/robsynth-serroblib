% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(2));
t32 = cos(qJ(2));
t42 = sin(pkin(10));
t44 = cos(pkin(6));
t37 = t44 * t42;
t43 = cos(pkin(10));
t11 = t31 * t43 + t32 * t37;
t38 = t44 * t43;
t9 = t31 * t42 - t32 * t38;
t55 = -g(1) * t11 - g(2) * t9;
t27 = sin(pkin(6));
t52 = g(3) * t27;
t23 = pkin(12) + qJ(6);
t19 = sin(t23);
t24 = pkin(11) + qJ(4);
t22 = cos(t24);
t51 = t19 * t22;
t21 = cos(t23);
t50 = t21 * t22;
t25 = sin(pkin(12));
t49 = t22 * t25;
t28 = cos(pkin(12));
t48 = t22 * t28;
t47 = t22 * t32;
t46 = t27 * t31;
t45 = t27 * t32;
t41 = t27 * t43;
t40 = t27 * t42;
t10 = t31 * t38 + t32 * t42;
t12 = -t31 * t37 + t32 * t43;
t39 = g(1) * t12 + g(2) * t10;
t20 = sin(t24);
t3 = t10 * t20 + t22 * t41;
t5 = t12 * t20 - t22 * t40;
t7 = t20 * t46 - t22 * t44;
t35 = g(1) * t5 + g(2) * t3 + g(3) * t7;
t4 = t10 * t22 - t20 * t41;
t6 = t12 * t22 + t20 * t40;
t8 = t20 * t44 + t22 * t46;
t34 = g(1) * t6 + g(2) * t4 + g(3) * t8;
t2 = g(3) * t45 + t55;
t33 = g(3) * t46 + t39;
t29 = cos(pkin(11));
t1 = t2 * t20;
t13 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t2, t33, -t2 * t29, t2 * sin(pkin(11)) -t33, -g(1) * (-pkin(2) * t11 + qJ(3) * t12) - g(2) * (-pkin(2) * t9 + qJ(3) * t10) - (pkin(2) * t32 + qJ(3) * t31) * t52, 0, 0, 0, 0, 0, -t2 * t22, t1, -g(1) * (-t11 * t48 + t12 * t25) - g(2) * (t10 * t25 - t48 * t9) - (t25 * t31 + t28 * t47) * t52, -g(1) * (t11 * t49 + t12 * t28) - g(2) * (t10 * t28 + t49 * t9) - (-t25 * t47 + t28 * t31) * t52, -t1 (t31 * t52 + t39) * (-pkin(8) - qJ(3)) + (-t32 * t52 - t55) * (t29 * pkin(3) + pkin(4) * t22 + qJ(5) * t20 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t11 * t50 + t12 * t19) - g(2) * (t10 * t19 - t50 * t9) - (t19 * t31 + t21 * t47) * t52, -g(1) * (t11 * t51 + t12 * t21) - g(2) * (t10 * t21 + t51 * t9) - (-t19 * t47 + t21 * t31) * t52; 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, t35 * t28, -t35 * t25, -t34, -g(1) * (-pkin(4) * t5 + qJ(5) * t6) - g(2) * (-pkin(4) * t3 + qJ(5) * t4) - g(3) * (-pkin(4) * t7 + qJ(5) * t8) 0, 0, 0, 0, 0, t35 * t21, -t35 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t21 - t19 * t6) - g(2) * (-t19 * t4 + t21 * t9) - g(3) * (-t19 * t8 - t21 * t45) -g(1) * (-t11 * t19 - t21 * t6) - g(2) * (-t19 * t9 - t21 * t4) - g(3) * (t19 * t45 - t21 * t8);];
taug_reg  = t13;
