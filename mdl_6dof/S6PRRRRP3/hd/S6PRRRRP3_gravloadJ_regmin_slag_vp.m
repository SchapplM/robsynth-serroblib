% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t43 = sin(pkin(11));
t45 = cos(pkin(6));
t39 = t45 * t43;
t44 = cos(pkin(11));
t10 = t31 * t44 + t34 * t39;
t56 = g(1) * t10;
t11 = -t31 * t39 + t34 * t44;
t55 = g(1) * t11;
t28 = sin(pkin(6));
t54 = g(3) * t28;
t27 = qJ(4) + qJ(5);
t23 = sin(t27);
t29 = sin(qJ(4));
t15 = pkin(4) * t29 + pkin(5) * t23;
t53 = pkin(8) + t15;
t33 = cos(qJ(3));
t52 = t23 * t33;
t24 = cos(t27);
t51 = t24 * t33;
t50 = t28 * t31;
t49 = t28 * t34;
t48 = t29 * t33;
t32 = cos(qJ(4));
t47 = t32 * t33;
t46 = t33 * t34;
t16 = pkin(4) * t32 + pkin(5) * t24;
t42 = t28 * t44;
t41 = t28 * t43;
t40 = t45 * t44;
t14 = pkin(3) + t16;
t26 = -qJ(6) - pkin(10) - pkin(9);
t30 = sin(qJ(3));
t38 = t14 * t33 - t26 * t30 + pkin(2);
t12 = t30 * t50 - t33 * t45;
t9 = t31 * t40 + t34 * t43;
t4 = t30 * t9 + t33 * t42;
t6 = t11 * t30 - t33 * t41;
t37 = g(1) * t6 + g(2) * t4 + g(3) * t12;
t13 = t30 * t45 + t33 * t50;
t5 = -t30 * t42 + t33 * t9;
t7 = t11 * t33 + t30 * t41;
t36 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t8 = t31 * t43 - t34 * t40;
t35 = -g(2) * t8 + g(3) * t49 - t56;
t1 = -g(1) * (t10 * t24 - t23 * t7) - g(2) * (-t23 * t5 + t24 * t8) - g(3) * (-t13 * t23 - t24 * t49);
t3 = t35 * t30;
t2 = -g(1) * (-t10 * t23 - t24 * t7) - g(2) * (-t23 * t8 - t24 * t5) - g(3) * (-t13 * t24 + t23 * t49);
t17 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t35, g(2) * t9 + g(3) * t50 + t55, 0, 0, 0, 0, 0, -t35 * t33, t3, 0, 0, 0, 0, 0, -g(1) * (-t10 * t47 + t11 * t29) - g(2) * (t29 * t9 - t47 * t8) - (t29 * t31 + t32 * t46) * t54, -g(1) * (t10 * t48 + t11 * t32) - g(2) * (t32 * t9 + t48 * t8) - (-t29 * t46 + t31 * t32) * t54, 0, 0, 0, 0, 0, -g(1) * (-t10 * t51 + t11 * t23) - g(2) * (t23 * t9 - t51 * t8) - (t23 * t31 + t24 * t46) * t54, -g(1) * (t10 * t52 + t11 * t24) - g(2) * (t24 * t9 + t52 * t8) - (-t23 * t46 + t24 * t31) * t54, -t3, -g(2) * (-t38 * t8 + t53 * t9) - t53 * t55 + t38 * t56 - (t31 * t53 + t34 * t38) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t36, 0, 0, 0, 0, 0, t37 * t32, -t37 * t29, 0, 0, 0, 0, 0, t37 * t24, -t37 * t23, -t36, -g(1) * (-t14 * t6 - t26 * t7) - g(2) * (-t14 * t4 - t26 * t5) - g(3) * (-t12 * t14 - t13 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t32 - t29 * t7) - g(2) * (-t29 * t5 + t32 * t8) - g(3) * (-t13 * t29 - t32 * t49) -g(1) * (-t10 * t29 - t32 * t7) - g(2) * (-t29 * t8 - t32 * t5) - g(3) * (-t13 * t32 + t29 * t49) 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t10 * t16 - t15 * t7) - g(2) * (-t15 * t5 + t16 * t8) - g(3) * (-t13 * t15 - t16 * t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37;];
taug_reg  = t17;
