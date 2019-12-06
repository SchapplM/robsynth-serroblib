% Calculate minimal parameter regressor of gravitation load for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t25 = cos(pkin(10));
t23 = sin(pkin(11));
t46 = sin(pkin(10));
t42 = t46 * t23;
t48 = cos(pkin(11));
t43 = t25 * t48;
t50 = cos(pkin(5));
t33 = -t50 * t43 + t42;
t24 = sin(pkin(5));
t47 = sin(pkin(6));
t45 = t24 * t47;
t49 = cos(pkin(6));
t56 = t25 * t45 + t33 * t49;
t38 = t46 * t48;
t53 = t25 * t23;
t34 = t50 * t38 + t53;
t44 = t24 * t46;
t55 = t34 * t49 - t47 * t44;
t54 = cos(qJ(3));
t26 = sin(qJ(5));
t30 = cos(qJ(4));
t52 = t26 * t30;
t29 = cos(qJ(5));
t51 = t29 * t30;
t40 = t50 * t47;
t39 = t49 * t48;
t28 = sin(qJ(3));
t12 = t28 * t40 + (t54 * t23 + t28 * t39) * t24;
t13 = -t25 * t24 * t49 + t33 * t47;
t14 = t34 * t47 + t49 * t44;
t17 = -t48 * t45 + t50 * t49;
t27 = sin(qJ(4));
t18 = t50 * t53 + t38;
t6 = t18 * t54 - t56 * t28;
t19 = -t50 * t42 + t43;
t8 = t19 * t54 - t55 * t28;
t36 = g(1) * (t14 * t30 - t8 * t27) + g(2) * (t13 * t30 - t6 * t27) + g(3) * (-t12 * t27 + t17 * t30);
t11 = -t54 * t40 + (t23 * t28 - t39 * t54) * t24;
t5 = t18 * t28 + t56 * t54;
t7 = t19 * t28 + t55 * t54;
t35 = g(1) * t7 + g(2) * t5 + g(3) * t11;
t10 = t12 * t30 + t17 * t27;
t4 = t14 * t27 + t8 * t30;
t2 = t13 * t27 + t6 * t30;
t1 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(3) * t50 + (-t46 * g(1) + g(2) * t25) * t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t35, g(1) * t8 + g(2) * t6 + g(3) * t12, 0, 0, 0, 0, 0, t35 * t30, -t35 * t27, 0, 0, 0, 0, 0, -g(1) * (t8 * t26 - t7 * t51) - g(2) * (t6 * t26 - t5 * t51) - g(3) * (-t11 * t51 + t12 * t26), -g(1) * (t8 * t29 + t7 * t52) - g(2) * (t6 * t29 + t5 * t52) - g(3) * (t11 * t52 + t12 * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, g(1) * t4 + g(2) * t2 + g(3) * t10, 0, 0, 0, 0, 0, -t36 * t29, t36 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t26 + t7 * t29) - g(2) * (-t2 * t26 + t5 * t29) - g(3) * (-t10 * t26 + t11 * t29), -g(1) * (-t7 * t26 - t4 * t29) - g(2) * (-t2 * t29 - t5 * t26) - g(3) * (-t10 * t29 - t11 * t26);];
taug_reg = t1;
