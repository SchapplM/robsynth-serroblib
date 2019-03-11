% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t12 = g(1) * t31 + g(2) * t28;
t26 = sin(qJ(4));
t23 = qJ(2) + pkin(10);
t18 = sin(t23);
t48 = g(3) * t18;
t19 = cos(t23);
t29 = cos(qJ(4));
t39 = t31 * t29;
t44 = t28 * t26;
t7 = t19 * t44 + t39;
t40 = t31 * t26;
t43 = t28 * t29;
t9 = -t19 * t40 + t43;
t53 = -g(1) * t9 + g(2) * t7 + t26 * t48;
t34 = -g(3) * t19 + t12 * t18;
t20 = qJ(4) + pkin(11) + qJ(6);
t14 = sin(t20);
t46 = t28 * t14;
t15 = cos(t20);
t45 = t28 * t15;
t42 = t31 * t14;
t41 = t31 * t15;
t25 = -qJ(3) - pkin(7);
t37 = pkin(4) * t26 - t25;
t11 = g(1) * t28 - g(2) * t31;
t16 = t29 * pkin(4) + pkin(3);
t24 = -qJ(5) - pkin(8);
t36 = t19 * t16 - t18 * t24;
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t32 = -g(3) * t30 + t12 * t27;
t21 = t30 * pkin(2);
t17 = t21 + pkin(1);
t13 = t31 * t17;
t10 = t19 * t39 + t44;
t8 = -t19 * t43 + t40;
t6 = t19 * t41 + t46;
t5 = -t19 * t42 + t45;
t4 = -t19 * t45 + t42;
t3 = t19 * t46 + t41;
t2 = g(1) * t6 - g(2) * t4 + t15 * t48;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t48;
t22 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t30, -t11 * t27, -t12, -g(1) * (-t28 * t17 - t31 * t25) - g(2) * (-t28 * t25 + t13) 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t11 * t18, -g(2) * t13 + (-g(1) * t37 - g(2) * t36) * t31 + (-g(1) * (-t17 - t36) - g(2) * t37) * t28, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t12 * t30, 0, t32 * pkin(2), 0, 0, 0, 0, 0, t34 * t29, -t34 * t26, -t12 * t19 - t48, -g(3) * (t21 + t36) + t12 * (pkin(2) * t27 + t16 * t18 + t19 * t24) 0, 0, 0, 0, 0, t34 * t15, -t34 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(1) * t10 - g(2) * t8 + t29 * t48, 0, t53 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t22;
