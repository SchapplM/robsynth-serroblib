% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:32:46
% EndTime: 2019-05-05 04:32:47
% DurationCPUTime: 0.33s
% Computational Cost: add. (312->82), mult. (526->152), div. (0->0), fcn. (656->14), ass. (0->48)
t22 = sin(pkin(11));
t24 = cos(pkin(11));
t31 = cos(qJ(2));
t28 = sin(qJ(2));
t40 = cos(pkin(6));
t39 = t28 * t40;
t10 = t22 * t31 + t24 * t39;
t12 = -t22 * t39 + t24 * t31;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t23 = sin(pkin(6));
t44 = t23 * t30;
t45 = t23 * t28;
t55 = -g(1) * (-t12 * t27 + t22 * t44) - g(2) * (-t10 * t27 - t24 * t44) - g(3) * (-t27 * t45 + t40 * t30);
t54 = g(3) * t23;
t20 = qJ(3) + pkin(12);
t17 = cos(t20);
t21 = qJ(5) + qJ(6);
t18 = sin(t21);
t53 = t17 * t18;
t19 = cos(t21);
t52 = t17 * t19;
t26 = sin(qJ(5));
t51 = t17 * t26;
t29 = cos(qJ(5));
t50 = t17 * t29;
t49 = t17 * t31;
t48 = t22 * t23;
t47 = t23 * t24;
t46 = t23 * t27;
t43 = t23 * t31;
t42 = t26 * t31;
t41 = t29 * t31;
t38 = t31 * t40;
t16 = sin(t20);
t37 = g(1) * (-t12 * t16 + t17 * t48) + g(2) * (-t10 * t16 - t17 * t47) + g(3) * (-t16 * t45 + t40 * t17);
t11 = t22 * t38 + t24 * t28;
t9 = t22 * t28 - t24 * t38;
t34 = -g(1) * t11 - g(2) * t9 + g(3) * t43;
t33 = g(1) * t12 + g(2) * t10 + g(3) * t45;
t25 = -qJ(4) - pkin(8);
t15 = t30 * pkin(3) + pkin(2);
t8 = t40 * t16 + t17 * t45;
t6 = t12 * t17 + t16 * t48;
t4 = t10 * t17 - t16 * t47;
t2 = -g(1) * (-t11 * t18 - t6 * t19) - g(2) * (-t9 * t18 - t4 * t19) - g(3) * (t18 * t43 - t8 * t19);
t1 = -g(1) * (t11 * t19 - t6 * t18) - g(2) * (-t4 * t18 + t9 * t19) - g(3) * (-t8 * t18 - t19 * t43);
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t34, t33, 0, 0, 0, 0, 0, -t34 * t30, t34 * t27, -t33, -g(1) * (-t11 * t15 - t12 * t25) - g(2) * (-t10 * t25 - t9 * t15) - (t15 * t31 - t25 * t28) * t54, 0, 0, 0, 0, 0, -g(1) * (-t11 * t50 + t12 * t26) - g(2) * (t10 * t26 - t9 * t50) - (t17 * t41 + t26 * t28) * t54, -g(1) * (t11 * t51 + t12 * t29) - g(2) * (t10 * t29 + t9 * t51) - (-t17 * t42 + t28 * t29) * t54, 0, 0, 0, 0, 0, -g(1) * (-t11 * t52 + t12 * t18) - g(2) * (t10 * t18 - t9 * t52) - (t18 * t28 + t19 * t49) * t54, -g(1) * (t11 * t53 + t12 * t19) - g(2) * (t10 * t19 + t9 * t53) - (-t18 * t49 + t19 * t28) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -g(1) * (-t12 * t30 - t22 * t46) - g(2) * (-t10 * t30 + t24 * t46) - g(3) * (-t40 * t27 - t28 * t44) 0, t55 * pkin(3), 0, 0, 0, 0, 0, -t37 * t29, t37 * t26, 0, 0, 0, 0, 0, -t37 * t19, t37 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t29 - t6 * t26) - g(2) * (-t4 * t26 + t9 * t29) - g(3) * (-t23 * t41 - t8 * t26) -g(1) * (-t11 * t26 - t6 * t29) - g(2) * (-t9 * t26 - t4 * t29) - g(3) * (t23 * t42 - t8 * t29) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
