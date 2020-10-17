% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:20:04
% EndTime: 2019-05-06 01:20:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (396->63), mult. (369->96), div. (0->0), fcn. (402->10), ass. (0->48)
t27 = qJ(1) + pkin(10);
t23 = sin(t27);
t24 = cos(t27);
t43 = g(1) * t24 + g(2) * t23;
t32 = cos(qJ(4));
t29 = sin(qJ(4));
t33 = cos(qJ(3));
t49 = t29 * t33;
t14 = t23 * t49 + t24 * t32;
t16 = t23 * t32 - t24 * t49;
t30 = sin(qJ(3));
t55 = g(3) * t30;
t60 = -g(1) * t16 + g(2) * t14 + t29 * t55;
t38 = -g(3) * t33 + t43 * t30;
t28 = qJ(4) + qJ(5);
t25 = sin(t28);
t53 = t25 * t30;
t52 = t25 * t33;
t26 = cos(t28);
t51 = t26 * t30;
t50 = t26 * t33;
t35 = -pkin(9) - pkin(8);
t48 = t30 * t35;
t47 = t32 * t33;
t45 = pkin(4) * t29 + pkin(7);
t11 = -t23 * t26 + t24 * t52;
t9 = t23 * t52 + t24 * t26;
t44 = g(1) * t9 - g(2) * t11;
t42 = g(1) * t23 - g(2) * t24;
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t41 = g(1) * t31 - g(2) * t34;
t22 = t32 * pkin(4) + pkin(3);
t40 = t33 * t22 + pkin(2) - t48;
t39 = pkin(5) * t26 + qJ(6) * t25 + t22;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t53;
t10 = t23 * t50 - t24 * t25;
t12 = t23 * t25 + t24 * t50;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t51;
t36 = -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-pkin(5) * t53 + qJ(6) * t51);
t18 = t42 * t30;
t17 = t23 * t29 + t24 * t47;
t15 = -t23 * t47 + t24 * t29;
t13 = t43 * t33 + t55;
t6 = t38 * t26;
t5 = t38 * t25;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, t41, g(1) * t34 + g(2) * t31, t41 * pkin(1), 0, 0, 0, 0, 0, t42 * t33, -t18, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, 0, 0, 0, 0, 0, t4, -t44, t4, t18, t44, -g(1) * (-t31 * pkin(1) - t10 * pkin(5) - t9 * qJ(6)) - g(2) * (t34 * pkin(1) + t12 * pkin(5) + t11 * qJ(6)) + (-g(1) * t45 - g(2) * t40) * t24 + (g(1) * t40 - g(2) * t45) * t23; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t13, 0, 0, 0, 0, 0, t38 * t32, -t38 * t29, 0, 0, 0, 0, 0, t6, -t5, t6, -t13, t5, -g(3) * (t39 * t33 - t48) + t43 * (t39 * t30 + t33 * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(1) * t17 - g(2) * t15 + t32 * t55, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t60 * pkin(4) + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
