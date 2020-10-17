% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:44:48
% EndTime: 2019-05-04 23:44:49
% DurationCPUTime: 0.24s
% Computational Cost: add. (292->72), mult. (513->116), div. (0->0), fcn. (616->12), ass. (0->45)
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t42 = sin(pkin(10));
t44 = cos(pkin(6));
t37 = t44 * t42;
t43 = cos(pkin(10));
t11 = t43 * t29 + t31 * t37;
t52 = g(1) * t11;
t12 = -t29 * t37 + t43 * t31;
t51 = g(1) * t12;
t24 = sin(pkin(6));
t50 = g(3) * t24;
t22 = pkin(11) + qJ(4);
t21 = cos(t22);
t28 = sin(qJ(5));
t49 = t21 * t28;
t30 = cos(qJ(5));
t48 = t21 * t30;
t47 = t24 * t29;
t46 = t28 * t31;
t45 = t30 * t31;
t41 = pkin(5) * t28 + pkin(8) + qJ(3);
t40 = t24 * t43;
t39 = t24 * t42;
t38 = t44 * t43;
t19 = t30 * pkin(5) + pkin(4);
t20 = sin(t22);
t25 = cos(pkin(11));
t26 = -qJ(6) - pkin(9);
t36 = t25 * pkin(3) + t19 * t21 - t20 * t26 + pkin(2);
t10 = t29 * t38 + t42 * t31;
t3 = t10 * t20 + t21 * t40;
t5 = t12 * t20 - t21 * t39;
t7 = t20 * t47 - t44 * t21;
t35 = g(1) * t5 + g(2) * t3 + g(3) * t7;
t4 = t10 * t21 - t20 * t40;
t6 = t12 * t21 + t20 * t39;
t8 = t44 * t20 + t21 * t47;
t34 = g(1) * t6 + g(2) * t4 + g(3) * t8;
t9 = t42 * t29 - t31 * t38;
t2 = -g(2) * t9 + t31 * t50 - t52;
t33 = g(2) * t10 + g(3) * t47 + t51;
t32 = -g(1) * (t11 * t30 - t6 * t28) - g(2) * (-t4 * t28 + t9 * t30) - g(3) * (-t24 * t45 - t8 * t28);
t1 = t2 * t20;
t13 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t2, t33, -t2 * t25, t2 * sin(pkin(11)) -t33, -g(1) * (-t11 * pkin(2) + t12 * qJ(3)) - g(2) * (-t9 * pkin(2) + t10 * qJ(3)) - (pkin(2) * t31 + qJ(3) * t29) * t50, 0, 0, 0, 0, 0, -t2 * t21, t1, 0, 0, 0, 0, 0, -g(1) * (-t11 * t48 + t12 * t28) - g(2) * (t10 * t28 - t9 * t48) - (t21 * t45 + t28 * t29) * t50, -g(1) * (t11 * t49 + t12 * t30) - g(2) * (t10 * t30 + t9 * t49) - (-t21 * t46 + t29 * t30) * t50, -t1, -g(2) * (t41 * t10 - t36 * t9) - t41 * t51 + t36 * t52 - (t41 * t29 + t36 * t31) * t50; 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, 0, 0, 0, 0, 0, t35 * t30, -t35 * t28, -t34, -g(1) * (-t5 * t19 - t6 * t26) - g(2) * (-t3 * t19 - t4 * t26) - g(3) * (-t7 * t19 - t8 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -g(1) * (-t11 * t28 - t6 * t30) - g(2) * (-t9 * t28 - t4 * t30) - g(3) * (t24 * t46 - t8 * t30) 0, t32 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35;];
taug_reg  = t13;
