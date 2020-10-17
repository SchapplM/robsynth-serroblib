% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:46:25
% EndTime: 2019-05-05 07:46:26
% DurationCPUTime: 0.29s
% Computational Cost: add. (333->83), mult. (649->144), div. (0->0), fcn. (805->12), ass. (0->49)
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t41 = sin(pkin(11));
t43 = cos(pkin(6));
t36 = t43 * t41;
t42 = cos(pkin(11));
t10 = t42 * t27 + t30 * t36;
t53 = g(1) * t10;
t11 = -t27 * t36 + t42 * t30;
t52 = g(1) * t11;
t23 = sin(pkin(6));
t51 = g(3) * t23;
t22 = qJ(4) + pkin(12) + qJ(6);
t19 = sin(t22);
t29 = cos(qJ(3));
t50 = t19 * t29;
t20 = cos(t22);
t49 = t20 * t29;
t48 = t23 * t27;
t47 = t23 * t30;
t25 = sin(qJ(4));
t46 = t25 * t29;
t28 = cos(qJ(4));
t45 = t28 * t29;
t44 = t29 * t30;
t40 = pkin(4) * t25 + pkin(8);
t39 = t23 * t42;
t38 = t23 * t41;
t37 = t43 * t42;
t21 = t28 * pkin(4) + pkin(3);
t24 = -qJ(5) - pkin(9);
t26 = sin(qJ(3));
t35 = t21 * t29 - t24 * t26 + pkin(2);
t12 = t26 * t48 - t43 * t29;
t9 = t27 * t37 + t41 * t30;
t4 = t9 * t26 + t29 * t39;
t6 = t11 * t26 - t29 * t38;
t34 = g(1) * t6 + g(2) * t4 + g(3) * t12;
t13 = t43 * t26 + t29 * t48;
t5 = -t26 * t39 + t9 * t29;
t7 = t11 * t29 + t26 * t38;
t33 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t8 = t41 * t27 - t30 * t37;
t32 = -g(2) * t8 + g(3) * t47 - t53;
t31 = -g(1) * (t10 * t28 - t7 * t25) - g(2) * (-t5 * t25 + t8 * t28) - g(3) * (-t13 * t25 - t28 * t47);
t3 = t32 * t26;
t2 = -g(1) * (-t10 * t19 - t7 * t20) - g(2) * (-t8 * t19 - t5 * t20) - g(3) * (-t13 * t20 + t19 * t47);
t1 = -g(1) * (t10 * t20 - t7 * t19) - g(2) * (-t5 * t19 + t8 * t20) - g(3) * (-t13 * t19 - t20 * t47);
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t32, g(2) * t9 + g(3) * t48 + t52, 0, 0, 0, 0, 0, -t32 * t29, t3, 0, 0, 0, 0, 0, -g(1) * (-t10 * t45 + t11 * t25) - g(2) * (t9 * t25 - t8 * t45) - (t25 * t27 + t28 * t44) * t51, -g(1) * (t10 * t46 + t11 * t28) - g(2) * (t9 * t28 + t8 * t46) - (-t25 * t44 + t27 * t28) * t51, -t3, -g(2) * (-t35 * t8 + t40 * t9) - t40 * t52 + t35 * t53 - (t40 * t27 + t35 * t30) * t51, 0, 0, 0, 0, 0, -g(1) * (-t10 * t49 + t11 * t19) - g(2) * (t9 * t19 - t8 * t49) - (t19 * t27 + t20 * t44) * t51, -g(1) * (t10 * t50 + t11 * t20) - g(2) * (t9 * t20 + t8 * t50) - (-t19 * t44 + t20 * t27) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, 0, 0, 0, 0, 0, t34 * t28, -t34 * t25, -t33, -g(1) * (-t6 * t21 - t7 * t24) - g(2) * (-t4 * t21 - t5 * t24) - g(3) * (-t12 * t21 - t13 * t24) 0, 0, 0, 0, 0, t34 * t20, -t34 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -g(1) * (-t10 * t25 - t7 * t28) - g(2) * (-t8 * t25 - t5 * t28) - g(3) * (-t13 * t28 + t25 * t47) 0, t31 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t14;
