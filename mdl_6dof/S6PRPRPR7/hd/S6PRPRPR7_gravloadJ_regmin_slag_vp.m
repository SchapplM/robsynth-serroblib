% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:21:39
% EndTime: 2019-05-04 23:21:40
% DurationCPUTime: 0.25s
% Computational Cost: add. (186->70), mult. (489->109), div. (0->0), fcn. (595->10), ass. (0->41)
t25 = sin(pkin(10));
t27 = cos(pkin(10));
t33 = cos(qJ(2));
t30 = sin(qJ(2));
t42 = cos(pkin(6));
t40 = t30 * t42;
t15 = t25 * t33 + t27 * t40;
t17 = -t25 * t40 + t27 * t33;
t52 = -g(1) * t17 - g(2) * t15;
t26 = sin(pkin(6));
t49 = g(3) * t26;
t29 = sin(qJ(4));
t48 = t26 * t29;
t47 = t26 * t30;
t32 = cos(qJ(4));
t46 = t26 * t32;
t45 = t26 * t33;
t28 = sin(qJ(6));
t44 = t28 * t32;
t31 = cos(qJ(6));
t43 = t31 * t32;
t41 = g(3) * (pkin(2) * t45 + qJ(3) * t47);
t39 = t33 * t42;
t38 = pkin(4) * t29 - qJ(5) * t32;
t18 = t42 * t29 + t32 * t45;
t16 = t25 * t39 + t27 * t30;
t6 = -t16 * t32 + t25 * t48;
t14 = t25 * t30 - t27 * t39;
t8 = t14 * t32 + t27 * t48;
t36 = g(1) * t6 - g(2) * t8 + g(3) * t18;
t19 = -t29 * t45 + t42 * t32;
t7 = t16 * t29 + t25 * t46;
t9 = -t14 * t29 + t27 * t46;
t35 = g(1) * t7 - g(2) * t9 + g(3) * t19;
t4 = -g(1) * t16 - g(2) * t14 + g(3) * t45;
t34 = g(3) * t47 - t52;
t13 = t16 * pkin(2);
t12 = t14 * pkin(2);
t3 = t34 * t32;
t2 = t34 * t29;
t1 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t4, t34, t4, -t34, -g(1) * (t17 * qJ(3) - t13) - g(2) * (t15 * qJ(3) - t12) - t41, 0, 0, 0, 0, 0, -t2, -t3, -t4, t2, t3, -g(1) * (-t16 * pkin(8) - t13) - g(2) * (-t14 * pkin(8) - t12) - t41 - (pkin(8) * t33 + t38 * t30) * t49 + t52 * (qJ(3) + t38) 0, 0, 0, 0, 0, -g(1) * (-t16 * t31 - t17 * t44) - g(2) * (-t14 * t31 - t15 * t44) - (-t30 * t44 + t31 * t33) * t49, -g(1) * (t16 * t28 - t17 * t43) - g(2) * (t14 * t28 - t15 * t43) - (-t28 * t33 - t30 * t43) * t49; 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35, 0, -t36, -t35, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t8 * pkin(4) - t9 * qJ(5)) - g(3) * (-t18 * pkin(4) + t19 * qJ(5)) 0, 0, 0, 0, 0, -t35 * t28, -t35 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t28 + t6 * t31) - g(2) * (-t15 * t28 - t8 * t31) - g(3) * (t18 * t31 - t28 * t47) -g(1) * (-t17 * t31 - t6 * t28) - g(2) * (-t15 * t31 + t8 * t28) - g(3) * (-t18 * t28 - t31 * t47);];
taug_reg  = t1;
