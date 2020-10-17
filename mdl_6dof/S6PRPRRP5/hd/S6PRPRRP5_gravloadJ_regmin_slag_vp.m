% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:57:01
% EndTime: 2019-05-04 23:57:02
% DurationCPUTime: 0.26s
% Computational Cost: add. (200->72), mult. (504->116), div. (0->0), fcn. (607->10), ass. (0->47)
t23 = sin(pkin(10));
t25 = cos(pkin(10));
t32 = cos(qJ(2));
t29 = sin(qJ(2));
t42 = cos(pkin(6));
t40 = t29 * t42;
t12 = t23 * t32 + t25 * t40;
t14 = -t23 * t40 + t25 * t32;
t56 = -g(1) * t14 - g(2) * t12;
t24 = sin(pkin(6));
t53 = g(3) * t24;
t39 = t32 * t42;
t11 = t23 * t29 - t25 * t39;
t27 = sin(qJ(5));
t52 = t11 * t27;
t13 = t23 * t39 + t25 * t29;
t51 = t13 * t27;
t28 = sin(qJ(4));
t50 = t24 * t28;
t49 = t24 * t29;
t31 = cos(qJ(4));
t48 = t24 * t31;
t47 = t24 * t32;
t46 = t27 * t28;
t45 = t27 * t29;
t30 = cos(qJ(5));
t44 = t28 * t30;
t43 = t29 * t30;
t41 = g(3) * (pkin(2) * t47 + qJ(3) * t49);
t22 = t30 * pkin(5) + pkin(4);
t26 = -qJ(6) - pkin(9);
t38 = t22 * t28 + t26 * t31;
t15 = t42 * t28 + t31 * t47;
t3 = -t13 * t31 + t23 * t50;
t5 = t11 * t31 + t25 * t50;
t36 = g(1) * t3 - g(2) * t5 + g(3) * t15;
t16 = -t28 * t47 + t42 * t31;
t4 = t13 * t28 + t23 * t48;
t6 = -t11 * t28 + t25 * t48;
t35 = g(1) * t4 - g(2) * t6 + g(3) * t16;
t2 = -g(1) * t13 - g(2) * t11 + g(3) * t47;
t34 = g(3) * t49 - t56;
t33 = -g(1) * (t14 * t30 - t4 * t27) - g(2) * (t12 * t30 + t6 * t27) - g(3) * (-t16 * t27 + t24 * t43);
t10 = t13 * pkin(2);
t9 = t11 * pkin(2);
t1 = t34 * t31;
t7 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t2, t34, t2, -t34, -g(1) * (t14 * qJ(3) - t10) - g(2) * (t12 * qJ(3) - t9) - t41, 0, 0, 0, 0, 0, -t34 * t28, -t1, 0, 0, 0, 0, 0, -g(1) * (t14 * t44 - t51) - g(2) * (t12 * t44 - t52) - (t27 * t32 + t28 * t43) * t53, -g(1) * (-t13 * t30 - t14 * t46) - g(2) * (-t11 * t30 - t12 * t46) - (-t28 * t45 + t30 * t32) * t53, t1, -g(1) * (-pkin(5) * t51 - t13 * pkin(8) - t10) - g(2) * (-pkin(5) * t52 - t11 * pkin(8) - t9) - t41 - ((pkin(5) * t27 + pkin(8)) * t32 + t38 * t29) * t53 + t56 * (qJ(3) + t38); 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35, 0, 0, 0, 0, 0, t36 * t30, -t36 * t27, -t35, -g(1) * (-t3 * t22 - t4 * t26) - g(2) * (t5 * t22 + t6 * t26) - g(3) * (-t15 * t22 - t16 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -g(1) * (-t14 * t27 - t4 * t30) - g(2) * (-t12 * t27 + t6 * t30) - g(3) * (-t16 * t30 - t24 * t45) 0, t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36;];
taug_reg  = t7;
