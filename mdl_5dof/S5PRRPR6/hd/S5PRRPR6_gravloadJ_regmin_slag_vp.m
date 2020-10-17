% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:47
% EndTime: 2019-12-05 16:32:48
% DurationCPUTime: 0.27s
% Computational Cost: add. (202->68), mult. (469->119), div. (0->0), fcn. (579->12), ass. (0->42)
t24 = sin(qJ(2));
t26 = cos(qJ(2));
t37 = cos(pkin(9));
t38 = cos(pkin(5));
t32 = t38 * t37;
t36 = sin(pkin(9));
t6 = t36 * t24 - t26 * t32;
t31 = t38 * t36;
t8 = t37 * t24 + t26 * t31;
t49 = -g(1) * t8 - g(2) * t6;
t21 = sin(pkin(5));
t46 = g(3) * t21;
t19 = pkin(10) + qJ(5);
t17 = sin(t19);
t25 = cos(qJ(3));
t45 = t17 * t25;
t18 = cos(t19);
t44 = t18 * t25;
t20 = sin(pkin(10));
t43 = t20 * t25;
t42 = t21 * t24;
t41 = t21 * t26;
t22 = cos(pkin(10));
t40 = t22 * t25;
t39 = t25 * t26;
t35 = t21 * t37;
t34 = t21 * t36;
t7 = t24 * t32 + t36 * t26;
t9 = -t24 * t31 + t37 * t26;
t33 = -g(1) * t9 - g(2) * t7;
t23 = sin(qJ(3));
t10 = t23 * t42 - t38 * t25;
t2 = t7 * t23 + t25 * t35;
t4 = t9 * t23 - t25 * t34;
t29 = g(1) * t4 + g(2) * t2 + g(3) * t10;
t11 = t38 * t23 + t25 * t42;
t3 = -t23 * t35 + t7 * t25;
t5 = t23 * t34 + t9 * t25;
t28 = g(1) * t5 + g(2) * t3 + g(3) * t11;
t27 = g(3) * t41 + t49;
t1 = t27 * t23;
t12 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t27, g(3) * t42 - t33, 0, 0, 0, 0, 0, -t27 * t25, t1, -g(1) * (t9 * t20 - t8 * t40) - g(2) * (t7 * t20 - t6 * t40) - (t20 * t24 + t22 * t39) * t46, -g(1) * (t9 * t22 + t8 * t43) - g(2) * (t7 * t22 + t6 * t43) - (-t20 * t39 + t22 * t24) * t46, -t1, (-t24 * t46 + t33) * pkin(7) + (-t26 * t46 - t49) * (pkin(3) * t25 + qJ(4) * t23 + pkin(2)), 0, 0, 0, 0, 0, -g(1) * (t9 * t17 - t8 * t44) - g(2) * (t7 * t17 - t6 * t44) - (t17 * t24 + t18 * t39) * t46, -g(1) * (t9 * t18 + t8 * t45) - g(2) * (t7 * t18 + t6 * t45) - (-t17 * t39 + t18 * t24) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t28, t29 * t22, -t29 * t20, -t28, -g(1) * (-t4 * pkin(3) + t5 * qJ(4)) - g(2) * (-t2 * pkin(3) + t3 * qJ(4)) - g(3) * (-t10 * pkin(3) + t11 * qJ(4)), 0, 0, 0, 0, 0, t29 * t18, -t29 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t5 * t17 + t8 * t18) - g(2) * (-t3 * t17 + t6 * t18) - g(3) * (-t11 * t17 - t18 * t41), -g(1) * (-t8 * t17 - t5 * t18) - g(2) * (-t6 * t17 - t3 * t18) - g(3) * (-t11 * t18 + t17 * t41);];
taug_reg = t12;
