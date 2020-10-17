% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:04:54
% EndTime: 2019-05-05 11:04:54
% DurationCPUTime: 0.30s
% Computational Cost: add. (421->87), mult. (736->163), div. (0->0), fcn. (946->14), ass. (0->44)
t23 = sin(pkin(6));
t47 = g(3) * t23;
t21 = qJ(4) + qJ(5);
t20 = qJ(6) + t21;
t16 = sin(t20);
t28 = cos(qJ(3));
t46 = t16 * t28;
t17 = cos(t20);
t45 = t17 * t28;
t18 = sin(t21);
t44 = t18 * t28;
t19 = cos(t21);
t43 = t19 * t28;
t26 = sin(qJ(2));
t42 = t23 * t26;
t41 = t23 * t28;
t29 = cos(qJ(2));
t40 = t23 * t29;
t24 = sin(qJ(4));
t39 = t24 * t28;
t27 = cos(qJ(4));
t38 = t27 * t28;
t37 = t28 * t29;
t36 = cos(pkin(6));
t35 = cos(pkin(12));
t22 = sin(pkin(12));
t34 = t22 * t36;
t33 = t23 * t35;
t32 = t36 * t35;
t10 = t22 * t29 + t26 * t32;
t12 = -t26 * t34 + t29 * t35;
t25 = sin(qJ(3));
t31 = g(1) * (-t12 * t25 + t22 * t41) + g(2) * (-t10 * t25 - t28 * t33) + g(3) * (-t25 * t42 + t28 * t36);
t11 = t26 * t35 + t29 * t34;
t9 = t22 * t26 - t29 * t32;
t30 = -g(1) * t11 - g(2) * t9 + g(3) * t40;
t14 = t25 * t36 + t26 * t41;
t8 = t22 * t23 * t25 + t12 * t28;
t6 = t10 * t28 - t25 * t33;
t4 = -g(1) * (-t11 * t18 - t19 * t8) - g(2) * (-t18 * t9 - t19 * t6) - g(3) * (-t14 * t19 + t18 * t40);
t3 = -g(1) * (t11 * t19 - t18 * t8) - g(2) * (-t18 * t6 + t19 * t9) - g(3) * (-t14 * t18 - t19 * t40);
t2 = -g(1) * (-t11 * t16 - t17 * t8) - g(2) * (-t16 * t9 - t17 * t6) - g(3) * (-t14 * t17 + t16 * t40);
t1 = -g(1) * (t11 * t17 - t16 * t8) - g(2) * (-t16 * t6 + t17 * t9) - g(3) * (-t14 * t16 - t17 * t40);
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t30, g(1) * t12 + g(2) * t10 + g(3) * t42, 0, 0, 0, 0, 0, -t30 * t28, t30 * t25, 0, 0, 0, 0, 0, -g(1) * (-t11 * t38 + t12 * t24) - g(2) * (t10 * t24 - t38 * t9) - (t24 * t26 + t27 * t37) * t47, -g(1) * (t11 * t39 + t12 * t27) - g(2) * (t10 * t27 + t39 * t9) - (-t24 * t37 + t26 * t27) * t47, 0, 0, 0, 0, 0, -g(1) * (-t11 * t43 + t12 * t18) - g(2) * (t10 * t18 - t43 * t9) - (t18 * t26 + t19 * t37) * t47, -g(1) * (t11 * t44 + t12 * t19) - g(2) * (t10 * t19 + t44 * t9) - (-t18 * t37 + t19 * t26) * t47, 0, 0, 0, 0, 0, -g(1) * (-t11 * t45 + t12 * t16) - g(2) * (t10 * t16 - t45 * t9) - (t16 * t26 + t17 * t37) * t47, -g(1) * (t11 * t46 + t12 * t17) - g(2) * (t10 * t17 + t46 * t9) - (-t16 * t37 + t17 * t26) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, g(1) * t8 + g(2) * t6 + g(3) * t14, 0, 0, 0, 0, 0, -t31 * t27, t31 * t24, 0, 0, 0, 0, 0, -t31 * t19, t31 * t18, 0, 0, 0, 0, 0, -t31 * t17, t31 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t27 - t24 * t8) - g(2) * (-t24 * t6 + t27 * t9) - g(3) * (-t14 * t24 - t27 * t40) -g(1) * (-t11 * t24 - t27 * t8) - g(2) * (-t24 * t9 - t27 * t6) - g(3) * (-t14 * t27 + t24 * t40) 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t5;
