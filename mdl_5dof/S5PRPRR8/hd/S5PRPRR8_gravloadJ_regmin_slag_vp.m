% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(pkin(5));
t36 = g(3) * t16;
t20 = sin(qJ(4));
t35 = t16 * t20;
t23 = cos(qJ(4));
t34 = t16 * t23;
t24 = cos(qJ(2));
t33 = t16 * t24;
t18 = cos(pkin(5));
t21 = sin(qJ(2));
t32 = t18 * t21;
t31 = t18 * t24;
t19 = sin(qJ(5));
t30 = t19 * t20;
t29 = t19 * t21;
t22 = cos(qJ(5));
t28 = t20 * t22;
t27 = t21 * t22;
t15 = sin(pkin(9));
t17 = cos(pkin(9));
t7 = t15 * t21 - t17 * t31;
t9 = t15 * t31 + t17 * t21;
t26 = g(1) * (-t15 * t35 + t9 * t23) + g(2) * (t17 * t35 + t7 * t23) + g(3) * (-t18 * t20 - t23 * t33);
t1 = -g(1) * t9 - g(2) * t7 + g(3) * t33;
t10 = -t15 * t32 + t17 * t24;
t8 = t15 * t24 + t17 * t32;
t25 = g(1) * t10 + g(2) * t8 + t21 * t36;
t12 = t18 * t23 - t20 * t33;
t5 = t17 * t34 - t7 * t20;
t3 = t15 * t34 + t9 * t20;
t2 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t1, t25, t1, -t25, -g(1) * (-t9 * pkin(2) + t10 * qJ(3)) - g(2) * (-t7 * pkin(2) + t8 * qJ(3)) - (pkin(2) * t24 + qJ(3) * t21) * t36, 0, 0, 0, 0, 0, -t25 * t20, -t25 * t23, 0, 0, 0, 0, 0, -g(1) * (t10 * t28 - t9 * t19) - g(2) * (-t7 * t19 + t8 * t28) - (t19 * t24 + t20 * t27) * t36, -g(1) * (-t10 * t30 - t9 * t22) - g(2) * (-t7 * t22 - t8 * t30) - (-t20 * t29 + t22 * t24) * t36; 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, g(1) * t3 - g(2) * t5 + g(3) * t12, 0, 0, 0, 0, 0, -t26 * t22, t26 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t22 - t3 * t19) - g(2) * (t5 * t19 + t8 * t22) - g(3) * (-t12 * t19 + t16 * t27), -g(1) * (-t10 * t19 - t3 * t22) - g(2) * (-t8 * t19 + t5 * t22) - g(3) * (-t12 * t22 - t16 * t29);];
taug_reg = t2;
