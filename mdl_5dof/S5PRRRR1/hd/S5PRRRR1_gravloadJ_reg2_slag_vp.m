% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_gravloadJ_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = sin(qJ(2));
t22 = cos(qJ(2));
t12 = g(1) * t22 + g(3) * t19;
t16 = qJ(3) + qJ(4);
t14 = sin(t16);
t15 = cos(t16);
t3 = g(2) * t15 + t12 * t14;
t31 = g(2) * t14;
t17 = sin(qJ(5));
t29 = t19 * t17;
t20 = cos(qJ(5));
t28 = t19 * t20;
t27 = t22 * t17;
t26 = t22 * t20;
t25 = g(1) * t19 - g(3) * t22;
t21 = cos(qJ(3));
t24 = t25 * t21;
t18 = sin(qJ(3));
t23 = g(2) * t21 + t12 * t18;
t11 = pkin(2) * t24;
t10 = t15 * t26 + t29;
t9 = -t15 * t27 + t28;
t8 = -t15 * t28 + t27;
t7 = t15 * t29 + t26;
t6 = t25 * t14;
t5 = t23 * pkin(2);
t4 = t12 * t15 - t31;
t2 = t3 * t20;
t1 = t3 * t17;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t12, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25 * t18, -t12, 0, 0, 0, 0, 0, 0, 0, t25 * t15, -t6, -t12, t11, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(3) * t10, -g(1) * t7 - g(3) * t9, t6, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -g(2) * t18 + t12 * t21, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(3) * t7 - t17 * t31, g(1) * t10 - g(3) * t8 - t20 * t31, 0, 0;];
taug_reg  = t13;
