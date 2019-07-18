% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t15 = g(1) * t27 + g(2) * t24;
t20 = qJ(2) + qJ(4);
t17 = sin(t20);
t18 = cos(t20);
t3 = -g(3) * t18 + t15 * t17;
t40 = g(3) * t17;
t4 = t15 * t18 + t40;
t26 = cos(qJ(2));
t43 = pkin(1) * t26;
t38 = t17 * pkin(4);
t22 = sin(qJ(5));
t37 = t24 * t22;
t25 = cos(qJ(5));
t36 = t24 * t25;
t28 = pkin(2) + pkin(1);
t35 = t26 * t28;
t21 = -pkin(3) - qJ(3);
t34 = t27 * t21;
t33 = t27 * t22;
t32 = t27 * t25;
t31 = -t24 * t21 + t27 * t35;
t14 = g(1) * t24 - g(2) * t27;
t30 = t35 + t38;
t23 = sin(qJ(2));
t5 = -g(3) * t26 + t15 * t23;
t13 = t14 * t26;
t12 = t14 * t23;
t11 = t18 * t32 + t37;
t10 = -t18 * t33 + t36;
t9 = -t18 * t36 + t33;
t8 = t18 * t37 + t32;
t7 = t14 * t17;
t6 = g(3) * t23 + t15 * t26;
t2 = t3 * t25;
t1 = t3 * t22;
t16 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t15, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t15, -g(1) * (t27 * qJ(3) - t24 * t43) - g(2) * (t24 * qJ(3) + t27 * t43), 0, 0, 0, 0, 0, 0, t14 * t18, -t7, -t15, -g(1) * (-t24 * t35 - t34) - g(2) * t31, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(1) * (-t30 * t24 - t34) - g(2) * (t27 * t38 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t5 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * t28, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t30 - t15 * (pkin(4) * t18 - t23 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t22 * t40, g(1) * t11 - g(2) * t9 + t25 * t40, 0, 0;];
taug_reg  = t16;
