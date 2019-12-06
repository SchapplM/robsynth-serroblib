% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(2));
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t37 = -g(1) * t14 - g(2) * t13;
t38 = t37 * t16;
t18 = cos(qJ(2));
t4 = g(3) * t16 - t18 * t37;
t32 = g(3) * t18;
t31 = t18 * pkin(2) + t16 * qJ(3);
t30 = t13 * t16;
t29 = t14 * t16;
t15 = sin(qJ(4));
t28 = t15 * t16;
t27 = t15 * t18;
t17 = cos(qJ(4));
t26 = t16 * t17;
t25 = qJ(3) * t18;
t5 = t13 * t25;
t6 = t14 * t25;
t24 = -g(1) * t6 - g(2) * t5;
t20 = -g(1) * (-t13 * t15 + t14 * t26) - g(2) * (t13 * t26 + t14 * t15) + t17 * t32;
t19 = -pkin(7) - pkin(6);
t12 = qJ(4) + qJ(5);
t8 = cos(t12);
t7 = sin(t12);
t3 = -t32 - t38;
t2 = -g(1) * (-t13 * t8 - t7 * t29) - g(2) * (t14 * t8 - t7 * t30) - t7 * t32;
t1 = -g(1) * (-t13 * t7 + t8 * t29) - g(2) * (t14 * t7 + t8 * t30) + t8 * t32;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (-pkin(2) * t29 + t6) - g(2) * (-pkin(2) * t30 + t5) - g(3) * t31, 0, 0, 0, 0, 0, 0, -t4 * t15, -t4 * t17, t3, -g(3) * (t18 * pkin(6) + t31) + t24 - (pkin(2) + pkin(6)) * t38, 0, 0, 0, 0, 0, 0, -t4 * t7, -t4 * t8, t3, -g(3) * (pkin(4) * t28 - t18 * t19 + t31) + t24 + t37 * (pkin(4) * t27 + (-pkin(2) + t19) * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -g(1) * (-t13 * t17 - t14 * t28) - g(2) * (-t13 * t28 + t14 * t17) - g(3) * t27, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
