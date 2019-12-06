% Calculate inertial parameters regressor of gravitation load for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = sin(pkin(7));
t21 = cos(pkin(7));
t30 = g(1) * t21 + g(2) * t19;
t23 = sin(qJ(2));
t24 = cos(qJ(2));
t8 = -g(3) * t24 + t30 * t23;
t36 = g(3) * t23;
t34 = t19 * t24;
t33 = t21 * t24;
t22 = -pkin(6) - qJ(3);
t32 = t22 * t24;
t20 = cos(pkin(8));
t13 = t20 * pkin(3) + pkin(2);
t31 = t24 * t13 - t23 * t22;
t17 = pkin(8) + qJ(4);
t14 = sin(t17);
t15 = cos(t17);
t28 = pkin(4) * t15 + qJ(5) * t14;
t4 = t14 * t34 + t21 * t15;
t6 = t14 * t33 - t19 * t15;
t1 = g(1) * t6 + g(2) * t4 + t14 * t36;
t5 = -t21 * t14 + t15 * t34;
t7 = t19 * t14 + t15 * t33;
t26 = g(1) * t7 + g(2) * t5 + t15 * t36;
t9 = t30 * t24 + t36;
t3 = t8 * t15;
t2 = t8 * t14;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t20, -t8 * sin(pkin(8)), -t9, -g(3) * (t24 * pkin(2) + t23 * qJ(3)) + t30 * (pkin(2) * t23 - qJ(3) * t24), 0, 0, 0, 0, 0, 0, t3, -t2, -t9, -g(3) * t31 + t30 * (t13 * t23 + t32), 0, 0, 0, 0, 0, 0, t3, -t9, t2, -g(3) * (t28 * t24 + t31) + t30 * (t32 - (-t13 - t28) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t26, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t26, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - (-pkin(4) * t14 + qJ(5) * t15) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
