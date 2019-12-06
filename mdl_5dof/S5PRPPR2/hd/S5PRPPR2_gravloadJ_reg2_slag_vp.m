% Calculate inertial parameters regressor of gravitation load for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t14 = sin(pkin(7));
t16 = cos(pkin(7));
t24 = g(1) * t16 + g(2) * t14;
t12 = qJ(2) + pkin(8);
t7 = sin(t12);
t9 = cos(t12);
t1 = -g(3) * t9 + t24 * t7;
t31 = g(3) * t7;
t18 = sin(qJ(2));
t29 = pkin(2) * t18;
t26 = t14 * t9;
t25 = t16 * t9;
t19 = cos(qJ(2));
t20 = -g(3) * t19 + t24 * t18;
t17 = -pkin(6) - qJ(4);
t15 = cos(pkin(9));
t11 = pkin(9) + qJ(5);
t10 = t19 * pkin(2);
t8 = cos(t11);
t6 = sin(t11);
t5 = t15 * pkin(4) + pkin(3);
t3 = -g(1) * t14 + g(2) * t16;
t2 = t24 * t9 + t31;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, g(3) * t18 + t24 * t19, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t20 * pkin(2), 0, 0, 0, 0, 0, 0, t1 * t15, -t1 * sin(pkin(9)), -t2, -g(3) * (t9 * pkin(3) + t7 * qJ(4) + t10) + t24 * (pkin(3) * t7 - qJ(4) * t9 + t29), 0, 0, 0, 0, 0, 0, t1 * t8, -t1 * t6, -t2, -g(3) * (-t7 * t17 + t9 * t5 + t10) + t24 * (t17 * t9 + t5 * t7 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t8 - t6 * t25) - g(2) * (-t16 * t8 - t6 * t26) + t6 * t31, -g(1) * (-t14 * t6 - t8 * t25) - g(2) * (t16 * t6 - t8 * t26) + t8 * t31, 0, 0;];
taug_reg = t4;
