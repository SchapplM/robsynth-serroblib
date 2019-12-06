% Calculate inertial parameters regressor of gravitation load for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t15 = pkin(7) + qJ(2);
t10 = sin(t15);
t12 = cos(t15);
t6 = g(1) * t12 + g(2) * t10;
t5 = g(1) * t10 - g(2) * t12;
t16 = qJ(3) + pkin(8);
t11 = sin(t16);
t13 = cos(t16);
t22 = t13 * pkin(4) + t11 * qJ(5);
t18 = sin(qJ(3));
t19 = cos(qJ(3));
t20 = -g(3) * t19 + t6 * t18;
t17 = -qJ(4) - pkin(6);
t14 = t19 * pkin(3);
t9 = t14 + pkin(2);
t7 = t12 * t9;
t4 = t5 * t13;
t3 = t5 * t11;
t2 = g(3) * t11 + t6 * t13;
t1 = -g(3) * t13 + t6 * t11;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t19, -t5 * t18, -t6, -g(1) * (-t10 * pkin(2) + t12 * pkin(6)) - g(2) * (t12 * pkin(2) + t10 * pkin(6)), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t10 * t9 - t12 * t17) - g(2) * (-t10 * t17 + t7), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(2) * t7 + (g(1) * t17 - g(2) * t22) * t12 + (-g(1) * (-t22 - t9) + g(2) * t17) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, g(3) * t18 + t6 * t19, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t20 * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * (t14 + t22) + t6 * (pkin(3) * t18 + pkin(4) * t11 - qJ(5) * t13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t8;
