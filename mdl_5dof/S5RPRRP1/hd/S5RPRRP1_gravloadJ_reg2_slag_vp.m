% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = -pkin(7) - pkin(6);
t17 = sin(qJ(3));
t24 = t17 * pkin(3);
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t23 = t20 * pkin(1) + t18 * qJ(2);
t8 = g(1) * t20 + g(2) * t18;
t7 = g(1) * t18 - g(2) * t20;
t16 = qJ(3) + qJ(4);
t10 = cos(t16);
t9 = sin(t16);
t2 = g(3) * t9 - t7 * t10;
t19 = cos(qJ(3));
t22 = g(3) * t17 - t7 * t19;
t15 = -qJ(5) + t21;
t12 = t20 * qJ(2);
t5 = pkin(4) * t9 + t24;
t4 = t8 * t10;
t3 = t8 * t9;
t1 = g(3) * t10 + t7 * t9;
t6 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t18 * pkin(1) + t12) - g(2) * t23, 0, 0, 0, 0, 0, 0, -t8 * t17, -t8 * t19, t7, -g(1) * (t12 + (-pkin(1) - pkin(6)) * t18) - g(2) * (t20 * pkin(6) + t23), 0, 0, 0, 0, 0, 0, -t3, -t4, t7, -g(1) * (t20 * t24 + t12 + (-pkin(1) + t21) * t18) - g(2) * (t18 * t24 - t20 * t21 + t23), 0, 0, 0, 0, 0, 0, -t3, -t4, t7, -g(1) * (t20 * t5 + t12 + (-pkin(1) + t15) * t18) - g(2) * (-t20 * t15 + t18 * t5 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t19 + t7 * t17, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t22 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, g(3) * t5 - t7 * (t19 * pkin(3) + pkin(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t6;
