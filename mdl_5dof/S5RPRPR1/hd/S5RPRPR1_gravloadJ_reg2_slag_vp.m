% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = sin(qJ(3));
t25 = t19 * pkin(3);
t18 = -qJ(4) - pkin(6);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t24 = t22 * pkin(1) + t20 * qJ(2);
t17 = qJ(3) + pkin(8);
t6 = g(1) * t22 + g(2) * t20;
t5 = g(1) * t20 - g(2) * t22;
t21 = cos(qJ(3));
t23 = g(3) * t19 - t5 * t21;
t16 = -pkin(7) + t18;
t13 = t22 * qJ(2);
t11 = qJ(5) + t17;
t10 = cos(t17);
t9 = sin(t17);
t8 = cos(t11);
t7 = sin(t11);
t3 = pkin(4) * t9 + t25;
t2 = g(3) * t7 - t5 * t8;
t1 = g(3) * t8 + t5 * t7;
t4 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -g(1) * (-t20 * pkin(1) + t13) - g(2) * t24, 0, 0, 0, 0, 0, 0, -t6 * t19, -t6 * t21, t5, -g(1) * (t13 + (-pkin(1) - pkin(6)) * t20) - g(2) * (t22 * pkin(6) + t24), 0, 0, 0, 0, 0, 0, -t6 * t9, -t6 * t10, t5, -g(1) * (t22 * t25 + t13 + (-pkin(1) + t18) * t20) - g(2) * (-t22 * t18 + t20 * t25 + t24), 0, 0, 0, 0, 0, 0, -t6 * t7, -t6 * t8, t5, -g(1) * (t22 * t3 + t13 + (-pkin(1) + t16) * t20) - g(2) * (-t22 * t16 + t20 * t3 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t21 + t5 * t19, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t9 - t5 * t10, g(3) * t10 + t5 * t9, 0, t23 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, g(3) * t3 - t5 * (t21 * pkin(3) + pkin(4) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t4;
