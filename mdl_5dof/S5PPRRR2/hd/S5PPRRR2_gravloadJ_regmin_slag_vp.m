% Calculate minimal parameter regressor of gravitation load for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t10 = cos(pkin(8));
t9 = sin(pkin(8));
t15 = g(1) * t10 + g(2) * t9;
t7 = pkin(9) + qJ(3);
t3 = sin(t7);
t4 = cos(t7);
t14 = -g(3) * t4 + t15 * t3;
t25 = g(3) * t3;
t8 = qJ(4) + qJ(5);
t5 = sin(t8);
t23 = t9 * t5;
t6 = cos(t8);
t22 = t9 * t6;
t21 = t10 * t5;
t20 = t10 * t6;
t11 = sin(qJ(4));
t19 = t9 * t11;
t12 = cos(qJ(4));
t18 = t9 * t12;
t17 = t10 * t11;
t16 = t10 * t12;
t2 = -g(1) * (-t4 * t20 - t23) - g(2) * (-t4 * t22 + t21) + t6 * t25;
t1 = -g(1) * (-t4 * t21 + t22) - g(2) * (-t4 * t23 - t20) + t5 * t25;
t13 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t9 + g(2) * t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t14, t15 * t4 + t25, 0, 0, 0, 0, 0, t14 * t12, -t14 * t11, 0, 0, 0, 0, 0, t14 * t6, -t14 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t17 + t18) - g(2) * (-t4 * t19 - t16) + t11 * t25, -g(1) * (-t4 * t16 - t19) - g(2) * (-t4 * t18 + t17) + t12 * t25, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t13;
