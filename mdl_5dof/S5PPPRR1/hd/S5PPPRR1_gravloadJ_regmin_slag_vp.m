% Calculate minimal parameter regressor of gravitation load for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% taug_reg [5x13]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t9 = sin(pkin(8));
t20 = g(3) * t9;
t13 = sin(qJ(5));
t19 = t13 * t9;
t14 = cos(qJ(5));
t18 = t14 * t9;
t10 = sin(pkin(7));
t11 = cos(pkin(8));
t17 = t10 * t11;
t12 = cos(pkin(7));
t16 = t11 * t12;
t8 = pkin(9) + qJ(4);
t6 = sin(t8);
t7 = cos(t8);
t15 = -g(1) * (t10 * t7 - t6 * t16) - g(2) * (-t12 * t7 - t6 * t17) + t6 * t20;
t5 = -g(1) * t10 + g(2) * t12;
t4 = t10 * t6 + t7 * t16;
t2 = -t12 * t6 + t7 * t17;
t1 = [-g(3), -g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t5, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(3) * t11 + (-g(1) * t12 - g(2) * t10) * t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t15, g(1) * t4 + g(2) * t2 + t7 * t20, 0, 0, 0, 0, 0, t15 * t14, -t15 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t18 - t4 * t13) - g(2) * (t10 * t18 - t2 * t13) - g(3) * (-t11 * t14 - t7 * t19), -g(1) * (-t12 * t19 - t4 * t14) - g(2) * (-t10 * t19 - t2 * t14) - g(3) * (t11 * t13 - t7 * t18);];
taug_reg = t1;
