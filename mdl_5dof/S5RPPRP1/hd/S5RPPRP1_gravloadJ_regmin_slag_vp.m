% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t11 = qJ(1) + pkin(7);
t10 = cos(t11);
t17 = cos(qJ(4));
t13 = cos(pkin(8));
t15 = sin(qJ(4));
t26 = t13 * t15;
t9 = sin(t11);
t1 = t10 * t17 + t9 * t26;
t12 = sin(pkin(8));
t28 = g(1) * t12;
t3 = t10 * t26 - t9 * t17;
t32 = -g(2) * t1 + g(3) * t3 + t15 * t28;
t29 = pkin(4) * t15;
t18 = cos(qJ(1));
t27 = t18 * pkin(1);
t25 = t13 * t17;
t16 = sin(qJ(1));
t23 = -t16 * pkin(1) + t10 * qJ(3);
t22 = g(2) * t9 - g(3) * t10;
t21 = g(2) * t10 + g(3) * t9;
t20 = g(2) * t18 + g(3) * t16;
t19 = t12 * (-qJ(5) - pkin(6)) - t13 * (t17 * pkin(4) + pkin(3)) - pkin(2);
t5 = t21 * t12;
t4 = -t10 * t25 - t9 * t15;
t2 = -t10 * t15 + t9 * t25;
t6 = [0, t20, -g(2) * t16 + g(3) * t18, t20 * pkin(1), t21 * t13, -t5, t22, -g(2) * (-t10 * pkin(2) - t9 * qJ(3) - t27) - g(3) * (-t9 * pkin(2) + t23), 0, 0, 0, 0, 0, -g(2) * t4 + g(3) * t2, -g(2) * t3 - g(3) * t1, t5, g(2) * t27 - g(3) * t23 + (-g(2) * t19 - g(3) * t29) * t10 + (-g(2) * (-qJ(3) - t29) - g(3) * t19) * t9; 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -g(2) * t2 - g(3) * t4 + t17 * t28, 0, t32 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t13 + t22 * t12;];
taug_reg = t6;
