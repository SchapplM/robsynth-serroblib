% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t28 = pkin(2) * t14;
t15 = cos(t17);
t27 = pkin(2) * t15;
t20 = sin(qJ(1));
t26 = t20 * pkin(1);
t21 = cos(qJ(1));
t25 = t21 * pkin(1);
t13 = pkin(8) + t17;
t10 = cos(t13);
t9 = sin(t13);
t5 = g(2) * t9 - g(3) * t10;
t24 = g(2) * t10 + g(3) * t9;
t7 = g(2) * t15 + g(3) * t14;
t23 = -t9 * pkin(3) + t10 * qJ(4) - t28;
t22 = -t10 * pkin(3) - t9 * qJ(4) - t27;
t16 = pkin(9) + qJ(5);
t12 = cos(t16);
t11 = sin(t16);
t6 = -g(2) * t14 + g(3) * t15;
t4 = t24 * cos(pkin(9));
t3 = t24 * sin(pkin(9));
t2 = t24 * t12;
t1 = t24 * t11;
t8 = [0, g(2) * t21 + g(3) * t20, -g(2) * t20 + g(3) * t21, 0, t7, t6, -g(2) * (-t25 - t27) - g(3) * (-t26 - t28), t4, -t3, t5, -g(2) * (t22 - t25) - g(3) * (t23 - t26), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t7, t6, t7 * pkin(2), t4, -t3, t5, -g(2) * t22 - g(3) * t23, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 - t5 * t11, g(1) * t11 - t5 * t12;];
taug_reg = t8;
