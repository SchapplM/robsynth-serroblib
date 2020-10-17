% Calculate minimal parameter regressor of gravitation load for
% S5PRPPR1
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
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:09
% EndTime: 2019-12-05 15:22:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (115->36), mult. (99->55), div. (0->0), fcn. (107->8), ass. (0->27)
t17 = pkin(7) + qJ(2);
t13 = sin(t17);
t29 = g(1) * t13;
t19 = sin(pkin(8));
t28 = g(3) * t19;
t15 = cos(t17);
t27 = t15 * pkin(2) + t13 * qJ(3);
t21 = cos(pkin(8));
t26 = t13 * t21;
t25 = t15 * t21;
t18 = sin(pkin(9));
t24 = t18 * t21;
t20 = cos(pkin(9));
t23 = t20 * t21;
t7 = g(1) * t15 + g(2) * t13;
t6 = -g(2) * t15 + t29;
t22 = pkin(3) * t21 + qJ(4) * t19;
t16 = pkin(9) + qJ(5);
t14 = cos(t16);
t12 = sin(t16);
t9 = t15 * qJ(3);
t5 = t6 * t19;
t4 = t13 * t12 + t14 * t25;
t3 = -t12 * t25 + t13 * t14;
t2 = t15 * t12 - t14 * t26;
t1 = t12 * t26 + t15 * t14;
t8 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t6, t7, t6 * t21, -t5, -t7, -g(1) * (-t13 * pkin(2) + t9) - g(2) * t27, -g(1) * (-t13 * t23 + t15 * t18) - g(2) * (t13 * t18 + t15 * t23), -g(1) * (t13 * t24 + t15 * t20) - g(2) * (t13 * t20 - t15 * t24), t5, -g(1) * t9 - g(2) * (t22 * t15 + t27) - (-pkin(2) - t22) * t29, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t21 - t7 * t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t12 * t28, g(1) * t4 - g(2) * t2 + t14 * t28;];
taug_reg = t8;
