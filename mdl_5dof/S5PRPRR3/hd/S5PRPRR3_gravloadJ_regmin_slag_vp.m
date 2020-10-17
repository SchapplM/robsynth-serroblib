% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:32
% EndTime: 2019-12-05 15:47:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (97->27), mult. (112->44), div. (0->0), fcn. (121->10), ass. (0->27)
t10 = cos(pkin(8));
t9 = sin(pkin(8));
t18 = g(1) * t10 + g(2) * t9;
t7 = qJ(2) + pkin(9);
t3 = sin(t7);
t4 = cos(t7);
t17 = -g(3) * t4 + t18 * t3;
t28 = g(3) * t3;
t8 = qJ(4) + qJ(5);
t5 = sin(t8);
t26 = t9 * t5;
t6 = cos(t8);
t25 = t9 * t6;
t24 = t10 * t5;
t23 = t10 * t6;
t11 = sin(qJ(4));
t22 = t9 * t11;
t13 = cos(qJ(4));
t21 = t9 * t13;
t20 = t10 * t11;
t19 = t10 * t13;
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t15 = -g(3) * t14 + t18 * t12;
t2 = -g(1) * (-t4 * t23 - t26) - g(2) * (-t4 * t25 + t24) + t6 * t28;
t1 = -g(1) * (-t4 * t24 + t25) - g(2) * (-t4 * t26 - t23) + t5 * t28;
t16 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t15, g(3) * t12 + t18 * t14, t15 * pkin(2), 0, 0, 0, 0, 0, t17 * t13, -t17 * t11, 0, 0, 0, 0, 0, t17 * t6, -t17 * t5; 0, 0, 0, 0, -g(1) * t9 + g(2) * t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t20 + t21) - g(2) * (-t4 * t22 - t19) + t11 * t28, -g(1) * (-t4 * t19 - t22) - g(2) * (-t4 * t21 + t20) + t13 * t28, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t16;
