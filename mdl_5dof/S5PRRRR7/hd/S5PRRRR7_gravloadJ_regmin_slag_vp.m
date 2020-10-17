% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:12:59
% DurationCPUTime: 0.14s
% Computational Cost: add. (163->33), mult. (178->58), div. (0->0), fcn. (202->10), ass. (0->24)
t14 = sin(qJ(2));
t16 = cos(qJ(2));
t11 = sin(pkin(9));
t12 = cos(pkin(9));
t19 = g(1) * t12 + g(2) * t11;
t18 = -g(3) * t16 + t19 * t14;
t25 = g(3) * t14;
t23 = t11 * t16;
t22 = t12 * t16;
t13 = sin(qJ(3));
t21 = t13 * t16;
t15 = cos(qJ(3));
t20 = t15 * t16;
t10 = qJ(3) + qJ(4);
t9 = qJ(5) + t10;
t8 = cos(t10);
t7 = sin(t10);
t6 = cos(t9);
t5 = sin(t9);
t4 = -g(1) * (-t11 * t7 - t8 * t22) - g(2) * (t12 * t7 - t8 * t23) + t8 * t25;
t3 = -g(1) * (t11 * t8 - t7 * t22) - g(2) * (-t12 * t8 - t7 * t23) + t7 * t25;
t2 = -g(1) * (-t11 * t5 - t6 * t22) - g(2) * (t12 * t5 - t6 * t23) + t6 * t25;
t1 = -g(1) * (t11 * t6 - t5 * t22) - g(2) * (-t12 * t6 - t5 * t23) + t5 * t25;
t17 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t18, t19 * t16 + t25, 0, 0, 0, 0, 0, t18 * t15, -t18 * t13, 0, 0, 0, 0, 0, t18 * t8, -t18 * t7, 0, 0, 0, 0, 0, t18 * t6, -t18 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t15 - t12 * t21) - g(2) * (-t11 * t21 - t12 * t15) + t13 * t25, -g(1) * (-t11 * t13 - t12 * t20) - g(2) * (-t11 * t20 + t12 * t13) + t15 * t25, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t17;
