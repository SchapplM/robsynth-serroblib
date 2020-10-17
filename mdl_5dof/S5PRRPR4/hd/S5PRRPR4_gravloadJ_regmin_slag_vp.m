% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:22
% EndTime: 2019-12-05 16:23:23
% DurationCPUTime: 0.13s
% Computational Cost: add. (109->31), mult. (143->50), div. (0->0), fcn. (149->8), ass. (0->23)
t10 = cos(pkin(8));
t9 = sin(pkin(8));
t19 = g(1) * t10 + g(2) * t9;
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t3 = -g(3) * t15 + t19 * t13;
t25 = g(3) * t13;
t23 = t15 * t9;
t22 = t10 * t15;
t12 = sin(qJ(3));
t21 = t12 * t15;
t14 = cos(qJ(3));
t20 = t14 * t15;
t16 = -g(1) * (-t10 * t21 + t9 * t14) - g(2) * (-t10 * t14 - t9 * t21) + t12 * t25;
t11 = -qJ(4) - pkin(6);
t8 = qJ(3) + pkin(9) + qJ(5);
t7 = t14 * pkin(3) + pkin(2);
t6 = cos(t8);
t5 = sin(t8);
t4 = t19 * t15 + t25;
t2 = -g(1) * (-t6 * t22 - t9 * t5) - g(2) * (t10 * t5 - t6 * t23) + t6 * t25;
t1 = -g(1) * (-t5 * t22 + t9 * t6) - g(2) * (-t10 * t6 - t5 * t23) + t5 * t25;
t17 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t14, -t3 * t12, -t4, -g(3) * (-t13 * t11 + t15 * t7) + t19 * (t11 * t15 + t13 * t7), 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -g(1) * (-t10 * t20 - t9 * t12) - g(2) * (t10 * t12 - t9 * t20) + t14 * t25, 0, t16 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t17;
