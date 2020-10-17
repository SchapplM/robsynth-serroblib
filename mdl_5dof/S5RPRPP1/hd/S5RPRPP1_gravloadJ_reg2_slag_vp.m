% Calculate inertial parameters regressor of gravitation load for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:14
% EndTime: 2019-12-31 18:09:14
% DurationCPUTime: 0.17s
% Computational Cost: add. (168->45), mult. (143->50), div. (0->0), fcn. (129->8), ass. (0->27)
t16 = qJ(3) + pkin(8);
t10 = sin(t16);
t12 = cos(t16);
t25 = t12 * pkin(4) + t10 * qJ(5);
t17 = qJ(1) + pkin(7);
t11 = sin(t17);
t13 = cos(t17);
t6 = g(1) * t13 + g(2) * t11;
t20 = sin(qJ(1));
t30 = t20 * pkin(1);
t22 = cos(qJ(1));
t15 = t22 * pkin(1);
t21 = cos(qJ(3));
t14 = t21 * pkin(3);
t9 = t14 + pkin(2);
t29 = t13 * t9 + t15;
t5 = g(1) * t11 - g(2) * t13;
t27 = g(1) * t20 - g(2) * t22;
t18 = -qJ(4) - pkin(6);
t26 = -t13 * t18 - t30;
t19 = sin(qJ(3));
t23 = -g(3) * t21 + t6 * t19;
t4 = t5 * t12;
t3 = t5 * t10;
t2 = g(3) * t10 + t6 * t12;
t1 = -g(3) * t12 + t6 * t10;
t7 = [0, 0, 0, 0, 0, 0, t27, g(1) * t22 + g(2) * t20, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, t5 * t21, -t5 * t19, -t6, -g(1) * (-t11 * pkin(2) + t13 * pkin(6) - t30) - g(2) * (t13 * pkin(2) + t11 * pkin(6) + t15), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t11 * t9 + t26) - g(2) * (-t11 * t18 + t29), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t26 - g(2) * (t25 * t13 + t29) + (-g(1) * (-t25 - t9) + g(2) * t18) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t19 + t6 * t21, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t23 * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * (t14 + t25) + t6 * (pkin(3) * t19 + pkin(4) * t10 - qJ(5) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
