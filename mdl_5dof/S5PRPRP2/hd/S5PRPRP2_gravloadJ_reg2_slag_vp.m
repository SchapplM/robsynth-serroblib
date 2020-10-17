% Calculate inertial parameters regressor of gravitation load for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:04
% EndTime: 2019-12-05 15:31:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (152->41), mult. (160->53), div. (0->0), fcn. (167->6), ass. (0->28)
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t24 = cos(qJ(4));
t37 = -t20 * (-qJ(5) - pkin(6)) + (t24 * pkin(4) + pkin(3)) * t21;
t23 = sin(qJ(4));
t33 = g(3) * t20;
t19 = pkin(7) + qJ(2);
t17 = sin(t19);
t18 = cos(t19);
t29 = t21 * t23;
t5 = t17 * t29 + t18 * t24;
t7 = t17 * t24 - t18 * t29;
t1 = -g(1) * t7 + g(2) * t5 + t23 * t33;
t34 = g(1) * t17;
t31 = t18 * t23;
t28 = t21 * t24;
t27 = t18 * pkin(2) + t17 * qJ(3);
t25 = pkin(3) * t21 + pkin(6) * t20;
t11 = g(1) * t18 + g(2) * t17;
t10 = -g(2) * t18 + t34;
t13 = t18 * qJ(3);
t9 = t10 * t20;
t8 = t17 * t23 + t18 * t28;
t6 = -t17 * t28 + t31;
t4 = -g(1) * t6 - g(2) * t8;
t3 = -g(1) * t5 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t6 + t24 * t33;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t21, -t9, -t11, -g(1) * (-t17 * pkin(2) + t13) - g(2) * t27, 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * t13 - g(2) * (t25 * t18 + t27) - (-pkin(2) - t25) * t34, 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * (pkin(4) * t31 + t13) - g(2) * (t37 * t18 + t27) + (-g(1) * (-pkin(2) - t37) - g(2) * pkin(4) * t23) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t21 - t11 * t20;];
taug_reg = t12;
