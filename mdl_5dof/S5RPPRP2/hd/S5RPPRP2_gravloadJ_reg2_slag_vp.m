% Calculate inertial parameters regressor of gravitation load for
% S5RPPRP2
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
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:24
% DurationCPUTime: 0.14s
% Computational Cost: add. (154->42), mult. (122->44), div. (0->0), fcn. (111->8), ass. (0->24)
t15 = pkin(8) + qJ(4);
t10 = sin(t15);
t12 = cos(t15);
t23 = t12 * pkin(4) + t10 * qJ(5);
t16 = qJ(1) + pkin(7);
t11 = sin(t16);
t13 = cos(t16);
t6 = g(1) * t13 + g(2) * t11;
t20 = sin(qJ(1));
t28 = t20 * pkin(1);
t21 = cos(qJ(1));
t14 = t21 * pkin(1);
t18 = cos(pkin(8));
t9 = t18 * pkin(3) + pkin(2);
t27 = t13 * t9 + t14;
t5 = g(1) * t11 - g(2) * t13;
t25 = g(1) * t20 - g(2) * t21;
t19 = -pkin(6) - qJ(3);
t24 = -t13 * t19 - t28;
t4 = t5 * t12;
t3 = t5 * t10;
t2 = g(3) * t10 + t6 * t12;
t1 = -g(3) * t12 + t6 * t10;
t7 = [0, 0, 0, 0, 0, 0, t25, g(1) * t21 + g(2) * t20, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t25 * pkin(1), 0, 0, 0, 0, 0, 0, t5 * t18, -t5 * sin(pkin(8)), -t6, -g(1) * (-t11 * pkin(2) + t13 * qJ(3) - t28) - g(2) * (t13 * pkin(2) + t11 * qJ(3) + t14), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t11 * t9 + t24) - g(2) * (-t11 * t19 + t27), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t24 - g(2) * (t23 * t13 + t27) + (-g(1) * (-t23 - t9) + g(2) * t19) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t23 + t6 * (pkin(4) * t10 - qJ(5) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
