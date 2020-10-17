% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (213->51), mult. (188->63), div. (0->0), fcn. (171->8), ass. (0->29)
t23 = pkin(8) + qJ(3);
t19 = qJ(4) + t23;
t14 = sin(t19);
t15 = cos(t19);
t35 = t15 * pkin(4) + t14 * qJ(5);
t34 = pkin(4) * t14;
t25 = cos(pkin(8));
t16 = t25 * pkin(2) + pkin(1);
t26 = -pkin(6) - qJ(2);
t32 = qJ(5) * t15;
t17 = sin(t23);
t31 = -pkin(3) * t17 - t34;
t27 = sin(qJ(1));
t28 = cos(qJ(1));
t10 = g(1) * t28 + g(2) * t27;
t9 = g(1) * t27 - g(2) * t28;
t18 = cos(t23);
t29 = -g(3) * t18 + t10 * t17;
t22 = -pkin(7) + t26;
t13 = pkin(3) * t18;
t8 = t28 * t32;
t7 = t27 * t32;
t6 = t13 + t16;
t5 = t28 * t6;
t4 = t9 * t15;
t3 = t9 * t14;
t2 = g(3) * t14 + t10 * t15;
t1 = -g(3) * t15 + t10 * t14;
t11 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t25, -t9 * sin(pkin(8)), -t10, -g(1) * (-t27 * pkin(1) + t28 * qJ(2)) - g(2) * (t28 * pkin(1) + t27 * qJ(2)), 0, 0, 0, 0, 0, 0, t9 * t18, -t9 * t17, -t10, -g(1) * (-t27 * t16 - t28 * t26) - g(2) * (t28 * t16 - t27 * t26), 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(1) * (-t28 * t22 - t27 * t6) - g(2) * (-t27 * t22 + t5), 0, 0, 0, 0, 0, 0, t4, -t10, t3, -g(2) * t5 + (g(1) * t22 - g(2) * t35) * t28 + (-g(1) * (-t35 - t6) + g(2) * t22) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(3) * t17 + t10 * t18, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t29 * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t31 * t28 + t8) - g(2) * (t31 * t27 + t7) - g(3) * (t13 + t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t28 * t34 + t8) - g(2) * (-t27 * t34 + t7) - g(3) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
