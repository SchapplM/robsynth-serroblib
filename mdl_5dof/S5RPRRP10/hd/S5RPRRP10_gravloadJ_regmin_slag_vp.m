% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP10
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
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:48
% EndTime: 2021-01-15 19:14:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->38), mult. (209->53), div. (0->0), fcn. (215->7), ass. (0->34)
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t15 = g(1) * t28 + g(2) * t26;
t21 = pkin(8) + qJ(3);
t19 = cos(t21);
t27 = cos(qJ(4));
t34 = t28 * t27;
t25 = sin(qJ(4));
t37 = t26 * t25;
t10 = t19 * t37 + t34;
t35 = t28 * t25;
t36 = t26 * t27;
t12 = -t19 * t35 + t36;
t18 = sin(t21);
t39 = g(3) * t18;
t1 = -g(1) * t12 + g(2) * t10 + t25 * t39;
t7 = -g(3) * t19 + t15 * t18;
t32 = pkin(4) * t25 + pkin(6) + qJ(2);
t14 = g(1) * t26 - g(2) * t28;
t17 = t27 * pkin(4) + pkin(3);
t23 = -qJ(5) - pkin(7);
t31 = t19 * t17 - t18 * t23;
t22 = cos(pkin(8));
t29 = t22 * pkin(2) + pkin(1) + t31;
t13 = t19 * t34 + t37;
t11 = -t19 * t36 + t35;
t9 = t14 * t18;
t8 = t15 * t19 + t39;
t6 = t7 * t27;
t5 = t7 * t25;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t27 * t39;
t16 = [0, t14, t15, t14 * t22, -t15, -g(1) * (-t26 * pkin(1) + t28 * qJ(2)) - g(2) * (t28 * pkin(1) + t26 * qJ(2)), 0, 0, 0, 0, 0, t14 * t19, -t9, 0, 0, 0, 0, 0, t4, t3, t4, t3, t9, (-g(1) * t32 - g(2) * t29) * t28 + (g(1) * t29 - g(2) * t32) * t26; 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t8, -g(3) * t31 + t15 * (t17 * t18 + t19 * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t16;
