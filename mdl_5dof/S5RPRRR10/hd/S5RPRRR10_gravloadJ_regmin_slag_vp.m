% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (141->35), mult. (164->58), div. (0->0), fcn. (182->10), ass. (0->33)
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t12 = g(1) * t24 + g(2) * t22;
t17 = pkin(9) + qJ(3);
t13 = sin(t17);
t14 = cos(t17);
t26 = -g(3) * t14 + t12 * t13;
t36 = g(3) * t13;
t18 = qJ(4) + qJ(5);
t15 = sin(t18);
t34 = t22 * t15;
t16 = cos(t18);
t33 = t22 * t16;
t21 = sin(qJ(4));
t32 = t22 * t21;
t23 = cos(qJ(4));
t31 = t22 * t23;
t30 = t24 * t15;
t29 = t24 * t16;
t28 = t24 * t21;
t27 = t24 * t23;
t11 = g(1) * t22 - g(2) * t24;
t10 = t14 * t27 + t32;
t9 = -t14 * t28 + t31;
t8 = -t14 * t31 + t28;
t7 = t14 * t32 + t27;
t6 = t14 * t29 + t34;
t5 = -t14 * t30 + t33;
t4 = -t14 * t33 + t30;
t3 = t14 * t34 + t29;
t2 = g(1) * t6 - g(2) * t4 + t16 * t36;
t1 = -g(1) * t5 + g(2) * t3 + t15 * t36;
t19 = [0, t11, t12, t11 * cos(pkin(9)), -t11 * sin(pkin(9)), -t12, -g(1) * (-t22 * pkin(1) + t24 * qJ(2)) - g(2) * (t24 * pkin(1) + t22 * qJ(2)), 0, 0, 0, 0, 0, t11 * t14, -t11 * t13, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t12 * t14 + t36, 0, 0, 0, 0, 0, t26 * t23, -t26 * t21, 0, 0, 0, 0, 0, t26 * t16, -t26 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t21 * t36, g(1) * t10 - g(2) * t8 + t23 * t36, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t19;
