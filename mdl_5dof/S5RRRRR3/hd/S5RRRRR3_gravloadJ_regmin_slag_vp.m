% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:29
% EndTime: 2019-12-05 18:56:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (176->34), mult. (206->56), div. (0->0), fcn. (226->10), ass. (0->40)
t22 = qJ(2) + qJ(3);
t18 = sin(t22);
t20 = cos(t22);
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t31 = g(1) * t28 + g(2) * t25;
t7 = -g(3) * t20 + t31 * t18;
t41 = g(3) * t18;
t21 = qJ(4) + qJ(5);
t17 = sin(t21);
t39 = t25 * t17;
t19 = cos(t21);
t38 = t25 * t19;
t23 = sin(qJ(4));
t37 = t25 * t23;
t26 = cos(qJ(4));
t36 = t25 * t26;
t35 = t28 * t17;
t34 = t28 * t19;
t33 = t28 * t23;
t32 = t28 * t26;
t30 = g(1) * t25 - g(2) * t28;
t27 = cos(qJ(2));
t24 = sin(qJ(2));
t16 = t20 * t32 + t37;
t15 = -t20 * t33 + t36;
t14 = -t20 * t36 + t33;
t13 = t20 * t37 + t32;
t12 = t20 * t34 + t39;
t11 = -t20 * t35 + t38;
t10 = -t20 * t38 + t35;
t9 = t20 * t39 + t34;
t8 = t31 * t20 + t41;
t6 = t7 * t26;
t5 = t7 * t23;
t4 = t7 * t19;
t3 = t7 * t17;
t2 = g(1) * t12 - g(2) * t10 + t19 * t41;
t1 = -g(1) * t11 + g(2) * t9 + t17 * t41;
t29 = [0, t30, t31, 0, 0, 0, 0, 0, t30 * t27, -t30 * t24, 0, 0, 0, 0, 0, t30 * t20, -t30 * t18, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t27 + t31 * t24, g(3) * t24 + t31 * t27, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t23 * t41, g(1) * t16 - g(2) * t14 + t26 * t41, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t29;
