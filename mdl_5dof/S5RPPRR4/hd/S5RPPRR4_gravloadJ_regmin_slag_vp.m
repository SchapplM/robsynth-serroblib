% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:58
% EndTime: 2022-01-23 09:16:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (134->42), mult. (143->68), div. (0->0), fcn. (163->10), ass. (0->41)
t24 = sin(pkin(8));
t41 = g(3) * t24;
t22 = pkin(9) + qJ(4);
t18 = qJ(5) + t22;
t14 = sin(t18);
t27 = sin(qJ(1));
t40 = t27 * t14;
t15 = cos(t18);
t39 = t27 * t15;
t16 = sin(t22);
t38 = t27 * t16;
t17 = cos(t22);
t37 = t27 * t17;
t23 = sin(pkin(9));
t36 = t27 * t23;
t25 = cos(pkin(9));
t35 = t27 * t25;
t28 = cos(qJ(1));
t34 = t28 * t14;
t33 = t28 * t15;
t32 = t28 * t16;
t31 = t28 * t17;
t30 = t28 * t23;
t29 = t28 * t25;
t13 = g(1) * t28 + g(2) * t27;
t12 = g(1) * t27 - g(2) * t28;
t26 = cos(pkin(8));
t20 = t28 * qJ(2);
t19 = t27 * qJ(2);
t11 = pkin(2) * t26 + t24 * qJ(3) + pkin(1);
t10 = t26 * t31 + t38;
t9 = -t26 * t32 + t37;
t8 = -t26 * t37 + t32;
t7 = t26 * t38 + t31;
t6 = t26 * t33 + t40;
t5 = -t26 * t34 + t39;
t4 = -t26 * t39 + t34;
t3 = t26 * t40 + t33;
t2 = g(1) * t6 - g(2) * t4 + t15 * t41;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t41;
t21 = [0, t12, t13, t12 * t26, -t13, -g(1) * (-t27 * pkin(1) + t20) - g(2) * (t28 * pkin(1) + t19), -g(1) * (-t26 * t35 + t30) - g(2) * (t26 * t29 + t36), -g(1) * (t26 * t36 + t29) - g(2) * (-t26 * t30 + t35), -g(1) * (-t11 * t27 + t20) - g(2) * (t11 * t28 + t19), 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, -t12, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t26 - t13 * t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t16 * t41, g(1) * t10 - g(2) * t8 + t17 * t41, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t21;
