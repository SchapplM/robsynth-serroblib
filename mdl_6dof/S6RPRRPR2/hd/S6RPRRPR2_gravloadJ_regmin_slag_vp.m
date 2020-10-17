% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:06:28
% EndTime: 2019-05-05 22:06:28
% DurationCPUTime: 0.20s
% Computational Cost: add. (238->46), mult. (217->74), div. (0->0), fcn. (227->10), ass. (0->41)
t20 = qJ(1) + pkin(10);
t17 = sin(t20);
t18 = cos(t20);
t35 = g(1) * t18 + g(2) * t17;
t25 = cos(qJ(4));
t22 = sin(qJ(4));
t26 = cos(qJ(3));
t39 = t22 * t26;
t11 = t17 * t25 - t18 * t39;
t23 = sin(qJ(3));
t43 = g(3) * t23;
t9 = t17 * t39 + t18 * t25;
t48 = -g(1) * t11 + g(2) * t9 + t22 * t43;
t7 = -g(3) * t26 + t35 * t23;
t41 = t17 * t26;
t40 = t18 * t26;
t38 = t25 * t26;
t36 = pkin(4) * t22 + pkin(7);
t34 = g(1) * t17 - g(2) * t18;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t33 = g(1) * t24 - g(2) * t27;
t16 = t25 * pkin(4) + pkin(3);
t21 = -qJ(5) - pkin(8);
t32 = t26 * t16 - t23 * t21;
t30 = pkin(2) + t32;
t29 = t33 * pkin(1);
t19 = qJ(4) + pkin(11) + qJ(6);
t15 = cos(t19);
t14 = sin(t19);
t13 = t34 * t23;
t12 = t17 * t22 + t18 * t38;
t10 = -t17 * t38 + t18 * t22;
t8 = t35 * t26 + t43;
t6 = t17 * t14 + t15 * t40;
t5 = -t14 * t40 + t17 * t15;
t4 = t18 * t14 - t15 * t41;
t3 = t14 * t41 + t18 * t15;
t2 = g(1) * t6 - g(2) * t4 + t15 * t43;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t43;
t28 = [0, t33, g(1) * t27 + g(2) * t24, t29, 0, 0, 0, 0, 0, t34 * t26, -t13, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, t29 + (-g(1) * t36 - g(2) * t30) * t18 + (g(1) * t30 - g(2) * t36) * t17, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t7 * t25, -t7 * t22, -t8, -g(3) * t32 + t35 * (t16 * t23 + t21 * t26) 0, 0, 0, 0, 0, t7 * t15, -t7 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(1) * t12 - g(2) * t10 + t25 * t43, 0, t48 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t28;
