% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:17:21
% EndTime: 2019-05-05 04:17:22
% DurationCPUTime: 0.25s
% Computational Cost: add. (312->66), mult. (448->115), div. (0->0), fcn. (548->12), ass. (0->42)
t21 = sin(pkin(11));
t23 = cos(pkin(11));
t30 = cos(qJ(2));
t27 = sin(qJ(2));
t39 = cos(pkin(6));
t38 = t27 * t39;
t12 = t21 * t30 + t23 * t38;
t14 = -t21 * t38 + t23 * t30;
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t22 = sin(pkin(6));
t42 = t22 * t29;
t43 = t22 * t27;
t50 = -g(1) * (-t14 * t26 + t21 * t42) - g(2) * (-t12 * t26 - t23 * t42) - g(3) * (-t26 * t43 + t39 * t29);
t49 = g(3) * t22;
t20 = qJ(3) + pkin(12) + qJ(5);
t18 = cos(t20);
t25 = sin(qJ(6));
t48 = t18 * t25;
t28 = cos(qJ(6));
t47 = t18 * t28;
t46 = t21 * t22;
t45 = t22 * t23;
t44 = t22 * t26;
t41 = t25 * t30;
t40 = t28 * t30;
t37 = t30 * t39;
t17 = sin(t20);
t36 = g(1) * (-t14 * t17 + t18 * t46) + g(2) * (-t12 * t17 - t18 * t45) + g(3) * (-t17 * t43 + t39 * t18);
t11 = t21 * t27 - t23 * t37;
t13 = t21 * t37 + t23 * t27;
t33 = -g(1) * t13 - g(2) * t11 + t30 * t49;
t32 = g(1) * t14 + g(2) * t12 + g(3) * t43;
t24 = -qJ(4) - pkin(8);
t19 = t29 * pkin(3) + pkin(2);
t10 = t39 * t17 + t18 * t43;
t8 = t14 * t18 + t17 * t46;
t6 = t12 * t18 - t17 * t45;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t10;
t2 = t36 * t28;
t1 = t36 * t25;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t33, t32, 0, 0, 0, 0, 0, -t33 * t29, t33 * t26, -t32, -g(1) * (-t13 * t19 - t14 * t24) - g(2) * (-t11 * t19 - t12 * t24) - (t19 * t30 - t24 * t27) * t49, 0, 0, 0, 0, 0, -t33 * t18, t33 * t17, 0, 0, 0, 0, 0, -g(1) * (-t13 * t47 + t14 * t25) - g(2) * (-t11 * t47 + t12 * t25) - (t18 * t40 + t25 * t27) * t49, -g(1) * (t13 * t48 + t14 * t28) - g(2) * (t11 * t48 + t12 * t28) - (-t18 * t41 + t27 * t28) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -g(1) * (-t14 * t29 - t21 * t44) - g(2) * (-t12 * t29 + t23 * t44) - g(3) * (-t39 * t26 - t27 * t42) 0, t50 * pkin(3), 0, 0, 0, 0, 0, -t36, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t28 - t8 * t25) - g(2) * (t11 * t28 - t6 * t25) - g(3) * (-t10 * t25 - t22 * t40) -g(1) * (-t13 * t25 - t8 * t28) - g(2) * (-t11 * t25 - t6 * t28) - g(3) * (-t10 * t28 + t22 * t41);];
taug_reg  = t3;
