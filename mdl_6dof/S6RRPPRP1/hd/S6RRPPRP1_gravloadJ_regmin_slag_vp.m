% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:04:52
% EndTime: 2019-05-06 09:04:53
% DurationCPUTime: 0.30s
% Computational Cost: add. (328->76), mult. (352->103), div. (0->0), fcn. (364->10), ass. (0->50)
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t13 = g(1) * t34 + g(2) * t32;
t26 = qJ(2) + pkin(9);
t20 = sin(t26);
t22 = cos(t26);
t37 = -g(3) * t22 + t13 * t20;
t31 = sin(qJ(2));
t58 = pkin(2) * t31;
t55 = g(3) * t20;
t30 = -pkin(8) - qJ(4);
t53 = t20 * t30;
t25 = pkin(10) + qJ(5);
t19 = sin(t25);
t52 = t32 * t19;
t21 = cos(t25);
t51 = t32 * t21;
t27 = sin(pkin(10));
t50 = t32 * t27;
t28 = cos(pkin(10));
t49 = t32 * t28;
t48 = t34 * t19;
t47 = t34 * t21;
t46 = t34 * t27;
t45 = t34 * t28;
t29 = -qJ(3) - pkin(7);
t44 = pkin(4) * t27 - t29;
t7 = t22 * t52 + t47;
t9 = t22 * t48 - t51;
t43 = g(1) * t7 - g(2) * t9;
t12 = g(1) * t32 - g(2) * t34;
t42 = t22 * pkin(3) + t20 * qJ(4);
t17 = t28 * pkin(4) + pkin(3);
t41 = t22 * t17 - t53;
t40 = pkin(5) * t21 + qJ(6) * t19 + t17;
t1 = g(1) * t9 + g(2) * t7 + t19 * t55;
t10 = t22 * t47 + t52;
t8 = t22 * t51 - t48;
t38 = g(1) * t10 + g(2) * t8 + t21 * t55;
t33 = cos(qJ(2));
t36 = -g(3) * t33 + t13 * t31;
t23 = t33 * pkin(2);
t18 = t23 + pkin(1);
t14 = t34 * t18;
t11 = t12 * t20;
t6 = -t13 * t22 - t55;
t4 = t37 * t21;
t3 = t37 * t19;
t2 = g(1) * t8 - g(2) * t10;
t5 = [0, t12, t13, 0, 0, 0, 0, 0, t12 * t33, -t12 * t31, -t13, -g(1) * (-t32 * t18 - t34 * t29) - g(2) * (-t32 * t29 + t14) -g(1) * (-t22 * t49 + t46) - g(2) * (t22 * t45 + t50) -g(1) * (t22 * t50 + t45) - g(2) * (-t22 * t46 + t49) t11, -g(2) * t14 + (g(1) * t29 - g(2) * t42) * t34 + (-g(1) * (-t18 - t42) + g(2) * t29) * t32, 0, 0, 0, 0, 0, t2, -t43, t2, t11, t43, -g(1) * (-t8 * pkin(5) - t7 * qJ(6)) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t14) + (-g(1) * t44 - g(2) * t41) * t34 + (-g(1) * (-t18 - t41) - g(2) * t44) * t32; 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t31 + t13 * t33, 0, t36 * pkin(2), t37 * t28, -t37 * t27, t6, -g(3) * (t23 + t42) + t13 * (pkin(3) * t20 - qJ(4) * t22 + t58) 0, 0, 0, 0, 0, t4, -t3, t4, t6, t3, -g(3) * (t40 * t22 + t23 - t53) + t13 * (t40 * t20 + t22 * t30 + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t38, t1, 0, -t38, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t19 + qJ(6) * t21) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
