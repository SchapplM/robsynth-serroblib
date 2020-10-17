% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:55:25
% EndTime: 2019-05-05 13:55:26
% DurationCPUTime: 0.31s
% Computational Cost: add. (326->76), mult. (251->95), div. (0->0), fcn. (247->12), ass. (0->44)
t23 = pkin(10) + qJ(4);
t16 = sin(t23);
t19 = cos(t23);
t37 = t19 * pkin(4) + t16 * qJ(5);
t24 = qJ(1) + pkin(9);
t17 = sin(t24);
t20 = cos(t24);
t9 = g(1) * t20 + g(2) * t17;
t1 = -g(3) * t19 + t9 * t16;
t53 = g(3) * t16;
t31 = sin(qJ(1));
t50 = t31 * pkin(1);
t49 = t17 * t19;
t25 = sin(pkin(11));
t48 = t17 * t25;
t27 = cos(pkin(11));
t47 = t17 * t27;
t22 = pkin(11) + qJ(6);
t15 = sin(t22);
t46 = t20 * t15;
t18 = cos(t22);
t45 = t20 * t18;
t44 = t20 * t25;
t43 = t20 * t27;
t28 = cos(pkin(10));
t14 = t28 * pkin(3) + pkin(2);
t32 = cos(qJ(1));
t21 = t32 * pkin(1);
t42 = t20 * t14 + t21;
t30 = -pkin(7) - qJ(3);
t40 = pkin(5) * t25 - t30;
t8 = g(1) * t17 - g(2) * t20;
t39 = g(1) * t31 - g(2) * t32;
t38 = -t20 * t30 - t50;
t13 = t27 * pkin(5) + pkin(4);
t29 = -pkin(8) - qJ(5);
t35 = t19 * t13 - t16 * t29;
t7 = t8 * t16;
t6 = t17 * t15 + t19 * t45;
t5 = t17 * t18 - t19 * t46;
t4 = -t18 * t49 + t46;
t3 = t15 * t49 + t45;
t2 = t9 * t19 + t53;
t10 = [0, 0, 0, 0, 0, 0, t39, g(1) * t32 + g(2) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t39 * pkin(1), 0, 0, 0, 0, 0, 0, t8 * t28, -t8 * sin(pkin(10)) -t9, -g(1) * (-t17 * pkin(2) + t20 * qJ(3) - t50) - g(2) * (t20 * pkin(2) + t17 * qJ(3) + t21) 0, 0, 0, 0, 0, 0, t8 * t19, -t7, -t9, -g(1) * (-t17 * t14 + t38) - g(2) * (-t17 * t30 + t42) 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t47 + t44) - g(2) * (t19 * t43 + t48) -g(1) * (t19 * t48 + t43) - g(2) * (-t19 * t44 + t47) t7, -g(1) * t38 - g(2) * (t37 * t20 + t42) + (-g(1) * (-t14 - t37) + g(2) * t30) * t17, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t7, g(1) * t50 - g(2) * t42 + (-g(1) * t40 - g(2) * t35) * t20 + (-g(1) * (-t14 - t35) - g(2) * t40) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t27, -t1 * t25, -t2, -g(3) * t37 + t9 * (pkin(4) * t16 - qJ(5) * t19) 0, 0, 0, 0, 0, 0, t1 * t18, -t1 * t15, -t2, -g(3) * t35 + t9 * (t13 * t16 + t19 * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t15 * t53, g(1) * t6 - g(2) * t4 + t18 * t53, 0, 0;];
taug_reg  = t10;
