% Calculate inertial parameters regressor of gravitation load for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:41:16
% EndTime: 2019-05-05 13:41:17
% DurationCPUTime: 0.26s
% Computational Cost: add. (208->61), mult. (334->75), div. (0->0), fcn. (399->10), ass. (0->36)
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t49 = sin(qJ(1));
t50 = cos(qJ(1));
t11 = -t49 * t44 - t50 * t45;
t12 = t50 * t44 - t49 * t45;
t37 = g(1) * t11 + g(2) * t12;
t26 = pkin(10) + qJ(5);
t19 = sin(t26);
t20 = cos(t26);
t55 = -g(3) * t20 + t37 * t19;
t52 = g(3) * t19;
t30 = sin(qJ(6));
t48 = t20 * t30;
t31 = cos(qJ(6));
t47 = t20 * t31;
t46 = t50 * pkin(1) + t49 * qJ(2);
t43 = t50 * pkin(2) + t46;
t28 = cos(pkin(10));
t18 = t28 * pkin(4) + pkin(3);
t29 = -pkin(7) - qJ(4);
t42 = -t11 * t18 - t12 * t29 + t43;
t41 = -t49 * pkin(1) + t50 * qJ(2);
t40 = t20 * pkin(5) + t19 * pkin(8);
t38 = g(1) * t12 - g(2) * t11;
t36 = t11 * t30 + t12 * t47;
t35 = -t11 * t31 + t12 * t48;
t34 = -t49 * pkin(2) + t41;
t32 = -t11 * t29 + t12 * t18 + t34;
t14 = g(1) * t50 + g(2) * t49;
t13 = g(1) * t49 - g(2) * t50;
t4 = -t11 * t47 + t12 * t30;
t3 = t11 * t48 + t12 * t31;
t2 = t38 * t19;
t1 = -t37 * t20 - t52;
t5 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t14, -g(1) * t41 - g(2) * t46, 0, 0, 0, 0, 0, 0, -t38, t37, 0, -g(1) * t34 - g(2) * t43, 0, 0, 0, 0, 0, 0, -t38 * t28, t38 * sin(pkin(10)) -t37, -g(1) * (t12 * pkin(3) + t11 * qJ(4) + t34) - g(2) * (-t11 * pkin(3) + t12 * qJ(4) + t43) 0, 0, 0, 0, 0, 0, -t38 * t20, t2, -t37, -g(1) * t32 - g(2) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t4, g(1) * t35 - g(2) * t3, -t2, -g(1) * (t40 * t12 + t32) - g(2) * (-t40 * t11 + t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t31, t55 * t30, -t1, g(3) * t40 - t37 * (pkin(5) * t19 - pkin(8) * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t35 - t30 * t52, g(1) * t4 - g(2) * t36 - t31 * t52, 0, 0;];
taug_reg  = t5;
