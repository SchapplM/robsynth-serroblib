% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:48:52
% EndTime: 2019-05-05 17:48:53
% DurationCPUTime: 0.26s
% Computational Cost: add. (328->73), mult. (340->100), div. (0->0), fcn. (355->10), ass. (0->45)
t33 = sin(qJ(1));
t34 = cos(qJ(1));
t14 = g(1) * t34 + g(2) * t33;
t26 = pkin(9) + qJ(3);
t21 = sin(t26);
t23 = cos(t26);
t5 = -g(3) * t23 + t14 * t21;
t51 = g(3) * t21;
t25 = pkin(10) + qJ(5);
t20 = sin(t25);
t49 = t33 * t20;
t22 = cos(t25);
t48 = t33 * t22;
t27 = sin(pkin(10));
t47 = t33 * t27;
t29 = cos(pkin(10));
t46 = t33 * t29;
t45 = t34 * t20;
t44 = t34 * t22;
t43 = t34 * t27;
t42 = t34 * t29;
t32 = -pkin(7) - qJ(2);
t41 = pkin(4) * t27 - t32;
t7 = t23 * t49 + t44;
t9 = t23 * t45 - t48;
t40 = g(1) * t7 - g(2) * t9;
t13 = g(1) * t33 - g(2) * t34;
t39 = t23 * pkin(3) + t21 * qJ(4);
t18 = t29 * pkin(4) + pkin(3);
t31 = -pkin(8) - qJ(4);
t37 = t23 * t18 - t21 * t31;
t36 = pkin(5) * t22 + qJ(6) * t20 + t18;
t1 = g(1) * t9 + g(2) * t7 + t20 * t51;
t10 = t23 * t44 + t49;
t8 = t23 * t48 - t45;
t35 = g(1) * t10 + g(2) * t8 + t22 * t51;
t30 = cos(pkin(9));
t19 = t30 * pkin(2) + pkin(1);
t12 = t34 * t19;
t11 = t13 * t21;
t6 = t14 * t23 + t51;
t4 = t5 * t22;
t3 = t5 * t20;
t2 = g(1) * t8 - g(2) * t10;
t15 = [0, t13, t14, t13 * t30, -t13 * sin(pkin(9)) -t14, -g(1) * (-t33 * pkin(1) + t34 * qJ(2)) - g(2) * (t34 * pkin(1) + t33 * qJ(2)) 0, 0, 0, 0, 0, t13 * t23, -t11, -g(1) * (-t23 * t46 + t43) - g(2) * (t23 * t42 + t47) -g(1) * (t23 * t47 + t42) - g(2) * (-t23 * t43 + t46) t11, -g(2) * t12 + (g(1) * t32 - g(2) * t39) * t34 + (-g(1) * (-t19 - t39) + g(2) * t32) * t33, 0, 0, 0, 0, 0, t2, -t40, t2, t11, t40, -g(1) * (-t8 * pkin(5) - t7 * qJ(6)) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t12) + (-g(1) * t41 - g(2) * t37) * t34 + (-g(1) * (-t19 - t37) - g(2) * t41) * t33; 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t5 * t29, -t5 * t27, -t6, -g(3) * t39 + t14 * (pkin(3) * t21 - qJ(4) * t23) 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3 (-g(3) * t36 + t14 * t31) * t23 + (g(3) * t31 + t14 * t36) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t35, t1, 0, -t35, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t20 + qJ(6) * t22) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t15;
