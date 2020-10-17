% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:55
% EndTime: 2020-01-03 11:22:56
% DurationCPUTime: 0.21s
% Computational Cost: add. (103->47), mult. (256->81), div. (0->0), fcn. (306->10), ass. (0->39)
t25 = sin(pkin(7));
t27 = cos(pkin(7));
t51 = pkin(2) * t27 + qJ(3) * t25;
t29 = sin(qJ(1));
t38 = cos(pkin(8));
t37 = t29 * t38;
t24 = sin(pkin(8));
t31 = cos(qJ(1));
t41 = t31 * t24;
t10 = t27 * t41 - t37;
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t36 = t31 * t38;
t42 = t29 * t24;
t11 = t27 * t36 + t42;
t23 = sin(pkin(9));
t26 = cos(pkin(9));
t43 = t25 * t31;
t4 = t11 * t26 + t23 * t43;
t50 = -t10 * t30 + t4 * t28;
t49 = t10 * t28 + t4 * t30;
t45 = t24 * t25;
t44 = t25 * t29;
t40 = t31 * pkin(1) + t29 * qJ(2);
t35 = t29 * pkin(1) - t31 * qJ(2);
t34 = t51 * t31 + t40;
t8 = t27 * t42 + t36;
t33 = g(2) * t10 + g(3) * t8;
t14 = g(2) * t31 + g(3) * t29;
t13 = g(2) * t29 - g(3) * t31;
t32 = t51 * t29 + t35;
t12 = t14 * t25;
t9 = t27 * t37 - t41;
t7 = t25 * t38 * t26 - t27 * t23;
t6 = g(1) * t27 - t13 * t25;
t3 = t23 * t44 + t9 * t26;
t2 = t8 * t28 + t3 * t30;
t1 = -t3 * t28 + t8 * t30;
t5 = [0, -t14, t13, -t14 * t27, t12, -t13, -g(2) * t40 - g(3) * t35, -g(2) * t11 - g(3) * t9, t33, -t12, -g(2) * t34 - g(3) * t32, -g(2) * t4 - g(3) * t3, -g(2) * (-t11 * t23 + t26 * t43) - g(3) * (-t9 * t23 + t26 * t44), -t33, -g(2) * (t11 * pkin(3) + t10 * qJ(4) + t34) - g(3) * (t9 * pkin(3) + t8 * qJ(4) + t32), 0, 0, 0, 0, 0, -g(2) * t49 - g(3) * t2, g(2) * t50 - g(3) * t1; 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, t14, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t8 + g(3) * t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t7 * t28 + t30 * t45) - g(2) * t1 - g(3) * t50, -g(1) * (-t28 * t45 - t7 * t30) + g(2) * t2 - g(3) * t49;];
taug_reg = t5;
