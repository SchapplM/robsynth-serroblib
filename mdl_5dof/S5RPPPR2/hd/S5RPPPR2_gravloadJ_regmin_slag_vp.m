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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:59:46
% EndTime: 2022-01-23 08:59:47
% DurationCPUTime: 0.27s
% Computational Cost: add. (103->61), mult. (240->106), div. (0->0), fcn. (286->10), ass. (0->43)
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t29 = cos(pkin(7));
t32 = cos(qJ(5));
t25 = sin(pkin(8));
t30 = sin(qJ(5));
t44 = t25 * t30;
t24 = sin(pkin(9));
t26 = sin(pkin(7));
t27 = cos(pkin(9));
t28 = cos(pkin(8));
t40 = t27 * t28;
t7 = t26 * t24 + t29 * t40;
t34 = t29 * t44 + t7 * t32;
t43 = t25 * t32;
t9 = t27 * t43 - t30 * t28;
t49 = t34 * t31 - t33 * t9;
t1 = t29 * t43 - t7 * t30;
t8 = t27 * t44 + t32 * t28;
t48 = -t1 * t33 + t31 * t8;
t45 = t26 * qJ(3) + pkin(1);
t42 = t26 * t31;
t41 = t26 * t33;
t39 = t31 * t25;
t38 = t31 * t28;
t37 = t33 * t25;
t36 = t33 * t28;
t10 = t29 * t39 + t36;
t12 = t29 * t37 - t38;
t35 = g(1) * t10 - g(2) * t12;
t19 = g(1) * t33 + g(2) * t31;
t18 = g(1) * t31 - g(2) * t33;
t22 = t33 * qJ(2);
t21 = t31 * qJ(2);
t16 = pkin(2) * t29 + t45;
t15 = -t25 * pkin(3) + qJ(4) * t28 - qJ(2);
t14 = t18 * t26;
t13 = t29 * t36 + t39;
t11 = -t29 * t38 + t37;
t6 = -t29 * t24 + t26 * t40;
t4 = g(3) * t29 - t19 * t26;
t3 = (t28 * pkin(3) + t25 * qJ(4) + pkin(2)) * t29 + t45;
t2 = [0, t18, t19, t18 * t29, -t14, -t19, -g(1) * (-t31 * pkin(1) + t22) - g(2) * (t33 * pkin(1) + t21), -g(1) * t11 - g(2) * t13, -t35, t14, -g(1) * (-t16 * t31 + t22) - g(2) * (t16 * t33 + t21), -g(1) * (t11 * t27 - t24 * t42) - g(2) * (t13 * t27 + t24 * t41), -g(1) * (-t11 * t24 - t27 * t42) - g(2) * (-t13 * t24 + t27 * t41), t35, -g(1) * (-t15 * t33 - t3 * t31) - g(2) * (-t15 * t31 + t3 * t33), 0, 0, 0, 0, 0, g(1) * t49 - g(2) * ((t27 * t39 + t7 * t33) * t32 + t12 * t30), -g(1) * (-t1 * t31 - t33 * t8) + g(2) * t48; 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, -t18, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 * t25 - g(1) * t12 - g(2) * t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * (-(-t27 * t37 + t7 * t31) * t30 + t10 * t32) - g(3) * (t26 * t43 - t30 * t6), -g(1) * (-t31 * t9 - t33 * t34) + g(2) * t49 - g(3) * (-t26 * t44 - t6 * t32);];
taug_reg = t2;
