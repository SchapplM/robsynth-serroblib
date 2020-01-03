% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t37 = cos(pkin(9));
t21 = t37 * pkin(3) + pkin(2);
t34 = pkin(9) + qJ(4);
t26 = cos(t34);
t13 = pkin(4) * t26 + t21;
t35 = qJ(1) + qJ(2);
t28 = sin(t35);
t29 = cos(t35);
t38 = -pkin(7) - qJ(3);
t33 = -pkin(8) + t38;
t48 = t28 * t13 + t29 * t33;
t47 = t28 * t21 + t29 * t38;
t46 = t29 * pkin(2) + t28 * qJ(3);
t45 = t29 * t13 - t28 * t33;
t44 = t29 * t21 - t28 * t38;
t43 = t28 * pkin(2) - t29 * qJ(3);
t12 = g(2) * t29 + g(3) * t28;
t11 = g(2) * t28 - g(3) * t29;
t39 = sin(qJ(1));
t40 = cos(qJ(1));
t42 = -g(2) * t40 - g(3) * t39;
t25 = sin(t34);
t41 = -g(1) * t26 + t11 * t25;
t32 = t40 * pkin(1);
t31 = t39 * pkin(1);
t27 = qJ(5) + t34;
t20 = cos(t27);
t19 = sin(t27);
t8 = t12 * t37;
t7 = t12 * sin(pkin(9));
t6 = t12 * t26;
t5 = t12 * t25;
t4 = t12 * t20;
t3 = t12 * t19;
t2 = -g(1) * t20 + t11 * t19;
t1 = g(1) * t19 + t11 * t20;
t9 = [0, 0, 0, 0, 0, 0, t42, g(2) * t39 - g(3) * t40, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, t42 * pkin(1), 0, 0, 0, 0, 0, 0, -t8, t7, -t11, -g(2) * (t32 + t46) - g(3) * (t31 + t43), 0, 0, 0, 0, 0, 0, -t6, t5, -t11, -g(2) * (t32 + t44) - g(3) * (t31 + t47), 0, 0, 0, 0, 0, 0, -t4, t3, -t11, -g(2) * (t32 + t45) - g(3) * (t31 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -t11, -g(2) * t46 - g(3) * t43, 0, 0, 0, 0, 0, 0, -t6, t5, -t11, -g(2) * t44 - g(3) * t47, 0, 0, 0, 0, 0, 0, -t4, t3, -t11, -g(2) * t45 - g(3) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(1) * t25 + t11 * t26, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t41 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t9;
