% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = pkin(9) + qJ(3);
t22 = sin(t27);
t24 = cos(t27);
t44 = t24 * pkin(3) + t22 * qJ(4);
t33 = -pkin(7) - qJ(2);
t56 = pkin(4) - t33;
t34 = sin(qJ(1));
t35 = cos(qJ(1));
t14 = g(1) * t35 + g(2) * t34;
t64 = t14 * t22;
t2 = g(3) * t22 + t14 * t24;
t63 = pkin(3) * t22;
t28 = sin(pkin(10));
t62 = pkin(5) * t28;
t58 = g(3) * t24;
t32 = -pkin(8) - qJ(5);
t55 = t24 * t32;
t26 = pkin(10) + qJ(6);
t21 = sin(t26);
t54 = t34 * t21;
t23 = cos(t26);
t53 = t34 * t23;
t52 = t34 * t28;
t30 = cos(pkin(10));
t51 = t34 * t30;
t50 = t35 * t21;
t49 = t35 * t23;
t48 = t35 * t28;
t47 = t35 * t30;
t46 = t35 * t33;
t43 = t30 * pkin(5) + t56;
t42 = qJ(4) * t24;
t41 = t24 * qJ(5);
t40 = t24 * t62;
t31 = cos(pkin(9));
t20 = t31 * pkin(2) + pkin(1);
t12 = t35 * t20;
t39 = g(2) * (t44 * t35 + t12);
t13 = g(1) * t34 - g(2) * t35;
t38 = t22 * t62 - t55;
t37 = -t20 - t44;
t11 = t35 * t42;
t9 = t34 * t42;
t8 = t13 * t24;
t7 = t13 * t22;
t6 = -t22 * t54 + t49;
t5 = t22 * t53 + t50;
t4 = t22 * t50 + t53;
t3 = t22 * t49 - t54;
t1 = -t58 + t64;
t10 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t31, -t13 * sin(pkin(9)) -t14, -g(1) * (-t34 * pkin(1) + t35 * qJ(2)) - g(2) * (t35 * pkin(1) + t34 * qJ(2)) 0, 0, 0, 0, 0, 0, t8, -t7, -t14, -g(1) * (-t34 * t20 - t46) - g(2) * (-t34 * t33 + t12) 0, 0, 0, 0, 0, 0, -t14, -t8, t7, g(1) * t46 - t39 + (-g(1) * t37 + g(2) * t33) * t34, 0, 0, 0, 0, 0, 0, -g(1) * (-t22 * t52 + t47) - g(2) * (t22 * t48 + t51) -g(1) * (-t22 * t51 - t48) - g(2) * (t22 * t47 - t52) t8, -t39 + (-g(1) * t56 - g(2) * t41) * t35 + (-g(1) * (t37 - t41) - g(2) * t56) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t8, -t39 + (-g(1) * t43 - g(2) * t38) * t35 + (-g(1) * (t37 - t38) - g(2) * t43) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t35 * t63 + t11) - g(2) * (-t34 * t63 + t9) - g(3) * t44, 0, 0, 0, 0, 0, 0, -t2 * t28, -t2 * t30, t1, -g(1) * t11 - g(2) * t9 - g(3) * (t41 + t44) + (pkin(3) + qJ(5)) * t64, 0, 0, 0, 0, 0, 0, -t2 * t21, -t2 * t23, t1, -g(1) * (t35 * t40 + t11) - g(2) * (t34 * t40 + t9) - g(3) * (t44 - t55) + (-g(3) * t62 + t14 * (pkin(3) - t32)) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t23 * t58, g(1) * t4 - g(2) * t6 - t21 * t58, 0, 0;];
taug_reg  = t10;
