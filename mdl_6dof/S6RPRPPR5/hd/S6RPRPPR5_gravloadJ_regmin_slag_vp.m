% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t30 = cos(pkin(9));
t26 = pkin(9) + qJ(3);
t21 = sin(t26);
t23 = cos(t26);
t39 = t23 * pkin(3) + t21 * qJ(4);
t57 = t30 * pkin(2) + pkin(1) + t39;
t33 = cos(qJ(1));
t56 = g(2) * t33;
t32 = sin(qJ(1));
t53 = g(1) * t33;
t14 = g(2) * t32 + t53;
t55 = t14 * t21;
t2 = g(3) * t21 + t14 * t23;
t54 = pkin(3) * t21;
t50 = g(3) * t23;
t31 = -pkin(7) - qJ(2);
t49 = pkin(4) - t31;
t25 = pkin(10) + qJ(6);
t20 = sin(t25);
t48 = t32 * t20;
t22 = cos(t25);
t47 = t32 * t22;
t27 = sin(pkin(10));
t46 = t32 * t27;
t29 = cos(pkin(10));
t45 = t32 * t29;
t44 = t33 * t20;
t43 = t33 * t22;
t42 = t33 * t27;
t41 = t33 * t29;
t38 = qJ(4) * t23;
t37 = t23 * qJ(5);
t36 = t57 * t56;
t13 = g(1) * t32 - t56;
t11 = t33 * t38;
t9 = t32 * t38;
t8 = t13 * t23;
t7 = t13 * t21;
t6 = -t21 * t48 + t43;
t5 = t21 * t47 + t44;
t4 = t21 * t44 + t47;
t3 = t21 * t43 - t48;
t1 = -t50 + t55;
t10 = [0, t13, t14, t13 * t30, -t13 * sin(pkin(9)) -t14, -g(1) * (-t32 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t32 * qJ(2)) 0, 0, 0, 0, 0, t8, -t7, -t14, -t8, t7, t31 * t53 - t36 + (g(1) * t57 + g(2) * t31) * t32, -g(1) * (-t21 * t46 + t41) - g(2) * (t21 * t42 + t45) -g(1) * (-t21 * t45 - t42) - g(2) * (t21 * t41 - t46) t8, -t36 + (-g(1) * t49 - g(2) * t37) * t33 + (-g(1) * (-t57 - t37) - g(2) * t49) * t32, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(1) * (-t33 * t54 + t11) - g(2) * (-t32 * t54 + t9) - g(3) * t39, -t2 * t27, -t2 * t29, t1, -g(1) * t11 - g(2) * t9 - g(3) * (t37 + t39) + (pkin(3) + qJ(5)) * t55, 0, 0, 0, 0, 0, -t2 * t20, -t2 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t22 * t50, g(1) * t4 - g(2) * t6 - t20 * t50;];
taug_reg  = t10;
