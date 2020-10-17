% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:42:45
% EndTime: 2019-05-05 16:42:46
% DurationCPUTime: 0.30s
% Computational Cost: add. (259->76), mult. (287->94), div. (0->0), fcn. (276->8), ass. (0->48)
t26 = qJ(1) + pkin(9);
t20 = cos(t26);
t19 = sin(t26);
t53 = g(2) * t19;
t10 = g(1) * t20 + t53;
t28 = sin(qJ(3));
t57 = t10 * t28;
t31 = cos(qJ(3));
t2 = g(3) * t28 + t10 * t31;
t56 = pkin(3) + pkin(4);
t55 = g(1) * t19;
t51 = g(3) * t31;
t23 = t31 * pkin(3);
t22 = t31 * pkin(4);
t50 = g(2) * qJ(5);
t49 = t20 * t28;
t48 = t20 * t31;
t27 = sin(qJ(6));
t47 = t27 * t28;
t30 = cos(qJ(6));
t46 = t28 * t30;
t21 = t28 * qJ(4);
t45 = t23 + t21;
t44 = qJ(4) * t31;
t43 = pkin(8) + t56;
t32 = cos(qJ(1));
t42 = t32 * pkin(1) + t20 * pkin(2) + t19 * pkin(7);
t41 = t22 + t45;
t29 = sin(qJ(1));
t40 = -t29 * pkin(1) + t20 * pkin(7);
t39 = -pkin(2) - t21;
t38 = g(1) * t43;
t37 = pkin(3) * t48 + t20 * t21 + t42;
t9 = -g(2) * t20 + t55;
t36 = g(1) * t29 - g(2) * t32;
t35 = pkin(4) * t48 + t37;
t34 = t39 - t23;
t33 = g(1) * (-t20 * qJ(5) + t40);
t13 = t20 * t44;
t11 = t19 * t44;
t8 = t9 * t31;
t7 = t9 * t28;
t6 = -t19 * t27 + t20 * t46;
t5 = -t19 * t30 - t20 * t47;
t4 = -t19 * t46 - t20 * t27;
t3 = t19 * t47 - t20 * t30;
t1 = -t51 + t57;
t12 = [0, 0, 0, 0, 0, 0, t36, g(1) * t32 + g(2) * t29, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t19 * pkin(2) + t40) - g(2) * t42, 0, 0, 0, 0, 0, 0, t8, -t10, t7, -g(1) * t40 - g(2) * t37 - t34 * t55, 0, 0, 0, 0, 0, 0, t7, -t8, t10, -t33 - g(2) * t35 + (-g(1) * (t34 - t22) + t50) * t19, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t8, -t33 - g(2) * (pkin(5) * t49 + pkin(8) * t48 + t35) + (-g(1) * (-t28 * pkin(5) + t39) + t50 + t31 * t38) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-pkin(3) * t49 + t13) - g(2) * (-t19 * t28 * pkin(3) + t11) - g(3) * t45, 0, 0, 0, 0, 0, 0, -t2, -t1, 0, -g(1) * t13 - g(2) * t11 - g(3) * t41 + t56 * t57, 0, 0, 0, 0, 0, 0, -t2 * t30, t2 * t27, t1, -g(1) * (pkin(5) * t48 + t13) - g(2) * (t19 * t31 * pkin(5) + t11) - g(3) * (t31 * pkin(8) + t41) + (-g(3) * pkin(5) + t20 * t38 + t43 * t53) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 - t27 * t51, g(1) * t6 - g(2) * t4 - t30 * t51, 0, 0;];
taug_reg  = t12;
