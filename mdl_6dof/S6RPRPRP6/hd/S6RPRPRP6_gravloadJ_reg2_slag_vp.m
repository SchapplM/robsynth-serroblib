% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:54:19
% EndTime: 2019-05-05 17:54:20
% DurationCPUTime: 0.35s
% Computational Cost: add. (284->81), mult. (360->94), div. (0->0), fcn. (354->8), ass. (0->51)
t30 = pkin(9) + qJ(3);
t27 = sin(t30);
t28 = cos(t30);
t47 = t28 * pkin(3) + t27 * qJ(4);
t34 = -pkin(7) - qJ(2);
t55 = pkin(4) - t34;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t20 = g(1) * t38 + g(2) * t36;
t67 = t20 * t27;
t37 = cos(qJ(5));
t48 = t38 * t37;
t35 = sin(qJ(5));
t52 = t36 * t35;
t11 = t27 * t48 - t52;
t49 = t38 * t35;
t51 = t36 * t37;
t13 = t27 * t51 + t49;
t58 = g(3) * t28;
t1 = -g(1) * t11 - g(2) * t13 + t37 * t58;
t8 = g(3) * t27 + t20 * t28;
t65 = pkin(3) * t27;
t64 = pkin(5) * t35;
t57 = t28 * pkin(8);
t33 = -qJ(6) - pkin(8);
t54 = t28 * t33;
t53 = t28 * t35;
t50 = t38 * t34;
t46 = t37 * pkin(5) + t55;
t45 = qJ(4) * t28;
t44 = pkin(5) * t53;
t32 = cos(pkin(9));
t25 = t32 * pkin(2) + pkin(1);
t18 = t38 * t25;
t42 = g(2) * (t47 * t38 + t18);
t19 = g(1) * t36 - g(2) * t38;
t41 = t27 * t64 - t54;
t40 = -t25 - t47;
t17 = t38 * t45;
t15 = t36 * t45;
t14 = -t27 * t52 + t48;
t12 = t27 * t49 + t51;
t10 = t19 * t28;
t9 = t19 * t27;
t7 = -t58 + t67;
t6 = t8 * t37;
t5 = t8 * t35;
t4 = -g(1) * t14 - g(2) * t12;
t3 = g(1) * t13 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t14 - g(3) * t53;
t16 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t32, -t19 * sin(pkin(9)) -t20, -g(1) * (-t36 * pkin(1) + t38 * qJ(2)) - g(2) * (t38 * pkin(1) + t36 * qJ(2)) 0, 0, 0, 0, 0, 0, t10, -t9, -t20, -g(1) * (-t36 * t25 - t50) - g(2) * (-t36 * t34 + t18) 0, 0, 0, 0, 0, 0, -t20, -t10, t9, g(1) * t50 - t42 + (-g(1) * t40 + g(2) * t34) * t36, 0, 0, 0, 0, 0, 0, t4, t3, t10, -t42 + (-g(1) * t55 - g(2) * t57) * t38 + (-g(1) * (t40 - t57) - g(2) * t55) * t36, 0, 0, 0, 0, 0, 0, t4, t3, t10, -t42 + (-g(1) * t46 - g(2) * t41) * t38 + (-g(1) * (t40 - t41) - g(2) * t46) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t38 * t65 + t17) - g(2) * (-t36 * t65 + t15) - g(3) * t47, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * t17 - g(2) * t15 - g(3) * (t47 + t57) + (pkin(3) + pkin(8)) * t67, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * (t38 * t44 + t17) - g(2) * (t36 * t44 + t15) - g(3) * (t47 - t54) + (-g(3) * t64 + t20 * (pkin(3) - t33)) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg  = t16;
