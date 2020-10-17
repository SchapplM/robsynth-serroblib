% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:13:30
% EndTime: 2019-05-05 19:13:31
% DurationCPUTime: 0.31s
% Computational Cost: add. (278->79), mult. (291->90), div. (0->0), fcn. (277->10), ass. (0->46)
t31 = qJ(3) + pkin(10);
t25 = qJ(5) + t31;
t21 = sin(t25);
t22 = cos(t25);
t44 = -t21 * pkin(5) + t22 * pkin(9);
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t57 = g(2) * t38;
t17 = g(1) * t35 - t57;
t59 = -g(3) * t21 + t17 * t22;
t58 = pkin(5) * t22;
t55 = g(3) * t22;
t34 = sin(qJ(3));
t53 = t34 * pkin(3);
t52 = t21 * t35;
t33 = sin(qJ(6));
t51 = t35 * t33;
t36 = cos(qJ(6));
t50 = t35 * t36;
t49 = t38 * t33;
t48 = t38 * t36;
t32 = -qJ(4) - pkin(7);
t47 = pkin(9) * t52 + t35 * t58;
t46 = t38 * pkin(1) + t35 * qJ(2);
t27 = t38 * qJ(2);
t45 = -t35 * pkin(1) + t27;
t23 = sin(t31);
t12 = pkin(4) * t23 + t53;
t43 = -pkin(9) * t21 - t58;
t18 = g(1) * t38 + g(2) * t35;
t30 = -pkin(8) + t32;
t41 = t38 * t12 + t35 * t30 + t45;
t40 = t35 * t12 - t38 * t30 + t46;
t37 = cos(qJ(3));
t39 = g(3) * t34 - t17 * t37;
t24 = cos(t31);
t13 = t37 * pkin(3) + pkin(4) * t24;
t9 = t21 * t48 - t51;
t8 = t21 * t49 + t50;
t7 = t21 * t50 + t49;
t6 = -t21 * t51 + t48;
t5 = t18 * t22;
t3 = g(1) * t52 - t21 * t57 + t55;
t2 = t59 * t36;
t1 = t59 * t33;
t4 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, -g(1) * t45 - g(2) * t46, 0, 0, 0, 0, 0, 0, -t18 * t34, -t18 * t37, t17, -g(1) * (t27 + (-pkin(1) - pkin(7)) * t35) - g(2) * (t38 * pkin(7) + t46) 0, 0, 0, 0, 0, 0, -t18 * t23, -t18 * t24, t17, -g(1) * (t38 * t53 + t27 + (-pkin(1) + t32) * t35) - g(2) * (-t38 * t32 + t35 * t53 + t46) 0, 0, 0, 0, 0, 0, -t18 * t21, -t5, t17, -g(1) * t41 - g(2) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7, g(1) * t8 - g(2) * t6, t5, -g(1) * (-t38 * t44 + t41) - g(2) * (-t35 * t44 + t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, g(3) * t37 + t17 * t34, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t23 - t17 * t24, g(3) * t24 + t17 * t23, 0, t39 * pkin(3), 0, 0, 0, 0, 0, 0, -t59, t3, 0, g(3) * t12 - t17 * t13, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t35 * t13 + t47) - g(3) * (-t12 + t44) - (-t13 + t43) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * t47 - g(3) * t44 - t43 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8 + t33 * t55, g(1) * t7 - g(2) * t9 + t36 * t55, 0, 0;];
taug_reg  = t4;
