% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:53:44
% EndTime: 2019-05-05 18:53:45
% DurationCPUTime: 0.36s
% Computational Cost: add. (410->92), mult. (515->113), div. (0->0), fcn. (563->10), ass. (0->58)
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t66 = g(1) * t42;
t22 = g(2) * t40 + t66;
t34 = pkin(10) + qJ(3);
t31 = sin(t34);
t71 = t22 * t31;
t38 = sin(qJ(6));
t32 = cos(t34);
t39 = sin(qJ(5));
t62 = cos(qJ(5));
t13 = t31 * t39 + t32 * t62;
t54 = t31 * t62;
t14 = -t32 * t39 + t54;
t7 = t14 * t40;
t60 = t32 * t42;
t9 = t39 * t60 - t42 * t54;
t45 = g(1) * t9 - g(2) * t7 + g(3) * t13;
t70 = t45 * t38;
t41 = cos(qJ(6));
t69 = t45 * t41;
t26 = t31 * qJ(4);
t58 = t32 * pkin(3) + t26;
t67 = pkin(3) * t31;
t64 = g(3) * t14;
t27 = t32 * pkin(4);
t37 = -pkin(7) - qJ(2);
t63 = pkin(8) + t37;
t59 = t42 * t37;
t57 = qJ(4) * t32;
t36 = cos(pkin(10));
t30 = t36 * pkin(2) + pkin(1);
t20 = t42 * t30;
t56 = pkin(3) * t60 + t42 * t26 + t20;
t55 = t27 + t58;
t53 = pkin(4) * t60 + t56;
t8 = t13 * t40;
t52 = t7 * pkin(5) + t8 * pkin(9);
t51 = -g(1) * t7 - g(2) * t9;
t10 = t13 * t42;
t50 = -t9 * pkin(5) + t10 * pkin(9);
t49 = -t13 * pkin(5) + t14 * pkin(9);
t21 = g(1) * t40 - g(2) * t42;
t48 = t8 * t38 - t42 * t41;
t47 = t42 * t38 + t8 * t41;
t46 = -t30 - t58;
t2 = g(1) * t10 + g(2) * t8 + t64;
t44 = (pkin(3) + pkin(4)) * t71;
t43 = (-g(1) * (t46 - t27) + g(2) * t63) * t40;
t19 = t42 * t57;
t17 = t40 * t57;
t12 = t21 * t32;
t11 = t21 * t31;
t6 = g(3) * t31 + t22 * t32;
t5 = -g(3) * t32 + t71;
t4 = t10 * t41 - t40 * t38;
t3 = -t10 * t38 - t40 * t41;
t1 = [0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t36, -t21 * sin(pkin(10)) -t22, -g(1) * (-t40 * pkin(1) + t42 * qJ(2)) - g(2) * (t42 * pkin(1) + t40 * qJ(2)) 0, 0, 0, 0, 0, 0, t12, -t11, -t22, -g(1) * (-t40 * t30 - t59) - g(2) * (-t40 * t37 + t20) 0, 0, 0, 0, 0, 0, t12, -t22, t11, g(1) * t59 - g(2) * t56 + (-g(1) * t46 + g(2) * t37) * t40, 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t10, -t51, t22, -g(2) * t53 + t63 * t66 + t43, 0, 0, 0, 0, 0, 0, g(1) * t47 - g(2) * t4, -g(1) * t48 - g(2) * t3, t51, -g(1) * (-t8 * pkin(5) - t42 * pkin(8) + t7 * pkin(9) - t59) - g(2) * (t10 * pkin(5) + t9 * pkin(9) + t53) + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (-t42 * t67 + t19) - g(2) * (-t40 * t67 + t17) - g(3) * t58, 0, 0, 0, 0, 0, 0, -t45, -t2, 0, -g(1) * t19 - g(2) * t17 - g(3) * t55 + t44, 0, 0, 0, 0, 0, 0, -t69, t70, t2, -g(1) * (t19 - t50) - g(2) * (t17 - t52) - g(3) * (-t49 + t55) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t2, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t70, -t2, -g(1) * t50 - g(2) * t52 - g(3) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t48 + t38 * t64, g(1) * t4 + g(2) * t47 + t41 * t64, 0, 0;];
taug_reg  = t1;
