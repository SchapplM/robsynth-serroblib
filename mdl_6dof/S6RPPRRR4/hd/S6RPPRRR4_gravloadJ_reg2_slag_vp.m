% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:41:33
% EndTime: 2019-05-05 15:41:34
% DurationCPUTime: 0.34s
% Computational Cost: add. (275->73), mult. (508->105), div. (0->0), fcn. (613->10), ass. (0->50)
t28 = sin(qJ(5));
t54 = sin(pkin(10));
t55 = cos(pkin(10));
t61 = sin(qJ(1));
t62 = cos(qJ(1));
t11 = -t61 * t54 - t62 * t55;
t12 = t62 * t54 - t61 * t55;
t30 = cos(qJ(5));
t31 = cos(qJ(4));
t58 = t28 * t31;
t40 = -t11 * t30 + t12 * t58;
t6 = t11 * t58 + t12 * t30;
t29 = sin(qJ(4));
t64 = g(3) * t29;
t68 = -g(1) * t6 - g(2) * t40 - t28 * t64;
t46 = g(1) * t11 + g(2) * t12;
t37 = -g(3) * t31 + t46 * t29;
t27 = qJ(5) + qJ(6);
t19 = sin(t27);
t60 = t19 * t31;
t20 = cos(t27);
t59 = t20 * t31;
t57 = t30 * t31;
t56 = t62 * pkin(1) + t61 * qJ(2);
t53 = t62 * pkin(2) + t56;
t52 = pkin(5) * t28 + pkin(7);
t51 = -t11 * pkin(3) + t53;
t50 = -t61 * pkin(1) + t62 * qJ(2);
t49 = t31 * pkin(4) + t29 * pkin(8);
t47 = g(1) * t12 - g(2) * t11;
t18 = t30 * pkin(5) + pkin(4);
t32 = -pkin(9) - pkin(8);
t45 = t31 * t18 - t29 * t32;
t43 = t11 * t19 + t12 * t59;
t42 = -t11 * t20 + t12 * t60;
t41 = t11 * t28 + t12 * t57;
t39 = t12 * pkin(7) + t51;
t38 = -t61 * pkin(2) + t50;
t35 = t12 * pkin(3) + t38;
t33 = t11 * pkin(7) + t35;
t14 = g(1) * t62 + g(2) * t61;
t13 = g(1) * t61 - g(2) * t62;
t8 = t47 * t29;
t7 = -t11 * t57 + t12 * t28;
t5 = -t46 * t31 - t64;
t4 = -t11 * t59 + t12 * t19;
t3 = t11 * t60 + t12 * t20;
t2 = g(1) * t4 - g(2) * t43 - t20 * t64;
t1 = -g(1) * t3 - g(2) * t42 - t19 * t64;
t9 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t14, -g(1) * t50 - g(2) * t56, 0, 0, 0, 0, 0, 0, -t47, t46, 0, -g(1) * t38 - g(2) * t53, 0, 0, 0, 0, 0, 0, -t47 * t31, t8, -t46, -g(1) * t33 - g(2) * t39, 0, 0, 0, 0, 0, 0, -g(1) * t41 - g(2) * t7, g(1) * t40 - g(2) * t6, -t8, -g(1) * (t49 * t12 + t33) - g(2) * (-t49 * t11 + t39) 0, 0, 0, 0, 0, 0, -g(1) * t43 - g(2) * t4, g(1) * t42 - g(2) * t3, -t8, -g(1) * t35 - g(2) * t51 + (-g(1) * t45 - g(2) * t52) * t12 + (-g(1) * t52 + g(2) * t45) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t30, t37 * t28, -t5, g(3) * t49 - t46 * (pkin(4) * t29 - pkin(8) * t31) 0, 0, 0, 0, 0, 0, -t37 * t20, t37 * t19, -t5, g(3) * t45 - t46 * (t18 * t29 + t31 * t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, g(1) * t7 - g(2) * t41 - t30 * t64, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t68 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t9;
