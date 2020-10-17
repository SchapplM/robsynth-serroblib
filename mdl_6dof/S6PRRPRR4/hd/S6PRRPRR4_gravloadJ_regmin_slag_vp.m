% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:10:48
% EndTime: 2019-05-05 05:10:49
% DurationCPUTime: 0.34s
% Computational Cost: add. (305->76), mult. (830->126), div. (0->0), fcn. (1059->12), ass. (0->50)
t32 = sin(qJ(6));
t35 = sin(qJ(2));
t39 = cos(qJ(2));
t31 = cos(pkin(11));
t54 = cos(pkin(6));
t50 = t31 * t54;
t53 = sin(pkin(11));
t21 = t35 * t50 + t53 * t39;
t34 = sin(qJ(3));
t30 = sin(pkin(6));
t38 = cos(qJ(3));
t56 = t30 * t38;
t15 = t21 * t34 + t31 * t56;
t16 = -t31 * t30 * t34 + t21 * t38;
t46 = t54 * t53;
t23 = t31 * t39 - t35 * t46;
t51 = t30 * t53;
t17 = t23 * t34 - t38 * t51;
t18 = t23 * t38 + t34 * t51;
t57 = t30 * t35;
t24 = t34 * t57 - t54 * t38;
t25 = t54 * t34 + t35 * t56;
t33 = sin(qJ(5));
t37 = cos(qJ(5));
t43 = g(1) * (t17 * t37 - t18 * t33) + g(2) * (t15 * t37 - t16 * t33) + g(3) * (t24 * t37 - t25 * t33);
t65 = t43 * t32;
t36 = cos(qJ(6));
t64 = t43 * t36;
t55 = t30 * t39;
t20 = -t53 * t35 + t39 * t50;
t59 = g(2) * t20;
t22 = -t31 * t35 - t39 * t46;
t60 = g(1) * t22;
t41 = g(3) * t55 + t59 + t60;
t14 = t24 * t33 + t25 * t37;
t4 = t15 * t33 + t16 * t37;
t7 = t17 * t33 + t18 * t37;
t63 = g(1) * t7 + g(2) * t4 + g(3) * t14;
t49 = -g(1) * t23 - g(2) * t21;
t47 = t33 * t34 + t37 * t38;
t45 = pkin(3) * t38 + qJ(4) * t34 + pkin(2);
t1 = g(1) * t17 + g(2) * t15 + g(3) * t24;
t42 = g(1) * t18 + g(2) * t16 + g(3) * t25;
t40 = g(3) * t57 - t49;
t19 = t47 * t55;
t11 = t47 * t22;
t10 = t47 * t20;
t9 = t41 * t38;
t8 = t41 * t34;
t2 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t41, t40, 0, 0, 0, 0, 0, -t9, t8, -t9, -t40, -t8, t49 * pkin(8) - t45 * t60 - t45 * t59 - g(3) * (pkin(8) * t35 + t45 * t39) * t30, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t10 - g(3) * t19, t41 * (t33 * t38 - t34 * t37) 0, 0, 0, 0, 0, -g(1) * (t11 * t36 - t23 * t32) - g(2) * (t10 * t36 - t21 * t32) - g(3) * (t19 * t36 - t32 * t57) -g(1) * (-t11 * t32 - t23 * t36) - g(2) * (-t10 * t32 - t21 * t36) - g(3) * (-t19 * t32 - t36 * t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t42, t1, 0, -t42, -g(1) * (-t17 * pkin(3) + t18 * qJ(4)) - g(2) * (-t15 * pkin(3) + t16 * qJ(4)) - g(3) * (-t24 * pkin(3) + t25 * qJ(4)) 0, 0, 0, 0, 0, t43, -t63, 0, 0, 0, 0, 0, t64, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t63, 0, 0, 0, 0, 0, -t64, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t36 - t7 * t32) - g(2) * (t20 * t36 - t4 * t32) - g(3) * (-t14 * t32 + t36 * t55) -g(1) * (-t22 * t32 - t7 * t36) - g(2) * (-t20 * t32 - t4 * t36) - g(3) * (-t14 * t36 - t32 * t55);];
taug_reg  = t2;
