% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [6x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:40:28
% EndTime: 2019-05-05 03:40:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (287->90), mult. (535->148), div. (0->0), fcn. (635->12), ass. (0->59)
t29 = sin(pkin(6));
t66 = g(3) * t29;
t28 = sin(pkin(10));
t30 = cos(pkin(10));
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t50 = cos(pkin(6));
t45 = t38 * t50;
t12 = t28 * t35 - t30 * t45;
t46 = t35 * t50;
t13 = t28 * t38 + t30 * t46;
t37 = cos(qJ(3));
t24 = t37 * pkin(3) + pkin(2);
t32 = -qJ(4) - pkin(8);
t65 = -t12 * t24 - t13 * t32;
t14 = t28 * t45 + t30 * t35;
t15 = -t28 * t46 + t30 * t38;
t64 = -t14 * t24 - t15 * t32;
t33 = sin(qJ(5));
t63 = t13 * t33;
t62 = t15 * t33;
t34 = sin(qJ(3));
t61 = t15 * t34;
t27 = qJ(3) + pkin(11);
t26 = cos(t27);
t60 = t26 * t33;
t36 = cos(qJ(5));
t59 = t26 * t36;
t58 = t28 * t29;
t57 = t29 * t30;
t56 = t29 * t34;
t55 = t29 * t35;
t54 = t29 * t37;
t53 = t29 * t38;
t52 = t33 * t38;
t51 = t36 * t38;
t49 = t28 * t54;
t48 = t34 * t55;
t47 = t30 * t54;
t44 = t50 * t37;
t23 = t36 * pkin(5) + pkin(4);
t25 = sin(t27);
t31 = -qJ(6) - pkin(9);
t43 = t23 * t26 - t25 * t31;
t42 = -t13 * t34 - t47;
t2 = t13 * t25 + t26 * t57;
t4 = t15 * t25 - t26 * t58;
t8 = t25 * t55 - t50 * t26;
t41 = g(1) * t4 + g(2) * t2 + g(3) * t8;
t1 = -g(1) * t14 - g(2) * t12 + g(3) * t53;
t40 = g(1) * t15 + g(2) * t13 + g(3) * t55;
t3 = t13 * t26 - t25 * t57;
t5 = t15 * t26 + t25 * t58;
t9 = t50 * t25 + t26 * t55;
t39 = -g(1) * (t14 * t36 - t5 * t33) - g(2) * (t12 * t36 - t3 * t33) - g(3) * (-t29 * t51 - t9 * t33);
t22 = pkin(3) * t44;
t18 = pkin(3) * t49;
t16 = t24 * t53;
t6 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t1, t40, 0, 0, 0, 0, 0, -t1 * t37, t1 * t34, -t40, -g(1) * t64 - g(2) * t65 - g(3) * (-t32 * t55 + t16) 0, 0, 0, 0, 0, -g(1) * (-t14 * t59 + t62) - g(2) * (-t12 * t59 + t63) - (t26 * t51 + t33 * t35) * t66, -g(1) * (t14 * t60 + t15 * t36) - g(2) * (t12 * t60 + t13 * t36) - (-t26 * t52 + t35 * t36) * t66, -t1 * t25, -g(1) * (pkin(5) * t62 - t43 * t14 + t64) - g(2) * (pkin(5) * t63 - t43 * t12 + t65) - g(3) * t16 - (t43 * t38 + (pkin(5) * t33 - t32) * t35) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t49 - t61) - g(2) * t42 - g(3) * (t44 - t48) -g(1) * (-t15 * t37 - t28 * t56) - g(2) * (-t13 * t37 + t30 * t56) - g(3) * (-t50 * t34 - t35 * t54) 0, -g(1) * t18 - g(3) * t22 + (g(2) * t47 + t34 * t40) * pkin(3), 0, 0, 0, 0, 0, t41 * t36, -t41 * t33, -g(1) * t5 - g(2) * t3 - g(3) * t9, -g(1) * (-pkin(3) * t61 - t4 * t23 - t5 * t31 + t18) - g(2) * (t42 * pkin(3) - t2 * t23 - t3 * t31) - g(3) * (-pkin(3) * t48 - t8 * t23 - t9 * t31 + t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -g(1) * (-t14 * t33 - t5 * t36) - g(2) * (-t12 * t33 - t3 * t36) - g(3) * (t29 * t52 - t9 * t36) 0, t39 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41;];
taug_reg  = t6;
