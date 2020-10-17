% Calculate inertial parameters regressor of gravitation load for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:24:30
% EndTime: 2019-05-06 08:24:31
% DurationCPUTime: 0.50s
% Computational Cost: add. (306->91), mult. (361->111), div. (0->0), fcn. (351->10), ass. (0->54)
t28 = qJ(2) + pkin(9);
t22 = sin(t28);
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t65 = -g(1) * t36 - g(2) * t34;
t68 = t22 * t65;
t24 = cos(t28);
t67 = t24 * pkin(3) + t22 * qJ(4);
t31 = -qJ(3) - pkin(7);
t58 = pkin(4) - t31;
t2 = g(3) * t22 - t24 * t65;
t33 = sin(qJ(2));
t64 = pkin(2) * t33;
t29 = sin(pkin(10));
t63 = pkin(5) * t29;
t59 = g(3) * t24;
t27 = pkin(10) + qJ(6);
t21 = sin(t27);
t57 = t21 * t36;
t23 = cos(t27);
t56 = t23 * t36;
t55 = t29 * t36;
t30 = cos(pkin(10));
t54 = t30 * t36;
t53 = t31 * t36;
t52 = t34 * t21;
t51 = t34 * t23;
t50 = t34 * t29;
t49 = t34 * t30;
t47 = pkin(5) * t30 + t58;
t46 = qJ(4) * t24;
t45 = qJ(5) * t24;
t35 = cos(qJ(2));
t25 = t35 * pkin(2);
t44 = t25 + t67;
t20 = t25 + pkin(1);
t15 = t36 * t20;
t43 = g(2) * (t67 * t36 + t15);
t42 = -pkin(3) * t22 - t64;
t12 = g(1) * t34 - g(2) * t36;
t32 = -pkin(8) - qJ(5);
t41 = t22 * t63 - t24 * t32;
t40 = -t20 - t67;
t38 = -g(3) * t35 - t33 * t65;
t11 = t36 * t46;
t9 = t34 * t46;
t8 = t12 * t24;
t7 = t12 * t22;
t6 = -t22 * t52 + t56;
t5 = t22 * t51 + t57;
t4 = t22 * t57 + t51;
t3 = t22 * t56 - t52;
t1 = -t59 - t68;
t10 = [0, 0, 0, 0, 0, 0, t12, -t65, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t35, -t12 * t33, t65, -g(1) * (-t34 * pkin(1) + pkin(7) * t36) - g(2) * (pkin(1) * t36 + t34 * pkin(7)) 0, 0, 0, 0, 0, 0, t8, -t7, t65, -g(1) * (-t34 * t20 - t53) - g(2) * (-t34 * t31 + t15) 0, 0, 0, 0, 0, 0, t65, -t8, t7, g(1) * t53 - t43 + (-g(1) * t40 + g(2) * t31) * t34, 0, 0, 0, 0, 0, 0, -g(1) * (-t22 * t50 + t54) - g(2) * (t22 * t55 + t49) -g(1) * (-t22 * t49 - t55) - g(2) * (t22 * t54 - t50) t8, -t43 + (-g(1) * t58 - g(2) * t45) * t36 + (-g(1) * (t40 - t45) - g(2) * t58) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t8, -t43 + (-g(1) * t47 - g(2) * t41) * t36 + (-g(1) * (t40 - t41) - g(2) * t47) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, g(3) * t33 - t35 * t65, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t38 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (t42 * t36 + t11) - g(2) * (t42 * t34 + t9) - g(3) * t44, 0, 0, 0, 0, 0, 0, -t2 * t29, -t2 * t30, t1, -g(1) * (-t36 * t64 + t11) - g(2) * (-t34 * t64 + t9) - g(3) * (t44 + t45) - (pkin(3) + qJ(5)) * t68, 0, 0, 0, 0, 0, 0, -t2 * t21, -t2 * t23, t1, -g(1) * t11 - g(2) * t9 - g(3) * (t41 + t44) + t65 * (t24 * t63 - t64 + (-pkin(3) + t32) * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t23 * t59, g(1) * t4 - g(2) * t6 - t21 * t59, 0, 0;];
taug_reg  = t10;
