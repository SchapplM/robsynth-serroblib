% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:21:50
% EndTime: 2019-05-07 04:21:51
% DurationCPUTime: 0.43s
% Computational Cost: add. (440->93), mult. (393->108), div. (0->0), fcn. (370->10), ass. (0->60)
t38 = qJ(2) + qJ(3);
t32 = pkin(10) + t38;
t28 = sin(t32);
t29 = cos(t32);
t73 = t29 * pkin(4) + t28 * qJ(5);
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t23 = g(1) * t44 + g(2) * t41;
t72 = t23 * t28;
t4 = g(3) * t28 + t23 * t29;
t45 = -pkin(8) - pkin(7);
t33 = sin(t38);
t69 = pkin(3) * t33;
t68 = pkin(4) * t28;
t64 = g(3) * t29;
t25 = t29 * pkin(9);
t37 = -qJ(4) + t45;
t63 = pkin(5) - t37;
t39 = sin(qJ(6));
t62 = t41 * t39;
t42 = cos(qJ(6));
t61 = t41 * t42;
t60 = t44 * t37;
t59 = t44 * t39;
t58 = t44 * t42;
t40 = sin(qJ(2));
t17 = -t40 * pkin(2) - t69;
t54 = qJ(5) * t29;
t18 = t41 * t54;
t57 = t41 * t17 + t18;
t20 = t44 * t54;
t56 = t44 * t17 + t20;
t34 = cos(t38);
t30 = pkin(3) * t34;
t43 = cos(qJ(2));
t35 = t43 * pkin(2);
t55 = t30 + t35;
t53 = t30 + t73;
t52 = t35 + t53;
t16 = pkin(1) + t55;
t13 = t44 * t16;
t51 = g(2) * (t73 * t44 + t13);
t50 = -t68 - t69;
t22 = g(1) * t41 - g(2) * t44;
t49 = -t16 - t73;
t5 = -g(3) * t34 + t23 * t33;
t47 = -g(3) * t43 + t23 * t40;
t46 = (pkin(4) + pkin(9)) * t72;
t31 = t35 + pkin(1);
t12 = -t28 * t62 + t58;
t11 = t28 * t61 + t59;
t10 = t28 * t59 + t61;
t9 = t28 * t58 - t62;
t8 = t22 * t29;
t7 = t22 * t28;
t6 = g(3) * t33 + t23 * t34;
t3 = -t64 + t72;
t2 = t4 * t42;
t1 = t4 * t39;
t14 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t43, -t22 * t40, -t23, -g(1) * (-t41 * pkin(1) + t44 * pkin(7)) - g(2) * (t44 * pkin(1) + t41 * pkin(7)) 0, 0, 0, 0, 0, 0, t22 * t34, -t22 * t33, -t23, -g(1) * (-t41 * t31 - t44 * t45) - g(2) * (t44 * t31 - t41 * t45) 0, 0, 0, 0, 0, 0, t8, -t7, -t23, -g(1) * (-t41 * t16 - t60) - g(2) * (-t41 * t37 + t13) 0, 0, 0, 0, 0, 0, -t23, -t8, t7, g(1) * t60 - t51 + (-g(1) * t49 + g(2) * t37) * t41, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t8, -t51 + (-g(1) * t63 - g(2) * t25) * t44 + (-g(1) * (t49 - t25) - g(2) * t63) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, g(3) * t40 + t23 * t43, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t47 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t55 - t23 * t17, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (-t44 * t68 + t56) - g(2) * (-t41 * t68 + t57) - g(3) * t52, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t56 - g(2) * t57 - g(3) * (t25 + t52) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t50 * t44 + t20) - g(2) * (t50 * t41 + t18) - g(3) * t53, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (-t44 * t69 + t20) - g(2) * (-t41 * t69 + t18) - g(3) * (t25 + t53) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t42 * t64, g(1) * t10 - g(2) * t12 - t39 * t64, 0, 0;];
taug_reg  = t14;
