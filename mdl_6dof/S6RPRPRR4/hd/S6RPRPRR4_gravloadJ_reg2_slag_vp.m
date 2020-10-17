% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR4
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:44:52
% EndTime: 2019-05-05 18:44:53
% DurationCPUTime: 0.43s
% Computational Cost: add. (337->86), mult. (357->113), div. (0->0), fcn. (355->10), ass. (0->57)
t35 = sin(qJ(3));
t32 = qJ(1) + pkin(10);
t24 = sin(t32);
t25 = cos(t32);
t73 = -g(1) * t25 - g(2) * t24;
t74 = t35 * t73;
t28 = t35 * qJ(4);
t38 = cos(qJ(3));
t54 = t38 * pkin(3) + t28;
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t56 = t35 * t37;
t11 = t24 * t56 + t25 * t34;
t63 = g(3) * t38;
t9 = -t24 * t34 + t25 * t56;
t72 = -g(1) * t9 - g(2) * t11 + t37 * t63;
t8 = g(3) * t35 - t38 * t73;
t69 = pkin(3) * t35;
t68 = g(1) * t24;
t62 = t38 * pkin(8);
t61 = t25 * t38;
t33 = qJ(5) + qJ(6);
t26 = sin(t33);
t60 = t26 * t35;
t27 = cos(t33);
t59 = t27 * t35;
t58 = t34 * t35;
t57 = t34 * t38;
t40 = -pkin(9) - pkin(8);
t55 = t38 * t40;
t53 = qJ(4) * t38;
t51 = t25 * t58;
t39 = cos(qJ(1));
t50 = t39 * pkin(1) + t25 * pkin(2) + t24 * pkin(7);
t36 = sin(qJ(1));
t49 = -t36 * pkin(1) + t25 * pkin(7);
t48 = pkin(3) * t61 + t25 * t28 + t50;
t16 = t24 * t53;
t18 = t25 * t53;
t47 = -g(1) * t18 - g(2) * t16;
t46 = -g(2) * t25 + t68;
t45 = g(1) * t36 - g(2) * t39;
t44 = pkin(5) * t58 - t55;
t43 = -pkin(2) - t54;
t23 = t37 * pkin(5) + pkin(4);
t14 = t46 * t38;
t13 = t46 * t35;
t12 = -t24 * t58 + t25 * t37;
t10 = t24 * t37 + t51;
t7 = -t63 - t74;
t6 = -t24 * t60 + t25 * t27;
t5 = t24 * t59 + t25 * t26;
t4 = t24 * t27 + t25 * t60;
t3 = -t24 * t26 + t25 * t59;
t2 = g(1) * t4 - g(2) * t6 - t26 * t63;
t1 = -g(1) * t3 - g(2) * t5 + t27 * t63;
t15 = [0, 0, 0, 0, 0, 0, t45, g(1) * t39 + g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t73, 0, t45 * pkin(1), 0, 0, 0, 0, 0, 0, t14, -t13, t73, -g(1) * (-t24 * pkin(2) + t49) - g(2) * t50, 0, 0, 0, 0, 0, 0, t73, -t14, t13, -g(1) * t49 - g(2) * t48 - t43 * t68, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t14, -g(1) * (t25 * pkin(4) + t49) - g(2) * (pkin(8) * t61 + t48) + (-g(1) * (t43 - t62) - g(2) * pkin(4)) * t24, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t14, -g(1) * (t25 * t23 + t49) - g(2) * (pkin(5) * t51 - t25 * t55 + t48) + (-g(1) * (t43 - t44) - g(2) * t23) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t25 * t69 + t18) - g(2) * (-t24 * t69 + t16) - g(3) * t54, 0, 0, 0, 0, 0, 0, -t8 * t34, -t8 * t37, t7, -g(3) * (t54 + t62) + t47 - (pkin(3) + pkin(8)) * t74, 0, 0, 0, 0, 0, 0, -t8 * t26, -t8 * t27, t7, -g(3) * (t44 + t54) + t47 + t73 * (pkin(5) * t57 + (-pkin(3) + t40) * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, g(1) * t10 - g(2) * t12 - g(3) * t57, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t72 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;
