% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:04:53
% EndTime: 2019-05-06 09:04:54
% DurationCPUTime: 0.38s
% Computational Cost: add. (403->89), mult. (437->112), div. (0->0), fcn. (442->10), ass. (0->54)
t34 = cos(pkin(10));
t23 = t34 * pkin(4) + pkin(3);
t32 = qJ(2) + pkin(9);
t26 = sin(t32);
t28 = cos(t32);
t36 = -pkin(8) - qJ(4);
t69 = t28 * t23 - t26 * t36;
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t16 = g(1) * t40 + g(2) * t38;
t5 = -g(3) * t28 + t16 * t26;
t37 = sin(qJ(2));
t68 = pkin(2) * t37;
t65 = g(3) * t26;
t31 = pkin(10) + qJ(5);
t25 = sin(t31);
t62 = t38 * t25;
t27 = cos(t31);
t61 = t38 * t27;
t33 = sin(pkin(10));
t60 = t38 * t33;
t59 = t38 * t34;
t58 = t40 * t25;
t57 = t40 * t27;
t56 = t40 * t33;
t55 = t40 * t34;
t35 = -qJ(3) - pkin(7);
t54 = t40 * t35;
t39 = cos(qJ(2));
t29 = t39 * pkin(2);
t24 = t29 + pkin(1);
t17 = t40 * t24;
t53 = -t38 * t35 + t17;
t7 = t28 * t62 + t57;
t9 = t28 * t58 - t61;
t52 = g(1) * t7 - g(2) * t9;
t51 = t29 + t69;
t15 = g(1) * t38 - g(2) * t40;
t50 = -t28 * t36 - t68;
t49 = t28 * pkin(3) + t26 * qJ(4);
t48 = pkin(5) * t27 + qJ(6) * t25;
t1 = g(1) * t9 + g(2) * t7 + t25 * t65;
t10 = t28 * t57 + t62;
t8 = t28 * t61 - t58;
t45 = g(1) * t10 + g(2) * t8 + t27 * t65;
t44 = pkin(4) * t60 + t69 * t40 + t53;
t43 = -g(3) * t39 + t16 * t37;
t42 = pkin(4) * t56 - t54 + (-t24 - t69) * t38;
t11 = t15 * t26;
t6 = t16 * t28 + t65;
t4 = t5 * t27;
t3 = t5 * t25;
t2 = g(1) * t8 - g(2) * t10;
t12 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t39, -t15 * t37, -t16, -g(1) * (-t38 * pkin(1) + t40 * pkin(7)) - g(2) * (t40 * pkin(1) + t38 * pkin(7)) 0, 0, 0, 0, 0, 0, t15 * t28, -t11, -t16, -g(1) * (-t38 * t24 - t54) - g(2) * t53, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t59 + t56) - g(2) * (t28 * t55 + t60) -g(1) * (t28 * t60 + t55) - g(2) * (-t28 * t56 + t59) t11, -g(2) * t17 + (g(1) * t35 - g(2) * t49) * t40 + (-g(1) * (-t24 - t49) + g(2) * t35) * t38, 0, 0, 0, 0, 0, 0, t2, -t52, t11, -g(1) * t42 - g(2) * t44, 0, 0, 0, 0, 0, 0, t2, t11, t52, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t42) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t37 + t16 * t39, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t43 * pkin(2), 0, 0, 0, 0, 0, 0, t5 * t34, -t5 * t33, -t6, -g(3) * (t29 + t49) + t16 * (pkin(3) * t26 - qJ(4) * t28 + t68) 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * t51 + t16 * (t23 * t26 - t50) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(3) * (t48 * t28 + t51) + t16 * (-(-t23 - t48) * t26 - t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t45, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t45, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t25 + qJ(6) * t27) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
