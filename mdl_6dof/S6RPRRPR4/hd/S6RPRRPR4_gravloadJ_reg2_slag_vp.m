% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:28:51
% EndTime: 2019-05-05 22:28:52
% DurationCPUTime: 0.41s
% Computational Cost: add. (446->89), mult. (370->113), div. (0->0), fcn. (357->12), ass. (0->57)
t36 = pkin(10) + qJ(3);
t31 = qJ(4) + t36;
t23 = sin(t31);
t24 = cos(t31);
t39 = cos(pkin(11));
t25 = t39 * pkin(5) + pkin(4);
t41 = -pkin(9) - qJ(5);
t73 = -t23 * t41 + t24 * t25;
t72 = t24 * pkin(4) + t23 * qJ(5);
t43 = sin(qJ(1));
t44 = cos(qJ(1));
t18 = g(1) * t44 + g(2) * t43;
t5 = -g(3) * t24 + t18 * t23;
t28 = sin(t36);
t71 = pkin(3) * t28;
t70 = pkin(4) * t23;
t30 = cos(t36);
t22 = pkin(3) * t30;
t40 = cos(pkin(10));
t26 = t40 * pkin(2) + pkin(1);
t14 = t22 + t26;
t12 = t44 * t14;
t68 = g(2) * t12;
t66 = g(3) * t23;
t35 = pkin(11) + qJ(6);
t27 = sin(t35);
t63 = t43 * t27;
t29 = cos(t35);
t62 = t43 * t29;
t37 = sin(pkin(11));
t61 = t43 * t37;
t60 = t43 * t39;
t59 = t44 * t27;
t58 = t44 * t29;
t57 = t44 * t37;
t56 = t44 * t39;
t42 = -pkin(7) - qJ(2);
t54 = qJ(5) * t24;
t34 = -pkin(8) + t42;
t53 = pkin(5) * t37 - t34;
t51 = -t70 - t71;
t17 = g(1) * t43 - g(2) * t44;
t48 = t23 * t25 + t24 * t41;
t45 = -g(3) * t30 + t18 * t28;
t16 = t44 * t54;
t15 = t43 * t54;
t11 = t17 * t23;
t10 = t24 * t58 + t63;
t9 = -t24 * t59 + t62;
t8 = -t24 * t62 + t59;
t7 = t24 * t63 + t58;
t6 = t18 * t24 + t66;
t4 = t5 * t39;
t3 = t5 * t37;
t2 = t5 * t29;
t1 = t5 * t27;
t13 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t40, -t17 * sin(pkin(10)) -t18, -g(1) * (-t43 * pkin(1) + t44 * qJ(2)) - g(2) * (t44 * pkin(1) + t43 * qJ(2)) 0, 0, 0, 0, 0, 0, t17 * t30, -t17 * t28, -t18, -g(1) * (-t43 * t26 - t44 * t42) - g(2) * (t44 * t26 - t43 * t42) 0, 0, 0, 0, 0, 0, t17 * t24, -t11, -t18, -g(1) * (-t43 * t14 - t44 * t34) - g(2) * (-t43 * t34 + t12) 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t60 + t57) - g(2) * (t24 * t56 + t61) -g(1) * (t24 * t61 + t56) - g(2) * (-t24 * t57 + t60) t11, -t68 + (g(1) * t34 - g(2) * t72) * t44 + (-g(1) * (-t14 - t72) + g(2) * t34) * t43, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t11, -t68 + (-g(1) * t53 - g(2) * t73) * t44 + (-g(1) * (-t14 - t73) - g(2) * t53) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t28 + t18 * t30, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t45 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t51 * t44 + t16) - g(2) * (t51 * t43 + t15) - g(3) * (t22 + t72) 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(3) * (t22 + t73) + t18 * (t48 + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t44 * t70 + t16) - g(2) * (-t43 * t70 + t15) - g(3) * t72, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(3) * t73 + t18 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t27 * t66, g(1) * t10 - g(2) * t8 + t29 * t66, 0, 0;];
taug_reg  = t13;
