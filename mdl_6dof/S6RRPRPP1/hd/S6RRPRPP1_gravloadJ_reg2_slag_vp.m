% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:19:12
% EndTime: 2019-05-06 12:19:13
% DurationCPUTime: 0.44s
% Computational Cost: add. (427->98), mult. (479->126), div. (0->0), fcn. (485->10), ass. (0->60)
t43 = cos(qJ(4));
t28 = t43 * pkin(4) + pkin(3);
t37 = qJ(2) + pkin(9);
t31 = sin(t37);
t33 = cos(t37);
t38 = -qJ(5) - pkin(8);
t75 = t33 * t28 - t31 * t38;
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t20 = g(1) * t45 + g(2) * t42;
t5 = -g(3) * t33 + t20 * t31;
t41 = sin(qJ(2));
t74 = pkin(2) * t41;
t71 = g(3) * t31;
t36 = qJ(4) + pkin(10);
t30 = sin(t36);
t68 = t42 * t30;
t32 = cos(t36);
t67 = t42 * t32;
t40 = sin(qJ(4));
t66 = t42 * t40;
t65 = t42 * t43;
t64 = t45 * t30;
t63 = t45 * t32;
t39 = -qJ(3) - pkin(7);
t62 = t45 * t39;
t61 = t45 * t40;
t60 = t45 * t43;
t59 = t33 * t61;
t44 = cos(qJ(2));
t34 = t44 * pkin(2);
t29 = t34 + pkin(1);
t21 = t45 * t29;
t58 = -t42 * t39 + t21;
t7 = t33 * t68 + t63;
t9 = t33 * t64 - t67;
t57 = g(1) * t7 - g(2) * t9;
t56 = t34 + t75;
t55 = t33 * pkin(3) + t31 * pkin(8);
t19 = g(1) * t42 - g(2) * t45;
t54 = -t33 * t38 - t74;
t53 = pkin(5) * t32 + qJ(6) * t30;
t12 = t33 * t66 + t60;
t1 = g(1) * t9 + g(2) * t7 + t30 * t71;
t10 = t33 * t63 + t68;
t8 = t33 * t67 - t64;
t50 = g(1) * t10 + g(2) * t8 + t32 * t71;
t49 = pkin(4) * t66 + t75 * t45 + t58;
t6 = t20 * t33 + t71;
t48 = -g(3) * t44 + t20 * t41;
t47 = pkin(4) * t61 - t62 + (-t29 - t75) * t42;
t25 = pkin(4) * t65;
t15 = t33 * t60 + t66;
t14 = -t59 + t65;
t13 = -t33 * t65 + t61;
t11 = t19 * t31;
t4 = t5 * t32;
t3 = t5 * t30;
t2 = g(1) * t8 - g(2) * t10;
t16 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t44, -t19 * t41, -t20, -g(1) * (-t42 * pkin(1) + t45 * pkin(7)) - g(2) * (t45 * pkin(1) + t42 * pkin(7)) 0, 0, 0, 0, 0, 0, t19 * t33, -t11, -t20, -g(1) * (-t42 * t29 - t62) - g(2) * t58, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, -g(1) * t12 - g(2) * t14, t11, -g(2) * t21 + (g(1) * t39 - g(2) * t55) * t45 + (-g(1) * (-t29 - t55) + g(2) * t39) * t42, 0, 0, 0, 0, 0, 0, t2, -t57, t11, -g(1) * t47 - g(2) * t49, 0, 0, 0, 0, 0, 0, t2, t11, t57, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t47) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(3) * t41 + t20 * t44, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t48 * pkin(2), 0, 0, 0, 0, 0, 0, t5 * t43, -t5 * t40, -t6, -g(3) * (t34 + t55) + t20 * (pkin(3) * t31 - pkin(8) * t33 + t74) 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * t56 + t20 * (t28 * t31 - t54) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(3) * (t33 * t53 + t56) + t20 * (-(-t28 - t53) * t31 - t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + g(2) * t12 + t40 * t71, g(1) * t15 - g(2) * t13 + t43 * t71, 0, 0, 0, 0, 0, 0, 0, 0, t1, t50, 0, -g(1) * t25 + (g(2) * t60 + t40 * t6) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t50, -g(1) * (-pkin(4) * t59 - t9 * pkin(5) + t10 * qJ(6) + t25) - g(2) * (-pkin(4) * t12 - t7 * pkin(5) + t8 * qJ(6)) - (-pkin(4) * t40 - pkin(5) * t30 + qJ(6) * t32) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t16;
