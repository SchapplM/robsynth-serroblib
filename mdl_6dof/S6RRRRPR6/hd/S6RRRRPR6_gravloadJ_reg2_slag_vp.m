% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:51:27
% EndTime: 2019-05-07 20:51:29
% DurationCPUTime: 0.68s
% Computational Cost: add. (594->123), mult. (615->180), div. (0->0), fcn. (626->12), ass. (0->83)
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t34 = g(1) * t59 + g(2) * t56;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t74 = t59 * t57;
t58 = cos(qJ(2));
t85 = t56 * t58;
t24 = t54 * t85 + t74;
t75 = t59 * t54;
t26 = t56 * t57 - t58 * t75;
t55 = sin(qJ(2));
t90 = g(3) * t55;
t99 = -g(1) * t26 + g(2) * t24 + t54 * t90;
t53 = qJ(3) + qJ(4);
t44 = sin(t53);
t45 = cos(t53);
t76 = t59 * t45;
t16 = t44 * t85 + t76;
t77 = t59 * t44;
t18 = t56 * t45 - t58 * t77;
t5 = -g(1) * t18 + g(2) * t16 + t44 * t90;
t22 = -g(3) * t58 + t34 * t55;
t60 = -pkin(9) - pkin(8);
t98 = pkin(4) * t44;
t95 = g(1) * t56;
t88 = t54 * pkin(3);
t87 = t55 * t59;
t86 = t55 * t60;
t84 = t58 * t59;
t43 = pkin(11) + t53;
t38 = sin(t43);
t28 = -pkin(5) * t38 - t98;
t20 = -t28 + t88;
t83 = t59 * t20;
t32 = t88 + t98;
t82 = t59 * t32;
t41 = qJ(6) + t43;
t36 = sin(t41);
t81 = t59 * t36;
t37 = cos(t41);
t80 = t59 * t37;
t79 = t59 * t38;
t39 = cos(t43);
t78 = t59 * t39;
t40 = pkin(4) * t45;
t48 = t57 * pkin(3);
t33 = t40 + t48;
t73 = t59 * pkin(1) + t56 * pkin(7);
t52 = -qJ(5) + t60;
t35 = pkin(5) * t39;
t21 = t35 + t33;
t70 = t58 * pkin(2) + t55 * pkin(8);
t68 = -g(2) * t59 + t95;
t15 = pkin(2) + t21;
t46 = -pkin(10) + t52;
t67 = t58 * t15 - t55 * t46;
t31 = pkin(2) + t33;
t65 = t58 * t31 - t55 * t52;
t42 = t48 + pkin(2);
t63 = t58 * t42 - t86;
t49 = t59 * pkin(7);
t30 = t68 * t55;
t29 = t35 + t40;
t27 = t56 * t54 + t58 * t74;
t25 = -t57 * t85 + t75;
t23 = t34 * t58 + t90;
t19 = t56 * t44 + t58 * t76;
t17 = -t45 * t85 + t77;
t14 = t56 * t38 + t58 * t78;
t13 = t56 * t39 - t58 * t79;
t12 = -t39 * t85 + t79;
t11 = t38 * t85 + t78;
t10 = t56 * t36 + t58 * t80;
t9 = t56 * t37 - t58 * t81;
t8 = -t37 * t85 + t81;
t7 = t36 * t85 + t80;
t6 = g(1) * t19 - g(2) * t17 + t45 * t90;
t4 = g(1) * t14 - g(2) * t12 + t39 * t90;
t3 = -g(1) * t13 + g(2) * t11 + t38 * t90;
t2 = g(1) * t10 - g(2) * t8 + t37 * t90;
t1 = -g(1) * t9 + g(2) * t7 + t36 * t90;
t47 = [0, 0, 0, 0, 0, 0, t68, t34, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t58, -t30, -t34, -g(1) * (-t56 * pkin(1) + t49) - g(2) * t73, 0, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t27, -g(1) * t24 - g(2) * t26, t30, -g(1) * t49 - g(2) * (t70 * t59 + t73) - (-pkin(1) - t70) * t95, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, -g(1) * t16 - g(2) * t18, t30, -g(1) * (pkin(3) * t75 + t49) - g(2) * (t42 * t84 - t59 * t86 + t73) + (-g(1) * (-pkin(1) - t63) - g(2) * t88) * t56, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t30, -g(1) * (t49 + t82) - g(2) * (t31 * t84 - t52 * t87 + t73) + (-g(1) * (-pkin(1) - t65) - g(2) * t32) * t56, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t30, -g(1) * (t49 + t83) - g(2) * (t15 * t84 - t46 * t87 + t73) + (-g(1) * (-pkin(1) - t67) - g(2) * t20) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t57, -t22 * t54, -t23, -g(3) * t70 + t34 * (pkin(2) * t55 - pkin(8) * t58) 0, 0, 0, 0, 0, 0, t22 * t45, -t22 * t44, -t23, -g(3) * t63 + t34 * (t42 * t55 + t58 * t60) 0, 0, 0, 0, 0, 0, t22 * t39, -t22 * t38, -t23, -g(3) * t65 + t34 * (t31 * t55 + t52 * t58) 0, 0, 0, 0, 0, 0, t22 * t37, -t22 * t36, -t23, -g(3) * t67 + t34 * (t15 * t55 + t46 * t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, g(1) * t27 - g(2) * t25 + t57 * t90, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t99 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (t56 * t33 - t58 * t82) - g(2) * (-t32 * t85 - t59 * t33) + t32 * t90, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t56 * t21 - t58 * t83) - g(2) * (-t20 * t85 - t59 * t21) + t20 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t28 * t84 + t56 * t29) - g(2) * (t28 * t85 - t59 * t29) - t28 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t47;
