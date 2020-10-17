% Calculate inertial parameters regressor of gravitation load for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:28:59
% EndTime: 2019-05-04 20:29:00
% DurationCPUTime: 0.47s
% Computational Cost: add. (779->98), mult. (2174->155), div. (0->0), fcn. (2814->14), ass. (0->70)
t72 = sin(pkin(12));
t73 = sin(pkin(11));
t57 = t73 * t72;
t76 = cos(pkin(12));
t77 = cos(pkin(11));
t64 = t77 * t76;
t79 = cos(pkin(6));
t48 = -t64 * t79 + t57;
t74 = sin(pkin(7));
t75 = sin(pkin(6));
t61 = t75 * t74;
t78 = cos(pkin(7));
t85 = t48 * t78 + t77 * t61;
t59 = t73 * t76;
t62 = t77 * t72;
t49 = t59 * t79 + t62;
t58 = t73 * t75;
t84 = t49 * t78 - t74 * t58;
t83 = t76 * t78 * t75 + t79 * t74;
t82 = cos(qJ(3));
t38 = sin(qJ(5));
t42 = cos(qJ(4));
t81 = t38 * t42;
t41 = cos(qJ(5));
t80 = t41 * t42;
t71 = pkin(5) * t38 + pkin(9);
t30 = t62 * t79 + t59;
t40 = sin(qJ(3));
t15 = t30 * t40 + t85 * t82;
t13 = t15 * pkin(3);
t16 = t30 * t82 - t85 * t40;
t70 = t16 * pkin(9) - t13;
t31 = -t57 * t79 + t64;
t17 = t31 * t40 + t84 * t82;
t14 = t17 * pkin(3);
t18 = t31 * t82 - t84 * t40;
t69 = t18 * pkin(9) - t14;
t60 = t75 * t72;
t24 = t40 * t60 - t83 * t82;
t23 = t24 * pkin(3);
t25 = t83 * t40 + t82 * t60;
t68 = t25 * pkin(9) - t23;
t39 = sin(qJ(4));
t67 = -pkin(4) * t42 - pkin(10) * t39;
t36 = t41 * pkin(5) + pkin(4);
t37 = -qJ(6) - pkin(10);
t66 = -t36 * t42 + t37 * t39;
t63 = t77 * t75;
t44 = t49 * t74 + t58 * t78;
t11 = t18 * t39 - t42 * t44;
t47 = -t61 * t76 + t78 * t79;
t19 = t25 * t39 - t42 * t47;
t43 = t48 * t74 - t63 * t78;
t9 = t16 * t39 - t42 * t43;
t56 = g(1) * t11 + g(2) * t9 + g(3) * t19;
t10 = t16 * t42 + t39 * t43;
t12 = t18 * t42 + t39 * t44;
t20 = t25 * t42 + t39 * t47;
t55 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t54 = g(1) * t17 + g(2) * t15 + g(3) * t24;
t53 = g(1) * t18 + g(2) * t16 + g(3) * t25;
t1 = -g(1) * (-t12 * t38 + t17 * t41) - g(2) * (-t10 * t38 + t15 * t41) - g(3) * (-t20 * t38 + t24 * t41);
t28 = -g(1) * t58 + g(2) * t63 - g(3) * t79;
t8 = t54 * t39;
t6 = t56 * t41;
t5 = t56 * t38;
t4 = -g(1) * (-t17 * t80 + t18 * t38) - g(2) * (-t15 * t80 + t16 * t38) - g(3) * (-t24 * t80 + t25 * t38);
t3 = -g(1) * (t17 * t81 + t18 * t41) - g(2) * (t15 * t81 + t16 * t41) - g(3) * (t24 * t81 + t25 * t41);
t2 = -g(1) * (-t12 * t41 - t17 * t38) - g(2) * (-t10 * t41 - t15 * t38) - g(3) * (-t20 * t41 - t24 * t38);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t42, -t8, -t53, -g(1) * t69 - g(2) * t70 - g(3) * t68, 0, 0, 0, 0, 0, 0, t4, t3, t8, -g(1) * (t17 * t67 + t69) - g(2) * (t15 * t67 + t70) - g(3) * (t24 * t67 + t68) 0, 0, 0, 0, 0, 0, t4, t3, t8, -g(1) * (t17 * t66 + t18 * t71 - t14) - g(2) * (t15 * t66 + t16 * t71 - t13) - g(3) * (t24 * t66 + t25 * t71 - t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t55, -g(1) * (-t11 * pkin(4) + t12 * pkin(10)) - g(2) * (-t9 * pkin(4) + t10 * pkin(10)) - g(3) * (-t19 * pkin(4) + t20 * pkin(10)) 0, 0, 0, 0, 0, 0, t6, -t5, -t55, -g(1) * (-t11 * t36 - t12 * t37) - g(2) * (-t10 * t37 - t9 * t36) - g(3) * (-t19 * t36 - t20 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56;];
taug_reg  = t7;
