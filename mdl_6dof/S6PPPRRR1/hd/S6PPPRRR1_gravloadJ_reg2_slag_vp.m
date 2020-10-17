% Calculate inertial parameters regressor of gravitation load for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPPRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:33:58
% EndTime: 2019-05-04 19:33:59
% DurationCPUTime: 0.48s
% Computational Cost: add. (1133->96), mult. (3232->165), div. (0->0), fcn. (4281->18), ass. (0->78)
t76 = sin(pkin(13));
t77 = sin(pkin(12));
t61 = t77 * t76;
t82 = cos(pkin(13));
t83 = cos(pkin(12));
t70 = t83 * t82;
t86 = cos(pkin(6));
t50 = t86 * t70 - t61;
t85 = cos(pkin(7));
t47 = t50 * t85;
t63 = t77 * t82;
t68 = t83 * t76;
t51 = t86 * t68 + t63;
t79 = sin(pkin(7));
t81 = cos(pkin(14));
t65 = t79 * t81;
t80 = sin(pkin(6));
t55 = t80 * t65;
t75 = sin(pkin(14));
t38 = -t81 * t47 + t51 * t75 + t83 * t55;
t69 = t83 * t80;
t44 = t50 * t79 + t85 * t69;
t78 = sin(pkin(8));
t84 = cos(pkin(8));
t92 = t38 * t84 + t44 * t78;
t52 = -t86 * t63 - t68;
t48 = t52 * t85;
t53 = -t86 * t61 + t70;
t39 = -t81 * t48 + t53 * t75 - t77 * t55;
t62 = t77 * t80;
t45 = t52 * t79 - t85 * t62;
t91 = t39 * t84 + t45 * t78;
t67 = t82 * t80;
t56 = t85 * t67;
t66 = t80 * t76;
t43 = -t81 * t56 - t86 * t65 + t75 * t66;
t49 = t79 * t67 - t86 * t85;
t90 = t43 * t84 + t49 * t78;
t89 = cos(qJ(4));
t31 = sin(qJ(6));
t35 = cos(qJ(5));
t88 = t31 * t35;
t34 = cos(qJ(6));
t87 = t34 * t35;
t64 = t79 * t75;
t54 = t80 * t64;
t22 = t75 * t47 + t51 * t81 - t83 * t54;
t33 = sin(qJ(4));
t10 = t22 * t33 + t92 * t89;
t11 = t22 * t89 - t92 * t33;
t74 = -t10 * pkin(4) + t11 * pkin(10);
t23 = t75 * t48 + t53 * t81 + t77 * t54;
t12 = t23 * t33 + t91 * t89;
t13 = t23 * t89 - t91 * t33;
t73 = -t12 * pkin(4) + t13 * pkin(10);
t28 = t75 * t56 + t86 * t64 + t81 * t66;
t15 = t28 * t33 + t90 * t89;
t16 = t28 * t89 - t90 * t33;
t72 = -t15 * pkin(4) + t16 * pkin(10);
t32 = sin(qJ(5));
t71 = -pkin(5) * t35 - pkin(11) * t32;
t17 = t38 * t78 - t44 * t84;
t2 = -t11 * t32 + t17 * t35;
t18 = t39 * t78 - t45 * t84;
t4 = -t13 * t32 + t18 * t35;
t24 = t43 * t78 - t49 * t84;
t8 = -t16 * t32 + t24 * t35;
t60 = g(1) * t4 + g(2) * t2 + g(3) * t8;
t3 = t11 * t35 + t17 * t32;
t5 = t13 * t35 + t18 * t32;
t9 = t16 * t35 + t24 * t32;
t59 = g(1) * t5 + g(2) * t3 + g(3) * t9;
t58 = g(1) * t12 + g(2) * t10 + g(3) * t15;
t57 = g(1) * t13 + g(2) * t11 + g(3) * t16;
t30 = -g(1) * t62 + g(2) * t69 - g(3) * t86;
t19 = g(1) * t45 + g(2) * t44 + g(3) * t49;
t1 = t58 * t32;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t35, -t1, -t57, -g(1) * t73 - g(2) * t74 - g(3) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t87 + t13 * t31) - g(2) * (-t10 * t87 + t11 * t31) - g(3) * (-t15 * t87 + t16 * t31) -g(1) * (t12 * t88 + t13 * t34) - g(2) * (t10 * t88 + t11 * t34) - g(3) * (t15 * t88 + t16 * t34) t1, -g(1) * (t71 * t12 + t73) - g(2) * (t71 * t10 + t74) - g(3) * (t71 * t15 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * t34, t60 * t31, -t59, -g(1) * (t4 * pkin(5) + t5 * pkin(11)) - g(2) * (t2 * pkin(5) + t3 * pkin(11)) - g(3) * (t8 * pkin(5) + t9 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t34 - t5 * t31) - g(2) * (t10 * t34 - t3 * t31) - g(3) * (t15 * t34 - t9 * t31) -g(1) * (-t12 * t31 - t5 * t34) - g(2) * (-t10 * t31 - t3 * t34) - g(3) * (-t15 * t31 - t9 * t34) 0, 0;];
taug_reg  = t6;
