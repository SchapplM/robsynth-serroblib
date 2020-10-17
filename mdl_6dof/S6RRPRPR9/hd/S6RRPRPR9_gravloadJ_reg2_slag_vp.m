% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:18:28
% EndTime: 2019-05-06 15:18:30
% DurationCPUTime: 0.93s
% Computational Cost: add. (703->155), mult. (1182->241), div. (0->0), fcn. (1424->14), ass. (0->72)
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t58 = cos(qJ(2));
t74 = cos(pkin(6));
t93 = cos(qJ(1));
t66 = t74 * t93;
t27 = t56 * t66 + t57 * t58;
t48 = pkin(11) + qJ(4);
t43 = sin(t48);
t45 = cos(t48);
t51 = sin(pkin(6));
t72 = t51 * t93;
t11 = t27 * t45 - t43 * t72;
t26 = t57 * t56 - t58 * t66;
t47 = pkin(12) + qJ(6);
t42 = sin(t47);
t44 = cos(t47);
t97 = t11 * t42 - t26 * t44;
t96 = t11 * t44 + t26 * t42;
t53 = cos(pkin(11));
t41 = t53 * pkin(3) + pkin(2);
t79 = t51 * t58;
t30 = t41 * t79;
t95 = g(3) * t30;
t94 = g(3) * t51;
t49 = sin(pkin(12));
t90 = t26 * t49;
t89 = t27 * t49;
t70 = t57 * t74;
t28 = t93 * t56 + t58 * t70;
t88 = t28 * t49;
t29 = -t56 * t70 + t93 * t58;
t87 = t29 * t49;
t86 = t42 * t45;
t85 = t44 * t45;
t84 = t45 * t49;
t52 = cos(pkin(12));
t83 = t45 * t52;
t82 = t45 * t58;
t81 = t51 * t56;
t80 = t51 * t57;
t55 = -pkin(9) - qJ(3);
t78 = t55 * t56;
t77 = -t26 * t41 - t27 * t55;
t76 = -t28 * t41 - t29 * t55;
t75 = t93 * pkin(1) + pkin(8) * t80;
t50 = sin(pkin(11));
t73 = t50 * t80;
t71 = -t57 * pkin(1) + pkin(8) * t72;
t69 = t50 * t72;
t68 = pkin(3) * t73 - t28 * t55 + t29 * t41 + t75;
t10 = t27 * t43 + t45 * t72;
t14 = t29 * t43 - t45 * t80;
t67 = -g(1) * t10 + g(2) * t14;
t9 = g(1) * t26 - g(2) * t28;
t65 = pkin(4) * t45 + qJ(5) * t43;
t40 = t52 * pkin(5) + pkin(4);
t54 = -pkin(10) - qJ(5);
t64 = t40 * t45 - t43 * t54;
t63 = g(1) * t93 + g(2) * t57;
t62 = pkin(3) * t69 + t26 * t55 - t27 * t41 + t71;
t20 = t43 * t81 - t74 * t45;
t61 = g(1) * t14 + g(2) * t10 + g(3) * t20;
t15 = t29 * t45 + t43 * t80;
t21 = t74 * t43 + t45 * t81;
t60 = g(1) * t15 + g(2) * t11 + g(3) * t21;
t7 = -g(1) * t28 - g(2) * t26 + g(3) * t79;
t59 = g(1) * t29 + g(2) * t27 + g(3) * t81;
t6 = t7 * t43;
t5 = t15 * t44 + t28 * t42;
t4 = -t15 * t42 + t28 * t44;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t93, t63, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t27 - g(2) * t29, -t9, -t63 * t51, -g(1) * t71 - g(2) * t75, 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t53 + t69) - g(2) * (t29 * t53 + t73) -g(1) * (t27 * t50 + t53 * t72) - g(2) * (-t29 * t50 + t53 * t80) t9, -g(1) * (-t27 * pkin(2) - t26 * qJ(3) + t71) - g(2) * (t29 * pkin(2) + t28 * qJ(3) + t75) 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t15, t67, t9, -g(1) * t62 - g(2) * t68, 0, 0, 0, 0, 0, 0, -g(1) * (-t11 * t52 - t90) - g(2) * (t15 * t52 + t88) -g(1) * (t11 * t49 - t26 * t52) - g(2) * (-t15 * t49 + t28 * t52) -t67, -g(1) * (-pkin(4) * t11 - qJ(5) * t10 + t62) - g(2) * (t15 * pkin(4) + t14 * qJ(5) + t68) 0, 0, 0, 0, 0, 0, g(1) * t96 - g(2) * t5, -g(1) * t97 - g(2) * t4, -t67, -g(1) * (-pkin(5) * t90 + t10 * t54 - t11 * t40 + t62) - g(2) * (pkin(5) * t88 - t14 * t54 + t15 * t40 + t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t53, t7 * t50, -t59, -g(1) * (-t28 * pkin(2) + t29 * qJ(3)) - g(2) * (-t26 * pkin(2) + t27 * qJ(3)) - (pkin(2) * t58 + qJ(3) * t56) * t94, 0, 0, 0, 0, 0, 0, -t7 * t45, t6, -t59, -g(1) * t76 - g(2) * t77 - g(3) * (-t51 * t78 + t30) 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t83 + t87) - g(2) * (-t26 * t83 + t89) - (t49 * t56 + t52 * t82) * t94, -g(1) * (t28 * t84 + t29 * t52) - g(2) * (t26 * t84 + t27 * t52) - (-t49 * t82 + t52 * t56) * t94, -t6, -g(1) * (-t65 * t28 + t76) - g(2) * (-t65 * t26 + t77) - t95 - (t65 * t58 - t78) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t85 + t29 * t42) - g(2) * (-t26 * t85 + t27 * t42) - (t42 * t56 + t44 * t82) * t94, -g(1) * (t28 * t86 + t29 * t44) - g(2) * (t26 * t86 + t27 * t44) - (-t42 * t82 + t44 * t56) * t94, -t6, -g(1) * (pkin(5) * t87 - t64 * t28 + t76) - g(2) * (pkin(5) * t89 - t64 * t26 + t77) - t95 - (t64 * t58 + (pkin(5) * t49 - t55) * t56) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t52, -t61 * t49, -t60, -g(1) * (-t14 * pkin(4) + t15 * qJ(5)) - g(2) * (-t10 * pkin(4) + t11 * qJ(5)) - g(3) * (-t20 * pkin(4) + t21 * qJ(5)) 0, 0, 0, 0, 0, 0, t61 * t44, -t61 * t42, -t60, -g(1) * (-t14 * t40 - t15 * t54) - g(2) * (-t10 * t40 - t11 * t54) - g(3) * (-t20 * t40 - t21 * t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t97 - g(3) * (-t21 * t42 - t44 * t79) g(1) * t5 + g(2) * t96 - g(3) * (-t21 * t44 + t42 * t79) 0, 0;];
taug_reg  = t1;
