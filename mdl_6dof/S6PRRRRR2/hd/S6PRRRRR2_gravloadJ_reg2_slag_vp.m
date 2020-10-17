% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:44:31
% EndTime: 2019-05-05 10:44:34
% DurationCPUTime: 0.82s
% Computational Cost: add. (742->146), mult. (1190->218), div. (0->0), fcn. (1437->14), ass. (0->79)
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t79 = cos(pkin(6));
t50 = sin(pkin(6));
t53 = sin(qJ(2));
t90 = t50 * t53;
t103 = -t52 * t90 + t55 * t79;
t56 = cos(qJ(2));
t49 = sin(pkin(12));
t72 = t49 * t79;
t78 = cos(pkin(12));
t34 = -t53 * t72 + t56 * t78;
t89 = t50 * t55;
t102 = -t34 * t52 + t49 * t89;
t42 = pkin(3) * t55 + pkin(2);
t88 = t50 * t56;
t35 = t42 * t88;
t101 = g(3) * t35;
t100 = g(3) * t50;
t65 = t79 * t78;
t32 = t49 * t56 + t53 * t65;
t51 = sin(qJ(5));
t99 = t32 * t51;
t98 = t34 * t51;
t47 = qJ(5) + qJ(6);
t43 = sin(t47);
t48 = qJ(3) + qJ(4);
t46 = cos(t48);
t96 = t43 * t46;
t45 = cos(t47);
t95 = t45 * t46;
t94 = t46 * t51;
t54 = cos(qJ(5));
t93 = t46 * t54;
t92 = t46 * t56;
t91 = t49 * t50;
t87 = t51 * t56;
t58 = -pkin(9) - pkin(8);
t86 = t53 * t58;
t85 = t54 * t56;
t44 = sin(t48);
t71 = t50 * t78;
t18 = -t32 * t44 - t46 * t71;
t19 = t32 * t46 - t44 * t71;
t41 = pkin(5) * t54 + pkin(4);
t57 = -pkin(11) - pkin(10);
t84 = t18 * t41 - t19 * t57;
t20 = -t34 * t44 + t46 * t91;
t21 = t34 * t46 + t44 * t91;
t83 = t20 * t41 - t21 * t57;
t27 = -t44 * t90 + t46 * t79;
t28 = t44 * t79 + t46 * t90;
t82 = t27 * t41 - t28 * t57;
t31 = t49 * t53 - t56 * t65;
t81 = -t31 * t42 - t32 * t58;
t33 = t53 * t78 + t56 * t72;
t80 = -t33 * t42 - t34 * t58;
t75 = pkin(4) * t18 + t19 * pkin(10);
t74 = pkin(4) * t20 + pkin(10) * t21;
t73 = pkin(4) * t27 + pkin(10) * t28;
t69 = t102 * pkin(3);
t68 = pkin(4) * t46 + pkin(10) * t44;
t67 = t41 * t46 - t44 * t57;
t66 = t103 * pkin(3);
t64 = g(1) * t20 + g(2) * t18 + g(3) * t27;
t9 = g(1) * t21 + g(2) * t19 + g(3) * t28;
t63 = -t32 * t52 - t55 * t71;
t62 = -g(1) * t33 - g(2) * t31 + g(3) * t88;
t61 = g(1) * t34 + g(2) * t32 + g(3) * t90;
t60 = t63 * pkin(3);
t59 = -g(1) * (-t21 * t51 + t33 * t54) - g(2) * (-t19 * t51 + t31 * t54) - g(3) * (-t28 * t51 - t50 * t85);
t10 = t62 * t44;
t6 = t64 * t54;
t5 = t64 * t51;
t4 = t64 * t45;
t3 = t64 * t43;
t2 = -g(1) * (-t21 * t45 - t33 * t43) - g(2) * (-t19 * t45 - t31 * t43) - g(3) * (-t28 * t45 + t43 * t88);
t1 = -g(1) * (-t21 * t43 + t33 * t45) - g(2) * (-t19 * t43 + t31 * t45) - g(3) * (-t28 * t43 - t45 * t88);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t61, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t55, t62 * t52, -t61, -g(1) * (-pkin(2) * t33 + pkin(8) * t34) - g(2) * (-pkin(2) * t31 + pkin(8) * t32) - (pkin(2) * t56 + pkin(8) * t53) * t100, 0, 0, 0, 0, 0, 0, -t62 * t46, t10, -t61, -g(1) * t80 - g(2) * t81 - g(3) * (-t50 * t86 + t35) 0, 0, 0, 0, 0, 0, -g(1) * (-t33 * t93 + t98) - g(2) * (-t31 * t93 + t99) - (t46 * t85 + t51 * t53) * t100, -g(1) * (t33 * t94 + t34 * t54) - g(2) * (t31 * t94 + t32 * t54) - (-t46 * t87 + t53 * t54) * t100, -t10, -g(1) * (-t33 * t68 + t80) - g(2) * (-t31 * t68 + t81) - t101 - (t56 * t68 - t86) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (-t33 * t95 + t34 * t43) - g(2) * (-t31 * t95 + t32 * t43) - (t43 * t53 + t45 * t92) * t100, -g(1) * (t33 * t96 + t34 * t45) - g(2) * (t31 * t96 + t32 * t45) - (-t43 * t92 + t45 * t53) * t100, -t10, -g(1) * (pkin(5) * t98 - t33 * t67 + t80) - g(2) * (pkin(5) * t99 - t31 * t67 + t81) - t101 - (t67 * t56 + (pkin(5) * t51 - t58) * t53) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t102 - g(2) * t63 - g(3) * t103, -g(1) * (-t34 * t55 - t52 * t91) - g(2) * (-t32 * t55 + t52 * t71) - g(3) * (-t52 * t79 - t53 * t89) 0, 0, 0, 0, 0, 0, 0, 0, -t64, t9, 0, -g(1) * t69 - g(2) * t60 - g(3) * t66, 0, 0, 0, 0, 0, 0, -t6, t5, -t9, -g(1) * (t69 + t74) - g(2) * (t60 + t75) - g(3) * (t66 + t73) 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(1) * (t69 + t83) - g(2) * (t60 + t84) - g(3) * (t66 + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t9, -g(1) * t74 - g(2) * t75 - g(3) * t73, 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(1) * t83 - g(2) * t84 - g(3) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -g(1) * (-t21 * t54 - t33 * t51) - g(2) * (-t19 * t54 - t31 * t51) - g(3) * (-t28 * t54 + t50 * t87) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t59 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t7;
