% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:16:23
% EndTime: 2019-05-05 07:16:25
% DurationCPUTime: 0.68s
% Computational Cost: add. (666->136), mult. (1065->196), div. (0->0), fcn. (1278->14), ass. (0->73)
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t77 = sin(pkin(11));
t79 = cos(pkin(6));
t63 = t79 * t77;
t78 = cos(pkin(11));
t33 = -t55 * t63 + t78 * t57;
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t51 = sin(pkin(6));
t74 = t51 * t77;
t99 = -t33 * t54 + t56 * t74;
t87 = t51 * t55;
t98 = -t54 * t87 + t79 * t56;
t43 = t56 * pkin(3) + pkin(2);
t86 = t51 * t57;
t34 = t43 * t86;
t97 = g(3) * t34;
t96 = g(3) * t51;
t64 = t79 * t78;
t31 = t55 * t64 + t77 * t57;
t50 = sin(pkin(12));
t95 = t31 * t50;
t94 = t33 * t50;
t48 = pkin(12) + qJ(6);
t44 = sin(t48);
t49 = qJ(3) + qJ(4);
t47 = cos(t49);
t92 = t44 * t47;
t45 = cos(t48);
t91 = t45 * t47;
t90 = t47 * t50;
t52 = cos(pkin(12));
t89 = t47 * t52;
t88 = t47 * t57;
t58 = -pkin(9) - pkin(8);
t85 = t55 * t58;
t46 = sin(t49);
t75 = t51 * t78;
t17 = t31 * t46 + t47 * t75;
t18 = t31 * t47 - t46 * t75;
t42 = t52 * pkin(5) + pkin(4);
t53 = -pkin(10) - qJ(5);
t84 = -t17 * t42 - t18 * t53;
t19 = t33 * t46 - t47 * t74;
t20 = t33 * t47 + t46 * t74;
t83 = -t19 * t42 - t20 * t53;
t26 = t46 * t87 - t79 * t47;
t27 = t79 * t46 + t47 * t87;
t82 = -t26 * t42 - t27 * t53;
t30 = t77 * t55 - t57 * t64;
t81 = -t30 * t43 - t31 * t58;
t32 = t78 * t55 + t57 * t63;
t80 = -t32 * t43 - t33 * t58;
t72 = -t17 * pkin(4) + t18 * qJ(5);
t71 = -t19 * pkin(4) + t20 * qJ(5);
t70 = -t26 * pkin(4) + t27 * qJ(5);
t69 = t99 * pkin(3);
t67 = pkin(4) * t47 + qJ(5) * t46;
t66 = t42 * t47 - t46 * t53;
t65 = t98 * pkin(3);
t6 = g(1) * t19 + g(2) * t17 + g(3) * t26;
t8 = g(1) * t20 + g(2) * t18 + g(3) * t27;
t62 = -t31 * t54 - t56 * t75;
t61 = -g(1) * t32 - g(2) * t30 + g(3) * t86;
t60 = g(1) * t33 + g(2) * t31 + g(3) * t87;
t59 = t62 * pkin(3);
t9 = t61 * t46;
t4 = t6 * t52;
t3 = t6 * t50;
t2 = t6 * t45;
t1 = t6 * t44;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t60, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t56, t61 * t54, -t60, -g(1) * (-t32 * pkin(2) + t33 * pkin(8)) - g(2) * (-t30 * pkin(2) + t31 * pkin(8)) - (pkin(2) * t57 + pkin(8) * t55) * t96, 0, 0, 0, 0, 0, 0, -t61 * t47, t9, -t60, -g(1) * t80 - g(2) * t81 - g(3) * (-t51 * t85 + t34) 0, 0, 0, 0, 0, 0, -g(1) * (-t32 * t89 + t94) - g(2) * (-t30 * t89 + t95) - (t50 * t55 + t52 * t88) * t96, -g(1) * (t32 * t90 + t33 * t52) - g(2) * (t30 * t90 + t31 * t52) - (-t50 * t88 + t52 * t55) * t96, -t9, -g(1) * (-t67 * t32 + t80) - g(2) * (-t67 * t30 + t81) - t97 - (t67 * t57 - t85) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (-t32 * t91 + t33 * t44) - g(2) * (-t30 * t91 + t31 * t44) - (t44 * t55 + t45 * t88) * t96, -g(1) * (t32 * t92 + t33 * t45) - g(2) * (t30 * t92 + t31 * t45) - (-t44 * t88 + t45 * t55) * t96, -t9, -g(1) * (pkin(5) * t94 - t66 * t32 + t80) - g(2) * (pkin(5) * t95 - t66 * t30 + t81) - t97 - (t66 * t57 + (pkin(5) * t50 - t58) * t55) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t62 - g(3) * t98, -g(1) * (-t33 * t56 - t54 * t74) - g(2) * (-t31 * t56 + t54 * t75) - g(3) * (-t79 * t54 - t56 * t87) 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, -g(1) * t69 - g(2) * t59 - g(3) * t65, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t69 + t71) - g(2) * (t59 + t72) - g(3) * (t65 + t70) 0, 0, 0, 0, 0, 0, t2, -t1, -t8, -g(1) * (t69 + t83) - g(2) * (t59 + t84) - g(3) * (t65 + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t71 - g(2) * t72 - g(3) * t70, 0, 0, 0, 0, 0, 0, t2, -t1, -t8, -g(1) * t83 - g(2) * t84 - g(3) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t20 * t44 + t32 * t45) - g(2) * (-t18 * t44 + t30 * t45) - g(3) * (-t27 * t44 - t45 * t86) -g(1) * (-t20 * t45 - t32 * t44) - g(2) * (-t18 * t45 - t30 * t44) - g(3) * (-t27 * t45 + t44 * t86) 0, 0;];
taug_reg  = t5;
