% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_inertiaJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:55:40
% EndTime: 2019-05-05 01:55:43
% DurationCPUTime: 1.03s
% Computational Cost: add. (1511->149), mult. (4353->306), div. (0->0), fcn. (5313->16), ass. (0->101)
t68 = cos(pkin(7));
t108 = pkin(2) * t68;
t62 = sin(pkin(14));
t66 = cos(pkin(14));
t64 = sin(pkin(7));
t81 = qJ(3) * t64;
t47 = t108 * t62 + t66 * t81;
t63 = sin(pkin(8));
t67 = cos(pkin(8));
t93 = t66 * t67;
t79 = t64 * t93;
t29 = (t63 * t68 + t79) * pkin(10) + t47;
t72 = sin(qJ(4));
t76 = cos(qJ(4));
t100 = t62 * t64;
t55 = t66 * t108;
t34 = t68 * pkin(3) + t55 + (-pkin(10) * t67 - qJ(3)) * t100;
t38 = (-pkin(10) * t62 * t63 - pkin(3) * t66 - pkin(2)) * t64;
t78 = t34 * t67 + t38 * t63;
t15 = -t29 * t72 + t76 * t78;
t99 = t63 * t72;
t33 = t68 * t99 + (t62 * t76 + t72 * t93) * t64;
t97 = t64 * t66;
t43 = t63 * t97 - t67 * t68;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t24 = t33 * t71 + t43 * t75;
t114 = -0.2e1 * t24;
t98 = t63 * t76;
t32 = t100 * t72 - t68 * t98 - t76 * t79;
t113 = 0.2e1 * t32;
t112 = -0.2e1 * t33;
t111 = -0.2e1 * t71;
t110 = 0.2e1 * t75;
t58 = t64 ^ 2;
t109 = pkin(2) * t58;
t74 = cos(qJ(6));
t107 = pkin(5) * t74;
t70 = sin(qJ(6));
t106 = pkin(11) * t70;
t21 = -t34 * t63 + t38 * t67;
t12 = pkin(4) * t32 - pkin(11) * t33 + t21;
t16 = t76 * t29 + t72 * t78;
t14 = -t43 * pkin(11) + t16;
t7 = t12 * t75 - t14 * t71;
t3 = -pkin(5) * t32 - t7;
t105 = t3 * t70;
t104 = t3 * t74;
t103 = t64 * pkin(2);
t25 = t33 * t75 - t43 * t71;
t20 = t25 * t74 + t32 * t70;
t102 = t20 * t70;
t65 = sin(pkin(6));
t73 = sin(qJ(2));
t77 = cos(qJ(2));
t92 = t68 * t77;
t69 = cos(pkin(6));
t96 = t64 * t69;
t30 = t66 * t96 + (-t62 * t73 + t66 * t92) * t65;
t101 = t30 * t67;
t95 = t65 * t73;
t94 = t65 * t77;
t91 = t70 * t24;
t90 = t70 * t71;
t89 = t70 * t74;
t88 = t70 * t75;
t87 = t71 * t24;
t86 = t71 * t32;
t85 = t74 * t24;
t84 = t74 * t71;
t83 = t74 * t75;
t82 = t75 * t32;
t80 = t71 * t110;
t8 = t12 * t71 + t14 * t75;
t13 = t43 * pkin(4) - t15;
t61 = t74 ^ 2;
t60 = t71 ^ 2;
t59 = t70 ^ 2;
t51 = -pkin(5) * t75 - pkin(12) * t71 - pkin(4);
t49 = t67 * t71 + t75 * t99;
t48 = -t67 * t75 + t71 * t99;
t46 = -t62 * t81 + t55;
t45 = -t64 * t94 + t68 * t69;
t41 = pkin(11) * t83 + t51 * t70;
t40 = -pkin(11) * t88 + t51 * t74;
t36 = t49 * t74 - t70 * t98;
t35 = -t49 * t70 - t74 * t98;
t31 = t66 * t95 + (t65 * t92 + t96) * t62;
t23 = -t30 * t63 + t45 * t67;
t19 = t25 * t70 - t32 * t74;
t18 = t31 * t76 + (t45 * t63 + t101) * t72;
t17 = -t101 * t76 + t31 * t72 - t45 * t98;
t11 = t18 * t75 + t23 * t71;
t10 = t18 * t71 - t23 * t75;
t9 = t24 * pkin(5) - t25 * pkin(12) + t13;
t6 = t11 * t74 + t17 * t70;
t5 = -t11 * t70 + t17 * t74;
t4 = pkin(12) * t32 + t8;
t2 = t4 * t74 + t70 * t9;
t1 = -t4 * t70 + t74 * t9;
t22 = [1, 0, 0, 0, 0, 0, 0, t30 ^ 2 + t31 ^ 2 + t45 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t94, -t95, t30 * t68 - t45 * t97, t100 * t45 - t31 * t68 (-t30 * t62 + t31 * t66) * t64, -t103 * t45 + t30 * t46 + t31 * t47, 0, 0, 0, 0, 0, t17 * t43 + t23 * t32, t18 * t43 + t23 * t33, 0, 0, 0, 0, 0, -t10 * t32 + t17 * t24, -t11 * t32 + t17 * t25, 0, 0, 0, 0, 0, t10 * t19 + t24 * t5, t10 * t20 - t24 * t6; 0, 1, 0, 0, 0.2e1 * t109 * t66 + 0.2e1 * t46 * t68, -0.2e1 * t109 * t62 - 0.2e1 * t47 * t68, 0.2e1 * (-t46 * t62 + t47 * t66) * t64, pkin(2) ^ 2 * t58 + t46 ^ 2 + t47 ^ 2, t33 ^ 2, t32 * t112, t43 * t112, t43 * t113, t43 ^ 2, -0.2e1 * t15 * t43 + 0.2e1 * t21 * t32, 0.2e1 * t16 * t43 + 0.2e1 * t21 * t33, t25 ^ 2, t25 * t114, t25 * t113, t32 * t114, t32 ^ 2, 0.2e1 * t13 * t24 + 0.2e1 * t32 * t7, 0.2e1 * t13 * t25 - 0.2e1 * t32 * t8, t20 ^ 2, -0.2e1 * t20 * t19, 0.2e1 * t20 * t24, t19 * t114, t24 ^ 2, 0.2e1 * t1 * t24 + 0.2e1 * t19 * t3, -0.2e1 * t2 * t24 + 0.2e1 * t20 * t3; 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t97, t100, 0, -t103, 0, 0, 0, 0, 0, t32 * t67 - t43 * t98, t33 * t67 + t43 * t99, 0, 0, 0, 0, 0, -t24 * t98 - t32 * t48, -t25 * t98 - t32 * t49, 0, 0, 0, 0, 0, t19 * t48 + t24 * t35, t20 * t48 - t24 * t36; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, 0, 0, 0, 0, -t17 * t75, t17 * t71, 0, 0, 0, 0, 0, t10 * t90 - t5 * t75, t10 * t84 + t6 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, -t43, t15, -t16, t25 * t71, t25 * t75 - t87, t86, t82, 0, -pkin(4) * t24 - pkin(11) * t86 - t13 * t75, -pkin(4) * t25 - pkin(11) * t82 + t13 * t71, t20 * t84 (-t19 * t74 - t102) * t71, -t20 * t75 + t24 * t84, t19 * t75 - t70 * t87, -t24 * t75, -t1 * t75 + t40 * t24 + (pkin(11) * t19 + t105) * t71, t2 * t75 - t41 * t24 + (pkin(11) * t20 + t104) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t99, 0, 0, 0, 0, 0, t75 * t98, -t71 * t98, 0, 0, 0, 0, 0, -t35 * t75 + t48 * t90, t36 * t75 + t48 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t60, t80, 0, 0, 0, pkin(4) * t110, pkin(4) * t111, t61 * t60, -0.2e1 * t60 * t89, t83 * t111, t70 * t80, t75 ^ 2, 0.2e1 * t106 * t60 - 0.2e1 * t40 * t75, 0.2e1 * pkin(11) * t60 * t74 + 0.2e1 * t41 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t74, t10 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, t32, t7, -t8, t102, -t19 * t70 + t20 * t74, t91, t85, 0, -pkin(5) * t19 - pkin(12) * t91 - t104, -pkin(5) * t20 - pkin(12) * t85 + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t49, 0, 0, 0, 0, 0, -t48 * t74, t48 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t75, 0, -t71 * pkin(11), -t75 * pkin(11), t70 * t84 (-t59 + t61) * t71, -t88, -t83, 0, -pkin(11) * t84 + (-pkin(5) * t71 + pkin(12) * t75) * t70, pkin(12) * t83 + (t106 - t107) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, 0.2e1 * t89, 0, 0, 0, 0.2e1 * t107, -0.2e1 * pkin(5) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t24, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t90, -t75, t40, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t74, 0, -t70 * pkin(12), -t74 * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t22;
