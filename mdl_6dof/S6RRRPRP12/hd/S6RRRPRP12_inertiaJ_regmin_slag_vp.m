% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:33:46
% EndTime: 2019-05-07 09:33:52
% DurationCPUTime: 1.30s
% Computational Cost: add. (1418->178), mult. (3112->322), div. (0->0), fcn. (3389->8), ass. (0->101)
t65 = sin(pkin(6));
t69 = sin(qJ(2));
t106 = t65 * t69;
t66 = cos(pkin(6));
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t37 = t71 * t106 + t66 * t68;
t117 = -0.2e1 * t37;
t116 = -0.2e1 * t68;
t115 = 0.2e1 * t71;
t114 = 2 * qJ(4);
t73 = -pkin(3) - pkin(10);
t113 = pkin(1) * t69;
t112 = t37 * pkin(5);
t111 = t68 * pkin(5);
t72 = cos(qJ(2));
t110 = t72 * pkin(1);
t51 = pkin(8) * t106;
t29 = t51 + (-pkin(2) - t110) * t66;
t79 = t68 * t106 - t66 * t71;
t13 = t79 * pkin(3) - t37 * qJ(4) + t29;
t11 = t79 * pkin(10) + t13;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t105 = t65 * t72;
t90 = pkin(8) * t105;
t30 = t90 + (pkin(9) + t113) * t66;
t31 = (-pkin(2) * t72 - pkin(9) * t69 - pkin(1)) * t65;
t16 = -t68 * t30 + t71 * t31;
t52 = pkin(3) * t105;
t15 = -t16 + t52;
t8 = t37 * pkin(4) + pkin(10) * t105 + t15;
t4 = t70 * t11 + t67 * t8;
t23 = t70 * t105 - t79 * t67;
t109 = t23 * t70;
t108 = t37 * t67;
t33 = t37 * t68;
t34 = t37 * t70;
t59 = t65 ^ 2;
t107 = t59 * t72;
t22 = t67 * t105 + t79 * t70;
t104 = t67 * t22;
t103 = t67 * t68;
t102 = t67 * t71;
t101 = t67 * t73;
t100 = t68 * t71;
t99 = t68 * t73;
t98 = t70 * t67;
t97 = t70 * t71;
t55 = t70 * t73;
t17 = t71 * t30 + t68 * t31;
t83 = -t68 * qJ(4) - pkin(2);
t40 = t73 * t71 + t83;
t56 = t68 * pkin(9);
t47 = t68 * pkin(4) + t56;
t21 = t70 * t40 + t67 * t47;
t57 = t71 * pkin(9);
t48 = t71 * pkin(4) + t57;
t60 = t67 ^ 2;
t62 = t70 ^ 2;
t50 = t60 + t62;
t61 = t68 ^ 2;
t63 = t71 ^ 2;
t96 = t61 + t63;
t95 = qJ(4) * t71;
t94 = t37 * qJ(6);
t93 = t68 * qJ(6);
t92 = 0.2e1 * t65 * t66;
t91 = -0.2e1 * t100;
t89 = t68 * t105;
t88 = t71 * t105;
t87 = t37 * t101;
t86 = t37 * t55;
t85 = qJ(4) * t105;
t84 = t67 * t11 - t70 * t8;
t82 = t67 * t40 - t70 * t47;
t81 = pkin(9) * t89;
t80 = pkin(9) * t88;
t1 = t94 + t4;
t2 = t84 - t112;
t78 = t1 * t67 - t2 * t70;
t77 = -pkin(3) * t68 + t95;
t46 = pkin(5) * t70 + t67 * qJ(6);
t76 = t67 * pkin(5) - t70 * qJ(6);
t14 = t85 - t17;
t75 = -t14 * t71 + t15 * t68;
t18 = t93 + t21;
t19 = t82 - t111;
t9 = t18 * t67 - t19 * t70;
t12 = -t79 * pkin(4) - t14;
t54 = t70 * t68;
t49 = t68 * t55;
t45 = -t71 * pkin(3) + t83;
t44 = qJ(4) + t76;
t43 = t50 * t73;
t39 = t66 * t113 + t90;
t38 = t66 * t110 - t51;
t36 = t37 ^ 2;
t25 = t46 * t71 + t48;
t5 = -t22 * pkin(5) + t23 * qJ(6) + t12;
t3 = [1, 0, 0, t59 * t69 ^ 2, 0.2e1 * t69 * t107, t69 * t92, t72 * t92, t66 ^ 2, 0.2e1 * pkin(1) * t107 + 0.2e1 * t38 * t66, -0.2e1 * t59 * t113 - 0.2e1 * t39 * t66, t36, t79 * t117, t105 * t117, 0.2e1 * t79 * t105, t59 * t72 ^ 2, -0.2e1 * t16 * t105 + 0.2e1 * t29 * t79, 0.2e1 * t17 * t105 + 0.2e1 * t29 * t37, 0.2e1 * t14 * t79 + 0.2e1 * t15 * t37, -0.2e1 * t15 * t105 - 0.2e1 * t13 * t79, 0.2e1 * t14 * t105 - 0.2e1 * t13 * t37, t13 ^ 2 + t14 ^ 2 + t15 ^ 2, t23 ^ 2, -0.2e1 * t23 * t22, t23 * t117, 0.2e1 * t22 * t37, t36, -0.2e1 * t12 * t22 - 0.2e1 * t37 * t84, -0.2e1 * t12 * t23 - 0.2e1 * t4 * t37, -0.2e1 * t2 * t37 - 0.2e1 * t5 * t22, 0.2e1 * t1 * t22 - 0.2e1 * t2 * t23, 0.2e1 * t1 * t37 + 0.2e1 * t5 * t23, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t106, t105, t66, t38, -t39, t33, t37 * t71 - t68 * t79, -t89, -t88, 0, -pkin(2) * t79 - t29 * t71 + t81, -pkin(2) * t37 + t29 * t68 + t80 (-t71 * t79 + t33) * pkin(9) + t75, t13 * t71 - t45 * t79 - t81, -t13 * t68 - t45 * t37 - t80, t75 * pkin(9) + t13 * t45, t23 * t102 (-t104 + t109) * t71, -t37 * t102 - t23 * t68, t22 * t68 - t37 * t97, t33, t12 * t97 - t48 * t22 - t37 * t82 - t68 * t84, -t102 * t12 - t21 * t37 - t48 * t23 - t4 * t68, -t19 * t37 - t2 * t68 - t25 * t22 + t5 * t97, t18 * t22 - t19 * t23 + (-t1 * t70 - t2 * t67) * t71, t1 * t68 + t102 * t5 + t18 * t37 + t25 * t23, t1 * t18 + t2 * t19 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t61, 0.2e1 * t100, 0, 0, 0, pkin(2) * t115, pkin(2) * t116, 0.2e1 * t96 * pkin(9), t45 * t115, t45 * t116, t96 * pkin(9) ^ 2 + t45 ^ 2, t60 * t63, 0.2e1 * t63 * t98, t67 * t91, t70 * t91, t61, 0.2e1 * t48 * t97 - 0.2e1 * t68 * t82, -0.2e1 * t102 * t48 - 0.2e1 * t21 * t68, -0.2e1 * t19 * t68 + 0.2e1 * t25 * t97 (-t18 * t70 - t19 * t67) * t115, 0.2e1 * t102 * t25 + 0.2e1 * t18 * t68, t18 ^ 2 + t19 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t79, -t105, t16, -t17, -t37 * pkin(3) - qJ(4) * t79, -t16 + 0.2e1 * t52, -0.2e1 * t85 + t17, -t15 * pkin(3) - t14 * qJ(4), -t109, t70 * t22 + t23 * t67, t34, -t108, 0, -qJ(4) * t22 + t12 * t67 + t86, -qJ(4) * t23 + t12 * t70 - t87, -t44 * t22 + t5 * t67 + t86 (t23 * t73 + t2) * t70 + (t22 * t73 - t1) * t67, t44 * t23 - t5 * t70 + t87, t5 * t44 + t73 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t71, 0, -t56, -t57, t77, t56, t57, t77 * pkin(9), -t67 * t97 (t60 - t62) * t71, t54, -t103, 0, t48 * t67 + t70 * t95 + t49, t48 * t70 + (-t95 - t99) * t67, t25 * t67 + t44 * t97 + t49, -t9, -t25 * t70 + (t44 * t71 + t99) * t67, t25 * t44 + t73 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t114, pkin(3) ^ 2 + (qJ(4) ^ 2) t62, -0.2e1 * t98, 0, 0, 0, t67 * t114, t70 * t114, 0.2e1 * t44 * t67, -0.2e1 * t43, -0.2e1 * t44 * t70, t50 * t73 ^ 2 + t44 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t105, 0, t15, 0, 0, 0, 0, 0, t34, -t108, t34, t104 + t109, t108, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, t56, 0, 0, 0, 0, 0, t54, -t103, t54, 0, t103, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, t37, -t84, -t4, -t84 + 0.2e1 * t112, t23 * pkin(5) + t22 * qJ(6), 0.2e1 * t94 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t97, t68, -t82, -t21, -t82 + 0.2e1 * t111, t76 * t71, 0.2e1 * t93 + t21, -t19 * pkin(5) + t18 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t67, 0, t55, -t101, t55, -t46, t101, t46 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t67, t70, 0, t67, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t23, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t102, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
