% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:44:31
% EndTime: 2019-05-05 10:44:34
% DurationCPUTime: 0.85s
% Computational Cost: add. (671->118), mult. (1482->190), div. (0->0), fcn. (1887->12), ass. (0->104)
t72 = sin(qJ(4));
t73 = sin(qJ(3));
t95 = cos(qJ(4));
t96 = cos(qJ(3));
t49 = t72 * t73 - t95 * t96;
t109 = -0.2e1 * t49;
t108 = 0.2e1 * t49;
t81 = t95 * pkin(3);
t62 = -t81 - pkin(4);
t76 = cos(qJ(5));
t99 = t76 * pkin(5);
t53 = t62 - t99;
t107 = 0.2e1 * t53;
t63 = -pkin(4) - t99;
t106 = 0.2e1 * t63;
t64 = -t96 * pkin(3) - pkin(2);
t105 = 0.2e1 * t64;
t104 = t49 * pkin(5);
t70 = sin(qJ(6));
t103 = t70 * pkin(5);
t102 = t72 * pkin(3);
t75 = cos(qJ(6));
t101 = t75 * pkin(5);
t51 = t72 * t96 + t95 * t73;
t27 = t49 * pkin(4) - t51 * pkin(10) + t64;
t71 = sin(qJ(5));
t55 = (-pkin(8) - pkin(9)) * t73;
t82 = t96 * pkin(8);
t57 = t96 * pkin(9) + t82;
t36 = t72 * t55 + t95 * t57;
t88 = t76 * t36;
t8 = t88 + (-pkin(11) * t51 + t27) * t71;
t100 = t75 * t8;
t98 = t76 * pkin(10);
t97 = pkin(4) - t62;
t69 = cos(pkin(6));
t68 = sin(pkin(6));
t74 = sin(qJ(2));
t92 = t68 * t74;
t39 = t69 * t96 - t73 * t92;
t80 = t68 * t96;
t40 = t69 * t73 + t74 * t80;
t20 = -t95 * t39 + t72 * t40;
t94 = t20 * t76;
t34 = -t95 * t55 + t72 * t57;
t93 = t34 * t76;
t77 = cos(qJ(2));
t91 = t68 * t77;
t90 = t71 * t51;
t89 = t71 * t76;
t87 = t76 * t51;
t61 = pkin(10) + t102;
t86 = t76 * t61;
t85 = t53 + t63;
t84 = 0.2e1 * t96;
t83 = t51 * t109;
t9 = t76 * t27 - t71 * t36;
t7 = -pkin(11) * t87 + t104 + t9;
t1 = t75 * t7 - t70 * t8;
t79 = -pkin(4) * t51 - pkin(10) * t49;
t78 = -t49 * t61 + t51 * t62;
t50 = t70 * t76 + t75 * t71;
t48 = t70 * t71 - t75 * t76;
t67 = t76 ^ 2;
t66 = t71 ^ 2;
t65 = t76 * pkin(11);
t58 = 0.2e1 * t89;
t56 = t65 + t98;
t54 = (-pkin(10) - pkin(11)) * t71;
t47 = t51 ^ 2;
t46 = t50 ^ 2;
t45 = t49 ^ 2;
t44 = t65 + t86;
t43 = (-pkin(11) - t61) * t71;
t42 = t76 * t49;
t41 = t71 * t49;
t38 = t71 * t87;
t35 = t70 * t54 + t75 * t56;
t33 = t75 * t54 - t70 * t56;
t32 = t50 * t49;
t31 = t48 * t49;
t30 = -0.2e1 * t50 * t48;
t29 = t34 * t71;
t28 = (-t66 + t67) * t51;
t26 = t70 * t43 + t75 * t44;
t25 = t75 * t43 - t70 * t44;
t23 = t48 * t51;
t22 = t50 * t51;
t21 = t72 * t39 + t95 * t40;
t19 = pkin(5) * t90 + t34;
t18 = t20 * t71;
t17 = t23 * t50;
t16 = t76 * t21 - t71 * t91;
t15 = -t71 * t21 - t76 * t91;
t14 = t20 * t50;
t13 = t20 * t48;
t12 = t19 * t50;
t11 = t19 * t48;
t10 = t71 * t27 + t88;
t5 = -t50 * t22 + t23 * t48;
t4 = t70 * t15 + t75 * t16;
t3 = t75 * t15 - t70 * t16;
t2 = t70 * t7 + t100;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t91, -t92, 0, 0, 0, 0, 0, t77 * t80, -t73 * t91, 0, 0, 0, 0, 0, -t49 * t91, -t51 * t91, 0, 0, 0, 0, 0, t15 * t49 + t20 * t90, -t16 * t49 + t20 * t87, 0, 0, 0, 0, 0, t20 * t22 + t3 * t49, -t20 * t23 - t4 * t49; 0, 1, 0, 0, t73 ^ 2, t73 * t84, 0, 0, 0, pkin(2) * t84, -0.2e1 * pkin(2) * t73, t47, t83, 0, 0, 0, t49 * t105, t51 * t105, t67 * t47, -0.2e1 * t47 * t89, t87 * t108, t71 * t83, t45, 0.2e1 * t34 * t90 + 0.2e1 * t9 * t49, -0.2e1 * t10 * t49 + 0.2e1 * t34 * t87, t23 ^ 2, 0.2e1 * t23 * t22, -t23 * t108, t22 * t109, t45, 0.2e1 * t1 * t49 + 0.2e1 * t19 * t22, -0.2e1 * t19 * t23 - 0.2e1 * t2 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, 0, 0, 0, 0, 0, -t20, -t21, 0, 0, 0, 0, 0, -t94, t18, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, t73, t96, 0, -t73 * pkin(8), -t82, 0, 0, t51, -t49, 0, -t34, -t36, t38, t28, t41, t42, 0, t78 * t71 - t93, t78 * t76 + t29, -t17, t5, t32, -t31, 0, t53 * t22 + t25 * t49 + t11, -t53 * t23 - t26 * t49 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t81, -0.2e1 * t102, t66, t58, 0, 0, 0, -0.2e1 * t62 * t76, 0.2e1 * t62 * t71, t46, t30, 0, 0, 0, t48 * t107, t50 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, 0, 0, 0, 0, -t94, t18, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, 0, -t34, -t36, t38, t28, t41, t42, 0, t79 * t71 - t93, t79 * t76 + t29, -t17, t5, t32, -t31, 0, t63 * t22 + t33 * t49 + t11, -t63 * t23 - t35 * t49 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t81, -t102, t66, t58, 0, 0, 0, t97 * t76, -t97 * t71, t46, t30, 0, 0, 0, t85 * t48, t85 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t66, t58, 0, 0, 0, 0.2e1 * pkin(4) * t76, -0.2e1 * pkin(4) * t71, t46, t30, 0, 0, 0, t48 * t106, t50 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t90, t49, t9, -t10, 0, 0, -t23, -t22, t49, t49 * t101 + t1, -t100 + (-t7 - t104) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t76, 0, -t71 * t61, -t86, 0, 0, t50, -t48, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t76, 0, -t71 * pkin(10), -t98, 0, 0, t50, -t48, 0, t33, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t101, -0.2e1 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22, t49, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t48, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t48, 0, t33, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t101, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
