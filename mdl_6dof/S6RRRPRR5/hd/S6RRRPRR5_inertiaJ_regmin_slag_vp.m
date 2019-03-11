% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t66 = sin(qJ(3));
t67 = sin(qJ(2));
t70 = cos(qJ(3));
t71 = cos(qJ(2));
t44 = t66 * t71 + t70 * t67;
t104 = -0.2e1 * t44;
t94 = t66 * pkin(2);
t54 = qJ(4) + t94;
t65 = sin(qJ(5));
t61 = t65 * pkin(5);
t48 = t54 + t61;
t103 = 0.2e1 * t48;
t102 = 0.2e1 * t54;
t56 = qJ(4) + t61;
t101 = 0.2e1 * t56;
t59 = -t71 * pkin(2) - pkin(1);
t100 = 0.2e1 * t59;
t99 = 0.2e1 * t71;
t73 = 0.2e1 * qJ(4);
t98 = pkin(3) + pkin(9);
t97 = -pkin(8) - pkin(7);
t96 = t44 * pkin(5);
t64 = sin(qJ(6));
t95 = t64 * pkin(5);
t68 = cos(qJ(6));
t93 = t68 * pkin(5);
t69 = cos(qJ(5));
t42 = t66 * t67 - t70 * t71;
t75 = -t44 * qJ(4) + t59;
t13 = t98 * t42 + t75;
t78 = pkin(10) * t42 + t13;
t50 = t97 * t67;
t51 = t97 * t71;
t30 = -t70 * t50 - t66 * t51;
t17 = t44 * pkin(4) + t30;
t87 = t65 * t17;
t5 = t78 * t69 + t87;
t92 = t68 * t5;
t91 = t69 * pkin(10);
t90 = t70 * pkin(2);
t41 = t64 * t69 + t68 * t65;
t28 = t41 * t44;
t89 = t44 * t42;
t88 = t54 * t42;
t86 = t65 * t42;
t85 = t65 * t44;
t84 = t69 * t42;
t83 = t69 * t65;
t82 = t48 + t56;
t81 = qJ(4) * t42;
t80 = qJ(4) + t54;
t79 = 0.2e1 * t89;
t58 = -pkin(3) - t90;
t14 = t69 * t17;
t4 = -t78 * t65 + t14 + t96;
t1 = t68 * t4 - t64 * t5;
t52 = -pkin(9) + t58;
t77 = -t44 * t52 + t88;
t31 = t66 * t50 - t70 * t51;
t76 = t44 * t98 + t81;
t74 = -0.2e1 * pkin(3);
t63 = t69 ^ 2;
t62 = t65 ^ 2;
t60 = t69 * t98;
t53 = -0.2e1 * t83;
t49 = t69 * t52;
t47 = -t60 - t91;
t46 = (-pkin(10) - t98) * t65;
t43 = -t64 * t65 + t68 * t69;
t40 = t44 ^ 2;
t39 = t43 ^ 2;
t38 = t42 ^ 2;
t37 = t49 - t91;
t36 = (-pkin(10) + t52) * t65;
t35 = t69 * t44;
t33 = t42 * t83;
t29 = t43 * t44;
t27 = -0.2e1 * t43 * t41;
t26 = t68 * t46 + t64 * t47;
t25 = -t64 * t46 + t68 * t47;
t24 = (-t62 + t63) * t42;
t23 = t42 * pkin(3) + t75;
t22 = t68 * t36 + t64 * t37;
t21 = -t64 * t36 + t68 * t37;
t20 = t41 * t42;
t19 = t64 * t86 - t68 * t84;
t18 = -t42 * pkin(4) + t31;
t16 = t18 * t69;
t15 = t18 * t65;
t12 = t20 * t43;
t11 = (-pkin(5) * t69 - pkin(4)) * t42 + t31;
t10 = t11 * t43;
t9 = t11 * t41;
t8 = t69 * t13 + t87;
t7 = -t65 * t13 + t14;
t6 = -t43 * t19 - t20 * t41;
t2 = t64 * t4 + t92;
t3 = [1, 0, 0, t67 ^ 2, t67 * t99, 0, 0, 0, pkin(1) * t99, -0.2e1 * pkin(1) * t67, t40, -0.2e1 * t89, 0, 0, 0, t42 * t100, t44 * t100, 0.2e1 * t30 * t44 - 0.2e1 * t31 * t42, -0.2e1 * t23 * t42, t23 * t104, t23 ^ 2 + t30 ^ 2 + t31 ^ 2, t62 * t38, 0.2e1 * t38 * t83, t65 * t79, t69 * t79, t40, -0.2e1 * t18 * t84 + 0.2e1 * t7 * t44, 0.2e1 * t18 * t86 - 0.2e1 * t8 * t44, t20 ^ 2, -0.2e1 * t20 * t19, 0.2e1 * t20 * t44, t19 * t104, t40, 0.2e1 * t1 * t44 + 0.2e1 * t11 * t19, 0.2e1 * t11 * t20 - 0.2e1 * t2 * t44; 0, 0, 0, 0, 0, t67, t71, 0, -t67 * pkin(7), -t71 * pkin(7), 0, 0, t44, -t42, 0, -t30, -t31, t58 * t44 - t88, t30, t31, t30 * t58 + t31 * t54, t33, t24, t35, -t85, 0, -t69 * t77 + t15, t65 * t77 + t16, t12, t6, t29, -t28, 0, t48 * t19 + t21 * t44 + t9, t48 * t20 - t22 * t44 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t90, -0.2e1 * t94, 0, 0.2e1 * t58, t102, t54 ^ 2 + t58 ^ 2, t63, t53, 0, 0, 0, t65 * t102, t69 * t102, t39, t27, 0, 0, 0, t41 * t103, t43 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, 0, -t30, -t31, -pkin(3) * t44 - t81, t30, t31, -t30 * pkin(3) + t31 * qJ(4), t33, t24, t35, -t85, 0, -t69 * t76 + t15, t65 * t76 + t16, t12, t6, t29, -t28, 0, t56 * t19 + t25 * t44 + t9, t56 * t20 - t26 * t44 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t90, -t94, 0, t74 - t90, t73 + t94, -t58 * pkin(3) + t54 * qJ(4), t63, t53, 0, 0, 0, t80 * t65, t80 * t69, t39, t27, 0, 0, 0, t82 * t41, t82 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t74, t73, pkin(3) ^ 2 + qJ(4) ^ 2, t63, t53, 0, 0, 0, t65 * t73, t69 * t73, t39, t27, 0, 0, 0, t41 * t101, t43 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, t30, 0, 0, 0, 0, 0, t35, -t85, 0, 0, 0, 0, 0, t29, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t84, t44, t7, -t8, 0, 0, t20, -t19, t44, t44 * t93 + t1, -t92 + (-t4 - t96) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t65, 0, t49, -t65 * t52, 0, 0, t43, -t41, 0, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t65, 0, -t60, t65 * t98, 0, 0, t43, -t41, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t65, 0, 0, 0, 0, 0, t43, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t93, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t44, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t93, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
