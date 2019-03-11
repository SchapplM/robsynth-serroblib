% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t75 = sin(qJ(3));
t76 = sin(qJ(2));
t79 = cos(qJ(3));
t80 = cos(qJ(2));
t54 = t75 * t80 + t79 * t76;
t71 = sin(pkin(11));
t72 = cos(pkin(11));
t83 = t75 * t76 - t79 * t80;
t35 = t71 * t54 + t72 * t83;
t109 = -0.2e1 * t35;
t108 = 0.2e1 * t35;
t68 = t79 * pkin(2);
t65 = t68 + pkin(3);
t98 = t75 * pkin(2);
t46 = t72 * t65 - t71 * t98;
t44 = -pkin(4) - t46;
t78 = cos(qJ(5));
t95 = t78 * pkin(5);
t42 = t44 - t95;
t107 = 0.2e1 * t42;
t63 = -t72 * pkin(3) - pkin(4);
t55 = t63 - t95;
t106 = 0.2e1 * t55;
t66 = -t80 * pkin(2) - pkin(1);
t105 = 0.2e1 * t66;
t74 = sin(qJ(5));
t104 = 0.2e1 * t74;
t103 = -0.2e1 * t78;
t102 = 0.2e1 * t80;
t101 = pkin(7) + pkin(8);
t100 = t35 * pkin(5);
t73 = sin(qJ(6));
t99 = t73 * pkin(5);
t77 = cos(qJ(6));
t97 = t77 * pkin(5);
t36 = t72 * t54 - t71 * t83;
t43 = t83 * pkin(3) + t66;
t15 = t35 * pkin(4) - t36 * pkin(9) + t43;
t57 = t101 * t76;
t58 = t101 * t80;
t41 = t75 * t57 - t79 * t58;
t27 = -t83 * qJ(4) - t41;
t40 = -t79 * t57 - t75 * t58;
t82 = -t54 * qJ(4) + t40;
t18 = t72 * t27 + t71 * t82;
t91 = t78 * t18;
t5 = t91 + (-pkin(10) * t36 + t15) * t74;
t96 = t77 * t5;
t16 = t71 * t27 - t72 * t82;
t94 = t16 * t78;
t53 = t73 * t78 + t77 * t74;
t25 = t53 * t35;
t29 = t74 * t35;
t93 = t74 * t36;
t92 = t74 * t78;
t90 = t78 * t36;
t47 = t71 * t65 + t72 * t98;
t45 = pkin(9) + t47;
t89 = t78 * t45;
t62 = t71 * pkin(3) + pkin(9);
t88 = t78 * t62;
t87 = t42 + t55;
t86 = t44 + t63;
t7 = t78 * t15 - t74 * t18;
t4 = -pkin(10) * t90 + t100 + t7;
t1 = t77 * t4 - t73 * t5;
t85 = -t35 * t45 + t36 * t44;
t84 = -t35 * t62 + t36 * t63;
t52 = t73 * t74 - t77 * t78;
t70 = t78 ^ 2;
t69 = t74 ^ 2;
t67 = t78 * pkin(10);
t61 = 0.2e1 * t92;
t51 = t53 ^ 2;
t50 = t67 + t88;
t49 = (-pkin(10) - t62) * t74;
t39 = t67 + t89;
t38 = (-pkin(10) - t45) * t74;
t37 = -0.2e1 * t53 * t52;
t34 = t36 ^ 2;
t33 = t35 ^ 2;
t32 = t73 * t49 + t77 * t50;
t31 = t77 * t49 - t73 * t50;
t30 = t78 * t35;
t28 = t74 * t90;
t24 = t52 * t35;
t23 = t73 * t38 + t77 * t39;
t22 = t77 * t38 - t73 * t39;
t21 = (-t69 + t70) * t36;
t20 = t52 * t36;
t19 = t53 * t36;
t14 = t16 * t74;
t13 = t20 * t53;
t11 = pkin(5) * t93 + t16;
t10 = t11 * t53;
t9 = t11 * t52;
t8 = t74 * t15 + t91;
t6 = -t53 * t19 + t20 * t52;
t2 = t73 * t4 + t96;
t3 = [1, 0, 0, t76 ^ 2, t76 * t102, 0, 0, 0, pkin(1) * t102, -0.2e1 * pkin(1) * t76, t54 ^ 2, -0.2e1 * t54 * t83, 0, 0, 0, t83 * t105, t54 * t105, 0.2e1 * t16 * t36 - 0.2e1 * t18 * t35, t16 ^ 2 + t18 ^ 2 + t43 ^ 2, t70 * t34, -0.2e1 * t34 * t92, t90 * t108, t93 * t109, t33, 0.2e1 * t16 * t93 + 0.2e1 * t7 * t35, 0.2e1 * t16 * t90 - 0.2e1 * t8 * t35, t20 ^ 2, 0.2e1 * t20 * t19, -t20 * t108, t19 * t109, t33, 0.2e1 * t1 * t35 + 0.2e1 * t11 * t19, -0.2e1 * t11 * t20 - 0.2e1 * t2 * t35; 0, 0, 0, 0, 0, t76, t80, 0, -t76 * pkin(7), -t80 * pkin(7), 0, 0, t54, -t83, 0, t40, t41, -t47 * t35 - t46 * t36, -t16 * t46 + t18 * t47, t28, t21, t29, t30, 0, t74 * t85 - t94, t78 * t85 + t14, -t13, t6, t25, -t24, 0, t42 * t19 + t22 * t35 + t9, -t42 * t20 - t23 * t35 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t98, 0, t46 ^ 2 + t47 ^ 2, t69, t61, 0, 0, 0, t44 * t103, t44 * t104, t51, t37, 0, 0, 0, t52 * t107, t53 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t83, 0, t40, t41 (-t35 * t71 - t36 * t72) * pkin(3) (-t16 * t72 + t18 * t71) * pkin(3), t28, t21, t29, t30, 0, t74 * t84 - t94, t78 * t84 + t14, -t13, t6, t25, -t24, 0, t55 * t19 + t31 * t35 + t9, -t55 * t20 - t32 * t35 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t68, -t98, 0 (t46 * t72 + t47 * t71) * pkin(3), t69, t61, 0, 0, 0, -t86 * t78, t86 * t74, t51, t37, 0, 0, 0, t87 * t52, t87 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t71 ^ 2 + t72 ^ 2) * pkin(3) ^ 2, t69, t61, 0, 0, 0, t63 * t103, t63 * t104, t51, t37, 0, 0, 0, t52 * t106, t53 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, -t24, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t93, t35, t7, -t8, 0, 0, -t20, -t19, t35, t35 * t97 + t1, -t96 + (-t4 - t100) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t78, 0, -t74 * t45, -t89, 0, 0, t53, -t52, 0, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t78, 0, -t74 * t62, -t88, 0, 0, t53, -t52, 0, t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t74, 0, 0, 0, 0, 0, -t52, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t97, -0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t19, t35, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, 0, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, 0, t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t97, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
