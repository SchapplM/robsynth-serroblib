% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t55 = sin(qJ(4));
t88 = t55 * pkin(3);
t40 = qJ(5) + t88;
t98 = t40 ^ 2;
t58 = cos(qJ(4));
t57 = sin(qJ(2));
t59 = cos(qJ(3));
t78 = t59 * t57;
t56 = sin(qJ(3));
t81 = t56 * t57;
t26 = -t55 * t81 + t58 * t78;
t97 = -0.2e1 * t26;
t30 = t55 * t59 + t58 * t56;
t96 = -0.2e1 * t30;
t45 = -t59 * pkin(3) - pkin(2);
t95 = 0.2e1 * t45;
t94 = -0.2e1 * t57;
t60 = cos(qJ(2));
t93 = 0.2e1 * t60;
t92 = -pkin(9) - pkin(8);
t91 = pkin(2) * t59;
t90 = pkin(7) * t56;
t25 = t30 * t57;
t89 = t25 * pkin(5);
t87 = t58 * pkin(3);
t86 = t60 * pkin(3);
t34 = t92 * t59;
t68 = t92 * t56;
t19 = -t55 * t34 - t58 * t68;
t85 = t19 * t60;
t20 = -t58 * t34 + t55 * t68;
t84 = t20 * t60;
t83 = t40 * t25;
t29 = t55 * t56 - t58 * t59;
t82 = t40 * t29;
t80 = t56 * t59;
t79 = t56 * t60;
t77 = t59 * t60;
t53 = pkin(4) + qJ(6);
t33 = -t60 * pkin(2) - t57 * pkin(8) - pkin(1);
t28 = t59 * t33;
t15 = -pkin(9) * t78 + t28 + (-pkin(3) - t90) * t60;
t70 = pkin(7) * t77;
t18 = t70 + (-pkin(9) * t57 + t33) * t56;
t76 = -t58 * t15 + t55 * t18;
t7 = t55 * t15 + t58 * t18;
t47 = t57 * pkin(7);
t32 = pkin(3) * t81 + t47;
t75 = qJ(5) * t25;
t74 = qJ(5) * t29;
t73 = t40 * qJ(5);
t72 = t60 * qJ(5);
t71 = t57 * t93;
t48 = t60 * pkin(4);
t4 = t48 + t76;
t69 = -0.2e1 * t72 + t7;
t44 = -pkin(4) - t87;
t67 = -t26 * pkin(5) - t4;
t3 = t72 - t7;
t66 = -t26 * qJ(5) + t32;
t65 = -t30 * qJ(5) + t45;
t64 = (-qJ(5) - t40) * t60 + t7;
t63 = -0.2e1 * pkin(4);
t62 = qJ(5) ^ 2;
t61 = 0.2e1 * qJ(5);
t52 = t60 ^ 2;
t51 = t59 ^ 2;
t50 = t57 ^ 2;
t49 = t56 ^ 2;
t42 = t61 + t88;
t37 = qJ(6) - t44;
t35 = 0.2e1 * t40;
t23 = t56 * t33 + t70;
t22 = -pkin(7) * t79 + t28;
t14 = t29 * pkin(4) + t65;
t11 = -t29 * pkin(5) + t20;
t10 = t30 * pkin(5) + t19;
t9 = t53 * t29 + t65;
t8 = t25 * pkin(4) + t66;
t5 = t53 * t25 + t66;
t2 = -t3 - t89;
t1 = t60 * qJ(6) - t67;
t6 = [1, 0, 0, t50, t71, 0, 0, 0, pkin(1) * t93, pkin(1) * t94, t51 * t50, -0.2e1 * t50 * t80, t77 * t94, t56 * t71, t52, -0.2e1 * t22 * t60 + 0.2e1 * t50 * t90, 0.2e1 * t50 * pkin(7) * t59 + 0.2e1 * t23 * t60, t26 ^ 2, t25 * t97, t60 * t97, t25 * t93, t52, 0.2e1 * t32 * t25 + 0.2e1 * t60 * t76, 0.2e1 * t32 * t26 + 0.2e1 * t7 * t60, 0.2e1 * t3 * t25 + 0.2e1 * t4 * t26, -0.2e1 * t8 * t25 - 0.2e1 * t4 * t60, -0.2e1 * t8 * t26 + 0.2e1 * t3 * t60, t3 ^ 2 + t4 ^ 2 + t8 ^ 2, 0.2e1 * t1 * t26 - 0.2e1 * t2 * t25, -0.2e1 * t2 * t60 - 0.2e1 * t5 * t26, 0.2e1 * t1 * t60 + 0.2e1 * t5 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t57, t60, 0, -t47, -t60 * pkin(7), t56 * t78 (-t49 + t51) * t57, -t79, -t77, 0, -pkin(7) * t78 + (-pkin(2) * t57 + pkin(8) * t60) * t56, pkin(8) * t77 + (t90 - t91) * t57, t26 * t30, -t30 * t25 - t26 * t29, -t30 * t60, t29 * t60, 0, t45 * t25 + t32 * t29 + t85, t45 * t26 + t32 * t30 + t84, t19 * t26 - t20 * t25 + t3 * t29 + t4 * t30, -t14 * t25 - t8 * t29 - t85, -t14 * t26 - t8 * t30 - t84, t8 * t14 + t4 * t19 - t3 * t20, t1 * t30 + t10 * t26 - t11 * t25 - t2 * t29, -t11 * t60 - t9 * t26 - t5 * t30, t10 * t60 + t9 * t25 + t5 * t29, t1 * t10 + t2 * t11 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t49, 0.2e1 * t80, 0, 0, 0, 0.2e1 * t91, -0.2e1 * pkin(2) * t56, t30 ^ 2, t29 * t96, 0, 0, 0, t29 * t95, t30 * t95, 0.2e1 * t19 * t30 - 0.2e1 * t20 * t29, -0.2e1 * t14 * t29, t14 * t96, t14 ^ 2 + t19 ^ 2 + t20 ^ 2, 0.2e1 * t10 * t30 - 0.2e1 * t11 * t29, t9 * t96, 0.2e1 * t9 * t29, t10 ^ 2 + t11 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t81, -t60, t22, -t23, 0, 0, t26, -t25, -t60, -t58 * t86 - t76, t55 * t86 - t7, t44 * t26 - t83, -t44 * t60 + t4, t64, -t3 * t40 + t4 * t44, -t37 * t26 - t83, t64 - t89 (-qJ(6) - t37) * t60 + t67, -t1 * t37 + t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t59, 0, -t56 * pkin(8), -t59 * pkin(8), 0, 0, t30, -t29, 0, -t19, -t20, t44 * t30 - t82, t19, t20, t19 * t44 + t20 * t40, -t37 * t30 - t82, t11, -t10, -t10 * t37 + t11 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t87, -0.2e1 * t88, 0, 0.2e1 * t44, t35, t44 ^ 2 + t98, 0, t35, 0.2e1 * t37, t37 ^ 2 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, -t60, -t76, -t7, -pkin(4) * t26 - t75, 0.2e1 * t48 + t76, t69, -t4 * pkin(4) - t3 * qJ(5), -t53 * t26 - t75, t69 - t89 (-qJ(6) - t53) * t60 + t67, t2 * qJ(5) - t1 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t19, -t20, -pkin(4) * t30 - t74, t19, t20, -t19 * pkin(4) + t20 * qJ(5), -t53 * t30 - t74, t11, -t10, t11 * qJ(5) - t10 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t87, -t88, 0, t63 - t87, t42, -t44 * pkin(4) + t73, 0, t42, 0.2e1 * pkin(4) + 0.2e1 * qJ(6) + t87, t37 * t53 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t63, t61, pkin(4) ^ 2 + t62, 0, t61, 0.2e1 * t53, t53 ^ 2 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t60, 0, t4, t26, 0, t60, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, t19, t30, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t44, 0, 0, -1, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t60, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
