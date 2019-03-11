% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t55 = cos(pkin(8));
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t51 = sin(pkin(8));
t60 = sin(qJ(4));
t86 = t51 * t60;
t37 = t59 * t55 + t63 * t86;
t58 = sin(qJ(6));
t62 = cos(qJ(6));
t64 = cos(qJ(4));
t85 = t51 * t64;
t25 = t58 * t37 + t62 * t85;
t97 = -0.2e1 * t25;
t96 = -0.2e1 * t37;
t95 = -0.2e1 * t59;
t94 = 0.2e1 * t63;
t93 = pkin(3) * t60;
t92 = pkin(5) * t62;
t91 = pkin(11) * t58;
t68 = pkin(10) * t85;
t32 = t68 + (pkin(11) + t93) * t55;
t33 = (-pkin(4) * t64 - pkin(11) * t60 - pkin(3)) * t51;
t18 = -t59 * t32 + t63 * t33;
t16 = pkin(5) * t85 - t18;
t90 = t16 * t58;
t89 = t16 * t62;
t26 = t62 * t37 - t58 * t85;
t88 = t26 * t58;
t45 = t51 ^ 2;
t87 = t45 * t64;
t52 = sin(pkin(7));
t61 = sin(qJ(3));
t84 = t52 * t61;
t65 = cos(qJ(3));
t83 = t52 * t65;
t54 = cos(pkin(14));
t56 = cos(pkin(7));
t82 = t54 * t56;
t81 = t55 * t60;
t80 = t55 * t64;
t79 = t55 * t65;
t36 = -t63 * t55 + t59 * t86;
t78 = t58 * t36;
t77 = t58 * t59;
t76 = t58 * t62;
t75 = t58 * t63;
t74 = t59 * t36;
t73 = t62 * t36;
t72 = t62 * t59;
t71 = t62 * t63;
t70 = 0.2e1 * t85;
t69 = t59 * t94;
t67 = t59 * t85;
t66 = t63 * t85;
t19 = t63 * t32 + t59 * t33;
t43 = pkin(10) * t86;
t31 = t43 + (-pkin(3) * t64 - pkin(4)) * t55;
t57 = cos(pkin(6));
t53 = sin(pkin(6));
t50 = sin(pkin(14));
t49 = t62 ^ 2;
t48 = t59 ^ 2;
t47 = t58 ^ 2;
t41 = -t63 * pkin(5) - t59 * pkin(12) - pkin(4);
t39 = pkin(3) * t81 + t68;
t38 = pkin(3) * t80 - t43;
t35 = -t51 * t83 + t55 * t56;
t34 = -t53 * t54 * t52 + t57 * t56;
t29 = pkin(11) * t71 + t58 * t41;
t28 = -pkin(11) * t75 + t62 * t41;
t24 = t56 * t86 + (t60 * t79 + t61 * t64) * t52;
t23 = -t64 * t52 * t79 - t56 * t85 + t60 * t84;
t22 = t57 * t84 + (t50 * t65 + t61 * t82) * t53;
t21 = t57 * t83 + (-t50 * t61 + t65 * t82) * t53;
t17 = -pkin(12) * t85 + t19;
t15 = t36 * pkin(5) - t37 * pkin(12) + t31;
t14 = t63 * t24 + t59 * t35;
t13 = t59 * t24 - t63 * t35;
t12 = -t21 * t51 + t34 * t55;
t10 = t62 * t14 + t58 * t23;
t9 = -t58 * t14 + t62 * t23;
t8 = t22 * t64 + (t21 * t55 + t34 * t51) * t60;
t7 = -t21 * t80 + t22 * t60 - t34 * t85;
t6 = t58 * t15 + t62 * t17;
t5 = t62 * t15 - t58 * t17;
t4 = t12 * t59 + t8 * t63;
t3 = -t12 * t63 + t8 * t59;
t2 = t4 * t62 + t7 * t58;
t1 = -t4 * t58 + t7 * t62;
t11 = [1, t57 ^ 2 + (t50 ^ 2 + t54 ^ 2) * t53 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t21, -t22, 0, 0, 0, 0, 0, -t12 * t85 - t7 * t55, t12 * t86 - t8 * t55, 0, 0, 0, 0, 0, t3 * t85 + t7 * t36, t7 * t37 + t4 * t85, 0, 0, 0, 0, 0, t1 * t36 + t3 * t25, -t2 * t36 + t3 * t26; 0, 0, 0, t83, -t84, 0, 0, 0, 0, 0, -t23 * t55 - t35 * t85, -t24 * t55 + t35 * t86, 0, 0, 0, 0, 0, t13 * t85 + t23 * t36, t14 * t85 + t23 * t37, 0, 0, 0, 0, 0, t13 * t25 + t9 * t36, -t10 * t36 + t13 * t26; 0, 0, 1, 0, 0, t45 * t60 ^ 2, 0.2e1 * t60 * t87, 0.2e1 * t51 * t81, t55 * t70, t55 ^ 2, 0.2e1 * pkin(3) * t87 + 0.2e1 * t38 * t55, -0.2e1 * t39 * t55 - 0.2e1 * t45 * t93, t37 ^ 2, t36 * t96, t85 * t96, t36 * t70, t45 * t64 ^ 2, -0.2e1 * t18 * t85 + 0.2e1 * t31 * t36, 0.2e1 * t19 * t85 + 0.2e1 * t31 * t37, t26 ^ 2, t26 * t97, 0.2e1 * t26 * t36, t36 * t97, t36 ^ 2, 0.2e1 * t16 * t25 + 0.2e1 * t5 * t36, 0.2e1 * t16 * t26 - 0.2e1 * t6 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, -t7 * t63, t7 * t59, 0, 0, 0, 0, 0, -t1 * t63 + t3 * t77, t2 * t63 + t3 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, 0, 0, 0, 0, -t23 * t63, t23 * t59, 0, 0, 0, 0, 0, t13 * t77 - t9 * t63, t10 * t63 + t13 * t72; 0, 0, 0, 0, 0, 0, 0, t86, t85, t55, t38, -t39, t37 * t59, t37 * t63 - t74, -t67, -t66, 0, -pkin(4) * t36 + pkin(11) * t67 - t31 * t63, -pkin(4) * t37 + pkin(11) * t66 + t31 * t59, t26 * t72 (-t25 * t62 - t88) * t59, -t26 * t63 + t36 * t72, t25 * t63 - t58 * t74, -t36 * t63, t28 * t36 - t5 * t63 + (pkin(11) * t25 + t90) * t59, -t29 * t36 + t6 * t63 + (pkin(11) * t26 + t89) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, t69, 0, 0, 0, pkin(4) * t94, pkin(4) * t95, t49 * t48, -0.2e1 * t48 * t76, t71 * t95, t58 * t69, t63 ^ 2, -0.2e1 * t28 * t63 + 0.2e1 * t48 * t91, 0.2e1 * t48 * pkin(11) * t62 + 0.2e1 * t29 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t3 * t62, t3 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, 0, 0, 0, 0, -t13 * t62, t13 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, -t85, t18, -t19, t88, -t58 * t25 + t26 * t62, t78, t73, 0, -pkin(5) * t25 - pkin(12) * t78 - t89, -pkin(5) * t26 - pkin(12) * t73 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t63, 0, -t59 * pkin(11), -t63 * pkin(11), t58 * t72 (-t47 + t49) * t59, -t75, -t71, 0, -pkin(11) * t72 + (-pkin(5) * t59 + pkin(12) * t63) * t58, pkin(12) * t71 + (t91 - t92) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t76, 0, 0, 0, 0.2e1 * t92, -0.2e1 * pkin(5) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t36, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t77, -t63, t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t62, 0, -t58 * pkin(12), -t62 * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
