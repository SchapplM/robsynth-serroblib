% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t94 = -t52 * pkin(5) + t55 * qJ(6);
t53 = sin(qJ(3));
t54 = sin(qJ(2));
t56 = cos(qJ(3));
t57 = cos(qJ(2));
t29 = t53 * t54 - t56 * t57;
t30 = t53 * t57 + t56 * t54;
t46 = -t57 * pkin(2) - pkin(1);
t61 = -t30 * qJ(4) + t46;
t17 = t29 * pkin(3) + t61;
t93 = -0.2e1 * t17;
t84 = t53 * pkin(2);
t42 = qJ(4) + t84;
t92 = 0.2e1 * t42;
t91 = 0.2e1 * t46;
t90 = 0.2e1 * t52;
t89 = -0.2e1 * t55;
t88 = 0.2e1 * t57;
t59 = 0.2e1 * qJ(4);
t87 = pkin(3) + pkin(9);
t86 = -pkin(8) - pkin(7);
t85 = t30 * pkin(5);
t83 = t56 * pkin(2);
t36 = t86 * t54;
t37 = t86 * t57;
t20 = t53 * t36 - t56 * t37;
t35 = pkin(5) * t55 + t52 * qJ(6);
t8 = (-pkin(4) - t35) * t29 + t20;
t82 = t8 * t55;
t10 = t87 * t29 + t61;
t19 = -t56 * t36 - t53 * t37;
t15 = t30 * pkin(4) + t19;
t6 = t55 * t10 + t52 * t15;
t81 = t30 * t29;
t45 = -pkin(3) - t83;
t40 = -pkin(9) + t45;
t80 = t30 * t40;
t79 = t30 * t87;
t78 = t42 * t29;
t24 = t52 * t29;
t77 = t52 * t30;
t76 = t52 * t40;
t75 = t52 * t87;
t74 = t55 * t29;
t25 = t55 * t30;
t73 = t55 * t52;
t33 = qJ(4) - t94;
t26 = t33 + t84;
t72 = t26 + t33;
t49 = t52 ^ 2;
t50 = t55 ^ 2;
t38 = t49 + t50;
t71 = qJ(4) * t29;
t70 = t30 * qJ(6);
t68 = qJ(4) + t42;
t67 = 0.2e1 * t81;
t32 = t38 * t87;
t66 = t52 * t10 - t55 * t15;
t3 = t70 + t6;
t4 = t66 - t85;
t1 = t3 * t52 - t4 * t55;
t65 = -t26 * t29 + t80;
t64 = -t29 * t33 - t79;
t63 = t78 - t80;
t62 = t71 + t79;
t60 = -0.2e1 * pkin(3);
t47 = t55 * t87;
t41 = -0.2e1 * t73;
t34 = t55 * t40;
t28 = t30 ^ 2;
t27 = t29 ^ 2;
t23 = t29 * t73;
t22 = t38 * t40;
t18 = (-t49 + t50) * t29;
t16 = -t29 * pkin(4) + t20;
t14 = t16 * t55;
t13 = t16 * t52;
t7 = t8 * t52;
t2 = [1, 0, 0, t54 ^ 2, t54 * t88, 0, 0, 0, pkin(1) * t88, -0.2e1 * pkin(1) * t54, t28, -0.2e1 * t81, 0, 0, 0, t29 * t91, t30 * t91, 0.2e1 * t19 * t30 - 0.2e1 * t20 * t29, t29 * t93, t30 * t93, t17 ^ 2 + t19 ^ 2 + t20 ^ 2, t49 * t27, 0.2e1 * t27 * t73, t52 * t67, t55 * t67, t28, -0.2e1 * t16 * t74 - 0.2e1 * t30 * t66, 0.2e1 * t16 * t24 - 0.2e1 * t6 * t30, -0.2e1 * t4 * t30 - 0.2e1 * t74 * t8, 0.2e1 * (t3 * t55 + t4 * t52) * t29, -0.2e1 * t24 * t8 + 0.2e1 * t3 * t30, t3 ^ 2 + t4 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, t54, t57, 0, -t54 * pkin(7), -t57 * pkin(7), 0, 0, t30, -t29, 0, -t19, -t20, t45 * t30 - t78, t19, t20, t19 * t45 + t20 * t42, t23, t18, t25, -t77, 0, -t55 * t63 + t13, t52 * t63 + t14, t55 * t65 + t7, -t1, t52 * t65 - t82, t1 * t40 + t8 * t26; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t83, -0.2e1 * t84, 0, 0.2e1 * t45, t92, t42 ^ 2 + t45 ^ 2, t50, t41, 0, 0, 0, t42 * t90, t55 * t92, t26 * t90, -0.2e1 * t22, t26 * t89, t38 * t40 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t19, -t20, -pkin(3) * t30 - t71, t19, t20, -t19 * pkin(3) + t20 * qJ(4), t23, t18, t25, -t77, 0, -t55 * t62 + t13, t52 * t62 + t14, t55 * t64 + t7, -t1, t52 * t64 - t82, -t1 * t87 + t8 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t83, -t84, 0, t60 - t83, t59 + t84, -t45 * pkin(3) + t42 * qJ(4), t50, t41, 0, 0, 0, t68 * t52, t68 * t55, t72 * t52 (-t40 + t87) * t38, -t72 * t55, t26 * t33 - t32 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t60, t59, pkin(3) ^ 2 + qJ(4) ^ 2, t50, t41, 0, 0, 0, t52 * t59, t55 * t59, t33 * t90, 0.2e1 * t32, t33 * t89, t38 * t87 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, t19, 0, 0, 0, 0, 0, t25, -t77, t25, 0, t77, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t74, t30, -t66, -t6, -t66 + 0.2e1 * t85, t94 * t29, 0.2e1 * t70 + t6, -t4 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t52, 0, t34, -t76, t34, -t35, t76, t35 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t52, 0, -t47, t75, -t47, -t35, -t75, -t35 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t52, t55, 0, t52, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t24, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
