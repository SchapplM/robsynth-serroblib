% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t53 = sin(pkin(10));
t55 = cos(pkin(10));
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t37 = t53 * t60 + t55 * t58;
t57 = sin(qJ(5));
t49 = t57 ^ 2;
t23 = t49 * t37;
t59 = cos(qJ(5));
t51 = t59 ^ 2;
t27 = t51 * t37;
t97 = t23 + t27;
t41 = t49 + t51;
t96 = -0.2e1 * t37;
t65 = pkin(5) * t57 - qJ(6) * t59;
t95 = t65 * t37;
t54 = sin(pkin(9));
t87 = t54 * pkin(1);
t46 = pkin(7) + t87;
t75 = qJ(4) + t46;
t31 = t75 * t60;
t71 = t75 * t58;
t11 = t53 * t31 + t55 * t71;
t94 = t11 ^ 2;
t35 = t53 * t58 - t55 * t60;
t93 = t35 ^ 2;
t56 = cos(pkin(9));
t85 = t56 * pkin(1);
t48 = -pkin(2) - t85;
t38 = -pkin(3) * t60 + t48;
t92 = 0.2e1 * t38;
t91 = 0.2e1 * t58;
t90 = -0.2e1 * t59;
t10 = pkin(4) * t35 - pkin(8) * t37 + t38;
t13 = t55 * t31 - t53 * t71;
t4 = t10 * t57 + t13 * t59;
t89 = t35 * pkin(5);
t88 = t53 * pkin(3);
t86 = t55 * pkin(3);
t84 = t11 * t35;
t45 = pkin(8) + t88;
t83 = t35 * t45;
t24 = t57 * t35;
t25 = t57 * t37;
t82 = t57 * t45;
t81 = t57 * t59;
t28 = t59 * t35;
t29 = t59 * t37;
t80 = t59 * t45;
t79 = t97 * t45;
t78 = t41 * t45 ^ 2;
t50 = t58 ^ 2;
t52 = t60 ^ 2;
t77 = t50 + t52;
t76 = t35 * qJ(6);
t74 = t35 * t25;
t34 = t37 ^ 2;
t73 = t34 * t81;
t47 = -pkin(4) - t86;
t72 = -t10 * t59 + t13 * t57;
t1 = t76 + t4;
t2 = t72 - t89;
t70 = t1 * t59 + t2 * t57;
t69 = t1 * t57 - t2 * t59;
t68 = t4 * t57 - t59 * t72;
t67 = t4 * t59 + t57 * t72;
t66 = pkin(5) * t59 + qJ(6) * t57;
t30 = t47 - t66;
t64 = t30 * t37 - t83;
t63 = t37 * t47 - t83;
t26 = t51 * t34;
t22 = t49 * t34;
t20 = t57 * t29;
t19 = 0.2e1 * t41 * t45;
t16 = 0.2e1 * t35 * t29;
t14 = t23 - t27;
t6 = t26 + t22 + t93;
t5 = t11 + t95;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t85, -0.2e1 * t87, 0 (t54 ^ 2 + t56 ^ 2) * pkin(1) ^ 2, t50, t60 * t91, 0, t52, 0, 0, -0.2e1 * t48 * t60, t48 * t91, 0.2e1 * t77 * t46, t46 ^ 2 * t77 + t48 ^ 2, t34, t35 * t96, 0, t93, 0, 0, t35 * t92, t37 * t92, 0.2e1 * t11 * t37 - 0.2e1 * t13 * t35, t13 ^ 2 + t38 ^ 2 + t94, t26, -0.2e1 * t73, t16, t22, -0.2e1 * t74, t93, 0.2e1 * t11 * t25 - 0.2e1 * t35 * t72, 0.2e1 * t11 * t29 - 0.2e1 * t35 * t4, t68 * t96, t4 ^ 2 + t72 ^ 2 + t94, t26, t16, 0.2e1 * t73, t93, 0.2e1 * t74, t22, -0.2e1 * t2 * t35 + 0.2e1 * t25 * t5, t69 * t96, 0.2e1 * t1 * t35 - 0.2e1 * t29 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t37 + t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t67 + t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t35 + t37 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 + t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t60, 0, -t58 * t46, -t60 * t46, 0, 0, 0, 0, t37, 0, -t35, 0, -t11, -t13 (-t35 * t53 - t37 * t55) * pkin(3) (-t11 * t55 + t13 * t53) * pkin(3), t20, -t14, t24, -t20, t28, 0, -t11 * t59 + t57 * t63, t11 * t57 + t59 * t63, t67, t11 * t47 + t45 * t67, t20, t24, t14, 0, -t28, -t20, -t5 * t59 + t57 * t64, t70, -t5 * t57 - t59 * t64, t5 * t30 + t45 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t37, 0 (-t35 * t55 + t37 * t53) * pkin(3), 0, 0, 0, 0, 0, 0, -t28, t24, t97, t35 * t47 + t79, 0, 0, 0, 0, 0, 0, -t28, t97, -t24, t30 * t35 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t86, -0.2e1 * t88, 0 (t53 ^ 2 + t55 ^ 2) * pkin(3) ^ 2, t49, 0.2e1 * t81, 0, t51, 0, 0, t47 * t90, 0.2e1 * t47 * t57, t19, t47 ^ 2 + t78, t49, 0, -0.2e1 * t81, 0, 0, t51, t30 * t90, t19, -0.2e1 * t30 * t57, t30 ^ 2 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t37, 0, t38, 0, 0, 0, 0, 0, 0, t28, -t24, -t97, t68, 0, 0, 0, 0, 0, 0, t28, -t97, t24, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t25, t35, -t72, -t4, 0, 0, 0, t29, 0, t35, t25, 0, -t72 + 0.2e1 * t89, -t66 * t37, 0.2e1 * t76 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t29, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, t29, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t59, 0, -t82, -t80, 0, 0, 0, t57, 0, 0, -t59, 0, -t82, -t65, t80, -t65 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t57, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t57, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t29, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;