% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t53 = sin(qJ(4));
t93 = t53 * pkin(3);
t43 = pkin(10) + t93;
t52 = sin(qJ(5));
t48 = t52 ^ 2;
t56 = cos(qJ(5));
t49 = t56 ^ 2;
t74 = t48 + t49;
t76 = t74 * t43;
t88 = cos(qJ(3));
t45 = -pkin(3) * t88 - pkin(2);
t99 = 0.2e1 * t45;
t98 = -0.2e1 * t52;
t97 = -0.2e1 * t56;
t54 = sin(qJ(3));
t87 = cos(qJ(4));
t31 = t53 * t54 - t87 * t88;
t96 = pkin(10) * t31;
t95 = t31 * pkin(5);
t94 = t52 * pkin(10);
t92 = t56 * pkin(10);
t36 = (-pkin(9) - pkin(8)) * t54;
t70 = t88 * pkin(8);
t37 = pkin(9) * t88 + t70;
t20 = -t36 * t87 + t53 * t37;
t32 = t53 * t88 + t54 * t87;
t62 = pkin(5) * t52 - qJ(6) * t56;
t7 = t32 * t62 + t20;
t91 = t7 * t52;
t90 = t7 * t56;
t69 = t87 * pkin(3);
t44 = -t69 - pkin(4);
t89 = pkin(4) - t44;
t51 = cos(pkin(6));
t50 = sin(pkin(6));
t55 = sin(qJ(2));
t83 = t50 * t55;
t24 = t51 * t88 - t54 * t83;
t67 = t50 * t88;
t25 = t51 * t54 + t55 * t67;
t12 = -t24 * t87 + t25 * t53;
t86 = t12 * t56;
t85 = t20 * t56;
t84 = t31 * t43;
t57 = cos(qJ(2));
t82 = t50 * t57;
t81 = t52 * t32;
t80 = t52 * t43;
t79 = t52 * t56;
t28 = t56 * t32;
t78 = t56 * t43;
t16 = pkin(4) * t31 - pkin(10) * t32 + t45;
t21 = t36 * t53 + t37 * t87;
t6 = t16 * t52 + t21 * t56;
t63 = -pkin(5) * t56 - qJ(6) * t52;
t34 = -pkin(4) + t63;
t29 = -t69 + t34;
t77 = -t29 - t34;
t75 = t74 * pkin(10);
t73 = t31 * qJ(6);
t72 = 0.2e1 * t88;
t71 = -0.2e1 * t32 * t31;
t13 = t24 * t53 + t25 * t87;
t9 = t13 * t52 + t56 * t82;
t68 = t12 * t81 - t31 * t9;
t66 = -t16 * t56 + t21 * t52;
t65 = -pkin(4) * t32 - t96;
t3 = t73 + t6;
t4 = t66 - t95;
t1 = t3 * t56 + t4 * t52;
t64 = -t32 * t34 + t96;
t10 = t13 * t56 - t52 * t82;
t2 = t10 * t56 + t52 * t9;
t61 = t29 * t32 - t84;
t60 = t32 * t44 - t84;
t59 = -t10 * t31 + t12 * t28;
t40 = 0.2e1 * t79;
t30 = t32 ^ 2;
t27 = t56 * t31;
t26 = t52 * t31;
t23 = t52 * t28;
t19 = t20 * t52;
t17 = (-t48 + t49) * t32;
t11 = t12 * t52;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t12 ^ 2 + t9 ^ 2; 0, 0, t82, -t83, 0, 0, 0, 0, 0, t57 * t67, -t54 * t82, 0, 0, 0, 0, 0, -t31 * t82, -t32 * t82, 0, 0, 0, 0, 0, t68, t59, t68 (-t10 * t52 + t56 * t9) * t32, -t59, t10 * t3 + t12 * t7 + t4 * t9; 0, 1, 0, 0, t54 ^ 2, t54 * t72, 0, 0, 0, pkin(2) * t72, -0.2e1 * pkin(2) * t54, t30, t71, 0, 0, 0, t31 * t99, t32 * t99, t49 * t30, -0.2e1 * t30 * t79, 0.2e1 * t31 * t28, t52 * t71, t31 ^ 2, 0.2e1 * t20 * t81 - 0.2e1 * t31 * t66, 0.2e1 * t20 * t28 - 0.2e1 * t31 * t6, -0.2e1 * t31 * t4 + 0.2e1 * t7 * t81, 0.2e1 * (-t3 * t52 + t4 * t56) * t32, -0.2e1 * t28 * t7 + 0.2e1 * t3 * t31, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t86, t11, -t86, t2, -t11, t12 * t29 + t2 * t43; 0, 0, 0, 0, 0, 0, t54, t88, 0, -t54 * pkin(8), -t70, 0, 0, t32, -t31, 0, -t20, -t21, t23, t17, t26, t27, 0, t52 * t60 - t85, t56 * t60 + t19, t52 * t61 - t90, t1, -t56 * t61 - t91, t1 * t43 + t7 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t93, t48, t40, 0, 0, 0, t44 * t97, 0.2e1 * t44 * t52, t29 * t97, 0.2e1 * t76, t29 * t98, t43 ^ 2 * t74 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t86, t11, -t86, t2, -t11, pkin(10) * t2 + t12 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, -t20, -t21, t23, t17, t26, t27, 0, t52 * t65 - t85, t56 * t65 + t19, -t52 * t64 - t90, t1, t56 * t64 - t91, pkin(10) * t1 + t7 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t93, t48, t40, 0, 0, 0, t89 * t56, -t89 * t52, t77 * t56, t75 + t76, t77 * t52, pkin(10) * t76 + t29 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, t40, 0, 0, 0, 0.2e1 * pkin(4) * t56, pkin(4) * t98, t34 * t97, 0.2e1 * t75, t34 * t98, pkin(10) ^ 2 * t74 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -t9, 0, t10, -pkin(5) * t9 + qJ(6) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t81, t31, -t66, -t6, -t66 + 0.2e1 * t95, t63 * t32, 0.2e1 * t73 + t6, -pkin(5) * t4 + qJ(6) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t56, 0, -t80, -t78, -t80, -t62, t78, -t62 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t56, 0, -t94, -t92, -t94, -t62, t92, -t62 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t28, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
