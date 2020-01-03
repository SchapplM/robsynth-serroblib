% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t62 = cos(pkin(5));
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t60 = sin(pkin(5));
t65 = sin(qJ(2));
t83 = t60 * t65;
t36 = t62 * t64 + t67 * t83;
t59 = sin(pkin(10));
t61 = cos(pkin(10));
t68 = cos(qJ(2));
t82 = t60 * t68;
t22 = t36 * t59 + t61 * t82;
t23 = t36 * t61 - t59 * t82;
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t10 = t66 * t22 + t63 * t23;
t96 = -0.2e1 * t10;
t41 = t63 * t59 - t66 * t61;
t34 = t41 * t64;
t95 = 0.2e1 * t34;
t94 = -0.2e1 * t36;
t53 = -t61 * pkin(4) - pkin(3);
t93 = 0.2e1 * t53;
t92 = 0.2e1 * t67;
t91 = pkin(1) * t65;
t90 = pkin(1) * t68;
t89 = pkin(8) * t59;
t54 = t64 * pkin(8);
t88 = t67 * pkin(8);
t75 = pkin(7) * t82;
t31 = t75 + (pkin(8) + t91) * t62;
t32 = (-pkin(2) * t68 - pkin(8) * t65 - pkin(1)) * t60;
t18 = -t64 * t31 + t67 * t32;
t17 = pkin(3) * t82 - t18;
t87 = t17 * t59;
t86 = t17 * t61;
t56 = t60 ^ 2;
t85 = t56 * t68;
t84 = t59 * t64;
t81 = t61 * t64;
t80 = t62 * t65;
t79 = pkin(9) + qJ(4);
t48 = pkin(7) * t83;
t30 = t48 + (-pkin(2) - t90) * t62;
t35 = -t62 * t67 + t64 * t83;
t15 = t35 * pkin(3) - t36 * qJ(4) + t30;
t19 = t67 * t31 + t64 * t32;
t16 = -qJ(4) * t82 + t19;
t6 = t59 * t15 + t61 * t16;
t44 = -t67 * pkin(3) - t64 * qJ(4) - pkin(2);
t29 = t59 * t44 + t61 * t88;
t78 = t59 ^ 2 + t61 ^ 2;
t77 = qJ(4) * t35;
t76 = 0.2e1 * t82;
t74 = t64 * t82;
t73 = t67 * t82;
t5 = t61 * t15 - t59 * t16;
t72 = -t5 * t59 + t6 * t61;
t71 = -pkin(3) * t64 + qJ(4) * t67;
t40 = t61 * t44;
t28 = -t59 * t88 + t40;
t70 = -t28 * t59 + t29 * t61;
t42 = t66 * t59 + t63 * t61;
t58 = t64 ^ 2;
t46 = t79 * t61;
t45 = t79 * t59;
t43 = pkin(4) * t84 + t54;
t38 = pkin(1) * t80 + t75;
t37 = t62 * t90 - t48;
t33 = t42 * t64;
t26 = -t63 * t45 + t66 * t46;
t25 = -t66 * t45 - t63 * t46;
t24 = -pkin(9) * t84 + t29;
t20 = -pkin(9) * t81 + t40 + (-pkin(4) - t89) * t67;
t11 = -t63 * t22 + t66 * t23;
t9 = t63 * t20 + t66 * t24;
t8 = t66 * t20 - t63 * t24;
t7 = t22 * pkin(4) + t17;
t4 = -t22 * pkin(9) + t6;
t3 = t35 * pkin(4) - t23 * pkin(9) + t5;
t2 = t63 * t3 + t66 * t4;
t1 = t66 * t3 - t63 * t4;
t12 = [1, 0, 0, t56 * t65 ^ 2, 0.2e1 * t65 * t85, 0.2e1 * t60 * t80, t62 * t76, t62 ^ 2, 0.2e1 * pkin(1) * t85 + 0.2e1 * t37 * t62, -0.2e1 * t38 * t62 - 0.2e1 * t56 * t91, t36 ^ 2, t35 * t94, t82 * t94, t35 * t76, t56 * t68 ^ 2, -0.2e1 * t18 * t82 + 0.2e1 * t30 * t35, 0.2e1 * t19 * t82 + 0.2e1 * t30 * t36, 0.2e1 * t17 * t22 + 0.2e1 * t5 * t35, 0.2e1 * t17 * t23 - 0.2e1 * t6 * t35, -0.2e1 * t6 * t22 - 0.2e1 * t5 * t23, t17 ^ 2 + t5 ^ 2 + t6 ^ 2, t11 ^ 2, t11 * t96, 0.2e1 * t11 * t35, t35 * t96, t35 ^ 2, 0.2e1 * t1 * t35 + 0.2e1 * t7 * t10, 0.2e1 * t7 * t11 - 0.2e1 * t2 * t35; 0, 0, 0, 0, 0, t83, t82, t62, t37, -t38, t36 * t64, -t64 * t35 + t36 * t67, -t74, -t73, 0, -pkin(2) * t35 + pkin(8) * t74 - t30 * t67, -pkin(2) * t36 + pkin(8) * t73 + t30 * t64, t28 * t35 - t5 * t67 + (pkin(8) * t22 + t87) * t64, -t29 * t35 + t6 * t67 + (pkin(8) * t23 + t86) * t64, -t29 * t22 - t28 * t23 + (-t5 * t61 - t59 * t6) * t64, t17 * t54 + t5 * t28 + t6 * t29, -t11 * t34, t34 * t10 - t11 * t33, -t11 * t67 - t34 * t35, t10 * t67 - t33 * t35, -t35 * t67, -t1 * t67 + t43 * t10 + t7 * t33 + t8 * t35, t43 * t11 + t2 * t67 - t7 * t34 - t9 * t35; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, t64 * t92, 0, 0, 0, pkin(2) * t92, -0.2e1 * pkin(2) * t64, -0.2e1 * t28 * t67 + 0.2e1 * t58 * t89, 0.2e1 * t58 * pkin(8) * t61 + 0.2e1 * t29 * t67, 0.2e1 * (-t28 * t61 - t29 * t59) * t64, t58 * pkin(8) ^ 2 + t28 ^ 2 + t29 ^ 2, t34 ^ 2, t33 * t95, t67 * t95, t33 * t92, t67 ^ 2, 0.2e1 * t43 * t33 - 0.2e1 * t8 * t67, -0.2e1 * t43 * t34 + 0.2e1 * t9 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, -t82, t18, -t19, -pkin(3) * t22 - t59 * t77 - t86, -pkin(3) * t23 - t61 * t77 + t87, (-t22 * t61 + t23 * t59) * qJ(4) + t72, -t17 * pkin(3) + qJ(4) * t72, t11 * t42, -t42 * t10 - t11 * t41, t42 * t35, -t41 * t35, 0, t53 * t10 + t25 * t35 + t7 * t41, t53 * t11 - t26 * t35 + t7 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t67, 0, -t54, -t88, -pkin(8) * t81 + t59 * t71, pkin(8) * t84 + t61 * t71, t70, -pkin(3) * t54 + qJ(4) * t70, -t34 * t42, -t42 * t33 + t34 * t41, -t42 * t67, t41 * t67, 0, -t25 * t67 + t53 * t33 + t43 * t41, t26 * t67 - t53 * t34 + t43 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t61, -0.2e1 * pkin(3) * t59, 0.2e1 * t78 * qJ(4), qJ(4) ^ 2 * t78 + pkin(3) ^ 2, t42 ^ 2, -0.2e1 * t42 * t41, 0, 0, 0, t41 * t93, t42 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t17, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t81, 0, t54, 0, 0, 0, 0, 0, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t59, 0, -pkin(3), 0, 0, 0, 0, 0, t41, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t35, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t33, -t67, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
