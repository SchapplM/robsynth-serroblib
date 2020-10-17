% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:42
% EndTime: 2019-12-31 22:43:45
% DurationCPUTime: 0.90s
% Computational Cost: add. (941->145), mult. (2336->287), div. (0->0), fcn. (2712->10), ass. (0->98)
t59 = sin(qJ(5));
t60 = sin(qJ(4));
t63 = cos(qJ(5));
t64 = cos(qJ(4));
t41 = t59 * t60 - t63 * t64;
t61 = sin(qJ(3));
t34 = t41 * t61;
t105 = 0.2e1 * t34;
t58 = cos(pkin(5));
t65 = cos(qJ(3));
t57 = sin(pkin(5));
t62 = sin(qJ(2));
t83 = t57 * t62;
t36 = -t58 * t65 + t61 * t83;
t104 = -0.2e1 * t36;
t103 = 0.2e1 * t36;
t37 = t58 * t61 + t65 * t83;
t102 = -0.2e1 * t37;
t51 = -t64 * pkin(4) - pkin(3);
t101 = 0.2e1 * t51;
t100 = -0.2e1 * t61;
t99 = 0.2e1 * t65;
t98 = pkin(9) + pkin(10);
t97 = pkin(1) * t62;
t66 = cos(qJ(2));
t96 = pkin(1) * t66;
t95 = pkin(3) * t64;
t94 = pkin(8) * t60;
t93 = t36 * pkin(4);
t92 = t59 * pkin(4);
t91 = t63 * pkin(4);
t82 = t57 * t66;
t22 = t37 * t60 + t64 * t82;
t47 = pkin(7) * t83;
t30 = t47 + (-pkin(2) - t96) * t58;
t14 = t36 * pkin(3) - t37 * pkin(9) + t30;
t69 = pkin(7) * t82;
t31 = t69 + (pkin(8) + t97) * t58;
t32 = (-pkin(2) * t66 - pkin(8) * t62 - pkin(1)) * t57;
t18 = t65 * t31 + t61 * t32;
t16 = -pkin(9) * t82 + t18;
t7 = t60 * t14 + t64 * t16;
t5 = -t22 * pkin(10) + t7;
t90 = t63 * t5;
t89 = t65 * pkin(4);
t17 = -t61 * t31 + t65 * t32;
t15 = pkin(3) * t82 - t17;
t88 = t15 * t60;
t87 = t15 * t64;
t23 = t37 * t64 - t60 * t82;
t86 = t23 * t60;
t85 = t36 * t65;
t52 = t57 ^ 2;
t84 = t52 * t66;
t81 = t58 * t62;
t80 = t60 * t36;
t79 = t60 * t64;
t78 = t60 * t65;
t77 = t61 * t36;
t44 = -t65 * pkin(3) - t61 * pkin(9) - pkin(2);
t73 = t64 * t65;
t70 = pkin(8) * t73;
t24 = t70 + (-pkin(10) * t61 + t44) * t60;
t76 = t63 * t24;
t75 = t64 * t36;
t74 = t64 * t61;
t72 = 0.2e1 * t82;
t71 = t61 * t99;
t68 = t65 * t82;
t67 = t61 * t82;
t6 = t64 * t14 - t60 * t16;
t4 = -t23 * pkin(10) + t6 + t93;
t1 = t63 * t4 - t59 * t5;
t40 = t64 * t44;
t21 = -pkin(10) * t74 + t40 + (-pkin(4) - t94) * t65;
t9 = t63 * t21 - t59 * t24;
t42 = t59 * t64 + t63 * t60;
t56 = t65 ^ 2;
t55 = t64 ^ 2;
t54 = t61 ^ 2;
t53 = t60 ^ 2;
t46 = t98 * t64;
t45 = t98 * t60;
t43 = (pkin(4) * t60 + pkin(8)) * t61;
t39 = pkin(1) * t81 + t69;
t38 = t58 * t96 - t47;
t35 = t36 ^ 2;
t33 = t42 * t61;
t29 = t60 * t44 + t70;
t28 = -pkin(8) * t78 + t40;
t26 = -t59 * t45 + t63 * t46;
t25 = -t63 * t45 - t59 * t46;
t12 = -t59 * t22 + t63 * t23;
t11 = t63 * t22 + t59 * t23;
t10 = t59 * t21 + t76;
t8 = t22 * pkin(4) + t15;
t2 = t59 * t4 + t90;
t3 = [1, 0, 0, t52 * t62 ^ 2, 0.2e1 * t62 * t84, 0.2e1 * t57 * t81, t58 * t72, t58 ^ 2, 0.2e1 * pkin(1) * t84 + 0.2e1 * t38 * t58, -0.2e1 * t39 * t58 - 0.2e1 * t52 * t97, t37 ^ 2, t36 * t102, t82 * t102, t36 * t72, t52 * t66 ^ 2, -0.2e1 * t17 * t82 + 0.2e1 * t30 * t36, 0.2e1 * t18 * t82 + 0.2e1 * t30 * t37, t23 ^ 2, -0.2e1 * t23 * t22, t23 * t103, t22 * t104, t35, 0.2e1 * t15 * t22 + 0.2e1 * t6 * t36, 0.2e1 * t15 * t23 - 0.2e1 * t7 * t36, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t103, t11 * t104, t35, 0.2e1 * t1 * t36 + 0.2e1 * t8 * t11, 0.2e1 * t8 * t12 - 0.2e1 * t2 * t36; 0, 0, 0, 0, 0, t83, t82, t58, t38, -t39, t37 * t61, t37 * t65 - t77, -t67, -t68, 0, -pkin(2) * t36 + pkin(8) * t67 - t30 * t65, -pkin(2) * t37 + pkin(8) * t68 + t30 * t61, t23 * t74, (-t22 * t64 - t86) * t61, -t23 * t65 + t36 * t74, t22 * t65 - t60 * t77, -t85, t28 * t36 - t6 * t65 + (pkin(8) * t22 + t88) * t61, -t29 * t36 + t7 * t65 + (pkin(8) * t23 + t87) * t61, -t12 * t34, t34 * t11 - t12 * t33, -t12 * t65 - t34 * t36, t11 * t65 - t33 * t36, -t85, -t1 * t65 + t43 * t11 + t8 * t33 + t9 * t36, -t10 * t36 + t43 * t12 + t2 * t65 - t8 * t34; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t54, t71, 0, 0, 0, pkin(2) * t99, pkin(2) * t100, t55 * t54, -0.2e1 * t54 * t79, t73 * t100, t60 * t71, t56, -0.2e1 * t28 * t65 + 0.2e1 * t54 * t94, 0.2e1 * t54 * pkin(8) * t64 + 0.2e1 * t29 * t65, t34 ^ 2, t33 * t105, t65 * t105, t33 * t99, t56, 0.2e1 * t43 * t33 - 0.2e1 * t9 * t65, 0.2e1 * t10 * t65 - 0.2e1 * t43 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, -t82, t17, -t18, t86, -t60 * t22 + t23 * t64, t80, t75, 0, -pkin(3) * t22 - pkin(9) * t80 - t87, -pkin(3) * t23 - pkin(9) * t75 + t88, t12 * t42, -t42 * t11 - t12 * t41, t42 * t36, -t41 * t36, 0, t51 * t11 + t25 * t36 + t8 * t41, t51 * t12 - t26 * t36 + t8 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t65, 0, -t61 * pkin(8), -t65 * pkin(8), t60 * t74, (-t53 + t55) * t61, -t78, -t73, 0, -pkin(8) * t74 + (-pkin(3) * t61 + pkin(9) * t65) * t60, pkin(9) * t73 + (t94 - t95) * t61, -t34 * t42, -t42 * t33 + t34 * t41, -t42 * t65, t41 * t65, 0, -t25 * t65 + t51 * t33 + t43 * t41, t26 * t65 - t51 * t34 + t43 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t53, 0.2e1 * t79, 0, 0, 0, 0.2e1 * t95, -0.2e1 * pkin(3) * t60, t42 ^ 2, -0.2e1 * t42 * t41, 0, 0, 0, t41 * t101, t42 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t36, t6, -t7, 0, 0, t12, -t11, t36, t36 * t91 + t1, -t90 + (-t4 - t93) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t60 * t61, -t65, t28, -t29, 0, 0, -t34, -t33, -t65, -t63 * t89 + t9, -t76 + (-t21 + t89) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * pkin(9), -t64 * pkin(9), 0, 0, t42, -t41, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t91, -0.2e1 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t36, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t33, -t65, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, 0, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t91, -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
