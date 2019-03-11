% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t32 = t53 * t58 + t54 * t56;
t30 = t32 ^ 2;
t34 = -t53 * t56 + t54 * t58;
t31 = t34 ^ 2;
t74 = -t31 - t30;
t55 = sin(qJ(5));
t48 = t55 ^ 2;
t57 = cos(qJ(5));
t50 = t57 ^ 2;
t40 = t48 + t50;
t99 = -0.2e1 * t34;
t98 = (t32 * t53 + t34 * t54) * pkin(3);
t97 = t74 * t57;
t96 = t74 * t55;
t59 = -pkin(1) - pkin(7);
t77 = -qJ(4) + t59;
t37 = t77 * t56;
t72 = t77 * t58;
t15 = t53 * t37 - t54 * t72;
t95 = t15 ^ 2;
t46 = pkin(3) * t56 + qJ(2);
t94 = 0.2e1 * t46;
t93 = -0.2e1 * t57;
t92 = 0.2e1 * qJ(2);
t91 = t32 * pkin(5);
t90 = t53 * pkin(3);
t89 = t54 * pkin(3);
t17 = t54 * t37 + t53 * t72;
t9 = pkin(4) * t32 - pkin(8) * t34 + t46;
t4 = t17 * t57 + t55 * t9;
t88 = t15 * t34;
t44 = pkin(8) + t90;
t87 = t32 * t44;
t45 = -pkin(4) - t89;
t67 = pkin(5) * t57 + qJ(6) * t55;
t28 = t45 - t67;
t86 = t34 * t28;
t85 = t34 * t45;
t24 = t55 * t32;
t84 = t55 * t34;
t83 = t55 * t44;
t82 = t55 * t57;
t26 = t57 * t32;
t27 = t57 * t34;
t81 = t57 * t44;
t80 = t40 * t87;
t79 = t40 * t44 ^ 2;
t49 = t56 ^ 2;
t51 = t58 ^ 2;
t41 = t49 + t51;
t78 = t32 * qJ(6);
t76 = t32 * t84;
t75 = t31 * t82;
t73 = t17 * t55 - t57 * t9;
t1 = t78 + t4;
t2 = t73 - t91;
t71 = t1 * t57 + t2 * t55;
t70 = t1 * t55 - t2 * t57;
t69 = t4 * t55 - t57 * t73;
t68 = t4 * t57 + t55 * t73;
t66 = pkin(5) * t55 - qJ(6) * t57;
t65 = -t17 * t32 + t88;
t64 = t86 - t87;
t63 = t85 - t87;
t60 = qJ(2) ^ 2;
t36 = t41 * t59;
t25 = t50 * t31;
t23 = t48 * t31;
t22 = t55 * t27;
t21 = 0.2e1 * t40 * t44;
t14 = 0.2e1 * t32 * t27;
t12 = t40 * t32;
t11 = t40 * t34;
t10 = (t48 - t50) * t34;
t6 = t30 * t40 + t31;
t5 = t34 * t66 + t15;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t92 (pkin(1) ^ 2) + t60, t51, -0.2e1 * t58 * t56, 0, t49, 0, 0, t56 * t92, t58 * t92, -0.2e1 * t36, t41 * t59 ^ 2 + t60, t31, t32 * t99, 0, t30, 0, 0, t32 * t94, t34 * t94, 0.2e1 * t65, t17 ^ 2 + t46 ^ 2 + t95, t25, -0.2e1 * t75, t14, t23, -0.2e1 * t76, t30, 0.2e1 * t15 * t84 - 0.2e1 * t32 * t73, 0.2e1 * t15 * t27 - 0.2e1 * t32 * t4, t69 * t99, t4 ^ 2 + t73 ^ 2 + t95, t25, t14, 0.2e1 * t75, t30, 0.2e1 * t76, t23, -0.2e1 * t2 * t32 + 0.2e1 * t5 * t84, t70 * t99, 0.2e1 * t1 * t32 - 0.2e1 * t27 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t41, t36, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t65, 0, 0, 0, 0, 0, 0, t96, t97, 0, t32 * t68 - t88, 0, 0, 0, 0, 0, 0, t96, 0, -t97, t32 * t71 - t5 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t56, 0, t58 * t59, -t56 * t59, 0, 0, 0, 0, t34, 0, -t32, 0, -t15, -t17, -t98 (-t15 * t54 + t17 * t53) * pkin(3), t22, -t10, t24, -t22, t26, 0, -t15 * t57 + t55 * t63, t15 * t55 + t57 * t63, t68, t15 * t45 + t44 * t68, t22, t24, t10, 0, -t26, -t22, -t5 * t57 + t55 * t64, t71, -t5 * t55 - t57 * t64, t5 * t28 + t44 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t56, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, t98, 0, 0, 0, 0, 0, 0, t27, -t84, t12, t80 - t85, 0, 0, 0, 0, 0, 0, t27, t12, t84, t80 - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t89, -0.2e1 * t90, 0 (t53 ^ 2 + t54 ^ 2) * pkin(3) ^ 2, t48, 0.2e1 * t82, 0, t50, 0, 0, t45 * t93, 0.2e1 * t45 * t55, t21, t45 ^ 2 + t79, t48, 0, -0.2e1 * t82, 0, 0, t50, t28 * t93, t21, -0.2e1 * t28 * t55, t28 ^ 2 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t34, 0, t46, 0, 0, 0, 0, 0, 0, t26, -t24, -t11, t69, 0, 0, 0, 0, 0, 0, t26, -t11, t24, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t84, t32, -t73, -t4, 0, 0, 0, t27, 0, t32, t84, 0, -t73 + 0.2e1 * t91, -t67 * t34, 0.2e1 * t78 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, t26, -t66 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t57, 0, -t83, -t81, 0, 0, 0, t55, 0, 0, -t57, 0, -t83, -t66, t81, -t66 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t55, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t55, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t27, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
