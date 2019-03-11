% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t55 = sin(pkin(9));
t60 = cos(qJ(4));
t57 = cos(pkin(9));
t86 = sin(qJ(4));
t69 = t86 * t57;
t35 = t55 * t60 + t69;
t54 = sin(pkin(10));
t56 = cos(pkin(10));
t59 = sin(qJ(6));
t87 = cos(qJ(6));
t63 = -t54 * t59 + t56 * t87;
t98 = t63 * t35;
t76 = t60 * t57;
t38 = -t55 * t86 + t76;
t15 = t63 * t38;
t97 = -0.2e1 * t38;
t28 = t35 ^ 2;
t94 = t38 ^ 2;
t96 = t28 + t94;
t58 = -pkin(1) - qJ(3);
t88 = -pkin(7) + t58;
t40 = t88 * t55;
t18 = t40 * t86 - t76 * t88;
t95 = t18 ^ 2;
t93 = 0.2e1 * t35;
t44 = pkin(3) * t55 + qJ(2);
t92 = 0.2e1 * t44;
t47 = -pkin(5) * t56 - pkin(4);
t91 = 0.2e1 * t47;
t90 = 0.2e1 * qJ(2);
t89 = t38 * pkin(4);
t16 = pkin(4) * t35 - qJ(5) * t38 + t44;
t20 = t60 * t40 + t69 * t88;
t6 = t16 * t54 + t20 * t56;
t85 = t15 * t63;
t84 = t18 * t38;
t37 = t54 * t87 + t56 * t59;
t83 = t37 * t35;
t82 = t38 * t37;
t81 = t54 * t35;
t80 = t54 * t38;
t79 = t54 * t56;
t78 = t56 * t38;
t75 = pkin(8) + qJ(5);
t49 = t54 ^ 2;
t51 = t56 ^ 2;
t74 = t49 + t51;
t50 = t55 ^ 2;
t52 = t57 ^ 2;
t43 = t50 + t52;
t73 = t35 * t97;
t72 = t54 * t78;
t5 = t16 * t56 - t54 * t20;
t68 = t74 * qJ(5);
t67 = t5 * t56 + t54 * t6;
t66 = -t5 * t54 + t56 * t6;
t65 = -qJ(5) * t35 - t89;
t64 = -t20 * t35 + t84;
t62 = qJ(2) ^ 2;
t42 = t75 * t56;
t41 = t75 * t54;
t32 = t43 * t58;
t29 = t37 ^ 2;
t27 = t63 ^ 2;
t25 = t56 * t35;
t22 = -t41 * t59 + t42 * t87;
t21 = -t41 * t87 - t42 * t59;
t8 = pkin(5) * t80 + t18;
t7 = t37 * t82;
t4 = -pkin(8) * t80 + t6;
t3 = pkin(5) * t35 - pkin(8) * t78 + t5;
t2 = t3 * t59 + t4 * t87;
t1 = t3 * t87 - t4 * t59;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t90 (pkin(1) ^ 2) + t62, t52, -0.2e1 * t57 * t55, 0, t50, 0, 0, t55 * t90, t57 * t90, -0.2e1 * t32, t43 * t58 ^ 2 + t62, t94, t73, 0, t28, 0, 0, t35 * t92, t38 * t92, 0.2e1 * t64, t20 ^ 2 + t44 ^ 2 + t95, t51 * t94, -0.2e1 * t94 * t79, t78 * t93, t49 * t94, t54 * t73, t28, 0.2e1 * t18 * t80 + 0.2e1 * t35 * t5, 0.2e1 * t18 * t78 - 0.2e1 * t35 * t6, t67 * t97, t5 ^ 2 + t6 ^ 2 + t95, t15 ^ 2, -0.2e1 * t15 * t82, t15 * t93, t82 ^ 2, -t82 * t93, t28, 0.2e1 * t1 * t35 + 0.2e1 * t8 * t82, 0.2e1 * t15 * t8 - 0.2e1 * t2 * t35, -0.2e1 * t1 * t15 - 0.2e1 * t2 * t82, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t43, t32, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t64, 0, 0, 0, 0, 0, 0, -t96 * t54, -t96 * t56, 0, t35 * t66 - t84, 0, 0, 0, 0, 0, 0, -t35 * t83 - t38 * t82, -t15 * t38 - t35 * t98, t15 * t83 - t82 * t98, -t1 * t83 + t2 * t98 - t38 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t74 + t94, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 ^ 2 + t98 ^ 2 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t57, 0, qJ(2), 0, 0, 0, 0, 0, 0, t35, t38, 0, t44, 0, 0, 0, 0, 0, 0, t25, -t81, -t74 * t38, t67, 0, 0, 0, 0, 0, 0, t98, -t83, -t7 - t85, t1 * t63 + t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t98 - t63 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t35, 0, -t18, -t20, 0, 0, t72 (-t49 + t51) * t38, t81, -t72, t25, 0, -t18 * t56 + t54 * t65, t18 * t54 + t56 * t65, t66, -t18 * pkin(4) + qJ(5) * t66, t15 * t37, -t7 + t85, t83, -t82 * t63, t98, 0, t21 * t35 + t47 * t82 - t63 * t8, t15 * t47 - t22 * t35 + t37 * t8, -t1 * t37 - t15 * t21 + t2 * t63 - t22 * t82, t1 * t21 + t2 * t22 + t47 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t80, t74 * t35, t35 * t68 + t89, 0, 0, 0, 0, 0, 0, t15, -t82, t37 * t83 + t63 * t98, -t21 * t83 + t22 * t98 - t38 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t63 + t22 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, 0.2e1 * t79, 0, t51, 0, 0, 0.2e1 * pkin(4) * t56, -0.2e1 * pkin(4) * t54, 0.2e1 * t68, qJ(5) ^ 2 * t74 + pkin(4) ^ 2, t29, 0.2e1 * t37 * t63, 0, t27, 0, 0, -t63 * t91, t37 * t91, -0.2e1 * t21 * t37 + 0.2e1 * t22 * t63, t21 ^ 2 + t22 ^ 2 + t47 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t78, 0, t18, 0, 0, 0, 0, 0, 0, t82, t15, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t54, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t63, t37, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t82, t35, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t98, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t63, 0, t21, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t9;
