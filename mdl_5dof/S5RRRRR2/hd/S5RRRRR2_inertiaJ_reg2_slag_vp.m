% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaJ_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:43
% EndTime: 2019-12-05 18:53:46
% DurationCPUTime: 0.80s
% Computational Cost: add. (271->58), mult. (980->121), div. (0->0), fcn. (1041->8), ass. (0->80)
t49 = sin(qJ(5));
t43 = t49 ^ 2;
t53 = cos(qJ(5));
t47 = t53 ^ 2;
t77 = t43 + t47;
t50 = sin(qJ(4));
t51 = sin(qJ(3));
t54 = cos(qJ(4));
t55 = cos(qJ(3));
t35 = t50 * t55 + t54 * t51;
t52 = sin(qJ(2));
t87 = t52 * pkin(1);
t15 = t35 * t87;
t93 = t15 ^ 2;
t33 = t50 * t51 - t54 * t55;
t92 = t33 ^ 2;
t56 = cos(qJ(2));
t84 = t56 * pkin(1);
t85 = t55 * pkin(2);
t37 = -t84 - t85;
t91 = 0.2e1 * t37;
t69 = t55 * t87;
t75 = t51 * t87;
t17 = -t50 * t75 + t54 * t69;
t5 = -t49 * t17 + t53 * t37;
t80 = t49 * t35;
t90 = t15 * t80 + t5 * t33;
t27 = t53 * t35;
t6 = t53 * t17 + t49 * t37;
t89 = t15 * t27 - t6 * t33;
t88 = t50 * pkin(2);
t86 = t54 * pkin(2);
t83 = t15 * t53;
t82 = t15 * t54;
t58 = pkin(1) ^ 2;
t81 = t52 ^ 2 * t58;
t79 = t49 * t53;
t70 = t35 * t85;
t78 = t77 * t70;
t45 = t51 ^ 2;
t48 = t55 ^ 2;
t76 = t45 + t48;
t13 = -0.2e1 * t35 * t33;
t74 = t53 * t88;
t73 = t49 * t86;
t72 = t53 * t86;
t71 = t33 * t85;
t68 = t51 * t84;
t67 = t55 * t84;
t66 = t37 - t85;
t57 = pkin(2) ^ 2;
t65 = t77 * t57;
t64 = t49 * t71;
t63 = t53 * t71;
t62 = -t49 * t6 - t5 * t53;
t1 = -t5 * t49 + t6 * t53;
t61 = -t33 * t74 - t35 * t72;
t60 = t62 * t35;
t9 = (-t33 * t50 - t35 * t54) * pkin(2);
t59 = t49 * t9;
t44 = t50 ^ 2;
t42 = t56 ^ 2 * t58;
t41 = t54 ^ 2 * t57;
t39 = 0.2e1 * t51 * t55;
t38 = 0.2e1 * t79;
t32 = t35 ^ 2;
t29 = t76 * t87;
t28 = t77 * t88;
t26 = t53 * t33;
t25 = t47 * t32;
t24 = t49 * t33;
t23 = t43 * t32;
t22 = t49 * t27;
t21 = -0.2e1 * t32 * t79;
t14 = t15 * t49;
t12 = 0.2e1 * t33 * t27;
t11 = t49 * t13;
t10 = (-t43 + t47) * t35;
t2 = t15 * t35 - t17 * t33;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t84, -0.2e1 * t87, 0, t42 + t81, t45, t39, 0, t48, 0, 0, 0.2e1 * t67, -0.2e1 * t68, 0.2e1 * t29, t76 * t81 + t42, t32, t13, 0, t92, 0, 0, t33 * t91, t35 * t91, 0.2e1 * t2, t17 ^ 2 + t37 ^ 2 + t93, t25, t21, t12, t23, t11, t92, 0.2e1 * t90, 0.2e1 * t89, 0.2e1 * t60, t5 ^ 2 + t6 ^ 2 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t84, -t87, 0, 0, t45, t39, 0, t48, 0, 0, t67, -t68, t29, 0, t32, t13, 0, t92, 0, 0, t66 * t33, t66 * t35, t2, -t37 * t85, t25, t21, t12, t23, t11, t92, -t63 + t90, t64 + t89, t60 + t78, t62 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, t39, 0, t48, 0, 0, 0, 0, 0, 0, t32, t13, 0, t92, 0, 0, -0.2e1 * t71, -0.2e1 * t70, 0, t48 * t57, t25, t21, t12, t23, t11, t92, -0.2e1 * t63, 0.2e1 * t64, 0.2e1 * t78, t48 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t55, 0, -t75, -t69, 0, 0, 0, 0, t35, 0, -t33, 0, -t15, -t17, t9, (t17 * t50 - t82) * pkin(2), t22, t10, t24, -t22, t26, 0, t59 - t83, t14 + t61, t1, (t1 * t50 - t82) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t55, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, 0, 0, 0, t9, 0, t22, t10, t24, -t22, t26, 0, t59, t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t86, -0.2e1 * t88, 0, t44 * t57 + t41, t43, t38, 0, t47, 0, 0, 0.2e1 * t72, -0.2e1 * t73, 0.2e1 * t28, t44 * t65 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, 0, -t15, -t17, 0, 0, t22, t10, t24, -t22, t26, 0, -t83, t14, t1, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, 0, 0, 0, 0, 0, t22, t10, t24, -t22, t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t86, -t88, 0, 0, t43, t38, 0, t47, 0, 0, t72, -t73, t28, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, t38, 0, t47, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t80, t33, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t80, t33, -t53 * t85, t49 * t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t53, 0, -t49 * t88, -t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t53, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
