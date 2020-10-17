% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:27
% EndTime: 2019-12-31 22:22:31
% DurationCPUTime: 1.02s
% Computational Cost: add. (1239->108), mult. (2385->206), div. (0->0), fcn. (2730->8), ass. (0->75)
t58 = cos(qJ(3));
t46 = t58 * pkin(2);
t41 = t46 + pkin(3);
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t54 = sin(qJ(3));
t85 = t54 * pkin(2);
t72 = t57 * t85;
t29 = t53 * t41 + t72;
t27 = pkin(9) + t29;
t52 = sin(qJ(5));
t48 = t52 ^ 2;
t56 = cos(qJ(5));
t50 = t56 ^ 2;
t75 = t48 + t50;
t78 = t75 * t27;
t86 = t53 * pkin(3);
t39 = pkin(9) + t86;
t94 = t75 * t39;
t55 = sin(qJ(2));
t88 = -pkin(7) - pkin(6);
t70 = t88 * t55;
t59 = cos(qJ(2));
t71 = t88 * t59;
t20 = t54 * t70 - t58 * t71;
t31 = -t54 * t55 + t58 * t59;
t11 = t31 * pkin(8) + t20;
t19 = t54 * t71 + t58 * t70;
t32 = t54 * t59 + t58 * t55;
t64 = -t32 * pkin(8) + t19;
t6 = t53 * t11 - t57 * t64;
t93 = t6 ^ 2;
t16 = -t57 * t31 + t53 * t32;
t92 = t16 ^ 2;
t42 = -t59 * pkin(2) - pkin(1);
t24 = -t31 * pkin(3) + t42;
t91 = 0.2e1 * t24;
t90 = 0.2e1 * t32;
t89 = 0.2e1 * t59;
t87 = pkin(4) * t52;
t84 = t6 * t56;
t68 = -t57 * t41 + t53 * t85;
t26 = -pkin(4) + t68;
t83 = t26 * t56;
t45 = t57 * pkin(3);
t40 = -t45 - pkin(4);
t82 = t40 * t56;
t18 = t53 * t31 + t57 * t32;
t81 = t52 * t18;
t80 = t52 * t56;
t79 = t56 * t18;
t76 = pkin(9) * t75;
t49 = t55 ^ 2;
t51 = t59 ^ 2;
t74 = t49 + t51;
t73 = -0.2e1 * t18 * t16;
t67 = -pkin(4) * t18 - pkin(9) * t16;
t5 = t16 * pkin(4) - t18 * pkin(9) + t24;
t8 = t57 * t11 + t53 * t64;
t2 = t56 * t5 - t52 * t8;
t3 = t52 * t5 + t56 * t8;
t1 = -t2 * t52 + t3 * t56;
t66 = -t16 * t27 + t18 * t26;
t65 = -t16 * t39 + t18 * t40;
t47 = pkin(4) * t56;
t37 = 0.2e1 * t80;
t36 = t40 * t52;
t23 = t26 * t52;
t15 = t18 ^ 2;
t14 = t56 * t16;
t13 = t52 * t16;
t12 = t52 * t79;
t9 = (-t48 + t50) * t18;
t4 = t6 * t52;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, t55 * t89, 0, t51, 0, 0, pkin(1) * t89, -0.2e1 * pkin(1) * t55, 0.2e1 * t74 * pkin(6), t74 * pkin(6) ^ 2 + pkin(1) ^ 2, t32 ^ 2, t31 * t90, 0, t31 ^ 2, 0, 0, -0.2e1 * t42 * t31, t42 * t90, -0.2e1 * t19 * t32 + 0.2e1 * t20 * t31, t19 ^ 2 + t20 ^ 2 + t42 ^ 2, t15, t73, 0, t92, 0, 0, t16 * t91, t18 * t91, -0.2e1 * t8 * t16 + 0.2e1 * t6 * t18, t24 ^ 2 + t8 ^ 2 + t93, t50 * t15, -0.2e1 * t15 * t80, 0.2e1 * t16 * t79, t48 * t15, t52 * t73, t92, 0.2e1 * t2 * t16 + 0.2e1 * t6 * t81, -0.2e1 * t3 * t16 + 0.2e1 * t6 * t79, 0.2e1 * (-t2 * t56 - t3 * t52) * t18, t2 ^ 2 + t3 ^ 2 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t59, 0, -t55 * pkin(6), -t59 * pkin(6), 0, 0, 0, 0, t32, 0, t31, 0, t19, -t20, (t31 * t54 - t32 * t58) * pkin(2), (t19 * t58 + t20 * t54) * pkin(2), 0, 0, t18, 0, -t16, 0, -t6, -t8, -t29 * t16 + t18 * t68, t8 * t29 + t6 * t68, t12, t9, t13, -t12, t14, 0, t52 * t66 - t84, t56 * t66 + t4, t1, t1 * t27 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t85, 0, (t54 ^ 2 + t58 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t68, -0.2e1 * t29, 0, t29 ^ 2 + t68 ^ 2, t48, t37, 0, t50, 0, 0, -0.2e1 * t83, 0.2e1 * t23, 0.2e1 * t78, t75 * t27 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t31, 0, t19, -t20, 0, 0, 0, 0, t18, 0, -t16, 0, -t6, -t8, (-t16 * t53 - t18 * t57) * pkin(3), (t53 * t8 - t57 * t6) * pkin(3), t12, t9, t13, -t12, t14, 0, t52 * t65 - t84, t56 * t65 + t4, t1, t1 * t39 + t6 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t85, 0, 0, 0, 0, 0, 0, 0, 1, t45 - t68, -t72 + (-pkin(3) - t41) * t53, 0, (t29 * t53 - t57 * t68) * pkin(3), t48, t37, 0, t50, 0, 0, (-t26 - t40) * t56, t36 + t23, t94 + t78, t26 * t40 + t27 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t86, 0, (t53 ^ 2 + t57 ^ 2) * pkin(3) ^ 2, t48, t37, 0, t50, 0, 0, -0.2e1 * t82, 0.2e1 * t36, 0.2e1 * t94, t75 * t39 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, -t6, -t8, 0, 0, t12, t9, t13, -t12, t14, 0, t52 * t67 - t84, t56 * t67 + t4, t1, -t6 * pkin(4) + pkin(9) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t68, -t29, 0, 0, t48, t37, 0, t50, 0, 0, t47 - t83, t23 - t87, t76 + t78, -t26 * pkin(4) + pkin(9) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, -t86, 0, 0, t48, t37, 0, t50, 0, 0, t47 - t82, t36 - t87, t76 + t94, -t40 * pkin(4) + pkin(9) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t48, t37, 0, t50, 0, 0, 0.2e1 * t47, -0.2e1 * t87, 0.2e1 * t76, t75 * pkin(9) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t81, t16, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t56, 0, -t52 * t27, -t56 * t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t56, 0, -t52 * t39, -t56 * t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t56, 0, -t52 * pkin(9), -t56 * pkin(9), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;
