% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR6
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:19
% EndTime: 2020-01-03 12:15:22
% DurationCPUTime: 0.87s
% Computational Cost: add. (941->101), mult. (1763->172), div. (0->0), fcn. (1952->8), ass. (0->77)
t69 = sin(qJ(2));
t90 = t69 * pkin(1);
t55 = pkin(7) + t90;
t68 = sin(qJ(3));
t64 = t68 ^ 2;
t72 = cos(qJ(3));
t65 = t72 ^ 2;
t79 = t64 + t65;
t81 = t79 * t55;
t67 = sin(qJ(4));
t71 = cos(qJ(4));
t42 = t67 * t68 - t71 * t72;
t58 = -t72 * pkin(3) - pkin(2);
t32 = t42 * pkin(4) + t58;
t73 = cos(qJ(2));
t88 = t73 * pkin(1);
t31 = t32 - t88;
t100 = 0.2e1 * t31;
t99 = 0.2e1 * t32;
t46 = t58 - t88;
t98 = 0.2e1 * t46;
t97 = 0.2e1 * t58;
t96 = 0.2e1 * t72;
t44 = t67 * t72 + t71 * t68;
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t25 = t70 * t42 + t66 * t44;
t27 = -t66 * t42 + t70 * t44;
t38 = (-pkin(8) - t55) * t68;
t63 = t72 * pkin(8);
t86 = t72 * t55;
t39 = t63 + t86;
t20 = t71 * t38 - t67 * t39;
t93 = t44 * pkin(9);
t12 = t20 - t93;
t21 = t67 * t38 + t71 * t39;
t37 = t42 * pkin(9);
t13 = t21 - t37;
t5 = t70 * t12 - t66 * t13;
t6 = t66 * t12 + t70 * t13;
t95 = -t6 * t25 - t5 * t27;
t47 = (-pkin(8) - pkin(7)) * t68;
t89 = t72 * pkin(7);
t48 = t63 + t89;
t29 = t71 * t47 - t67 * t48;
t16 = t29 - t93;
t30 = t67 * t47 + t71 * t48;
t17 = t30 - t37;
t8 = t70 * t16 - t66 * t17;
t9 = t66 * t16 + t70 * t17;
t94 = -t9 * t25 - t8 * t27;
t92 = t66 * pkin(4);
t91 = t67 * pkin(3);
t57 = -pkin(2) - t88;
t87 = pkin(2) - t57;
t85 = -t20 * t44 - t21 * t42;
t84 = -t29 * t44 - t30 * t42;
t83 = t31 + t32;
t82 = t46 + t58;
t80 = t79 * pkin(7);
t78 = t70 * t91;
t61 = t71 * pkin(3);
t56 = t61 + pkin(4);
t34 = t70 * t56 - t66 * t91;
t60 = t70 * pkin(4);
t52 = t68 * t96;
t41 = t44 ^ 2;
t40 = t42 ^ 2;
t35 = t66 * t56 + t78;
t28 = -0.2e1 * t44 * t42;
t24 = t27 ^ 2;
t23 = t25 ^ 2;
t22 = (-t42 * t67 - t44 * t71) * pkin(3);
t11 = -0.2e1 * t27 * t25;
t10 = (-t25 * t66 - t27 * t70) * pkin(4);
t7 = -t35 * t25 - t34 * t27;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t88, -0.2e1 * t90, 0, (t69 ^ 2 + t73 ^ 2) * pkin(1) ^ 2, t64, t52, 0, t65, 0, 0, -0.2e1 * t57 * t72, 0.2e1 * t57 * t68, 0.2e1 * t81, t79 * t55 ^ 2 + t57 ^ 2, t41, t28, 0, t40, 0, 0, t42 * t98, t44 * t98, 0.2e1 * t85, t20 ^ 2 + t21 ^ 2 + t46 ^ 2, t24, t11, 0, t23, 0, 0, t25 * t100, t27 * t100, 0.2e1 * t95, t31 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t88, -t90, 0, 0, t64, t52, 0, t65, 0, 0, t87 * t72, -t87 * t68, t80 + t81, -t57 * pkin(2) + pkin(7) * t81, t41, t28, 0, t40, 0, 0, t82 * t42, t82 * t44, t84 + t85, t20 * t29 + t21 * t30 + t46 * t58, t24, t11, 0, t23, 0, 0, t83 * t25, t83 * t27, t94 + t95, t31 * t32 + t5 * t8 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t64, t52, 0, t65, 0, 0, pkin(2) * t96, -0.2e1 * pkin(2) * t68, 0.2e1 * t80, t79 * pkin(7) ^ 2 + pkin(2) ^ 2, t41, t28, 0, t40, 0, 0, t42 * t97, t44 * t97, 0.2e1 * t84, t29 ^ 2 + t30 ^ 2 + t58 ^ 2, t24, t11, 0, t23, 0, 0, t25 * t99, t27 * t99, 0.2e1 * t94, t32 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t72, 0, -t68 * t55, -t86, 0, 0, 0, 0, t44, 0, -t42, 0, t20, -t21, t22, (t20 * t71 + t21 * t67) * pkin(3), 0, 0, t27, 0, -t25, 0, t5, -t6, t7, t5 * t34 + t6 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t72, 0, -t68 * pkin(7), -t89, 0, 0, 0, 0, t44, 0, -t42, 0, t29, -t30, t22, (t29 * t71 + t30 * t67) * pkin(3), 0, 0, t27, 0, -t25, 0, t8, -t9, t7, t8 * t34 + t9 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t91, 0, (t67 ^ 2 + t71 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35, 0, t34 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t42, 0, t20, -t21, 0, 0, 0, 0, t27, 0, -t25, 0, t5, -t6, t10, (t5 * t70 + t6 * t66) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t42, 0, t29, -t30, 0, 0, 0, 0, t27, 0, -t25, 0, t8, -t9, t10, (t66 * t9 + t70 * t8) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t61, -t91, 0, 0, 0, 0, 0, 0, 0, 1, t34 + t60, -t78 + (-pkin(4) - t56) * t66, 0, (t34 * t70 + t35 * t66) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t92, 0, (t66 ^ 2 + t70 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t60, -t92, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
