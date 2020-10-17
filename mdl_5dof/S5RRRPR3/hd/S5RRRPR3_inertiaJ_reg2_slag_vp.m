% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:29
% EndTime: 2020-01-03 12:09:32
% DurationCPUTime: 0.75s
% Computational Cost: add. (801->83), mult. (1495->151), div. (0->0), fcn. (1670->8), ass. (0->73)
t66 = sin(qJ(2));
t84 = t66 * pkin(1);
t54 = pkin(7) + t84;
t65 = sin(qJ(3));
t60 = t65 ^ 2;
t67 = cos(qJ(3));
t61 = t67 ^ 2;
t72 = t60 + t61;
t74 = t72 * t54;
t62 = sin(pkin(9));
t63 = cos(pkin(9));
t41 = t62 * t65 - t63 * t67;
t56 = -t67 * pkin(3) - pkin(2);
t31 = t41 * pkin(4) + t56;
t68 = cos(qJ(2));
t82 = t68 * pkin(1);
t30 = t31 - t82;
t94 = 0.2e1 * t30;
t93 = 0.2e1 * t31;
t45 = t56 - t82;
t92 = 0.2e1 * t45;
t91 = 0.2e1 * t56;
t90 = 0.2e1 * t67;
t43 = t62 * t67 + t63 * t65;
t64 = sin(qJ(5));
t80 = cos(qJ(5));
t24 = t80 * t41 + t64 * t43;
t26 = -t64 * t41 + t80 * t43;
t37 = (-qJ(4) - t54) * t65;
t57 = t67 * qJ(4);
t79 = t67 * t54;
t38 = t57 + t79;
t19 = t63 * t37 - t62 * t38;
t87 = t43 * pkin(8);
t11 = t19 - t87;
t20 = t62 * t37 + t63 * t38;
t36 = t41 * pkin(8);
t12 = t20 - t36;
t5 = t80 * t11 - t64 * t12;
t6 = t64 * t11 + t80 * t12;
t89 = -t6 * t24 - t5 * t26;
t46 = (-qJ(4) - pkin(7)) * t65;
t83 = t67 * pkin(7);
t47 = t57 + t83;
t28 = t63 * t46 - t62 * t47;
t15 = t28 - t87;
t29 = t62 * t46 + t63 * t47;
t16 = t29 - t36;
t8 = t80 * t15 - t64 * t16;
t9 = t64 * t15 + t80 * t16;
t88 = -t9 * t24 - t8 * t26;
t86 = t62 * pkin(3);
t85 = t63 * pkin(3);
t55 = -pkin(2) - t82;
t81 = pkin(2) - t55;
t78 = -t19 * t43 - t20 * t41;
t77 = -t28 * t43 - t29 * t41;
t76 = t30 + t31;
t75 = t45 + t56;
t73 = t72 * pkin(7);
t53 = pkin(4) + t85;
t50 = t65 * t90;
t40 = t43 ^ 2;
t39 = t41 ^ 2;
t34 = t64 * t53 + t80 * t86;
t33 = t80 * t53 - t64 * t86;
t27 = -0.2e1 * t43 * t41;
t23 = t26 ^ 2;
t22 = t24 ^ 2;
t21 = (-t41 * t62 - t43 * t63) * pkin(3);
t10 = -0.2e1 * t26 * t24;
t7 = -t34 * t24 - t33 * t26;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t82, -0.2e1 * t84, 0, (t66 ^ 2 + t68 ^ 2) * pkin(1) ^ 2, t60, t50, 0, t61, 0, 0, -0.2e1 * t55 * t67, 0.2e1 * t55 * t65, 0.2e1 * t74, t72 * t54 ^ 2 + t55 ^ 2, t40, t27, 0, t39, 0, 0, t41 * t92, t43 * t92, 0.2e1 * t78, t19 ^ 2 + t20 ^ 2 + t45 ^ 2, t23, t10, 0, t22, 0, 0, t24 * t94, t26 * t94, 0.2e1 * t89, t30 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t82, -t84, 0, 0, t60, t50, 0, t61, 0, 0, t81 * t67, -t81 * t65, t73 + t74, -t55 * pkin(2) + pkin(7) * t74, t40, t27, 0, t39, 0, 0, t75 * t41, t75 * t43, t77 + t78, t19 * t28 + t20 * t29 + t45 * t56, t23, t10, 0, t22, 0, 0, t76 * t24, t76 * t26, t88 + t89, t30 * t31 + t5 * t8 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t60, t50, 0, t61, 0, 0, pkin(2) * t90, -0.2e1 * pkin(2) * t65, 0.2e1 * t73, t72 * pkin(7) ^ 2 + pkin(2) ^ 2, t40, t27, 0, t39, 0, 0, t41 * t91, t43 * t91, 0.2e1 * t77, t28 ^ 2 + t29 ^ 2 + t56 ^ 2, t23, t10, 0, t22, 0, 0, t24 * t93, t26 * t93, 0.2e1 * t88, t31 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t67, 0, -t65 * t54, -t79, 0, 0, 0, 0, t43, 0, -t41, 0, t19, -t20, t21, (t19 * t63 + t20 * t62) * pkin(3), 0, 0, t26, 0, -t24, 0, t5, -t6, t7, t5 * t33 + t6 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t67, 0, -t65 * pkin(7), -t83, 0, 0, 0, 0, t43, 0, -t41, 0, t28, -t29, t21, (t28 * t63 + t29 * t62) * pkin(3), 0, 0, t26, 0, -t24, 0, t8, -t9, t7, t8 * t33 + t9 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t85, -0.2e1 * t86, 0, (t62 ^ 2 + t63 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t34, 0, t33 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t43, 0, t45, 0, 0, 0, 0, 0, 0, t24, t26, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t43, 0, t56, 0, 0, 0, 0, 0, 0, t24, t26, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
