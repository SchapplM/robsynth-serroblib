% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:54
% EndTime: 2020-01-03 12:05:56
% DurationCPUTime: 0.38s
% Computational Cost: add. (365->84), mult. (778->117), div. (0->0), fcn. (788->8), ass. (0->78)
t59 = sin(pkin(9));
t93 = -0.2e1 * t59;
t92 = 0.2e1 * t59;
t60 = cos(pkin(9));
t91 = -0.2e1 * t60;
t90 = 0.2e1 * t60;
t89 = pkin(8) * t59;
t61 = sin(qJ(5));
t88 = t61 * pkin(4);
t87 = sin(qJ(2)) * pkin(1);
t64 = cos(qJ(5));
t86 = t64 * pkin(4);
t85 = cos(qJ(2)) * pkin(1);
t54 = -pkin(2) - t85;
t84 = pkin(2) - t54;
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t35 = -t61 * t62 + t64 * t65;
t28 = t35 * t59;
t39 = -t60 * pkin(3) - t59 * pkin(7) - pkin(2);
t32 = t39 - t85;
t29 = t65 * t32;
t53 = qJ(3) + t87;
t52 = t65 * t59;
t70 = pkin(8) * t52;
t8 = -t70 + t29 + (-t53 * t62 - pkin(4)) * t60;
t78 = t65 * t60;
t69 = t53 * t78;
t9 = t69 + (t32 - t89) * t62;
t3 = t61 * t8 + t64 * t9;
t79 = t62 * t59;
t48 = pkin(4) * t79;
t30 = t59 * t53 + t48;
t83 = t30 * t28 + t3 * t60;
t34 = t59 * qJ(3) + t48;
t33 = t65 * t39;
t72 = qJ(3) * t62;
t15 = -t70 + t33 + (-pkin(4) - t72) * t60;
t68 = qJ(3) * t78;
t20 = t68 + (t39 - t89) * t62;
t6 = t61 * t15 + t64 * t20;
t82 = t34 * t28 + t6 * t60;
t81 = t35 * t60;
t57 = t59 ^ 2;
t42 = t57 * t53;
t80 = t57 * t65;
t50 = t62 * t60;
t19 = t62 * t32 + t69;
t77 = t19 * t60 + t53 * t80;
t58 = t60 ^ 2;
t76 = t58 * t53 + t42;
t25 = t62 * t39 + t68;
t55 = t57 * qJ(3);
t75 = t25 * t60 + t65 * t55;
t74 = t58 * qJ(3) + t55;
t73 = t57 + t58;
t71 = t60 * t86;
t2 = -t61 * t9 + t64 * t8;
t5 = t64 * t15 - t61 * t20;
t36 = t61 * t65 + t64 * t62;
t51 = t65 ^ 2 * t57;
t47 = t60 * t88;
t45 = t62 * t55;
t44 = -0.2e1 * t62 * t80;
t41 = t78 * t93;
t40 = t50 * t92;
t37 = t62 * t42;
t31 = t36 * t60;
t27 = t36 * t59;
t26 = t28 ^ 2;
t24 = -t60 * t72 + t33;
t23 = t28 * t91;
t22 = t27 * t90;
t18 = -t53 * t50 + t29;
t16 = t34 * t27;
t12 = t30 * t27;
t10 = -0.2e1 * t28 * t27;
t1 = [1, 0, 0, 1, 0.2e1 * t85, -0.2e1 * t87, t54 * t91, t54 * t92, 0.2e1 * t76, t73 * t53 ^ 2 + t54 ^ 2, t51, t44, t41, t40, t58, -0.2e1 * t18 * t60 + 0.2e1 * t37, 0.2e1 * t77, t26, t10, t23, t22, t58, -0.2e1 * t2 * t60 + 0.2e1 * t12, 0.2e1 * t83; 0, 0, 0, 1, t85, -t87, t84 * t60, -t84 * t59, t74 + t76, t73 * t53 * qJ(3) - t54 * pkin(2), t51, t44, t41, t40, t58, t37 + t45 + (-t18 - t24) * t60, t75 + t77, t26, t10, t23, t22, t58, t12 + t16 + (-t2 - t5) * t60, t82 + t83; 0, 0, 0, 1, 0, 0, pkin(2) * t90, pkin(2) * t93, 0.2e1 * t74, t73 * qJ(3) ^ 2 + pkin(2) ^ 2, t51, t44, t41, t40, t58, -0.2e1 * t24 * t60 + 0.2e1 * t45, 0.2e1 * t75, t26, t10, t23, t22, t58, -0.2e1 * t5 * t60 + 0.2e1 * t16, 0.2e1 * t82; 0, 0, 0, 0, 0, 0, -t60, t59, 0, t54, 0, 0, 0, 0, 0, -t78, t50, 0, 0, 0, 0, 0, -t81, t31; 0, 0, 0, 0, 0, 0, -t60, t59, 0, -pkin(2), 0, 0, 0, 0, 0, -t78, t50, 0, 0, 0, 0, 0, -t81, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t79, -t60, t18, -t19, 0, 0, t28, -t27, -t60, t2 - t71, -t3 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t79, -t60, t24, -t25, 0, 0, t28, -t27, -t60, t5 - t71, t47 - t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, 0, 0, 0, 0, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t86, -0.2e1 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, -t60, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, -t60, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t86, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
