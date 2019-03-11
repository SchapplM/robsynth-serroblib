% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t50 = sin(pkin(9));
t55 = cos(qJ(4));
t51 = cos(pkin(9));
t80 = sin(qJ(4));
t68 = t80 * t51;
t31 = t55 * t50 + t68;
t28 = t31 ^ 2;
t75 = t55 * t51;
t33 = -t80 * t50 + t75;
t29 = t33 ^ 2;
t69 = -t29 - t28;
t53 = sin(qJ(5));
t48 = t53 ^ 2;
t54 = cos(qJ(5));
t49 = t54 ^ 2;
t39 = t48 + t49;
t93 = -0.2e1 * t33;
t92 = t69 * t54;
t91 = t69 * t53;
t52 = -pkin(1) - qJ(3);
t81 = -pkin(7) + t52;
t35 = t81 * t50;
t14 = t80 * t35 - t81 * t75;
t90 = t14 ^ 2;
t40 = t50 * pkin(3) + qJ(2);
t89 = 0.2e1 * t40;
t88 = -0.2e1 * t53;
t87 = 0.2e1 * qJ(2);
t86 = pkin(8) * t31;
t85 = t31 * pkin(5);
t84 = t33 * pkin(4);
t83 = t53 * pkin(8);
t82 = t54 * pkin(8);
t16 = t55 * t35 + t81 * t68;
t9 = t31 * pkin(4) - t33 * pkin(8) + t40;
t4 = t54 * t16 + t53 * t9;
t79 = t14 * t33;
t60 = t54 * pkin(5) + t53 * qJ(6);
t37 = -pkin(4) - t60;
t78 = t33 * t37;
t23 = t53 * t31;
t77 = t53 * t33;
t76 = t53 * t54;
t25 = t54 * t31;
t26 = t54 * t33;
t74 = t39 * t86;
t73 = t39 * pkin(8) ^ 2;
t45 = t50 ^ 2;
t46 = t51 ^ 2;
t38 = t45 + t46;
t72 = t31 * qJ(6);
t71 = t31 * t77;
t70 = t29 * t76;
t67 = t53 * t16 - t54 * t9;
t66 = -t84 - t86;
t1 = t72 + t4;
t2 = t67 - t85;
t65 = t1 * t54 + t2 * t53;
t64 = t1 * t53 - t2 * t54;
t63 = t4 * t53 - t54 * t67;
t62 = t4 * t54 + t53 * t67;
t61 = -t78 + t86;
t59 = pkin(5) * t53 - t54 * qJ(6);
t58 = -t16 * t31 + t79;
t56 = qJ(2) ^ 2;
t36 = 0.2e1 * t39 * pkin(8);
t30 = t38 * t52;
t24 = t49 * t29;
t22 = t48 * t29;
t19 = t53 * t26;
t17 = 0.2e1 * t31 * t26;
t13 = t39 * t31;
t12 = t39 * t33;
t11 = (t48 - t49) * t33;
t6 = t39 * t28 + t29;
t5 = t59 * t33 + t14;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t87 (pkin(1) ^ 2) + t56, t46, -0.2e1 * t51 * t50, 0, t45, 0, 0, t50 * t87, t51 * t87, -0.2e1 * t30, t38 * t52 ^ 2 + t56, t29, t31 * t93, 0, t28, 0, 0, t31 * t89, t33 * t89, 0.2e1 * t58, t16 ^ 2 + t40 ^ 2 + t90, t24, -0.2e1 * t70, t17, t22, -0.2e1 * t71, t28, 0.2e1 * t14 * t77 - 0.2e1 * t31 * t67, 0.2e1 * t14 * t26 - 0.2e1 * t4 * t31, t63 * t93, t4 ^ 2 + t67 ^ 2 + t90, t24, t17, 0.2e1 * t70, t28, 0.2e1 * t71, t22, -0.2e1 * t2 * t31 + 0.2e1 * t5 * t77, t64 * t93, 0.2e1 * t1 * t31 - 0.2e1 * t26 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t38, t30, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t58, 0, 0, 0, 0, 0, 0, t91, t92, 0, t31 * t62 - t79, 0, 0, 0, 0, 0, 0, t91, 0, -t92, t31 * t65 - t5 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, 0, qJ(2), 0, 0, 0, 0, 0, 0, t31, t33, 0, t40, 0, 0, 0, 0, 0, 0, t25, -t23, -t12, t63, 0, 0, 0, 0, 0, 0, t25, -t12, t23, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, -t14, -t16, 0, 0, t19, -t11, t23, -t19, t25, 0, -t14 * t54 + t66 * t53, t14 * t53 + t54 * t66, t62, -t14 * pkin(4) + pkin(8) * t62, t19, t23, t11, 0, -t25, -t19, -t5 * t54 - t53 * t61, t65, -t5 * t53 + t54 * t61, pkin(8) * t65 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t77, t13, t74 + t84, 0, 0, 0, 0, 0, 0, t26, t13, t77, t74 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t48, 0.2e1 * t76, 0, t49, 0, 0, 0.2e1 * pkin(4) * t54, pkin(4) * t88, t36, pkin(4) ^ 2 + t73, t48, 0, -0.2e1 * t76, 0, 0, t49, -0.2e1 * t37 * t54, t36, t37 * t88, t37 ^ 2 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t77, t31, -t67, -t4, 0, 0, 0, t26, 0, t31, t77, 0, -t67 + 0.2e1 * t85, -t60 * t33, 0.2e1 * t72 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, t25, -t59 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t53, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t53, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t54, 0, -t83, -t82, 0, 0, 0, t53, 0, 0, -t54, 0, -t83, -t59, t82, -t59 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
