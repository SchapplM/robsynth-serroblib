% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = sin(pkin(10));
t54 = cos(pkin(10));
t43 = t52 ^ 2 + t54 ^ 2;
t53 = sin(pkin(9));
t55 = cos(pkin(9));
t57 = sin(qJ(2));
t78 = cos(qJ(2));
t36 = t53 * t78 + t55 * t57;
t84 = -0.2e1 * t36;
t56 = sin(qJ(6));
t58 = cos(qJ(6));
t37 = t58 * t52 - t56 * t54;
t68 = t78 * pkin(7);
t42 = t78 * qJ(3) + t68;
t67 = (-qJ(3) - pkin(7)) * t57;
t21 = t53 * t42 - t55 * t67;
t83 = t21 ^ 2;
t35 = t56 * t52 + t58 * t54;
t12 = t35 * t36;
t82 = -0.2e1 * t12;
t48 = -t55 * pkin(2) - pkin(3);
t60 = t52 * qJ(5) - t48;
t79 = pkin(4) + pkin(5);
t25 = t79 * t54 + t60;
t81 = 0.2e1 * t25;
t80 = -0.2e1 * t54;
t33 = t53 * t57 - t55 * t78;
t45 = t53 * pkin(2) + qJ(4);
t77 = t33 * t45;
t76 = t37 * t33;
t75 = t52 * t33;
t27 = t52 * t36;
t74 = t54 * t33;
t28 = t54 * t36;
t49 = -t78 * pkin(2) - pkin(1);
t16 = t33 * pkin(3) - t36 * qJ(4) + t49;
t23 = t55 * t42 + t53 * t67;
t9 = t52 * t16 + t54 * t23;
t71 = t43 * t45 ^ 2;
t70 = qJ(5) * t54;
t69 = 0.2e1 * t78;
t5 = t33 * qJ(5) + t9;
t18 = t52 * t23;
t8 = t54 * t16 - t18;
t6 = -t33 * pkin(4) - t8;
t66 = t5 * t54 + t6 * t52;
t65 = t5 * t52 - t6 * t54;
t64 = t9 * t52 + t8 * t54;
t63 = -t8 * t52 + t9 * t54;
t29 = -t54 * pkin(4) - t60;
t62 = t29 * t36 - t77;
t61 = t36 * t48 - t77;
t40 = t52 * t45;
t31 = (-pkin(8) + t45) * t54;
t30 = -t52 * pkin(8) + t40;
t24 = 0.2e1 * t43 * t45;
t20 = t35 * t33;
t17 = t43 * t36;
t15 = t56 * t30 + t58 * t31;
t14 = t58 * t30 - t56 * t31;
t11 = t37 * t36;
t10 = (pkin(4) * t52 - t70) * t36 + t21;
t7 = (-t79 * t52 + t70) * t36 - t21;
t4 = pkin(8) * t27 + t5;
t3 = t18 + (-pkin(8) * t36 - t16) * t54 - t79 * t33;
t2 = t56 * t3 + t58 * t4;
t1 = t58 * t3 - t56 * t4;
t13 = [1, 0, 0, t57 ^ 2, t57 * t69, 0, 0, 0, pkin(1) * t69, -0.2e1 * pkin(1) * t57, 0.2e1 * t21 * t36 - 0.2e1 * t23 * t33, t23 ^ 2 + t49 ^ 2 + t83, 0.2e1 * t21 * t27 + 0.2e1 * t8 * t33, 0.2e1 * t21 * t28 - 0.2e1 * t9 * t33, t64 * t84, t8 ^ 2 + t9 ^ 2 + t83, 0.2e1 * t10 * t27 - 0.2e1 * t6 * t33, t65 * t84, -0.2e1 * t10 * t28 + 0.2e1 * t5 * t33, t10 ^ 2 + t5 ^ 2 + t6 ^ 2, t12 ^ 2, -t11 * t82, t33 * t82, -0.2e1 * t11 * t33, t33 ^ 2, -0.2e1 * t1 * t33 - 0.2e1 * t7 * t11, 0.2e1 * t7 * t12 + 0.2e1 * t2 * t33; 0, 0, 0, 0, 0, t57, t78, 0, -t57 * pkin(7), -t68 (-t33 * t53 - t36 * t55) * pkin(2) (-t21 * t55 + t23 * t53) * pkin(2), -t21 * t54 + t61 * t52, t21 * t52 + t61 * t54, t63, t21 * t48 + t63 * t45, -t10 * t54 + t62 * t52, t66, -t10 * t52 - t62 * t54, t10 * t29 + t66 * t45, t12 * t37, t37 * t11 - t12 * t35, -t76, t20, 0, -t25 * t11 - t14 * t33 + t7 * t35, t25 * t12 + t15 * t33 + t7 * t37; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t53 ^ 2 + t55 ^ 2) * pkin(2) ^ 2, t48 * t80, 0.2e1 * t48 * t52, t24, t48 ^ 2 + t71, t29 * t80, t24, -0.2e1 * t29 * t52, t29 ^ 2 + t71, t37 ^ 2, -0.2e1 * t37 * t35, 0, 0, 0, t35 * t81, t37 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t74, -t75, -t17, t64, t74, -t17, t75, t65, 0, 0, 0, 0, 0, t20, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t43, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28, 0, t21, t27, 0, -t28, t10, 0, 0, 0, 0, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t52, 0, t48, -t54, 0, -t52, t29, 0, 0, 0, 0, 0, -t35, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t28, 0, t6, 0, 0, 0, 0, 0, -t58 * t33, t56 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, -t33, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
