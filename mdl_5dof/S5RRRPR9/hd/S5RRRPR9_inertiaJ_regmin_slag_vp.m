% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR9
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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:55
% EndTime: 2021-01-15 23:23:59
% DurationCPUTime: 0.48s
% Computational Cost: add. (650->98), mult. (1379->197), div. (0->0), fcn. (1525->8), ass. (0->68)
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t34 = t51 * t57 + t52 * t54;
t55 = sin(qJ(2));
t28 = t34 * t55;
t65 = t57 * t55;
t68 = t54 * t55;
t29 = -t51 * t68 + t52 * t65;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t11 = -t28 * t53 + t29 * t56;
t78 = -0.2e1 * t11;
t33 = t51 * t54 - t52 * t57;
t45 = -pkin(3) * t57 - pkin(2);
t27 = pkin(4) * t33 + t45;
t77 = 0.2e1 * t27;
t76 = 0.2e1 * t45;
t75 = -0.2e1 * t55;
t58 = cos(qJ(2));
t74 = 0.2e1 * t58;
t73 = pkin(2) * t57;
t72 = pkin(6) * t54;
t71 = t51 * pkin(3);
t70 = t52 * pkin(3);
t69 = t58 * pkin(3);
t67 = t54 * t57;
t66 = t54 * t58;
t64 = t57 * t58;
t63 = qJ(4) + pkin(7);
t38 = -pkin(2) * t58 - pkin(7) * t55 - pkin(1);
t35 = t57 * t38;
t62 = qJ(4) * t55;
t18 = -t57 * t62 + t35 + (-pkin(3) - t72) * t58;
t60 = pkin(6) * t64;
t21 = t60 + (t38 - t62) * t54;
t9 = t51 * t18 + t52 * t21;
t46 = t55 * pkin(6);
t37 = pkin(3) * t68 + t46;
t61 = t55 * t74;
t8 = t52 * t18 - t21 * t51;
t4 = -pkin(4) * t58 - pkin(8) * t29 + t8;
t7 = -pkin(8) * t28 + t9;
t1 = t56 * t4 - t53 * t7;
t39 = t63 * t54;
t40 = t63 * t57;
t22 = -t52 * t39 - t40 * t51;
t2 = t4 * t53 + t56 * t7;
t23 = -t39 * t51 + t40 * t52;
t50 = t58 ^ 2;
t49 = t57 ^ 2;
t48 = t55 ^ 2;
t47 = t54 ^ 2;
t44 = pkin(4) + t70;
t31 = t44 * t53 + t56 * t71;
t30 = t44 * t56 - t53 * t71;
t26 = t38 * t54 + t60;
t25 = -pkin(6) * t66 + t35;
t19 = pkin(4) * t28 + t37;
t17 = -t33 * t53 + t34 * t56;
t16 = t56 * t33 + t34 * t53;
t13 = -pkin(8) * t33 + t23;
t12 = -pkin(8) * t34 + t22;
t10 = t56 * t28 + t29 * t53;
t6 = t12 * t53 + t13 * t56;
t5 = t12 * t56 - t13 * t53;
t3 = [1, 0, 0, t48, t61, 0, 0, 0, pkin(1) * t74, pkin(1) * t75, t49 * t48, -0.2e1 * t48 * t67, t64 * t75, t54 * t61, t50, -0.2e1 * t25 * t58 + 0.2e1 * t48 * t72, 0.2e1 * pkin(6) * t48 * t57 + 0.2e1 * t26 * t58, 0.2e1 * t28 * t37 - 0.2e1 * t58 * t8, 0.2e1 * t29 * t37 + 0.2e1 * t58 * t9, -0.2e1 * t28 * t9 - 0.2e1 * t29 * t8, t37 ^ 2 + t8 ^ 2 + t9 ^ 2, t11 ^ 2, t10 * t78, t58 * t78, t10 * t74, t50, -0.2e1 * t1 * t58 + 0.2e1 * t10 * t19, 0.2e1 * t11 * t19 + 0.2e1 * t2 * t58; 0, 0, 0, 0, 0, t55, t58, 0, -t46, -t58 * pkin(6), t54 * t65, (-t47 + t49) * t55, -t66, -t64, 0, -pkin(6) * t65 + (-pkin(2) * t55 + pkin(7) * t58) * t54, pkin(7) * t64 + (t72 - t73) * t55, -t22 * t58 + t28 * t45 + t33 * t37, t23 * t58 + t29 * t45 + t34 * t37, -t22 * t29 - t23 * t28 - t33 * t9 - t34 * t8, t22 * t8 + t23 * t9 + t37 * t45, t11 * t17, -t10 * t17 - t11 * t16, -t17 * t58, t16 * t58, 0, t10 * t27 + t16 * t19 - t5 * t58, t11 * t27 + t17 * t19 + t58 * t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t67, 0, 0, 0, 0.2e1 * t73, -0.2e1 * pkin(2) * t54, t33 * t76, t34 * t76, -0.2e1 * t22 * t34 - 0.2e1 * t23 * t33, t22 ^ 2 + t23 ^ 2 + t45 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t77, t17 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t68, -t58, t25, -t26, -t52 * t69 + t8, t51 * t69 - t9, (-t28 * t51 - t29 * t52) * pkin(3), (t51 * t9 + t52 * t8) * pkin(3), 0, 0, t11, -t10, -t58, -t30 * t58 + t1, t31 * t58 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t57, 0, -t54 * pkin(7), -t57 * pkin(7), t22, -t23, (-t33 * t51 - t34 * t52) * pkin(3), (t22 * t52 + t23 * t51) * pkin(3), 0, 0, t17, -t16, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t70, -0.2e1 * t71, 0, (t51 ^ 2 + t52 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29, 0, t37, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t34, 0, t45, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, -t58, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
