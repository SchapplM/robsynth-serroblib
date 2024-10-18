% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 0.75s
% Computational Cost: add. (649->103), mult. (1774->177), div. (0->0), fcn. (1975->12), ass. (0->66)
t48 = sin(pkin(6));
t49 = sin(pkin(5));
t54 = sin(qJ(3));
t55 = sin(qJ(2));
t58 = cos(qJ(3));
t59 = cos(qJ(2));
t21 = (-t54 * t55 + t58 * t59) * t49;
t22 = (t54 * t59 + t55 * t58) * t49;
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t11 = t57 * t21 - t53 * t22;
t50 = cos(pkin(6));
t51 = cos(pkin(5));
t6 = -t48 * t11 + t51 * t50;
t75 = t48 * t6;
t74 = t53 * pkin(3);
t73 = t54 * pkin(2);
t43 = t58 * pkin(2);
t40 = t43 + pkin(3);
t24 = t57 * t40 - t53 * t73;
t23 = pkin(4) + t24;
t44 = t48 ^ 2;
t56 = cos(qJ(5));
t69 = t44 * t56;
t63 = t57 * t73;
t25 = t53 * t40 + t63;
t41 = t48 * pkin(10);
t19 = t25 + t41;
t52 = sin(qJ(5));
t67 = t50 * t56;
t9 = -t52 * t19 + t23 * t67;
t72 = t23 * t69 + t9 * t50;
t42 = t57 * pkin(3);
t39 = t42 + pkin(4);
t71 = t39 * t44;
t70 = t44 * t52;
t37 = t48 * t52;
t38 = t48 * t56;
t68 = t50 * t52;
t32 = t41 + t74;
t15 = -t52 * t32 + t39 * t67;
t66 = t15 * t50 + t39 * t69;
t26 = pkin(4) * t67 - pkin(10) * t37;
t65 = pkin(4) * t69 + t26 * t50;
t64 = 0.2e1 * t48 * t50;
t62 = t11 * t50 + t48 * t51;
t47 = t51 ^ 2;
t46 = t50 ^ 2;
t36 = t44 * t56 ^ 2;
t35 = t44 * t52 ^ 2;
t31 = 0.2e1 * t52 * t69;
t30 = t56 * t64;
t29 = t52 * t64;
t27 = pkin(4) * t68 + pkin(10) * t38;
t18 = t27 * t38;
t16 = t56 * t32 + t39 * t68;
t13 = t16 * t38;
t12 = t53 * t21 + t57 * t22;
t10 = t56 * t19 + t23 * t68;
t7 = t10 * t38;
t5 = t56 * t12 + t62 * t52;
t4 = -t52 * t12 + t62 * t56;
t3 = t6 * t37 - t5 * t50;
t2 = -t6 * t38 + t4 * t50;
t1 = (-t4 * t52 + t5 * t56) * t48;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 + (t55 ^ 2 + t59 ^ 2) * t49 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 ^ 2 + t22 ^ 2 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t12 ^ 2 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t59, -t49 * t55, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, (t21 * t58 + t22 * t54) * pkin(2), 0, 0, 0, 0, 0, 0, t11, -t12, 0, t11 * t24 + t12 * t25, 0, 0, 0, 0, 0, 0, t2, t3, t1, t5 * t10 - t23 * t75 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t73, 0, (t54 ^ 2 + t58 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, -0.2e1 * t25, 0, t24 ^ 2 + t25 ^ 2, t35, t31, t29, t36, t30, t46, 0.2e1 * t72, -0.2e1 * t10 * t50 - 0.2e1 * t23 * t70, -0.2e1 * t9 * t37 + 0.2e1 * t7, t44 * t23 ^ 2 + t10 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, 0, (t11 * t57 + t12 * t53) * pkin(3), 0, 0, 0, 0, 0, 0, t2, t3, t1, t4 * t15 + t5 * t16 - t39 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t73, 0, 0, 0, 0, 0, 0, 0, 1, t24 + t42, -t63 + (-pkin(3) - t40) * t53, 0, (t24 * t57 + t25 * t53) * pkin(3), t35, t31, t29, t36, t30, t46, t66 + t72, (-t10 - t16) * t50 + (-t23 - t39) * t70, t13 + t7 + (-t15 - t9) * t37, t10 * t16 + t9 * t15 + t23 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t74, 0, (t53 ^ 2 + t57 ^ 2) * pkin(3) ^ 2, t35, t31, t29, t36, t30, t46, 0.2e1 * t66, -0.2e1 * t16 * t50 - 0.2e1 * t39 * t70, -0.2e1 * t15 * t37 + 0.2e1 * t13, t44 * t39 ^ 2 + t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, 0, 0, 0, 0, 0, 0, 0, 0, t2, t3, t1, -pkin(4) * t75 + t4 * t26 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, -t25, 0, 0, t35, t31, t29, t36, t30, t46, t65 + t72, (-t10 - t27) * t50 + (-pkin(4) - t23) * t70, t18 + t7 + (-t26 - t9) * t37, t44 * t23 * pkin(4) + t10 * t27 + t9 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t42, -t74, 0, 0, t35, t31, t29, t36, t30, t46, t65 + t66, (-t16 - t27) * t50 + (-pkin(4) - t39) * t70, t13 + t18 + (-t15 - t26) * t37, pkin(4) * t71 + t15 * t26 + t16 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t35, t31, t29, t36, t30, t46, 0.2e1 * t65, -0.2e1 * pkin(4) * t70 - 0.2e1 * t27 * t50, -0.2e1 * t26 * t37 + 0.2e1 * t18, t44 * pkin(4) ^ 2 + t26 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t38, t50, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t38, t50, t15, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t38, t50, t26, -t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
