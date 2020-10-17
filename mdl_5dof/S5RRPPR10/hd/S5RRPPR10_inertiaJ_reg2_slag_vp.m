% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:41
% EndTime: 2019-12-31 19:44:44
% DurationCPUTime: 0.84s
% Computational Cost: add. (474->103), mult. (1004->200), div. (0->0), fcn. (1011->6), ass. (0->68)
t49 = sin(pkin(8));
t45 = t49 ^ 2;
t50 = cos(pkin(8));
t46 = t50 ^ 2;
t79 = t45 + t46;
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t22 = t53 * t49 - t51 * t50;
t52 = sin(qJ(2));
t14 = t22 * t52;
t78 = 0.2e1 * t14;
t59 = t49 * qJ(4) + pkin(2);
t73 = pkin(3) + pkin(4);
t17 = t73 * t50 + t59;
t77 = 0.2e1 * t17;
t76 = -0.2e1 * t50;
t75 = 0.2e1 * t52;
t54 = cos(qJ(2));
t74 = 0.2e1 * t54;
t72 = pkin(2) * t49;
t47 = t52 ^ 2;
t71 = t47 * pkin(6);
t70 = t52 * pkin(6);
t69 = t49 * t50;
t36 = t49 * t52;
t68 = t49 * t54;
t39 = t50 * t52;
t67 = t50 * t54;
t65 = t52 * t54;
t25 = -t54 * pkin(2) - t52 * qJ(3) - pkin(1);
t12 = pkin(6) * t67 + t49 * t25;
t63 = t79 * qJ(3) ^ 2;
t62 = qJ(3) * t54;
t41 = t49 * qJ(3);
t61 = t47 * t69;
t60 = t49 * t65;
t33 = pkin(6) * t68;
t11 = t50 * t25 - t33;
t9 = -t54 * qJ(4) + t12;
t44 = t54 * pkin(3);
t10 = -t11 + t44;
t58 = t10 * t49 + t9 * t50;
t57 = -t11 * t49 + t12 * t50;
t20 = t51 * t49 + t53 * t50;
t56 = pkin(6) ^ 2;
t48 = t54 ^ 2;
t43 = t47 * t56;
t38 = t46 * t47;
t35 = t45 * t47;
t32 = qJ(4) * t39;
t31 = t54 * t41;
t29 = t49 * t39;
t28 = t65 * t76;
t27 = (-pkin(7) + qJ(3)) * t50;
t26 = -t49 * pkin(7) + t41;
t24 = -t50 * pkin(3) - t59;
t23 = 0.2e1 * t79 * qJ(3);
t19 = (t45 - t46) * t52;
t16 = t20 * t52;
t13 = -t32 + (pkin(3) * t49 + pkin(6)) * t52;
t7 = -t32 - (-t73 * t49 - pkin(6)) * t52;
t6 = t51 * t26 + t53 * t27;
t5 = t53 * t26 - t51 * t27;
t4 = pkin(7) * t36 + t9;
t3 = t54 * pkin(4) + t33 + t44 + (-pkin(7) * t52 - t25) * t50;
t2 = t51 * t3 + t53 * t4;
t1 = t53 * t3 - t51 * t4;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t47, 0.2e1 * t65, 0, t48, 0, 0, pkin(1) * t74, -0.2e1 * pkin(1) * t52, 0.2e1 * (t47 + t48) * pkin(6), pkin(1) ^ 2 + t48 * t56 + t43, t38, -0.2e1 * t61, t28, t35, 0.2e1 * t60, t48, -0.2e1 * t11 * t54 + 0.2e1 * t49 * t71, 0.2e1 * t12 * t54 + 0.2e1 * t50 * t71, (-t11 * t50 - t12 * t49) * t75, t11 ^ 2 + t12 ^ 2 + t43, t38, t28, 0.2e1 * t61, t48, -0.2e1 * t60, t35, 0.2e1 * t10 * t54 + 0.2e1 * t13 * t36, (t10 * t50 - t49 * t9) * t75, -0.2e1 * t13 * t39 - 0.2e1 * t9 * t54, t10 ^ 2 + t13 ^ 2 + t9 ^ 2, t16 ^ 2, t16 * t78, t16 * t74, t14 ^ 2, t54 * t78, t48, 0.2e1 * t1 * t54 + 0.2e1 * t7 * t14, -0.2e1 * t7 * t16 - 0.2e1 * t2 * t54, -0.2e1 * t1 * t16 + 0.2e1 * t2 * t14, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t54, 0, -t70, -t54 * pkin(6), 0, 0, t29, -t19, -t68, -t29, -t67, 0, t31 + (-pkin(6) * t50 - t72) * t52, pkin(6) * t36 + (-pkin(2) * t52 + t62) * t50, t57, -pkin(2) * t70 + t57 * qJ(3), t29, -t68, t19, 0, t67, -t29, -t13 * t50 + t24 * t36 + t31, t58, -t13 * t49 + (-t24 * t52 - t62) * t50, t58 * qJ(3) + t13 * t24, t16 * t22, t22 * t14 - t16 * t20, t22 * t54, -t14 * t20, -t20 * t54, 0, -t17 * t14 - t7 * t20 + t5 * t54, t17 * t16 - t7 * t22 - t6 * t54, -t1 * t22 + t6 * t14 - t5 * t16 - t2 * t20, t1 * t5 - t7 * t17 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, 0.2e1 * t69, 0, t46, 0, 0, 0.2e1 * pkin(2) * t50, -0.2e1 * t72, t23, pkin(2) ^ 2 + t63, t45, 0, -0.2e1 * t69, 0, 0, t46, t24 * t76, t23, -0.2e1 * t24 * t49, t24 ^ 2 + t63, t22 ^ 2, -0.2e1 * t22 * t20, 0, t20 ^ 2, 0, 0, t20 * t77, t22 * t77, -0.2e1 * t6 * t20 - 0.2e1 * t5 * t22, t17 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, 0, t70, 0, 0, 0, 0, 0, 0, t36, 0, -t39, t13, 0, 0, 0, 0, 0, 0, t14, -t16, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t50, 0, -t49, t24, 0, 0, 0, 0, 0, 0, -t20, -t22, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t39, 0, t10, 0, 0, 0, 0, 0, 0, t53 * t54, -t51 * t54, t51 * t14 - t53 * t16, t1 * t53 + t2 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t20 - t53 * t22, t5 * t53 + t6 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 ^ 2 + t53 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t14, t54, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t20, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
