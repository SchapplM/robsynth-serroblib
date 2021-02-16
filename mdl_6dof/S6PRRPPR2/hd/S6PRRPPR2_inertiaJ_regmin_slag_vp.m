% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:18
% EndTime: 2021-01-16 02:21:20
% DurationCPUTime: 0.53s
% Computational Cost: add. (387->74), mult. (835->139), div. (0->0), fcn. (1024->10), ass. (0->60)
t34 = sin(pkin(11));
t36 = cos(pkin(11));
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t21 = t34 * t38 - t36 * t40;
t22 = t34 * t40 + t36 * t38;
t31 = -t40 * pkin(3) - pkin(2);
t45 = -t22 * qJ(5) + t31;
t12 = t21 * pkin(4) + t45;
t69 = -0.2e1 * t12;
t65 = t34 * pkin(3);
t26 = qJ(5) + t65;
t68 = 0.2e1 * t26;
t67 = 0.2e1 * t31;
t66 = 0.2e1 * t40;
t64 = t36 * pkin(3);
t63 = t26 * t21;
t35 = sin(pkin(6));
t41 = cos(qJ(2));
t62 = t35 * t41;
t37 = sin(qJ(6));
t61 = t37 * t21;
t60 = t37 * t22;
t39 = cos(qJ(6));
t59 = t39 * t21;
t18 = t39 * t22;
t58 = t39 * t37;
t57 = qJ(4) + pkin(8);
t56 = cos(pkin(6));
t55 = 0.2e1 * t21 * t22;
t54 = t21 * t62;
t53 = t22 * t62;
t24 = t57 * t40;
t50 = t57 * t38;
t13 = t34 * t24 + t36 * t50;
t15 = t36 * t24 - t34 * t50;
t52 = t13 ^ 2 + t15 ^ 2;
t30 = -pkin(4) - t64;
t51 = t35 * sin(qJ(2));
t19 = t56 * t38 + t40 * t51;
t43 = -t38 * t51 + t56 * t40;
t10 = t36 * t19 + t34 * t43;
t8 = t34 * t19 - t36 * t43;
t49 = t35 ^ 2 * t41 ^ 2 + t10 ^ 2 + t8 ^ 2;
t48 = t10 * t15 + t8 * t13;
t47 = -t10 * t21 + t8 * t22;
t25 = -pkin(9) + t30;
t46 = -t22 * t25 + t63;
t44 = 0.2e1 * t13 * t22 - 0.2e1 * t15 * t21;
t33 = t39 ^ 2;
t32 = t37 ^ 2;
t20 = t21 ^ 2;
t7 = -t21 * pkin(5) + t15;
t6 = t22 * pkin(5) + t13;
t5 = (pkin(4) + pkin(9)) * t21 + t45;
t4 = -t37 * t8 + t39 * t62;
t3 = t37 * t62 + t39 * t8;
t2 = t37 * t6 + t39 * t5;
t1 = -t37 * t5 + t39 * t6;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0; 0, 0, t62, -t51, 0, 0, 0, 0, 0, t40 * t62, -t38 * t62, -t54, -t53, t47, -t31 * t62 + t48, t47, t54, t53, -t12 * t62 + t48, 0, 0, 0, 0, 0, -t10 * t59 + t3 * t22, t10 * t61 + t4 * t22; 0, 1, 0, 0, t38 ^ 2, t38 * t66, 0, 0, 0, pkin(2) * t66, -0.2e1 * pkin(2) * t38, t21 * t67, t22 * t67, t44, t31 ^ 2 + t52, t44, t21 * t69, t22 * t69, t12 ^ 2 + t52, t32 * t20, 0.2e1 * t20 * t58, t37 * t55, t39 * t55, t22 ^ 2, 0.2e1 * t1 * t22 - 0.2e1 * t7 * t59, -0.2e1 * t2 * t22 + 0.2e1 * t7 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t19, -t8, -t10, 0, (t10 * t34 - t36 * t8) * pkin(3), 0, t8, t10, t10 * t26 + t8 * t30, 0, 0, 0, 0, 0, t10 * t37, t10 * t39; 0, 0, 0, 0, 0, 0, t38, t40, 0, -t38 * pkin(8), -t40 * pkin(8), -t13, -t15, (-t21 * t34 - t22 * t36) * pkin(3), (-t13 * t36 + t15 * t34) * pkin(3), t30 * t22 - t63, t13, t15, t13 * t30 + t15 * t26, t21 * t58, (-t32 + t33) * t21, t18, -t60, 0, t7 * t37 - t46 * t39, t46 * t37 + t7 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t64, -0.2e1 * t65, 0, (t34 ^ 2 + t36 ^ 2) * pkin(3) ^ 2, 0, 0.2e1 * t30, t68, t26 ^ 2 + t30 ^ 2, t33, -0.2e1 * t58, 0, 0, 0, t37 * t68, t39 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, t31, 0, -t21, -t22, t12, 0, 0, 0, 0, 0, -t60, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, t13, 0, 0, 0, 0, 0, t18, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t59, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, t39 * t25, -t37 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
