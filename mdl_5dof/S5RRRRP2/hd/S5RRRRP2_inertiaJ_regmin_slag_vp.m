% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:27
% EndTime: 2022-01-20 11:49:29
% DurationCPUTime: 0.32s
% Computational Cost: add. (367->56), mult. (667->93), div. (0->0), fcn. (733->6), ass. (0->53)
t38 = sin(qJ(4));
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t49 = cos(qJ(4));
t22 = t38 * t39 - t49 * t41;
t61 = 0.2e1 * t22;
t23 = t38 * t41 + t49 * t39;
t60 = 0.2e1 * t23;
t33 = -t41 * pkin(3) - pkin(2);
t51 = cos(qJ(2)) * pkin(1);
t25 = t33 - t51;
t59 = 0.2e1 * t25;
t58 = 0.2e1 * t33;
t57 = 0.2e1 * t41;
t53 = sin(qJ(2)) * pkin(1);
t30 = pkin(7) + t53;
t19 = (-pkin(8) - t30) * t39;
t36 = t41 * pkin(8);
t48 = t41 * t30;
t20 = t36 + t48;
t10 = t49 * t19 - t38 * t20;
t44 = t23 * qJ(5);
t5 = t10 - t44;
t11 = -t38 * t19 - t49 * t20;
t45 = t22 * qJ(5);
t6 = -t11 - t45;
t56 = -t6 * t22 - t5 * t23;
t26 = (-pkin(7) - pkin(8)) * t39;
t52 = t41 * pkin(7);
t27 = t36 + t52;
t13 = t49 * t26 - t38 * t27;
t7 = t13 - t44;
t14 = -t38 * t26 - t49 * t27;
t8 = -t14 - t45;
t55 = -t8 * t22 - t7 * t23;
t54 = t38 * pkin(3);
t32 = -pkin(2) - t51;
t50 = pkin(2) - t32;
t16 = t22 * pkin(4) + t33;
t15 = t16 - t51;
t47 = t15 + t16;
t46 = t25 + t33;
t43 = 0.2e1 * pkin(4);
t37 = t39 ^ 2;
t35 = t49 * pkin(3);
t34 = -0.2e1 * t54;
t31 = t35 + pkin(4);
t28 = t39 * t57;
t21 = t23 ^ 2;
t18 = t23 * pkin(4);
t12 = -0.2e1 * t23 * t22;
t9 = -t22 * t54 - t31 * t23;
t1 = [1, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t53, t37, t28, 0, 0, 0, -0.2e1 * t32 * t41, 0.2e1 * t32 * t39, t21, t12, 0, 0, 0, t22 * t59, t23 * t59, t15 * t61, t15 * t60, 0.2e1 * t56, t15 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 1, t51, -t53, t37, t28, 0, 0, 0, t50 * t41, -t50 * t39, t21, t12, 0, 0, 0, t46 * t22, t46 * t23, t47 * t22, t47 * t23, t55 + t56, t15 * t16 + t5 * t7 + t6 * t8; 0, 0, 0, 1, 0, 0, t37, t28, 0, 0, 0, pkin(2) * t57, -0.2e1 * pkin(2) * t39, t21, t12, 0, 0, 0, t22 * t58, t23 * t58, t16 * t61, t16 * t60, 0.2e1 * t55, t16 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, -t39 * t30, -t48, 0, 0, t23, -t22, 0, t10, t11, t5, -t6, t9, t5 * t31 + t6 * t54; 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, -t39 * pkin(7), -t52, 0, 0, t23, -t22, 0, t13, t14, t7, -t8, t9, t7 * t31 + t8 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, t34, 0.2e1 * t31, t34, 0, t38 ^ 2 * pkin(3) ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t10, t11, t5, -t6, -t18, t5 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t13, t14, t7, -t8, -t18, t7 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t54, t43 + t35, -t54, 0, t31 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t43, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
