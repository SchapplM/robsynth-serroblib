% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:00:31
% EndTime: 2019-05-05 03:00:33
% DurationCPUTime: 0.45s
% Computational Cost: add. (176->71), mult. (364->118), div. (0->0), fcn. (408->8), ass. (0->50)
t33 = sin(qJ(3));
t60 = -0.2e1 * t33;
t36 = cos(qJ(3));
t59 = -0.2e1 * t36;
t58 = 0.2e1 * t36;
t38 = -pkin(3) - pkin(4);
t57 = t36 * pkin(8);
t30 = sin(pkin(6));
t56 = t30 * sin(qJ(2));
t37 = cos(qJ(2));
t18 = t30 * t37;
t32 = sin(qJ(6));
t55 = t32 * t33;
t35 = cos(qJ(6));
t54 = t32 * t35;
t53 = t32 * t36;
t52 = t35 * t33;
t51 = t35 * t36;
t47 = cos(pkin(6));
t8 = t47 * t33 + t36 * t56;
t50 = t8 * qJ(4);
t27 = t33 ^ 2;
t29 = t36 ^ 2;
t49 = t27 + t29;
t48 = t36 * qJ(4);
t46 = t33 * t58;
t11 = -t36 * pkin(3) - t33 * qJ(4) - pkin(2);
t7 = t33 * t56 - t47 * t36;
t45 = t30 ^ 2 * t37 ^ 2 + t7 ^ 2 + t8 ^ 2;
t10 = t36 * pkin(4) - t11;
t44 = t7 * t33 + t8 * t36;
t43 = -t33 * pkin(3) + t48;
t25 = -pkin(9) + t38;
t31 = qJ(4) + pkin(5);
t42 = -t25 * t33 - t31 * t36;
t40 = qJ(4) ^ 2;
t39 = 0.2e1 * qJ(4);
t28 = t35 ^ 2;
t26 = t32 ^ 2;
t22 = t33 * pkin(8);
t16 = t36 * t18;
t15 = t33 * t18;
t13 = -t36 * qJ(5) + t57;
t12 = -t33 * qJ(5) + t22;
t5 = t33 * pkin(5) + t36 * pkin(9) + t10;
t4 = t32 * t18 + t7 * t35;
t3 = t35 * t18 - t7 * t32;
t2 = t35 * t12 + t32 * t5;
t1 = -t32 * t12 + t35 * t5;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, t18, -t56, 0, 0, 0, 0, 0, t16, -t15, t16, t44, t15, t44 * pkin(8) - t11 * t18, t15, -t16, -t44, t10 * t18 + t7 * t12 + t8 * t13, 0, 0, 0, 0, 0, t3 * t33 - t8 * t53, -t4 * t33 - t8 * t51; 0, 1, 0, 0, t27, t46, 0, 0, 0, pkin(2) * t58, pkin(2) * t60, t11 * t59, 0.2e1 * t49 * pkin(8), t11 * t60, t49 * pkin(8) ^ 2 + t11 ^ 2, 0.2e1 * t10 * t33, t10 * t59, -0.2e1 * t12 * t33 - 0.2e1 * t13 * t36, t10 ^ 2 + t12 ^ 2 + t13 ^ 2, t28 * t29, -0.2e1 * t29 * t54, t51 * t60, t32 * t46, t27, 0.2e1 * t1 * t33 - 0.2e1 * t13 * t53, -0.2e1 * t13 * t51 - 0.2e1 * t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7, 0, t8, -t7 * pkin(3) + t50, t8, t7, 0, t7 * t38 + t50, 0, 0, 0, 0, 0, t8 * t35, -t8 * t32; 0, 0, 0, 0, 0, 0, t33, t36, 0, -t22, -t57, -t22, t43, t57, t43 * pkin(8), t13, t12, -t38 * t33 - t48, t13 * qJ(4) + t12 * t38, t32 * t51 (-t26 + t28) * t36, -t55, -t52, 0, t13 * t35 + t42 * t32, -t13 * t32 + t42 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t39, pkin(3) ^ 2 + t40, t39, 0.2e1 * t38, 0, t38 ^ 2 + t40, t26, 0.2e1 * t54, 0, 0, 0, 0.2e1 * t31 * t35, -0.2e1 * t31 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t22, 0, 0, -t33, t12, 0, 0, 0, 0, 0, -t55, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 1, 0, t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t36, 0, t10, 0, 0, 0, 0, 0, t52, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t53, t33, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t35, 0, -t32 * t25, -t35 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
