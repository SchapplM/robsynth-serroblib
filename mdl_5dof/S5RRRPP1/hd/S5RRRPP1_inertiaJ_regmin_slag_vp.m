% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:47
% EndTime: 2021-01-15 22:14:50
% DurationCPUTime: 0.34s
% Computational Cost: add. (410->58), mult. (732->99), div. (0->0), fcn. (764->6), ass. (0->50)
t40 = cos(qJ(3));
t34 = t40 * qJ(4);
t57 = sin(qJ(2)) * pkin(1);
t47 = pkin(7) + t57;
t44 = t40 * t47;
t21 = t44 + t34;
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t38 = sin(qJ(3));
t43 = (-qJ(4) - t47) * t38;
t11 = t36 * t21 - t37 * t43;
t13 = t37 * t21 + t36 * t43;
t22 = t36 * t38 - t37 * t40;
t23 = t36 * t40 + t37 * t38;
t69 = t11 * t23 - t13 * t22;
t73 = 0.2e1 * t69;
t56 = t40 * pkin(7);
t26 = t34 + t56;
t45 = (-qJ(4) - pkin(7)) * t38;
t17 = t36 * t26 - t37 * t45;
t19 = t37 * t26 + t36 * t45;
t70 = t17 * t23 - t19 * t22;
t72 = 0.2e1 * t70;
t71 = t70 + t69;
t68 = 0.2e1 * t22;
t67 = -0.2e1 * t23;
t33 = -t40 * pkin(3) - pkin(2);
t55 = cos(qJ(2)) * pkin(1);
t25 = t33 - t55;
t66 = 0.2e1 * t25;
t65 = 0.2e1 * t33;
t64 = 0.2e1 * t40;
t59 = t36 * pkin(3);
t58 = t37 * pkin(3);
t32 = -pkin(2) - t55;
t54 = pkin(2) - t32;
t15 = t22 * pkin(4) - t23 * qJ(5) + t33;
t9 = t15 - t55;
t53 = t15 + t9;
t50 = t25 + t33;
t49 = t11 ^ 2 + t13 ^ 2;
t48 = t17 ^ 2 + t19 ^ 2;
t46 = t11 * t17 + t13 * t19;
t35 = t38 ^ 2;
t30 = pkin(4) + t58;
t28 = qJ(5) + t59;
t27 = t38 * t64;
t14 = (-t22 * t36 - t23 * t37) * pkin(3);
t5 = -t28 * t22 - t30 * t23;
t1 = [1, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t57, t35, t27, 0, 0, 0, -0.2e1 * t32 * t40, 0.2e1 * t32 * t38, t22 * t66, t23 * t66, t73, t25 ^ 2 + t49, t9 * t68, t73, t9 * t67, t9 ^ 2 + t49; 0, 0, 0, 1, t55, -t57, t35, t27, 0, 0, 0, t54 * t40, -t54 * t38, t50 * t22, t50 * t23, t71, t25 * t33 + t46, t53 * t22, t71, -t53 * t23, t9 * t15 + t46; 0, 0, 0, 1, 0, 0, t35, t27, 0, 0, 0, pkin(2) * t64, -0.2e1 * pkin(2) * t38, t22 * t65, t23 * t65, t72, t33 ^ 2 + t48, t15 * t68, t72, t15 * t67, t15 ^ 2 + t48; 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t38 * t47, -t44, -t11, -t13, t14, (-t11 * t37 + t13 * t36) * pkin(3), -t11, t5, t13, -t11 * t30 + t13 * t28; 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t38 * pkin(7), -t56, -t17, -t19, t14, (-t17 * t37 + t19 * t36) * pkin(3), -t17, t5, t19, -t17 * t30 + t19 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t58, -0.2e1 * t59, 0, (t36 ^ 2 + t37 ^ 2) * pkin(3) ^ 2, 0.2e1 * t30, 0, 0.2e1 * t28, t28 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t25, t22, 0, -t23, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t33, t22, 0, -t23, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
