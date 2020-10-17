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
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:49:58
% DurationCPUTime: 0.28s
% Computational Cost: add. (368->53), mult. (652->89), div. (0->0), fcn. (674->6), ass. (0->45)
t40 = cos(qJ(3));
t34 = t40 * qJ(4);
t55 = t40 * pkin(7);
t26 = t34 + t55;
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t38 = sin(qJ(3));
t45 = (-qJ(4) - pkin(7)) * t38;
t17 = t36 * t26 - t37 * t45;
t19 = t37 * t26 + t36 * t45;
t22 = t36 * t38 - t37 * t40;
t23 = t36 * t40 + t37 * t38;
t64 = t17 * t23 - t19 * t22;
t68 = 0.2e1 * t64;
t56 = sin(qJ(2)) * pkin(1);
t47 = pkin(7) + t56;
t44 = t40 * t47;
t21 = t44 + t34;
t43 = (-qJ(4) - t47) * t38;
t11 = t36 * t21 - t37 * t43;
t13 = t37 * t21 + t36 * t43;
t65 = t11 * t23 - t13 * t22;
t67 = 0.2e1 * t65;
t66 = t64 + t65;
t63 = 0.2e1 * t22;
t62 = -0.2e1 * t23;
t61 = 0.2e1 * t40;
t54 = cos(qJ(2)) * pkin(1);
t32 = -pkin(2) - t54;
t53 = pkin(2) - t32;
t33 = -t40 * pkin(3) - pkin(2);
t15 = t22 * pkin(4) - t23 * qJ(5) + t33;
t9 = t15 - t54;
t52 = t15 + t9;
t49 = t11 ^ 2 + t13 ^ 2;
t48 = t17 ^ 2 + t19 ^ 2;
t46 = t11 * t17 + t13 * t19;
t35 = t38 ^ 2;
t30 = t37 * pkin(3) + pkin(4);
t28 = t36 * pkin(3) + qJ(5);
t27 = t38 * t61;
t25 = t33 - t54;
t14 = (-t22 * t36 - t23 * t37) * pkin(3);
t5 = -t28 * t22 - t30 * t23;
t1 = [1, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t56, t35, t27, 0, 0, 0, -0.2e1 * t32 * t40, 0.2e1 * t32 * t38, t67, t25 ^ 2 + t49, t9 * t63, t67, t9 * t62, t9 ^ 2 + t49; 0, 0, 0, 1, t54, -t56, t35, t27, 0, 0, 0, t53 * t40, -t53 * t38, t66, t25 * t33 + t46, t52 * t22, t66, -t52 * t23, t9 * t15 + t46; 0, 0, 0, 1, 0, 0, t35, t27, 0, 0, 0, pkin(2) * t61, -0.2e1 * pkin(2) * t38, t68, t33 ^ 2 + t48, t15 * t63, t68, t15 * t62, t15 ^ 2 + t48; 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t38 * t47, -t44, t14, (-t11 * t37 + t13 * t36) * pkin(3), -t11, t5, t13, -t11 * t30 + t13 * t28; 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t38 * pkin(7), -t55, t14, (-t17 * t37 + t19 * t36) * pkin(3), -t17, t5, t19, -t17 * t30 + t19 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t36 ^ 2 + t37 ^ 2) * pkin(3) ^ 2, 0.2e1 * t30, 0, 0.2e1 * t28, t28 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t22, 0, -t23, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t22, 0, -t23, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
