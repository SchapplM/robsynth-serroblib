% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:32
% EndTime: 2021-01-15 19:47:34
% DurationCPUTime: 0.34s
% Computational Cost: add. (458->69), mult. (920->142), div. (0->0), fcn. (1062->8), ass. (0->52)
t38 = sin(pkin(8));
t40 = cos(pkin(8));
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t24 = t38 * t44 + t40 * t42;
t37 = sin(pkin(9));
t39 = cos(pkin(9));
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t25 = t43 * t37 + t41 * t39;
t8 = t25 * t24;
t62 = -0.2e1 * t8;
t51 = -qJ(3) - pkin(6);
t28 = t51 * t44;
t49 = t51 * t42;
t17 = -t38 * t28 - t40 * t49;
t61 = t17 ^ 2;
t56 = t40 * pkin(2);
t33 = -pkin(3) - t56;
t27 = -t39 * pkin(4) + t33;
t60 = 0.2e1 * t27;
t34 = -t44 * pkin(2) - pkin(1);
t59 = 0.2e1 * t34;
t58 = 0.2e1 * t44;
t57 = t38 * pkin(2);
t30 = qJ(4) + t57;
t55 = pkin(7) + t30;
t22 = t38 * t42 - t40 * t44;
t54 = t25 * t22;
t53 = t37 * t24;
t52 = t39 * t24;
t14 = t22 * pkin(3) - t24 * qJ(4) + t34;
t19 = -t40 * t28 + t38 * t49;
t6 = t37 * t14 + t39 * t19;
t50 = t37 ^ 2 + t39 ^ 2;
t5 = t39 * t14 - t37 * t19;
t48 = t6 * t37 + t5 * t39;
t47 = -t5 * t37 + t6 * t39;
t46 = -t22 * t30 + t24 * t33;
t23 = t41 * t37 - t43 * t39;
t21 = t55 * t39;
t20 = t55 * t37;
t16 = t23 * t22;
t13 = -t41 * t20 + t43 * t21;
t12 = -t43 * t20 - t41 * t21;
t9 = t23 * t24;
t7 = pkin(4) * t53 + t17;
t4 = -pkin(7) * t53 + t6;
t3 = t22 * pkin(4) - pkin(7) * t52 + t5;
t2 = t41 * t3 + t43 * t4;
t1 = t43 * t3 - t41 * t4;
t10 = [1, 0, 0, t42 ^ 2, t42 * t58, 0, 0, 0, pkin(1) * t58, -0.2e1 * pkin(1) * t42, t22 * t59, t24 * t59, 0.2e1 * t17 * t24 - 0.2e1 * t19 * t22, t19 ^ 2 + t34 ^ 2 + t61, 0.2e1 * t17 * t53 + 0.2e1 * t5 * t22, 0.2e1 * t17 * t52 - 0.2e1 * t6 * t22, -0.2e1 * t48 * t24, t5 ^ 2 + t6 ^ 2 + t61, t9 ^ 2, -t9 * t62, -0.2e1 * t9 * t22, t22 * t62, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t7 * t8, -0.2e1 * t2 * t22 - 0.2e1 * t7 * t9; 0, 0, 0, 0, 0, t42, t44, 0, -t42 * pkin(6), -t44 * pkin(6), -t17, -t19, (-t22 * t38 - t24 * t40) * pkin(2), (-t17 * t40 + t19 * t38) * pkin(2), -t17 * t39 + t46 * t37, t17 * t37 + t46 * t39, t47, t17 * t33 + t47 * t30, -t9 * t25, t9 * t23 - t25 * t8, t54, -t16, 0, t12 * t22 + t7 * t23 + t27 * t8, -t13 * t22 + t7 * t25 - t27 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t56, -0.2e1 * t57, 0, (t38 ^ 2 + t40 ^ 2) * pkin(2) ^ 2, -0.2e1 * t33 * t39, 0.2e1 * t33 * t37, 0.2e1 * t50 * t30, t50 * t30 ^ 2 + t33 ^ 2, t25 ^ 2, -0.2e1 * t25 * t23, 0, 0, 0, t23 * t60, t25 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t24, 0, t34, t39 * t22, -t37 * t22, -t50 * t24, t48, 0, 0, 0, 0, 0, -t16, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t52, 0, t17, 0, 0, 0, 0, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t37, 0, t33, 0, 0, 0, 0, 0, t23, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t8, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
