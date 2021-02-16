% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR1
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
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:57
% EndTime: 2021-01-15 22:48:59
% DurationCPUTime: 0.36s
% Computational Cost: add. (577->57), mult. (1114->113), div. (0->0), fcn. (1315->8), ass. (0->51)
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t47 = cos(qJ(3));
t48 = cos(qJ(2));
t28 = t44 * t45 - t47 * t48;
t29 = t44 * t48 + t47 * t45;
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t16 = t42 * t28 + t41 * t29;
t38 = -t48 * pkin(2) - pkin(1);
t21 = t28 * pkin(3) + t38;
t58 = 0.2e1 * t16 * pkin(4) + 0.2e1 * t21;
t57 = 0.2e1 * t21;
t56 = 0.2e1 * t38;
t55 = 0.2e1 * t48;
t54 = pkin(6) + pkin(7);
t53 = t41 * pkin(3);
t52 = t44 * pkin(2);
t51 = t42 * t52;
t40 = t47 * pkin(2);
t37 = t40 + pkin(3);
t25 = t41 * t37 + t51;
t50 = -t25 - t53;
t33 = t54 * t45;
t34 = t54 * t48;
t18 = -t47 * t33 - t44 * t34;
t13 = -t29 * qJ(4) + t18;
t19 = t44 * t33 - t47 * t34;
t14 = -t28 * qJ(4) - t19;
t5 = t42 * t13 - t41 * t14;
t23 = t42 * t37 - t41 * t52;
t6 = t41 * t13 + t42 * t14;
t46 = cos(qJ(5));
t43 = sin(qJ(5));
t39 = t42 * pkin(3);
t35 = t39 + pkin(4);
t32 = t46 * t35;
t26 = -t43 * t35 - t46 * t53;
t24 = -t43 * t53 + t32;
t22 = pkin(4) + t23;
t20 = t46 * t22;
t17 = -t41 * t28 + t42 * t29;
t12 = -t43 * t22 - t46 * t25;
t11 = -t43 * t25 + t20;
t8 = -t43 * t16 + t46 * t17;
t7 = t46 * t16 + t43 * t17;
t4 = -t16 * pkin(8) + t6;
t3 = -t17 * pkin(8) + t5;
t2 = -t43 * t3 - t46 * t4;
t1 = t46 * t3 - t43 * t4;
t9 = [1, 0, 0, t45 ^ 2, t45 * t55, 0, 0, 0, pkin(1) * t55, -0.2e1 * pkin(1) * t45, t29 ^ 2, -0.2e1 * t29 * t28, 0, 0, 0, t28 * t56, t29 * t56, t16 * t57, t17 * t57, -0.2e1 * t6 * t16 - 0.2e1 * t5 * t17, t21 ^ 2 + t5 ^ 2 + t6 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t58, t8 * t58; 0, 0, 0, 0, 0, t45, t48, 0, -t45 * pkin(6), -t48 * pkin(6), 0, 0, t29, -t28, 0, t18, t19, t5, -t6, -t25 * t16 - t23 * t17, t5 * t23 + t6 * t25, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t52, 0.2e1 * t23, -0.2e1 * t25, 0, t23 ^ 2 + t25 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t18, t19, t5, -t6, (-t16 * t41 - t17 * t42) * pkin(3), (t41 * t6 + t42 * t5) * pkin(3), 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t52, t23 + t39, -t51 + (-pkin(3) - t37) * t41, 0, (t23 * t42 + t25 * t41) * pkin(3), 0, 0, 0, 0, 1, t50 * t43 + t20 + t32, t50 * t46 + (-t22 - t35) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, -0.2e1 * t53, 0, (t41 ^ 2 + t42 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t24, 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, t21, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
