% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:57
% EndTime: 2021-01-15 20:29:59
% DurationCPUTime: 0.36s
% Computational Cost: add. (440->70), mult. (849->142), div. (0->0), fcn. (944->6), ass. (0->50)
t33 = cos(qJ(2));
t25 = -t33 * pkin(2) - pkin(1);
t54 = 0.2e1 * t25;
t30 = sin(qJ(4));
t53 = 0.2e1 * t30;
t32 = cos(qJ(4));
t52 = -0.2e1 * t32;
t51 = 0.2e1 * t33;
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(2));
t17 = t28 * t31 - t29 * t33;
t50 = t17 * pkin(4);
t49 = t28 * pkin(2);
t48 = t29 * pkin(2);
t47 = t32 * pkin(4);
t11 = t30 * t17;
t18 = t28 * t33 + t29 * t31;
t46 = t30 * t18;
t45 = t30 * t32;
t43 = -qJ(3) - pkin(6);
t21 = t43 * t33;
t39 = t43 * t31;
t10 = -t29 * t21 + t28 * t39;
t44 = t32 * t10;
t13 = t32 * t18;
t26 = t30 ^ 2;
t27 = t32 ^ 2;
t42 = t26 + t27;
t41 = qJ(5) * t18;
t23 = pkin(7) + t49;
t40 = qJ(5) + t23;
t24 = -pkin(3) - t48;
t7 = t17 * pkin(3) - t18 * pkin(7) + t25;
t3 = -t30 * t10 + t32 * t7;
t8 = -t28 * t21 - t29 * t39;
t35 = -t32 * t41 + t3;
t1 = t35 + t50;
t2 = t44 + (t7 - t41) * t30;
t38 = t1 * t32 + t2 * t30;
t14 = t40 * t30;
t15 = t40 * t32;
t37 = -t14 * t32 + t15 * t30;
t36 = -t17 * t23 + t18 * t24;
t20 = t24 - t47;
t16 = t18 ^ 2;
t12 = t32 * t17;
t5 = pkin(4) * t46 + t8;
t4 = t30 * t7 + t44;
t6 = [1, 0, 0, t31 ^ 2, t31 * t51, 0, 0, 0, pkin(1) * t51, -0.2e1 * pkin(1) * t31, t17 * t54, t18 * t54, -0.2e1 * t10 * t17 + 0.2e1 * t8 * t18, t10 ^ 2 + t25 ^ 2 + t8 ^ 2, t27 * t16, -0.2e1 * t16 * t45, 0.2e1 * t17 * t13, -0.2e1 * t17 * t46, t17 ^ 2, 0.2e1 * t3 * t17 + 0.2e1 * t8 * t46, 0.2e1 * t8 * t13 - 0.2e1 * t4 * t17, 0.2e1 * t1 * t17 + 0.2e1 * t5 * t46, 0.2e1 * t5 * t13 - 0.2e1 * t2 * t17, -0.2e1 * t38 * t18, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t31, t33, 0, -t31 * pkin(6), -t33 * pkin(6), -t8, -t10, (-t17 * t28 - t18 * t29) * pkin(2), (t10 * t28 - t29 * t8) * pkin(2), t30 * t13, (-t26 + t27) * t18, t11, t12, 0, t30 * t36 - t8 * t32, t8 * t30 + t32 * t36, -t14 * t17 + t20 * t46 - t5 * t32, t20 * t13 - t15 * t17 + t5 * t30, -t1 * t30 - t18 * t37 + t2 * t32, -t1 * t14 + t2 * t15 + t5 * t20; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t48, -0.2e1 * t49, 0, (t28 ^ 2 + t29 ^ 2) * pkin(2) ^ 2, t26, 0.2e1 * t45, 0, 0, 0, t24 * t52, t24 * t53, t20 * t52, t20 * t53, 0.2e1 * t14 * t30 + 0.2e1 * t15 * t32, t14 ^ 2 + t15 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, t25, 0, 0, 0, 0, 0, t12, -t11, t12, -t11, -t42 * t18, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t46, t17, t3, -t4, t35 + 0.2e1 * t50, -t2, -pkin(4) * t13, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, -t30 * t23, -t32 * t23, -t14, -t15, -t30 * pkin(4), -t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, t32, -t30, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t13, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
