% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:18
% EndTime: 2021-01-15 21:13:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (148->43), mult. (336->85), div. (0->0), fcn. (382->6), ass. (0->47)
t48 = 2 * pkin(1);
t29 = cos(qJ(2));
t30 = pkin(1) + pkin(2);
t16 = t30 * t29;
t47 = -0.2e1 * t16;
t27 = cos(qJ(5));
t39 = pkin(3) + qJ(3);
t15 = t39 * t29;
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t26 = sin(qJ(2));
t34 = t39 * t26;
t5 = t25 * t15 + t28 * t34;
t46 = t5 * t27;
t13 = t25 * t26 - t28 * t29;
t24 = sin(qJ(5));
t9 = t24 * t13;
t14 = t25 * t29 + t28 * t26;
t45 = t24 * t14;
t44 = t24 * t27;
t43 = t25 * t30;
t42 = t26 * t29;
t10 = t27 * t13;
t41 = t27 * t14;
t40 = t28 * t30;
t38 = t26 * qJ(3);
t37 = -0.2e1 * t14 * t13;
t36 = t24 * t40;
t35 = t27 * t40;
t17 = pkin(4) + t43;
t33 = -t13 * t17 - t14 * t40;
t32 = pkin(1) ^ 2;
t31 = qJ(3) ^ 2;
t23 = t29 ^ 2;
t22 = t27 ^ 2;
t21 = t26 ^ 2;
t20 = t24 ^ 2;
t18 = 0.2e1 * t44;
t11 = t14 ^ 2;
t8 = t24 * t41;
t7 = -t14 * pkin(4) - t16;
t6 = t28 * t15 - t25 * t34;
t4 = t5 * t24;
t3 = (-t20 + t22) * t14;
t2 = t24 * t7 + t27 * t6;
t1 = -t24 * t6 + t27 * t7;
t12 = [1, 0, 0, t21, 0.2e1 * t42, 0, 0, 0, 0, 0, t23 * t48, -0.2e1 * pkin(1) * t42, 0.2e1 * (t21 + t23) * qJ(3), t21 * t31 + (t31 + t32) * t23, t11, t37, 0, 0, 0, t13 * t47, t14 * t47, t22 * t11, -0.2e1 * t11 * t44, 0.2e1 * t13 * t41, t24 * t37, t13 ^ 2, 0.2e1 * t1 * t13 + 0.2e1 * t5 * t45, -0.2e1 * t2 * t13 + 0.2e1 * t5 * t41; 0, 0, 0, 0, 0, t26, t29, 0, 0, 0, -t38, -t29 * qJ(3), -t26 * pkin(1), -pkin(1) * t38, 0, 0, t14, -t13, 0, -t5, -t6, t8, t3, t9, t10, 0, t33 * t24 - t46, t33 * t27 + t4; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, 0, 0, t32, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t43, t20, t18, 0, 0, 0, 0.2e1 * t35, -0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t26, 0, -t29 * pkin(1), 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, -t5, -t6, t8, t3, t9, t10, 0, -pkin(4) * t9 - t46, -pkin(4) * t10 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t43, t20, t18, 0, 0, 0, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t20, t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t45, t13, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t27, 0, -t24 * t17, -t27 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t27, 0, -t24 * pkin(4), -t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
