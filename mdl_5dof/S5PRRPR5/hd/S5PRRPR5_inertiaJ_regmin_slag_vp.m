% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:43
% EndTime: 2021-01-15 16:04:45
% DurationCPUTime: 0.28s
% Computational Cost: add. (225->58), mult. (528->125), div. (0->0), fcn. (648->10), ass. (0->44)
t32 = cos(qJ(3));
t23 = -t32 * pkin(3) - pkin(2);
t49 = 0.2e1 * t23;
t48 = 0.2e1 * t32;
t26 = sin(pkin(10));
t47 = t26 * pkin(3);
t28 = cos(pkin(10));
t46 = t28 * pkin(3);
t27 = sin(pkin(5));
t33 = cos(qJ(2));
t45 = t27 * t33;
t30 = sin(qJ(3));
t16 = t26 * t30 - t28 * t32;
t29 = sin(qJ(5));
t44 = t29 * t16;
t17 = t26 * t32 + t28 * t30;
t43 = t29 * t17;
t31 = cos(qJ(5));
t42 = t29 * t31;
t41 = t31 * t17;
t40 = qJ(4) + pkin(7);
t39 = cos(pkin(5));
t38 = t27 * sin(qJ(2));
t37 = t40 * t30;
t21 = pkin(8) + t47;
t22 = -pkin(4) - t46;
t36 = -t16 * t21 + t17 * t22;
t35 = -t30 * t38 + t39 * t32;
t25 = t31 ^ 2;
t24 = t29 ^ 2;
t19 = t40 * t32;
t15 = t17 ^ 2;
t14 = t39 * t30 + t32 * t38;
t13 = t31 * t16;
t11 = t28 * t19 - t26 * t37;
t9 = t26 * t19 + t28 * t37;
t8 = t16 * pkin(4) - t17 * pkin(8) + t23;
t7 = t28 * t14 + t26 * t35;
t5 = t26 * t14 - t28 * t35;
t4 = -t29 * t45 + t31 * t7;
t3 = -t29 * t7 - t31 * t45;
t2 = t31 * t11 + t29 * t8;
t1 = -t29 * t11 + t31 * t8;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 ^ 2 * t33 ^ 2 + t5 ^ 2 + t7 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t45, -t38, 0, 0, 0, 0, 0, t32 * t45, -t30 * t45, -t16 * t45, -t17 * t45, -t7 * t16 + t5 * t17, t7 * t11 - t23 * t45 + t5 * t9, 0, 0, 0, 0, 0, t3 * t16 + t5 * t43, -t4 * t16 + t5 * t41; 0, 1, 0, 0, t30 ^ 2, t30 * t48, 0, 0, 0, pkin(2) * t48, -0.2e1 * pkin(2) * t30, t16 * t49, t17 * t49, -0.2e1 * t11 * t16 + 0.2e1 * t9 * t17, t11 ^ 2 + t23 ^ 2 + t9 ^ 2, t25 * t15, -0.2e1 * t15 * t42, 0.2e1 * t16 * t41, -0.2e1 * t16 * t43, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t9 * t43, -0.2e1 * t2 * t16 + 0.2e1 * t9 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t14, -t5, -t7, 0, (t26 * t7 - t28 * t5) * pkin(3), 0, 0, 0, 0, 0, -t5 * t31, t5 * t29; 0, 0, 0, 0, 0, 0, t30, t32, 0, -t30 * pkin(7), -t32 * pkin(7), -t9, -t11, (-t16 * t26 - t17 * t28) * pkin(3), (t11 * t26 - t28 * t9) * pkin(3), t29 * t41, (-t24 + t25) * t17, t44, t13, 0, t36 * t29 - t9 * t31, t9 * t29 + t36 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t46, -0.2e1 * t47, 0, (t26 ^ 2 + t28 ^ 2) * pkin(3) ^ 2, t24, 0.2e1 * t42, 0, 0, 0, -0.2e1 * t22 * t31, 0.2e1 * t22 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, t23, 0, 0, 0, 0, 0, t13, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t43, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, -t29 * t21, -t31 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
