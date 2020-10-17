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
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:27:56
% DurationCPUTime: 0.28s
% Computational Cost: add. (203->52), mult. (474->118), div. (0->0), fcn. (580->10), ass. (0->41)
t33 = cos(qJ(3));
t47 = 0.2e1 * t33;
t27 = sin(pkin(5));
t46 = t27 * sin(qJ(2));
t34 = cos(qJ(2));
t45 = t27 * t34;
t26 = sin(pkin(10));
t28 = cos(pkin(10));
t30 = sin(qJ(3));
t16 = t26 * t30 - t28 * t33;
t29 = sin(qJ(5));
t44 = t29 * t16;
t17 = t26 * t33 + t28 * t30;
t43 = t29 * t17;
t32 = cos(qJ(5));
t42 = t29 * t32;
t41 = t32 * t17;
t40 = -qJ(4) - pkin(7);
t39 = cos(pkin(5));
t23 = -t33 * pkin(3) - pkin(2);
t38 = t40 * t30;
t21 = t26 * pkin(3) + pkin(8);
t22 = -t28 * pkin(3) - pkin(4);
t37 = -t16 * t21 + t17 * t22;
t36 = -t30 * t46 + t39 * t33;
t25 = t32 ^ 2;
t24 = t29 ^ 2;
t19 = t40 * t33;
t15 = t17 ^ 2;
t14 = t39 * t30 + t33 * t46;
t13 = t32 * t16;
t11 = -t28 * t19 + t26 * t38;
t9 = -t26 * t19 - t28 * t38;
t8 = t16 * pkin(4) - t17 * pkin(8) + t23;
t7 = t28 * t14 + t26 * t36;
t5 = t26 * t14 - t28 * t36;
t4 = -t29 * t45 + t32 * t7;
t3 = -t29 * t7 - t32 * t45;
t2 = t32 * t11 + t29 * t8;
t1 = -t29 * t11 + t32 * t8;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 ^ 2 * t34 ^ 2 + t5 ^ 2 + t7 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t45, -t46, 0, 0, 0, 0, 0, t33 * t45, -t30 * t45, -t7 * t16 + t5 * t17, t7 * t11 - t23 * t45 + t5 * t9, 0, 0, 0, 0, 0, t3 * t16 + t5 * t43, -t4 * t16 + t5 * t41; 0, 1, 0, 0, t30 ^ 2, t30 * t47, 0, 0, 0, pkin(2) * t47, -0.2e1 * pkin(2) * t30, -0.2e1 * t11 * t16 + 0.2e1 * t9 * t17, t11 ^ 2 + t23 ^ 2 + t9 ^ 2, t25 * t15, -0.2e1 * t15 * t42, 0.2e1 * t16 * t41, -0.2e1 * t16 * t43, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t9 * t43, -0.2e1 * t2 * t16 + 0.2e1 * t9 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t14, 0, (t26 * t7 - t28 * t5) * pkin(3), 0, 0, 0, 0, 0, -t5 * t32, t5 * t29; 0, 0, 0, 0, 0, 0, t30, t33, 0, -t30 * pkin(7), -t33 * pkin(7), (-t16 * t26 - t17 * t28) * pkin(3), (t11 * t26 - t28 * t9) * pkin(3), t29 * t41, (-t24 + t25) * t17, t44, t13, 0, t37 * t29 - t9 * t32, t9 * t29 + t37 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t26 ^ 2 + t28 ^ 2) * pkin(3) ^ 2, t24, 0.2e1 * t42, 0, 0, 0, -0.2e1 * t22 * t32, 0.2e1 * t22 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t13, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t43, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t32, 0, -t29 * t21, -t32 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
