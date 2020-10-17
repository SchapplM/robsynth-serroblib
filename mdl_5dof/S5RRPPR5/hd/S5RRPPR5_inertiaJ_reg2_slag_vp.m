% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:57
% EndTime: 2019-12-31 19:29:59
% DurationCPUTime: 0.51s
% Computational Cost: add. (404->62), mult. (751->114), div. (0->0), fcn. (845->6), ass. (0->42)
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t36 = sin(qJ(2));
t38 = cos(qJ(2));
t21 = t33 * t36 - t34 * t38;
t23 = t33 * t38 + t34 * t36;
t30 = -t38 * pkin(2) - pkin(1);
t42 = t23 * qJ(4) - t30;
t4 = (-pkin(3) - pkin(4)) * t21 + t42;
t56 = 0.2e1 * t4;
t55 = t21 ^ 2;
t54 = 0.2e1 * t30;
t53 = 0.2e1 * t38;
t52 = t33 * pkin(2);
t51 = t34 * pkin(2);
t50 = t23 * t21;
t49 = -qJ(3) - pkin(6);
t31 = t36 ^ 2;
t32 = t38 ^ 2;
t48 = t31 + t32;
t45 = t49 * t36;
t46 = t49 * t38;
t13 = -t33 * t46 - t34 * t45;
t15 = t33 * t45 - t34 * t46;
t47 = t13 ^ 2 + t15 ^ 2;
t28 = pkin(3) + t51;
t44 = -pkin(4) - t28;
t43 = 0.2e1 * t13 * t23 - 0.2e1 * t15 * t21;
t41 = -t23 * pkin(7) + t13;
t37 = cos(qJ(5));
t35 = sin(qJ(5));
t26 = qJ(4) + t52;
t20 = t23 ^ 2;
t18 = t37 * t26 + t35 * t44;
t16 = t35 * t26 - t37 * t44;
t11 = t35 * t21 + t37 * t23;
t9 = -t37 * t21 + t35 * t23;
t8 = t21 * pkin(3) - t42;
t6 = t21 * pkin(7) + t15;
t3 = t35 * t41 + t37 * t6;
t1 = t35 * t6 - t37 * t41;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, t36 * t53, 0, t32, 0, 0, pkin(1) * t53, -0.2e1 * pkin(1) * t36, 0.2e1 * t48 * pkin(6), t48 * pkin(6) ^ 2 + pkin(1) ^ 2, t20, -0.2e1 * t50, 0, t55, 0, 0, t21 * t54, t23 * t54, t43, t30 ^ 2 + t47, t20, 0, 0.2e1 * t50, 0, 0, t55, 0.2e1 * t8 * t21, t43, -0.2e1 * t8 * t23, t8 ^ 2 + t47, t11 ^ 2, -0.2e1 * t11 * t9, 0, t9 ^ 2, 0, 0, t9 * t56, t11 * t56, 0.2e1 * t1 * t11 - 0.2e1 * t3 * t9, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t38, 0, -t36 * pkin(6), -t38 * pkin(6), 0, 0, 0, 0, t23, 0, -t21, 0, -t13, -t15, (-t21 * t33 - t23 * t34) * pkin(2), (-t13 * t34 + t15 * t33) * pkin(2), 0, t23, 0, 0, t21, 0, -t13, -t26 * t21 - t28 * t23, t15, -t13 * t28 + t15 * t26, 0, 0, -t11, 0, t9, 0, t1, t3, t16 * t11 - t18 * t9, t1 * t16 + t3 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t52, 0, (t33 ^ 2 + t34 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t28, 0, 0.2e1 * t26, t26 ^ 2 + t28 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t16, 0.2e1 * t18, 0, t16 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, t30, 0, 0, 0, 0, 0, 0, t21, 0, -t23, t8, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t11 - t35 * t9, -t1 * t37 + t3 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t28, 0, 0, 0, 0, 0, 0, -t37, t35, 0, -t16 * t37 + t18 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, 0, -t1, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t16, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
