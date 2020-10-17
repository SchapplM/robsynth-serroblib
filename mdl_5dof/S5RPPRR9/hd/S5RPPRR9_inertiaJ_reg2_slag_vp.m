% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:50
% DurationCPUTime: 0.55s
% Computational Cost: add. (298->72), mult. (468->127), div. (0->0), fcn. (460->6), ass. (0->54)
t32 = sin(qJ(4));
t58 = -0.2e1 * t32;
t34 = cos(qJ(4));
t57 = 0.2e1 * t34;
t33 = cos(qJ(5));
t56 = pkin(4) * t33;
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t35 = -pkin(1) - pkin(2);
t13 = t30 * qJ(2) + t29 * t35;
t10 = -pkin(6) + t13;
t55 = t10 * t29;
t26 = t32 ^ 2;
t31 = sin(qJ(5));
t54 = t26 * t31;
t53 = t26 * t33;
t27 = t33 ^ 2;
t52 = t27 * t32;
t51 = t31 * t32;
t50 = t31 * t33;
t49 = t31 * t34;
t48 = t32 * t10;
t47 = t32 * t29;
t46 = t33 * t32;
t20 = t33 * t34;
t45 = t34 * t10;
t44 = t34 * t29;
t25 = t31 ^ 2;
t43 = t25 + t27;
t28 = t34 ^ 2;
t42 = t26 + t28;
t41 = t32 * t57;
t40 = t31 * t46;
t39 = t43 * pkin(7);
t11 = qJ(2) * t29 - t30 * t35;
t9 = pkin(3) + t11;
t22 = t34 * pkin(4);
t3 = pkin(7) * t32 + t22 + t9;
t1 = t3 * t33 - t31 * t45;
t2 = t3 * t31 + t33 * t45;
t38 = -t1 * t31 + t2 * t33;
t5 = -t30 * t33 - t31 * t44;
t6 = -t30 * t31 + t33 * t44;
t37 = -t31 * t5 + t33 * t6;
t24 = t30 ^ 2;
t23 = t29 ^ 2;
t19 = t27 * t26;
t18 = t25 * t32;
t17 = t25 * t26;
t14 = t26 * t23;
t8 = t10 ^ 2;
t7 = t26 * t8;
t4 = t26 * t55;
t12 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t13, 0, t11 ^ 2 + t13 ^ 2, t26, t41, 0, t28, 0, 0, t9 * t57, t9 * t58, -0.2e1 * t42 * t10, t28 * t8 + t9 ^ 2 + t7, t19, -0.2e1 * t26 * t50, t20 * t58, t17, t31 * t41, t28, 0.2e1 * t1 * t34 - 0.2e1 * t10 * t54, -0.2e1 * t10 * t53 - 0.2e1 * t2 * t34, 0.2e1 * (t1 * t33 + t2 * t31) * t32, t1 ^ 2 + t2 ^ 2 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t30, t29, 0, -t11 * t30 + t13 * t29, 0, 0, 0, 0, 0, 0, -t30 * t34, t32 * t30, -t42 * t29, t28 * t55 - t30 * t9 + t4, 0, 0, 0, 0, 0, 0, -t29 * t54 + t34 * t5, -t29 * t53 - t34 * t6, (t31 * t6 + t33 * t5) * t32, t1 * t5 + t2 * t6 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t28 + t14 + t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t38 - t45) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t37 - t44) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 + t17 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t34, 0, -t48, -t45, 0, 0, -t40, t18 - t52, t49, t40, t20, 0, -t10 * t46 + (pkin(4) * t32 - pkin(7) * t34) * t31, -pkin(7) * t20 + (t10 * t31 + t56) * t32, t38, -pkin(4) * t48 + t38 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t44, 0, 0, 0, 0, 0, 0, 0, 0, -t29 * t46, t31 * t47, t37, -pkin(4) * t47 + t37 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t49, t18 + t52, t32 * t39 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, 0.2e1 * t50, 0, t27, 0, 0, 0.2e1 * t56, -0.2e1 * pkin(4) * t31, 0.2e1 * t39, t43 * pkin(7) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, t51, t34, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t33, 0, -t31 * pkin(7), -t33 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t12;
