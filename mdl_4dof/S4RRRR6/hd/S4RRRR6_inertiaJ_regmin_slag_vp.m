% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:30:52
% DurationCPUTime: 0.36s
% Computational Cost: add. (279->74), mult. (773->173), div. (0->0), fcn. (864->8), ass. (0->60)
t27 = cos(pkin(4));
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t26 = sin(pkin(4));
t30 = sin(qJ(2));
t48 = t26 * t30;
t16 = t27 * t29 + t32 * t48;
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t33 = cos(qJ(2));
t47 = t26 * t33;
t8 = t16 * t28 + t31 * t47;
t60 = -0.2e1 * t8;
t59 = -0.2e1 * t16;
t58 = -0.2e1 * t29;
t57 = 0.2e1 * t32;
t56 = pkin(1) * t30;
t55 = pkin(1) * t33;
t54 = pkin(3) * t31;
t53 = pkin(7) * t28;
t36 = pkin(6) * t47;
t13 = t36 + (pkin(7) + t56) * t27;
t14 = (-pkin(2) * t33 - pkin(7) * t30 - pkin(1)) * t26;
t6 = -t29 * t13 + t32 * t14;
t4 = pkin(3) * t47 - t6;
t52 = t4 * t28;
t51 = t4 * t31;
t9 = t16 * t31 - t28 * t47;
t50 = t9 * t28;
t22 = t26 ^ 2;
t49 = t22 * t33;
t46 = t27 * t30;
t15 = -t27 * t32 + t29 * t48;
t45 = t28 * t15;
t44 = t28 * t31;
t43 = t28 * t32;
t42 = t29 * t15;
t41 = t31 * t15;
t40 = t31 * t29;
t39 = t31 * t32;
t38 = 0.2e1 * t47;
t37 = t29 * t57;
t35 = t32 * t47;
t34 = t29 * t47;
t7 = t32 * t13 + t29 * t14;
t20 = pkin(6) * t48;
t12 = t20 + (-pkin(2) - t55) * t27;
t25 = t31 ^ 2;
t24 = t29 ^ 2;
t23 = t28 ^ 2;
t19 = -t32 * pkin(3) - t29 * pkin(8) - pkin(2);
t18 = pkin(1) * t46 + t36;
t17 = t27 * t55 - t20;
t11 = pkin(7) * t39 + t28 * t19;
t10 = -pkin(7) * t43 + t31 * t19;
t5 = -pkin(8) * t47 + t7;
t3 = t15 * pkin(3) - t16 * pkin(8) + t12;
t2 = t28 * t3 + t31 * t5;
t1 = -t28 * t5 + t31 * t3;
t21 = [1, 0, 0, t22 * t30 ^ 2, 0.2e1 * t30 * t49, 0.2e1 * t26 * t46, t27 * t38, t27 ^ 2, 0.2e1 * pkin(1) * t49 + 0.2e1 * t17 * t27, -0.2e1 * t18 * t27 - 0.2e1 * t22 * t56, t16 ^ 2, t15 * t59, t47 * t59, t15 * t38, t22 * t33 ^ 2, 0.2e1 * t12 * t15 - 0.2e1 * t6 * t47, 0.2e1 * t12 * t16 + 0.2e1 * t7 * t47, t9 ^ 2, t9 * t60, 0.2e1 * t9 * t15, t15 * t60, t15 ^ 2, 0.2e1 * t1 * t15 + 0.2e1 * t4 * t8, -0.2e1 * t2 * t15 + 0.2e1 * t4 * t9; 0, 0, 0, 0, 0, t48, t47, t27, t17, -t18, t16 * t29, t16 * t32 - t42, -t34, -t35, 0, -pkin(2) * t15 + pkin(7) * t34 - t12 * t32, -pkin(2) * t16 + pkin(7) * t35 + t12 * t29, t9 * t40, (-t31 * t8 - t50) * t29, t15 * t40 - t9 * t32, -t28 * t42 + t8 * t32, -t15 * t32, -t1 * t32 + t10 * t15 + (pkin(7) * t8 + t52) * t29, -t11 * t15 + t2 * t32 + (pkin(7) * t9 + t51) * t29; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, t37, 0, 0, 0, pkin(2) * t57, pkin(2) * t58, t25 * t24, -0.2e1 * t24 * t44, t39 * t58, t28 * t37, t32 ^ 2, -0.2e1 * t10 * t32 + 0.2e1 * t24 * t53, 0.2e1 * t24 * pkin(7) * t31 + 0.2e1 * t11 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t47, t6, -t7, t50, -t28 * t8 + t9 * t31, t45, t41, 0, -pkin(3) * t8 - pkin(8) * t45 - t51, -pkin(3) * t9 - pkin(8) * t41 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t32, 0, -t29 * pkin(7), -t32 * pkin(7), t28 * t40, (-t23 + t25) * t29, -t43, -t39, 0, -pkin(7) * t40 + (-pkin(3) * t29 + pkin(8) * t32) * t28, pkin(8) * t39 + (t53 - t54) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, 0.2e1 * t44, 0, 0, 0, 0.2e1 * t54, -0.2e1 * pkin(3) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, t15, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t28 * t29, -t32, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(8), -t31 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t21;
