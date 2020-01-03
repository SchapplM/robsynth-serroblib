% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t46 = sin(qJ(2));
t41 = t46 ^ 2;
t47 = cos(qJ(2));
t42 = t47 ^ 2;
t64 = t41 + t42;
t32 = -t47 * pkin(2) - t46 * qJ(3) - pkin(1);
t18 = t47 * pkin(3) - t32;
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t21 = t46 * t43 + t47 * t44;
t12 = t21 * pkin(4) + t18;
t63 = 0.2e1 * t12;
t62 = 0.2e1 * t18;
t61 = -0.2e1 * t46;
t60 = -pkin(2) - pkin(3);
t59 = t46 * pkin(6);
t58 = t47 * pkin(6);
t57 = cos(qJ(5));
t56 = t46 * t47;
t55 = pkin(6) - qJ(4);
t52 = t55 * t47;
t53 = t55 * t46;
t15 = t43 * t53 + t44 * t52;
t54 = t64 * pkin(6) ^ 2;
t28 = t43 * qJ(3) - t44 * t60;
t51 = -pkin(4) - t28;
t50 = -t46 * pkin(2) + t47 * qJ(3);
t13 = t43 * t52 - t44 * t53;
t23 = -t47 * t43 + t46 * t44;
t49 = -t23 * pkin(7) - t13;
t45 = sin(qJ(5));
t31 = 0.2e1 * t64 * pkin(6);
t30 = t44 * qJ(3) + t43 * t60;
t24 = t57 * t43 + t45 * t44;
t19 = t45 * t43 - t57 * t44;
t11 = t57 * t30 + t45 * t51;
t9 = t45 * t30 - t57 * t51;
t8 = -t45 * t21 + t57 * t23;
t6 = t57 * t21 + t45 * t23;
t5 = -t21 * pkin(7) + t15;
t3 = t45 * t49 + t57 * t5;
t1 = t45 * t5 - t57 * t49;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, 0.2e1 * t56, 0, t42, 0, 0, 0.2e1 * pkin(1) * t47, pkin(1) * t61, t31, pkin(1) ^ 2 + t54, t41, 0, -0.2e1 * t56, 0, 0, t42, -0.2e1 * t32 * t47, t31, t32 * t61, t32 ^ 2 + t54, t23 ^ 2, -0.2e1 * t23 * t21, 0, t21 ^ 2, 0, 0, t21 * t62, t23 * t62, 0.2e1 * t13 * t23 - 0.2e1 * t15 * t21, t13 ^ 2 + t15 ^ 2 + t18 ^ 2, t8 ^ 2, -0.2e1 * t8 * t6, 0, t6 ^ 2, 0, 0, t6 * t63, t8 * t63, 0.2e1 * t1 * t8 - 0.2e1 * t3 * t6, t1 ^ 2 + t12 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t47, 0, -t59, -t58, 0, 0, 0, t46, 0, 0, -t47, 0, -t59, t50, t58, t50 * pkin(6), 0, 0, -t23, 0, t21, 0, t13, t15, -t30 * t21 + t28 * t23, t13 * t28 + t15 * t30, 0, 0, -t8, 0, t6, 0, t1, t3, -t11 * t6 + t9 * t8, t1 * t9 + t3 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, 0.2e1 * t30, 0, t28 ^ 2 + t30 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t9, 0.2e1 * t11, 0, t11 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t21 - t44 * t23, -t13 * t44 + t15 * t43, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t8 - t24 * t6, t1 * t19 + t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t44, t43, 0, -t28 * t44 + t30 * t43, 0, 0, 0, 0, 0, 0, t19, t24, 0, t11 * t24 + t9 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 ^ 2 + t44 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, t18, 0, 0, 0, 0, 0, 0, t6, t8, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, -t6, 0, -t1, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t9, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
