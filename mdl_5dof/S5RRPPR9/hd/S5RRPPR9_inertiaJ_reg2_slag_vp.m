% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t30 = sin(qJ(2));
t25 = t30 ^ 2;
t32 = cos(qJ(2));
t27 = t32 ^ 2;
t56 = t25 + t27;
t51 = t32 * pkin(6);
t10 = -t32 * qJ(4) + t51;
t55 = t10 ^ 2;
t54 = -0.2e1 * t30;
t53 = -0.2e1 * t32;
t52 = 0.2e1 * t32;
t33 = -pkin(2) - pkin(3);
t29 = sin(qJ(5));
t24 = t29 ^ 2;
t50 = t24 * t32;
t49 = t29 * t30;
t31 = cos(qJ(5));
t48 = t29 * t31;
t47 = t29 * t32;
t46 = t30 * t32;
t45 = t31 * t30;
t44 = t31 * t32;
t43 = t56 * pkin(6) ^ 2;
t26 = t31 ^ 2;
t12 = t26 + t24;
t42 = t32 * qJ(3);
t41 = -0.2e1 * t46;
t13 = 0.2e1 * t46;
t8 = -t32 * pkin(2) - t30 * qJ(3) - pkin(1);
t40 = t29 * t44;
t5 = t32 * pkin(3) - t8;
t4 = t30 * pkin(4) + t32 * pkin(7) + t5;
t18 = t30 * pkin(6);
t9 = -t30 * qJ(4) + t18;
t2 = -t29 * t9 + t31 * t4;
t3 = t29 * t4 + t31 * t9;
t39 = t2 * t31 + t3 * t29;
t1 = -t2 * t29 + t3 * t31;
t38 = -t30 * pkin(2) + t42;
t23 = -pkin(7) + t33;
t28 = qJ(3) + pkin(4);
t37 = -t23 * t30 - t28 * t32;
t35 = qJ(3) ^ 2;
t34 = 0.2e1 * qJ(3);
t14 = t26 * t32;
t7 = 0.2e1 * t56 * pkin(6);
t6 = t12 * t23;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, t13, 0, t27, 0, 0, pkin(1) * t52, pkin(1) * t54, t7, pkin(1) ^ 2 + t43, t25, 0, t41, 0, 0, t27, t8 * t53, t7, t8 * t54, t8 ^ 2 + t43, t27, t13, 0, t25, 0, 0, 0.2e1 * t5 * t30, t5 * t53, -0.2e1 * t10 * t32 - 0.2e1 * t9 * t30, t5 ^ 2 + t9 ^ 2 + t55, t26 * t27, -0.2e1 * t27 * t48, t31 * t41, t24 * t27, t29 * t13, t25, -0.2e1 * t10 * t47 + 0.2e1 * t2 * t30, -0.2e1 * t10 * t44 - 0.2e1 * t3 * t30, t39 * t52, t2 ^ 2 + t3 ^ 2 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t32, 0, -t18, -t51, 0, 0, 0, t30, 0, 0, -t32, 0, -t18, t38, t51, t38 * pkin(6), 0, 0, t32, 0, t30, 0, t10, t9, -t33 * t30 - t42, t10 * qJ(3) + t9 * t33, t40, t14 - t50, -t49, -t40, -t45, 0, t10 * t31 + t29 * t37, -t10 * t29 + t31 * t37, -t1, t1 * t23 + t10 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t34, pkin(2) ^ 2 + t35, 0, 0, 0, 0, 0, 1, t34, 0.2e1 * t33, 0, t33 ^ 2 + t35, t24, 0.2e1 * t48, 0, t26, 0, 0, 0.2e1 * t28 * t31, -0.2e1 * t28 * t29, -0.2e1 * t6, t12 * t23 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t9, 0, 0, 0, 0, 0, 0, -t49, -t45, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 1, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t32, 0, t5, 0, 0, 0, 0, 0, 0, t45, -t49, t14 + t50, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, t47, t30, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t31, 0, -t29 * t23, -t31 * t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t11;
