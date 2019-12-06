% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t14 = t28 * t30 - t29 * t32;
t52 = t14 ^ 2;
t23 = -t32 * pkin(3) - pkin(2);
t51 = 0.2e1 * t23;
t50 = 0.2e1 * t32;
t49 = t28 * pkin(3);
t48 = t29 * pkin(3);
t16 = t28 * t32 + t29 * t30;
t47 = t16 * t14;
t33 = cos(qJ(2));
t46 = t33 * t14;
t45 = t33 * t16;
t44 = -qJ(4) - pkin(6);
t24 = t30 ^ 2;
t26 = t32 ^ 2;
t43 = t24 + t26;
t18 = t44 * t32;
t38 = t44 * t30;
t6 = -t28 * t18 - t29 * t38;
t8 = -t29 * t18 + t28 * t38;
t42 = t6 ^ 2 + t8 ^ 2;
t31 = sin(qJ(2));
t10 = t16 * t31;
t12 = t14 * t31;
t41 = t10 * t6 - t12 * t8;
t40 = t10 * t16 + t12 * t14;
t27 = t33 ^ 2;
t39 = t10 ^ 2 + t12 ^ 2 + t27;
t37 = t43 * t31;
t36 = -0.2e1 * t8 * t14 + 0.2e1 * t6 * t16;
t25 = t31 ^ 2;
t21 = pkin(4) + t48;
t19 = qJ(5) + t49;
t13 = t16 ^ 2;
t3 = t14 * pkin(4) - t16 * qJ(5) + t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t25 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t32, -t33 * t30, t37, t33 * pkin(2) + pkin(6) * t37, 0, 0, 0, 0, 0, 0, -t46, -t45, t40, -t33 * t23 + t41, 0, 0, 0, 0, 0, 0, -t46, t40, t45, -t33 * t3 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t24, t30 * t50, 0, t26, 0, 0, pkin(2) * t50, -0.2e1 * pkin(2) * t30, 0.2e1 * t43 * pkin(6), t43 * pkin(6) ^ 2 + pkin(2) ^ 2, t13, -0.2e1 * t47, 0, t52, 0, 0, t14 * t51, t16 * t51, t36, t23 ^ 2 + t42, t13, 0, 0.2e1 * t47, 0, 0, t52, 0.2e1 * t3 * t14, t36, -0.2e1 * t3 * t16, t3 ^ 2 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t31, -t32 * t31, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t12, 0, (-t10 * t29 - t12 * t28) * pkin(3), 0, 0, 0, 0, 0, 0, -t10, 0, -t12, -t10 * t21 - t12 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t32, 0, -t30 * pkin(6), -t32 * pkin(6), 0, 0, 0, 0, t16, 0, -t14, 0, -t6, -t8, (-t14 * t28 - t16 * t29) * pkin(3), (t28 * t8 - t29 * t6) * pkin(3), 0, t16, 0, 0, t14, 0, -t6, -t19 * t14 - t21 * t16, t8, t8 * t19 - t6 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t49, 0, (t28 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t21, 0, 0.2e1 * t19, t19 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, t23, 0, 0, 0, 0, 0, 0, t14, 0, -t16, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
