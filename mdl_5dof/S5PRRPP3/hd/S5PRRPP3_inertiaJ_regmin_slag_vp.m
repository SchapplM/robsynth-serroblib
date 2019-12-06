% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPP3
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
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t59 = t30 ^ 2 + t31 ^ 2;
t15 = -pkin(4) * t31 - qJ(5) * t30 - pkin(3);
t58 = -0.2e1 * t15;
t32 = sin(qJ(3));
t57 = 0.2e1 * t32;
t34 = cos(qJ(3));
t56 = 0.2e1 * t34;
t55 = pkin(3) * t30;
t54 = pkin(6) * t30;
t53 = pkin(6) * t31;
t26 = t32 * pkin(6);
t52 = t34 * pkin(6);
t35 = cos(qJ(2));
t33 = sin(qJ(2));
t49 = t34 * t33;
t12 = -t30 * t35 + t31 * t49;
t8 = t12 * t31;
t29 = t32 ^ 2;
t51 = t29 * t33;
t20 = t30 * t32;
t16 = -pkin(3) * t34 - qJ(4) * t32 - pkin(2);
t50 = t31 * t16;
t21 = t31 * t32;
t23 = t32 * t33;
t5 = t30 * t16 + t31 * t52;
t48 = t59 * qJ(4) ^ 2;
t47 = qJ(4) * t34;
t46 = t30 * qJ(4);
t45 = t30 * t23;
t44 = t31 * t23;
t10 = t30 * t49 + t31 * t35;
t43 = t10 * t30 + t8;
t42 = t29 * t33 ^ 2 + t10 ^ 2 + t12 ^ 2;
t41 = t10 * t34 + t30 * t51;
t40 = qJ(4) * t8 + t10 * t46;
t2 = -qJ(5) * t34 + t5;
t3 = -t50 + (pkin(4) + t54) * t34;
t39 = t2 * t31 + t3 * t30;
t4 = -t30 * t52 + t50;
t38 = -t30 * t4 + t31 * t5;
t37 = (t10 * t31 - t12 * t30) * t32;
t18 = t34 * t46;
t14 = 0.2e1 * t59 * qJ(4);
t7 = t26 + (pkin(4) * t30 - qJ(5) * t31) * t32;
t1 = t12 * t34 + t31 * t51;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, t42; 0, 0, t35, -t33, 0, 0, 0, 0, 0, t35 * t34, -t35 * t32, t41, t1, t37, pkin(6) * t51 - t10 * t4 + t12 * t5, t41, t37, -t1, t10 * t3 + t12 * t2 + t7 * t23; 0, 1, 0, 0, t29, t32 * t56, 0, 0, 0, pkin(2) * t56, -0.2e1 * pkin(2) * t32, 0.2e1 * t29 * t54 - 0.2e1 * t34 * t4, 0.2e1 * t29 * t53 + 0.2e1 * t34 * t5, (-t30 * t5 - t31 * t4) * t57, pkin(6) ^ 2 * t29 + t4 ^ 2 + t5 ^ 2, 0.2e1 * t7 * t20 + 0.2e1 * t3 * t34, (-t2 * t30 + t3 * t31) * t57, -0.2e1 * t2 * t34 - 0.2e1 * t7 * t21, t2 ^ 2 + t3 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t49, -t44, t45, t43, -pkin(3) * t23 + t40, -t44, t43, -t45, t15 * t23 + t40; 0, 0, 0, 0, 0, 0, t32, t34, 0, -t26, -t52, t18 + (-t53 - t55) * t32, pkin(6) * t20 + (-pkin(3) * t32 + t47) * t31, t38, -pkin(3) * t26 + t38 * qJ(4), t15 * t20 - t31 * t7 + t18, t39, -t7 * t30 + (-t15 * t32 - t47) * t31, t39 * qJ(4) + t7 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t31, -0.2e1 * t55, t14, pkin(3) ^ 2 + t48, t31 * t58, t14, t30 * t58, t15 ^ 2 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, t26, t20, 0, -t21, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t30, 0, -pkin(3), -t31, 0, -t30, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t21, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
