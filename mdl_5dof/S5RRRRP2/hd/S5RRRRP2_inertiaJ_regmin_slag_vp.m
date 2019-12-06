% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t41 = cos(qJ(3));
t34 = -t41 * pkin(3) - pkin(2);
t48 = cos(qJ(2)) * pkin(1);
t26 = t34 - t48;
t56 = 0.2e1 * t26;
t55 = 0.2e1 * t34;
t54 = 0.2e1 * t41;
t38 = sin(qJ(4));
t39 = sin(qJ(3));
t46 = cos(qJ(4));
t23 = t38 * t39 - t46 * t41;
t24 = t38 * t41 + t46 * t39;
t50 = sin(qJ(2)) * pkin(1);
t31 = pkin(7) + t50;
t20 = (-pkin(8) - t31) * t39;
t36 = t41 * pkin(8);
t45 = t41 * t31;
t21 = t36 + t45;
t10 = t46 * t20 - t38 * t21;
t43 = t24 * qJ(5);
t5 = t10 - t43;
t11 = -t38 * t20 - t46 * t21;
t18 = t23 * qJ(5);
t6 = -t18 - t11;
t53 = -t6 * t23 - t5 * t24;
t27 = (-pkin(7) - pkin(8)) * t39;
t49 = t41 * pkin(7);
t28 = t36 + t49;
t13 = t46 * t27 - t38 * t28;
t7 = t13 - t43;
t14 = -t38 * t27 - t46 * t28;
t8 = -t18 - t14;
t52 = -t8 * t23 - t7 * t24;
t51 = t38 * pkin(3);
t33 = -pkin(2) - t48;
t47 = pkin(2) - t33;
t44 = t26 + t34;
t16 = t23 * pkin(4) + t34;
t37 = t39 ^ 2;
t35 = t46 * pkin(3);
t32 = t35 + pkin(4);
t29 = t39 * t54;
t22 = t24 ^ 2;
t19 = pkin(4) * t24;
t15 = t16 - t48;
t12 = -0.2e1 * t24 * t23;
t9 = -t23 * t51 - t32 * t24;
t1 = [1, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t50, t37, t29, 0, 0, 0, -0.2e1 * t33 * t41, 0.2e1 * t33 * t39, t22, t12, 0, 0, 0, t23 * t56, t24 * t56, 0.2e1 * t53, t15 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 1, t48, -t50, t37, t29, 0, 0, 0, t47 * t41, -t47 * t39, t22, t12, 0, 0, 0, t44 * t23, t44 * t24, t52 + t53, t15 * t16 + t5 * t7 + t6 * t8; 0, 0, 0, 1, 0, 0, t37, t29, 0, 0, 0, pkin(2) * t54, -0.2e1 * pkin(2) * t39, t22, t12, 0, 0, 0, t23 * t55, t24 * t55, 0.2e1 * t52, t16 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, -t39 * t31, -t45, 0, 0, t24, -t23, 0, t10, t11, t9, t5 * t32 + t6 * t51; 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, -t39 * pkin(7), -t49, 0, 0, t24, -t23, 0, t13, t14, t9, t7 * t32 + t8 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t51, 0, t38 ^ 2 * pkin(3) ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t10, t11, -t19, t5 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t13, t14, -t19, t7 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t51, 0, t32 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
