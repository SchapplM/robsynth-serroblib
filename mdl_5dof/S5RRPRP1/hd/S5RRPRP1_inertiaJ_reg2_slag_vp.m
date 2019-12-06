% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t31 = sin(qJ(4));
t27 = t31 ^ 2;
t33 = cos(qJ(4));
t28 = t33 ^ 2;
t19 = t27 + t28;
t29 = sin(pkin(8));
t49 = t29 * pkin(2);
t23 = pkin(7) + t49;
t41 = t19 * t23;
t52 = 0.2e1 * t31;
t51 = -0.2e1 * t33;
t34 = cos(qJ(2));
t26 = t34 * pkin(1);
t25 = t26 + pkin(2);
t30 = cos(pkin(8));
t32 = sin(qJ(2));
t46 = t32 * pkin(1);
t38 = t30 * t46;
t12 = t25 * t29 + t38;
t9 = pkin(7) + t12;
t50 = t19 * t9;
t48 = t30 * pkin(2);
t47 = t31 * pkin(4);
t45 = t33 * pkin(4);
t37 = -pkin(3) - t45;
t15 = t37 - t48;
t40 = -t30 * t25 + t29 * t46;
t4 = t37 + t40;
t44 = t15 + t4;
t24 = -pkin(3) - t48;
t8 = -pkin(3) + t40;
t43 = t24 + t8;
t42 = qJ(5) + t9;
t39 = qJ(5) + t23;
t22 = t33 * t52;
t14 = t39 * t33;
t13 = t39 * t31;
t10 = t14 * t33;
t3 = t42 * t33;
t2 = t42 * t31;
t1 = t3 * t33;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t46, 0, (t32 ^ 2 + t34 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t40, -0.2e1 * t12, 0, t12 ^ 2 + t40 ^ 2, t27, t22, 0, t28, 0, 0, t8 * t51, t8 * t52, 0.2e1 * t50, t19 * t9 ^ 2 + t8 ^ 2, t27, t22, 0, t28, 0, 0, t4 * t51, t4 * t52, 0.2e1 * t2 * t31 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t46, 0, 0, 0, 0, 0, 0, 0, 1, -t40 + t48, -t38 + (-pkin(2) - t25) * t29, 0, (t12 * t29 - t30 * t40) * pkin(2), t27, t22, 0, t28, 0, 0, -t43 * t33, t43 * t31, t41 + t50, t24 * t8 + t9 * t41, t27, t22, 0, t28, 0, 0, -t44 * t33, t44 * t31, t1 + t10 + (t13 + t2) * t31, t13 * t2 + t14 * t3 + t15 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t49, 0, (t29 ^ 2 + t30 ^ 2) * pkin(2) ^ 2, t27, t22, 0, t28, 0, 0, t24 * t51, t24 * t52, 0.2e1 * t41, t19 * t23 ^ 2 + t24 ^ 2, t27, t22, 0, t28, 0, 0, t15 * t51, t15 * t52, 0.2e1 * t13 * t31 + 0.2e1 * t10, t13 ^ 2 + t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t33 + t3 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t33 + t14 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t33, 0, -t31 * t9, -t33 * t9, 0, 0, 0, 0, t31, 0, t33, 0, -t2, -t3, -t47, -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t33, 0, -t31 * t23, -t33 * t23, 0, 0, 0, 0, t31, 0, t33, 0, -t13, -t14, -t47, -t13 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t31, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t31, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
