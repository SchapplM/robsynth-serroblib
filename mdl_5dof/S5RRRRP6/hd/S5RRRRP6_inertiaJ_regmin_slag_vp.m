% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP6
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
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:35
% DurationCPUTime: 0.36s
% Computational Cost: add. (403->75), mult. (778->137), div. (0->0), fcn. (859->6), ass. (0->58)
t43 = cos(qJ(2));
t34 = -t43 * pkin(2) - pkin(1);
t61 = 0.2e1 * t34;
t60 = 0.2e1 * t43;
t59 = -pkin(7) - pkin(6);
t38 = sin(qJ(4));
t58 = t38 * pkin(4);
t39 = sin(qJ(3));
t57 = t39 * pkin(2);
t41 = cos(qJ(4));
t56 = t41 * pkin(8);
t42 = cos(qJ(3));
t55 = t42 * pkin(2);
t32 = -pkin(3) - t55;
t54 = pkin(3) - t32;
t40 = sin(qJ(2));
t27 = t59 * t40;
t28 = t59 * t43;
t11 = -t42 * t27 - t39 * t28;
t53 = t11 * t41;
t22 = t39 * t43 + t42 * t40;
t52 = t38 * t22;
t51 = t38 * t41;
t12 = t39 * t27 - t42 * t28;
t50 = t41 * t12;
t49 = t41 * t22;
t31 = pkin(8) + t57;
t48 = t41 * t31;
t35 = t41 * qJ(5);
t21 = t39 * t40 - t42 * t43;
t47 = -0.2e1 * t22 * t21;
t33 = -t41 * pkin(4) - pkin(3);
t8 = t21 * pkin(3) - t22 * pkin(8) + t34;
t4 = -t38 * t12 + t41 * t8;
t1 = t21 * pkin(4) - t22 * t35 + t4;
t3 = t50 + (-qJ(5) * t22 + t8) * t38;
t46 = -t1 * t38 + t3 * t41;
t45 = -pkin(3) * t22 - pkin(8) * t21;
t44 = -t21 * t31 + t22 * t32;
t37 = t41 ^ 2;
t36 = t38 ^ 2;
t29 = 0.2e1 * t51;
t26 = t35 + t56;
t25 = (-qJ(5) - pkin(8)) * t38;
t24 = t33 - t55;
t20 = t26 * t41;
t19 = t22 ^ 2;
t18 = t35 + t48;
t17 = (-qJ(5) - t31) * t38;
t16 = t41 * t21;
t15 = t38 * t21;
t14 = t18 * t41;
t13 = t38 * t49;
t10 = t11 * t38;
t9 = (-t36 + t37) * t22;
t6 = pkin(4) * t52 + t11;
t5 = t38 * t8 + t50;
t2 = [1, 0, 0, t40 ^ 2, t40 * t60, 0, 0, 0, pkin(1) * t60, -0.2e1 * pkin(1) * t40, t19, t47, 0, 0, 0, t21 * t61, t22 * t61, t37 * t19, -0.2e1 * t19 * t51, 0.2e1 * t21 * t49, t38 * t47, t21 ^ 2, 0.2e1 * t11 * t52 + 0.2e1 * t4 * t21, 0.2e1 * t11 * t49 - 0.2e1 * t5 * t21, 0.2e1 * (-t1 * t41 - t3 * t38) * t22, t1 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(6), -t43 * pkin(6), 0, 0, t22, -t21, 0, -t11, -t12, t13, t9, t15, t16, 0, t44 * t38 - t53, t44 * t41 + t10, (-t17 * t41 - t18 * t38) * t22 + t46, t1 * t17 + t3 * t18 + t6 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t57, t36, t29, 0, 0, 0, -0.2e1 * t32 * t41, 0.2e1 * t32 * t38, -0.2e1 * t17 * t38 + 0.2e1 * t14, t17 ^ 2 + t18 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, -t11, -t12, t13, t9, t15, t16, 0, t45 * t38 - t53, t45 * t41 + t10, (-t25 * t41 - t26 * t38) * t22 + t46, t1 * t25 + t3 * t26 + t6 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t57, t36, t29, 0, 0, 0, t54 * t41, -t54 * t38, t14 + t20 + (-t17 - t25) * t38, t17 * t25 + t18 * t26 + t24 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, t29, 0, 0, 0, 0.2e1 * pkin(3) * t41, -0.2e1 * pkin(3) * t38, -0.2e1 * t25 * t38 + 0.2e1 * t20, t25 ^ 2 + t26 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t52, t21, t4, -t5, -pkin(4) * t49, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * t31, -t48, -t58, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * pkin(8), -t56, -t58, t25 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
