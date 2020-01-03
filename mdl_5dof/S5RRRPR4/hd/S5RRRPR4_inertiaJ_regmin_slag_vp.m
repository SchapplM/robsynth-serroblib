% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t54 = sin(qJ(2)) * pkin(1);
t27 = pkin(7) + t54;
t40 = sin(qJ(3));
t37 = t40 ^ 2;
t43 = cos(qJ(3));
t48 = t43 ^ 2 + t37;
t50 = t48 * t27;
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t12 = t40 * t39 + t43 * t42;
t61 = 0.2e1 * t12;
t13 = -t43 * t39 + t40 * t42;
t60 = 0.2e1 * t13;
t59 = -0.2e1 * t40;
t58 = -0.2e1 * t43;
t57 = 0.2e1 * t43;
t34 = t43 * pkin(4);
t36 = cos(qJ(2)) * pkin(1);
t47 = t43 * pkin(3) + t40 * qJ(4) + pkin(2);
t7 = -t36 - t47;
t6 = t34 - t7;
t8 = t34 + t47;
t56 = t6 + t8;
t55 = t40 * pkin(8);
t53 = t43 * pkin(8);
t28 = -t36 - pkin(2);
t52 = pkin(2) - t28;
t51 = t47 - t7;
t49 = t48 * pkin(7);
t17 = -t40 * pkin(3) + t43 * qJ(4);
t45 = -pkin(3) - pkin(4);
t33 = t43 * pkin(7);
t31 = t40 * pkin(7);
t24 = t40 * t57;
t23 = t43 * t27;
t21 = t40 * t27;
t19 = t33 - t53;
t18 = t31 - t55;
t15 = t42 * qJ(4) + t39 * t45;
t14 = t39 * qJ(4) - t42 * t45;
t11 = t13 ^ 2;
t10 = t23 - t53;
t9 = t21 - t55;
t5 = t39 * t18 + t42 * t19;
t4 = -t42 * t18 + t39 * t19;
t3 = -0.2e1 * t13 * t12;
t2 = t42 * t10 + t39 * t9;
t1 = t39 * t10 - t42 * t9;
t16 = [1, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t54, t37, t24, 0, 0, 0, t28 * t58, 0.2e1 * t28 * t40, t7 * t58, 0.2e1 * t50, t7 * t59, t48 * t27 ^ 2 + t7 ^ 2, t11, t3, 0, 0, 0, t6 * t61, t6 * t60; 0, 0, 0, 1, t36, -t54, t37, t24, 0, 0, 0, t52 * t43, -t52 * t40, t51 * t43, t49 + t50, t51 * t40, pkin(7) * t50 - t47 * t7, t11, t3, 0, 0, 0, t56 * t12, t56 * t13; 0, 0, 0, 1, 0, 0, t37, t24, 0, 0, 0, pkin(2) * t57, pkin(2) * t59, -t47 * t58, 0.2e1 * t49, -t47 * t59, t48 * pkin(7) ^ 2 + t47 ^ 2, t11, t3, 0, 0, 0, t8 * t61, t8 * t60; 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t21, -t23, -t21, t17, t23, t17 * t27, 0, 0, -t13, t12, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t31, -t33, -t31, t17, t33, t17 * pkin(7), 0, 0, -t13, t12, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t14, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t42, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t16;
