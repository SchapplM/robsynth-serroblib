% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t43 = sin(qJ(2));
t63 = -0.2e1 * t43;
t45 = cos(qJ(2));
t62 = 0.2e1 * t45;
t44 = cos(qJ(3));
t61 = pkin(2) * t44;
t42 = sin(qJ(3));
t60 = pkin(6) * t42;
t27 = -t45 * pkin(2) - t43 * pkin(7) - pkin(1);
t24 = t44 * t27;
t53 = qJ(4) * t43;
t11 = -t44 * t53 + t24 + (-pkin(3) - t60) * t45;
t55 = t44 * t45;
t50 = pkin(6) * t55;
t13 = t50 + (t27 - t53) * t42;
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t4 = t40 * t11 + t41 * t13;
t59 = t42 * t43;
t58 = t42 * t44;
t57 = t42 * t45;
t56 = t44 * t43;
t54 = -qJ(4) - pkin(7);
t36 = t43 * pkin(6);
t26 = pkin(3) * t59 + t36;
t52 = t43 * t62;
t28 = t54 * t44;
t48 = t54 * t42;
t15 = -t40 * t28 - t41 * t48;
t17 = -t41 * t28 + t40 * t48;
t51 = t15 ^ 2 + t17 ^ 2;
t35 = -t44 * pkin(3) - pkin(2);
t23 = t40 * t44 + t41 * t42;
t20 = t23 * t43;
t21 = -t40 * t59 + t41 * t56;
t49 = t15 * t21 - t17 * t20;
t3 = t41 * t11 - t40 * t13;
t22 = t40 * t42 - t41 * t44;
t47 = 0.2e1 * t15 * t23 - 0.2e1 * t17 * t22;
t39 = t44 ^ 2;
t38 = t43 ^ 2;
t37 = t42 ^ 2;
t33 = t41 * pkin(3) + pkin(4);
t31 = t40 * pkin(3) + qJ(5);
t19 = t42 * t27 + t50;
t18 = -pkin(6) * t57 + t24;
t8 = t22 * pkin(4) - t23 * qJ(5) + t35;
t5 = t20 * pkin(4) - t21 * qJ(5) + t26;
t2 = t45 * pkin(4) - t3;
t1 = -t45 * qJ(5) + t4;
t6 = [1, 0, 0, t38, t52, 0, 0, 0, pkin(1) * t62, pkin(1) * t63, t39 * t38, -0.2e1 * t38 * t58, t55 * t63, t42 * t52, t45 ^ 2, -0.2e1 * t18 * t45 + 0.2e1 * t38 * t60, 0.2e1 * t38 * pkin(6) * t44 + 0.2e1 * t19 * t45, -0.2e1 * t4 * t20 - 0.2e1 * t3 * t21, t26 ^ 2 + t3 ^ 2 + t4 ^ 2, 0.2e1 * t2 * t45 + 0.2e1 * t5 * t20, -0.2e1 * t1 * t20 + 0.2e1 * t2 * t21, -0.2e1 * t1 * t45 - 0.2e1 * t5 * t21, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t43, t45, 0, -t36, -t45 * pkin(6), t42 * t56, (-t37 + t39) * t43, -t57, -t55, 0, -pkin(6) * t56 + (-pkin(2) * t43 + pkin(7) * t45) * t42, pkin(7) * t55 + (t60 - t61) * t43, -t4 * t22 - t3 * t23 + t49, -t3 * t15 + t4 * t17 + t26 * t35, t15 * t45 + t8 * t20 + t5 * t22, -t1 * t22 + t2 * t23 + t49, -t17 * t45 - t8 * t21 - t5 * t23, t1 * t17 + t2 * t15 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, 0.2e1 * t58, 0, 0, 0, 0.2e1 * t61, -0.2e1 * pkin(2) * t42, t47, t35 ^ 2 + t51, 0.2e1 * t8 * t22, t47, -0.2e1 * t8 * t23, t8 ^ 2 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t59, -t45, t18, -t19, (-t20 * t40 - t21 * t41) * pkin(3), (t3 * t41 + t4 * t40) * pkin(3), (-pkin(4) - t33) * t45 + t3, -t31 * t20 - t33 * t21, (-qJ(5) - t31) * t45 + t4, t1 * t31 - t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t44, 0, -t42 * pkin(7), -t44 * pkin(7), (-t22 * t40 - t23 * t41) * pkin(3), (-t15 * t41 + t17 * t40) * pkin(3), -t15, -t31 * t22 - t33 * t23, t17, -t15 * t33 + t17 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t40 ^ 2 + t41 ^ 2) * pkin(3) ^ 2, 0.2e1 * t33, 0, 0.2e1 * t31, t31 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t20, 0, -t21, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t22, 0, -t23, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
