% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR11
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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(2));
t43 = cos(qJ(3));
t30 = t43 * t41;
t39 = sin(qJ(5));
t40 = sin(qJ(3));
t42 = cos(qJ(5));
t57 = t42 * t40;
t13 = t39 * t30 - t41 * t57;
t70 = -0.2e1 * t13;
t53 = t40 * qJ(4);
t65 = pkin(3) + pkin(4);
t15 = t65 * t43 + pkin(2) + t53;
t69 = 0.2e1 * t15;
t47 = -t43 * pkin(3) - t53;
t23 = -pkin(2) + t47;
t68 = -0.2e1 * t23;
t67 = -0.2e1 * t41;
t44 = cos(qJ(2));
t66 = 0.2e1 * t44;
t64 = pkin(2) * t40;
t63 = pkin(2) * t43;
t62 = pkin(6) * t40;
t61 = pkin(6) * t43;
t60 = t40 * t41;
t59 = t40 * t43;
t58 = t40 * t44;
t56 = t43 * t44;
t24 = -t44 * pkin(2) - t41 * pkin(7) - pkin(1);
t55 = pkin(6) * t58 - t43 * t24;
t11 = pkin(6) * t56 + t40 * t24;
t35 = t40 ^ 2;
t37 = t43 ^ 2;
t54 = t35 + t37;
t52 = t43 * qJ(4);
t51 = t44 * qJ(4);
t50 = t41 * t66;
t34 = t44 * pkin(3);
t9 = t34 + t55;
t32 = t40 * pkin(7);
t49 = -t40 * pkin(8) + t32;
t8 = -t51 + t11;
t3 = t44 * pkin(4) - pkin(8) * t30 + t9;
t4 = pkin(8) * t60 + t8;
t1 = t42 * t3 - t39 * t4;
t2 = t39 * t3 + t42 * t4;
t48 = t9 * t40 + t8 * t43;
t46 = -pkin(3) * t40 + t52;
t18 = t39 * t40 + t42 * t43;
t38 = t44 ^ 2;
t36 = t41 ^ 2;
t33 = t43 * pkin(7);
t28 = pkin(7) * t58;
t25 = -t43 * pkin(8) + t33;
t22 = t42 * qJ(4) - t39 * t65;
t21 = t39 * qJ(4) + t42 * t65;
t19 = -t39 * t43 + t57;
t14 = t18 * t41;
t12 = (pkin(6) - t46) * t41;
t7 = t42 * t25 + t39 * t49;
t6 = t39 * t25 - t42 * t49;
t5 = (-t65 * t40 - pkin(6) + t52) * t41;
t10 = [1, 0, 0, t36, t50, 0, 0, 0, pkin(1) * t66, pkin(1) * t67, t37 * t36, -0.2e1 * t36 * t59, t56 * t67, t40 * t50, t38, 0.2e1 * t36 * t62 + 0.2e1 * t44 * t55, 0.2e1 * t11 * t44 + 0.2e1 * t36 * t61, 0.2e1 * t12 * t60 + 0.2e1 * t9 * t44, 0.2e1 * (-t40 * t8 + t43 * t9) * t41, -0.2e1 * t12 * t30 - 0.2e1 * t8 * t44, t12 ^ 2 + t8 ^ 2 + t9 ^ 2, t14 ^ 2, t14 * t70, t14 * t66, t44 * t70, t38, 0.2e1 * t1 * t44 + 0.2e1 * t5 * t13, 0.2e1 * t5 * t14 - 0.2e1 * t2 * t44; 0, 0, 0, 0, 0, t41, t44, 0, -t41 * pkin(6), -t44 * pkin(6), t40 * t30, (-t35 + t37) * t41, -t58, -t56, 0, t28 + (-t61 - t64) * t41, pkin(7) * t56 + (t62 - t63) * t41, -t12 * t43 + t23 * t60 + t28, t48, -t12 * t40 + (-pkin(7) * t44 - t23 * t41) * t43, t48 * pkin(7) + t12 * t23, t14 * t19, -t19 * t13 - t14 * t18, t19 * t44, -t18 * t44, 0, t15 * t13 + t5 * t18 - t6 * t44, t15 * t14 + t5 * t19 - t7 * t44; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, 0.2e1 * t59, 0, 0, 0, 0.2e1 * t63, -0.2e1 * t64, t43 * t68, 0.2e1 * t54 * pkin(7), t40 * t68, t54 * pkin(7) ^ 2 + t23 ^ 2, t19 ^ 2, -0.2e1 * t19 * t18, 0, 0, 0, t18 * t69, t19 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t60, -t44, -t55, -t11, -0.2e1 * t34 - t55, t47 * t41, -0.2e1 * t51 + t11, -t9 * pkin(3) + t8 * qJ(4), 0, 0, -t14, t13, -t44, -t21 * t44 - t1, -t22 * t44 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t32, -t33, -t32, t46, t33, t46 * pkin(7), 0, 0, -t19, t18, 0, t6, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t30, 0, t9, 0, 0, 0, 0, 0, t42 * t44, -t39 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t42, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t44, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
