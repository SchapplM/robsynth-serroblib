% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t27 = sin(qJ(3));
t31 = pkin(3) + pkin(4);
t29 = cos(qJ(3));
t44 = t29 * qJ(4);
t63 = t31 * t27 - t44;
t45 = t27 * qJ(4);
t62 = t29 * t31 + t45;
t9 = pkin(2) + t62;
t61 = 0.2e1 * t9;
t38 = -t29 * pkin(3) - t45;
t12 = -pkin(2) + t38;
t60 = -0.2e1 * t12;
t28 = sin(qJ(2));
t59 = -0.2e1 * t28;
t58 = 0.2e1 * t28;
t30 = cos(qJ(2));
t57 = 0.2e1 * t30;
t56 = pkin(2) * t27;
t55 = pkin(2) * t29;
t54 = pkin(6) * t27;
t53 = pkin(6) * t29;
t52 = t27 * pkin(7);
t51 = t27 * t28;
t50 = t27 * t29;
t49 = t27 * t30;
t20 = t29 * t28;
t48 = t29 * t30;
t13 = -t30 * pkin(2) - t28 * pkin(7) - pkin(1);
t47 = pkin(6) * t49 - t29 * t13;
t7 = pkin(6) * t48 + t27 * t13;
t24 = t27 ^ 2;
t26 = t29 ^ 2;
t46 = t24 + t26;
t43 = t29 * qJ(5);
t42 = t30 * qJ(4);
t41 = t28 * t57;
t23 = t30 * pkin(3);
t5 = t23 + t47;
t40 = -0.2e1 * t42 + t7;
t4 = -t42 + t7;
t39 = t5 * t27 + t4 * t29;
t37 = -pkin(3) * t27 + t44;
t36 = t28 * t43 - t5;
t34 = qJ(4) ^ 2;
t33 = 0.2e1 * qJ(4);
t25 = t28 ^ 2;
t22 = t29 * pkin(7);
t18 = pkin(7) * t49;
t16 = qJ(5) * t51;
t15 = t22 - t43;
t14 = (pkin(7) - qJ(5)) * t27;
t8 = (pkin(6) - t37) * t28;
t3 = (-pkin(6) - t63) * t28;
t2 = t16 + t4;
t1 = t30 * pkin(4) - t36;
t6 = [1, 0, 0, t25, t41, 0, 0, 0, pkin(1) * t57, pkin(1) * t59, t26 * t25, -0.2e1 * t25 * t50, t48 * t59, t27 * t41, t30 ^ 2, 0.2e1 * t25 * t54 + 0.2e1 * t30 * t47, 0.2e1 * t25 * t53 + 0.2e1 * t7 * t30, 0.2e1 * t5 * t30 + 0.2e1 * t8 * t51, (-t27 * t4 + t29 * t5) * t58, -0.2e1 * t8 * t20 - 0.2e1 * t4 * t30, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, 0.2e1 * t1 * t30 - 0.2e1 * t3 * t51, -0.2e1 * t2 * t30 + 0.2e1 * t3 * t20, (-t1 * t29 + t2 * t27) * t58, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t28, t30, 0, -t28 * pkin(6), -t30 * pkin(6), t27 * t20, (-t24 + t26) * t28, -t49, -t48, 0, t18 + (-t53 - t56) * t28, pkin(7) * t48 + (t54 - t55) * t28, t12 * t51 - t8 * t29 + t18, t39, -t8 * t27 + (-pkin(7) * t30 - t12 * t28) * t29, pkin(7) * t39 + t8 * t12, t14 * t30 + t3 * t29 - t9 * t51, -t15 * t30 + t9 * t20 + t3 * t27, (-t14 * t28 - t2) * t29 + (t15 * t28 - t1) * t27, t1 * t14 + t2 * t15 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, 0.2e1 * t50, 0, 0, 0, 0.2e1 * t55, -0.2e1 * t56, t29 * t60, 0.2e1 * t46 * pkin(7), t27 * t60, t46 * pkin(7) ^ 2 + t12 ^ 2, t29 * t61, t27 * t61, -0.2e1 * t14 * t27 - 0.2e1 * t15 * t29, t14 ^ 2 + t15 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t51, -t30, -t47, -t7, -0.2e1 * t23 - t47, t38 * t28, t40, -t5 * pkin(3) + t4 * qJ(4), (-pkin(4) - t31) * t30 + t36, t16 + t40, t62 * t28, t2 * qJ(4) - t1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, -t52, -t22, -t52, t37, t22, t37 * pkin(7), -t14, t15, t63, t15 * qJ(4) - t14 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t33, pkin(3) ^ 2 + t34, 0.2e1 * t31, t33, 0, t31 ^ 2 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t20, 0, t5, t30, 0, -t20, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t52, 0, 0, -t27, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -1, 0, 0, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t20, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t27, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
