% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t51 = sin(qJ(3));
t54 = cos(qJ(3));
t33 = -t48 * t51 + t49 * t54;
t45 = -t54 * pkin(3) - pkin(2);
t23 = -t33 * pkin(4) + t45;
t61 = cos(qJ(2)) * pkin(1);
t22 = t23 - t61;
t69 = 0.2e1 * t22;
t68 = 0.2e1 * t23;
t67 = 0.2e1 * t54;
t63 = sin(qJ(2)) * pkin(1);
t43 = pkin(7) + t63;
t31 = (-qJ(4) - t43) * t51;
t46 = t54 * qJ(4);
t59 = t54 * t43;
t32 = t46 + t59;
t14 = t49 * t31 - t48 * t32;
t15 = t48 * t31 + t49 * t32;
t34 = t48 * t54 + t49 * t51;
t66 = -t14 * t34 + t15 * t33;
t65 = pkin(3) * t48;
t64 = t34 * pkin(8);
t62 = t54 * pkin(7);
t44 = -pkin(2) - t61;
t60 = pkin(2) - t44;
t39 = (-qJ(4) - pkin(7)) * t51;
t40 = t46 + t62;
t20 = t49 * t39 - t48 * t40;
t21 = t48 * t39 + t49 * t40;
t58 = -t20 * t34 + t21 * t33;
t57 = t22 + t23;
t53 = cos(qJ(5));
t50 = sin(qJ(5));
t47 = t51 ^ 2;
t42 = t49 * pkin(3) + pkin(4);
t41 = t51 * t67;
t38 = t45 - t61;
t30 = t33 * pkin(8);
t28 = -t50 * t42 - t53 * t65;
t27 = t53 * t42 - t50 * t65;
t19 = t50 * t33 + t53 * t34;
t18 = -t53 * t33 + t50 * t34;
t17 = t19 ^ 2;
t16 = (t33 * t48 - t34 * t49) * pkin(3);
t11 = t30 + t21;
t10 = t20 - t64;
t7 = t30 + t15;
t6 = t14 - t64;
t5 = -0.2e1 * t19 * t18;
t4 = -t50 * t10 - t53 * t11;
t3 = t53 * t10 - t50 * t11;
t2 = -t50 * t6 - t53 * t7;
t1 = -t50 * t7 + t53 * t6;
t8 = [1, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t63, t47, t41, 0, 0, 0, -0.2e1 * t44 * t54, 0.2e1 * t44 * t51, 0.2e1 * t66, t14 ^ 2 + t15 ^ 2 + t38 ^ 2, t17, t5, 0, 0, 0, t18 * t69, t19 * t69; 0, 0, 0, 1, t61, -t63, t47, t41, 0, 0, 0, t60 * t54, -t60 * t51, t58 + t66, t14 * t20 + t15 * t21 + t38 * t45, t17, t5, 0, 0, 0, t57 * t18, t57 * t19; 0, 0, 0, 1, 0, 0, t47, t41, 0, 0, 0, pkin(2) * t67, -0.2e1 * pkin(2) * t51, 0.2e1 * t58, t20 ^ 2 + t21 ^ 2 + t45 ^ 2, t17, t5, 0, 0, 0, t18 * t68, t19 * t68; 0, 0, 0, 0, 0, 0, 0, 0, t51, t54, 0, -t51 * t43, -t59, t16, (t14 * t49 + t15 * t48) * pkin(3), 0, 0, t19, -t18, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t51, t54, 0, -t51 * pkin(7), -t62, t16, (t20 * t49 + t21 * t48) * pkin(3), 0, 0, t19, -t18, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t48 ^ 2 + t49 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t27, 0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
