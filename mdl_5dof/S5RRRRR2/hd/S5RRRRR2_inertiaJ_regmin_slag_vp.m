% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t34 = sin(qJ(4));
t35 = sin(qJ(3));
t38 = cos(qJ(4));
t39 = cos(qJ(3));
t24 = t34 * t39 + t38 * t35;
t65 = -0.2e1 * t24;
t57 = cos(qJ(2)) * pkin(1);
t58 = t39 * pkin(2);
t26 = -t57 - t58;
t64 = 0.2e1 * t26;
t60 = sin(qJ(2)) * pkin(1);
t12 = t24 * t60;
t23 = t34 * t35 - t38 * t39;
t48 = t39 * t60;
t53 = t35 * t60;
t13 = -t34 * t53 + t38 * t48;
t33 = sin(qJ(5));
t37 = cos(qJ(5));
t3 = -t33 * t13 + t37 * t26;
t55 = t33 * t24;
t63 = t12 * t55 + t3 * t23;
t20 = t37 * t24;
t4 = t37 * t13 + t33 * t26;
t62 = t12 * t20 - t4 * t23;
t61 = t34 * pkin(2);
t59 = t38 * pkin(2);
t56 = t12 * t37;
t54 = t33 * t37;
t10 = t23 * t65;
t52 = t37 * t61;
t51 = t33 * t59;
t50 = t37 * t59;
t49 = t23 * t58;
t47 = t35 * t57;
t46 = t39 * t57;
t45 = t26 - t58;
t44 = t33 * t49;
t43 = t37 * t49;
t42 = -t23 * t52 - t24 * t50;
t41 = (-t23 * t34 - t24 * t38) * t33 * pkin(2);
t32 = t37 ^ 2;
t31 = t35 ^ 2;
t30 = t33 ^ 2;
t28 = 0.2e1 * t35 * t39;
t27 = 0.2e1 * t54;
t22 = t24 ^ 2;
t21 = t23 ^ 2;
t19 = t37 * t23;
t18 = t32 * t22;
t17 = t33 * t23;
t16 = t33 * t20;
t15 = -0.2e1 * t22 * t54;
t11 = t12 * t33;
t9 = 0.2e1 * t23 * t20;
t8 = t33 * t10;
t7 = (-t30 + t32) * t24;
t1 = [1, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t60, t31, t28, 0, 0, 0, 0.2e1 * t46, -0.2e1 * t47, t22, t10, 0, 0, 0, t23 * t64, t24 * t64, t18, t15, t9, t8, t21, 0.2e1 * t63, 0.2e1 * t62; 0, 0, 0, 1, t57, -t60, t31, t28, 0, 0, 0, t46, -t47, t22, t10, 0, 0, 0, t45 * t23, t45 * t24, t18, t15, t9, t8, t21, -t43 + t63, t44 + t62; 0, 0, 0, 1, 0, 0, t31, t28, 0, 0, 0, 0, 0, t22, t10, 0, 0, 0, -0.2e1 * t49, t58 * t65, t18, t15, t9, t8, t21, -0.2e1 * t43, 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, t35, t39, 0, -t53, -t48, 0, 0, t24, -t23, 0, -t12, -t13, t16, t7, t17, t19, 0, t41 - t56, t11 + t42; 0, 0, 0, 0, 0, 0, 0, 0, t35, t39, 0, 0, 0, 0, 0, t24, -t23, 0, 0, 0, t16, t7, t17, t19, 0, t41, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t61, t30, t27, 0, 0, 0, 0.2e1 * t50, -0.2e1 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t12, -t13, t16, t7, t17, t19, 0, -t56, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, 0, 0, t16, t7, t17, t19, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t61, t30, t27, 0, 0, 0, t50, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t30, t27, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t55, t23, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t55, t23, -t37 * t58, t33 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t37, 0, -t33 * t61, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t37, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t1;
