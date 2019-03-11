% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t40 = sin(qJ(2));
t42 = cos(qJ(2));
t21 = t37 * t40 - t38 * t42;
t22 = t37 * t42 + t38 * t40;
t34 = -t42 * pkin(2) - pkin(1);
t45 = -t22 * qJ(4) + t34;
t11 = t21 * pkin(3) + t45;
t66 = -0.2e1 * t11;
t65 = 0.2e1 * t21;
t30 = t37 * pkin(2) + qJ(4);
t64 = 0.2e1 * t30;
t63 = 0.2e1 * t42;
t39 = sin(qJ(5));
t62 = t39 * pkin(5);
t53 = -qJ(3) - pkin(7);
t26 = t53 * t40;
t27 = t53 * t42;
t12 = -t38 * t26 - t37 * t27;
t8 = t22 * pkin(4) + t12;
t61 = t39 * t8;
t41 = cos(qJ(5));
t60 = t41 * pkin(5);
t59 = t30 * t21;
t35 = t39 ^ 2;
t58 = t35 * t21;
t57 = t39 * t21;
t56 = t39 * t22;
t55 = t41 * t21;
t17 = t41 * t22;
t54 = t41 * t39;
t52 = t22 * t65;
t14 = t37 * t26 - t38 * t27;
t51 = t12 ^ 2 + t14 ^ 2;
t33 = -t38 * pkin(2) - pkin(3);
t6 = (pkin(3) + pkin(8)) * t21 + t45;
t50 = qJ(6) * t21 + t6;
t7 = t41 * t8;
t1 = t22 * pkin(5) - t50 * t39 + t7;
t2 = t50 * t41 + t61;
t49 = t1 * t41 + t2 * t39;
t48 = -t1 * t39 + t2 * t41;
t29 = -pkin(8) + t33;
t18 = (-qJ(6) + t29) * t39;
t25 = t41 * t29;
t19 = -t41 * qJ(6) + t25;
t47 = t18 * t41 - t19 * t39;
t46 = -t22 * t29 + t59;
t44 = 0.2e1 * t12 * t22 - 0.2e1 * t14 * t21;
t36 = t41 ^ 2;
t28 = -t35 - t36;
t24 = t30 + t62;
t20 = t21 ^ 2;
t16 = t36 * t21;
t10 = t18 * t39 + t19 * t41;
t9 = -t21 * pkin(4) + t14;
t5 = (-pkin(4) - t60) * t21 + t14;
t4 = t41 * t6 + t61;
t3 = -t39 * t6 + t7;
t13 = [1, 0, 0, t40 ^ 2, t40 * t63, 0, 0, 0, pkin(1) * t63, -0.2e1 * pkin(1) * t40, t44, t34 ^ 2 + t51, t44, t21 * t66, t22 * t66, t11 ^ 2 + t51, t35 * t20, 0.2e1 * t20 * t54, t39 * t52, t41 * t52, t22 ^ 2, 0.2e1 * t3 * t22 - 0.2e1 * t9 * t55, -0.2e1 * t4 * t22 + 0.2e1 * t9 * t57, t48 * t65, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t40, t42, 0, -t40 * pkin(7), -t42 * pkin(7) (-t21 * t37 - t22 * t38) * pkin(2) (-t12 * t38 + t14 * t37) * pkin(2), t33 * t22 - t59, t12, t14, t12 * t33 + t14 * t30, t21 * t54, t16 - t58, t17, -t56, 0, t9 * t39 - t46 * t41, t46 * t39 + t9 * t41, t47 * t21 - t49, t1 * t19 + t2 * t18 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t37 ^ 2 + t38 ^ 2) * pkin(2) ^ 2, 0, 0.2e1 * t33, t64, t30 ^ 2 + t33 ^ 2, t36, -0.2e1 * t54, 0, 0, 0, t39 * t64, t41 * t64, -0.2e1 * t10, t18 ^ 2 + t19 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t21, -t22, t11, 0, 0, 0, 0, 0, -t56, -t17, t16 + t58, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, t12, 0, 0, 0, 0, 0, t17, -t56, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t33, 0, 0, 0, 0, 0, 0, 0, t28, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t55, t22, t3, -t4, -pkin(5) * t57, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, t25, -t39 * t29, -t60, t19 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t41, 0, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t13;
