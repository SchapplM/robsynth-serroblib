% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t35 = sin(qJ(4));
t64 = -0.2e1 * t35;
t38 = cos(qJ(4));
t63 = 0.2e1 * t38;
t62 = 2 * qJ(5);
t40 = -pkin(4) - pkin(10);
t28 = sin(pkin(12));
t30 = sin(pkin(6));
t33 = cos(pkin(6));
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t31 = cos(pkin(12));
t32 = cos(pkin(7));
t57 = t31 * t32;
t29 = sin(pkin(7));
t58 = t29 * t39;
t6 = -t33 * t58 + (t28 * t36 - t39 * t57) * t30;
t61 = t6 * t35;
t60 = t6 * t38;
t59 = t29 * t36;
t34 = sin(qJ(6));
t56 = t34 * t35;
t55 = t34 * t38;
t54 = t35 * t38;
t37 = cos(qJ(6));
t53 = t37 * t34;
t52 = t37 * t38;
t25 = t35 ^ 2;
t27 = t38 ^ 2;
t51 = t25 + t27;
t50 = qJ(5) * t38;
t49 = -0.2e1 * t54;
t48 = t35 * t58;
t47 = t38 * t58;
t46 = -t35 * qJ(5) - pkin(3);
t12 = -t30 * t31 * t29 + t33 * t32;
t7 = t33 * t59 + (t28 * t39 + t36 * t57) * t30;
t3 = -t12 * t38 + t7 * t35;
t4 = t12 * t35 + t7 * t38;
t45 = t3 * t35 + t4 * t38;
t44 = -pkin(4) * t35 + t50;
t13 = -t38 * t32 + t35 * t59;
t14 = t35 * t32 + t38 * t59;
t43 = t13 * t35 + t14 * t38;
t42 = t35 * t40 + t50;
t26 = t37 ^ 2;
t24 = t34 ^ 2;
t22 = t38 * pkin(9);
t21 = t35 * pkin(9);
t20 = t37 * t35;
t19 = t38 * pkin(5) + t22;
t18 = t35 * pkin(5) + t21;
t17 = -t38 * pkin(4) + t46;
t16 = t40 * t38 + t46;
t11 = -t34 * t13 + t37 * t58;
t10 = t37 * t13 + t34 * t58;
t9 = t37 * t16 + t34 * t18;
t8 = -t34 * t16 + t37 * t18;
t2 = t3 * t34 + t6 * t37;
t1 = t3 * t37 - t6 * t34;
t5 = [1, t33 ^ 2 + (t28 ^ 2 + t31 ^ 2) * t30 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t6 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t13 + t4 * t14 - t6 * t58, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 ^ 2 * t39 ^ 2 + t13 ^ 2 + t14 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t60, t61, t45, t60, -t61, t45 * pkin(9) + t6 * t17, 0, 0, 0, 0, 0, t1 * t35 + t4 * t52, -t2 * t35 - t4 * t55; 0, 0, 0, t58, -t59, 0, 0, 0, 0, 0, t47, -t48, t43, -t47, t48, t43 * pkin(9) - t17 * t58, 0, 0, 0, 0, 0, t10 * t35 + t14 * t52, t11 * t35 - t14 * t55; 0, 0, 1, 0, 0, t25, 0.2e1 * t54, 0, 0, 0, pkin(3) * t63, pkin(3) * t64, 0.2e1 * t51 * pkin(9), t17 * t63, t17 * t64, t51 * pkin(9) ^ 2 + t17 ^ 2, t24 * t27, 0.2e1 * t27 * t53, t34 * t49, t37 * t49, t25, 0.2e1 * t19 * t52 + 0.2e1 * t8 * t35, -0.2e1 * t19 * t55 - 0.2e1 * t9 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t3, t4, -t3 * pkin(4) + t4 * qJ(5), 0, 0, 0, 0, 0, t4 * t34, t4 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, t13, t14, -t13 * pkin(4) + t14 * qJ(5), 0, 0, 0, 0, 0, t14 * t34, t14 * t37; 0, 0, 0, 0, 0, 0, 0, t35, t38, 0, -t21, -t22, t44, t21, t22, t44 * pkin(9), -t34 * t52 (t24 - t26) * t38, t20, -t56, 0, t19 * t34 + t42 * t37, t19 * t37 - t42 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t62, pkin(4) ^ 2 + (qJ(5) ^ 2) t26, -0.2e1 * t53, 0, 0, 0, t34 * t62, t37 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, t21, 0, 0, 0, 0, 0, t20, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t52, t35, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t34, 0, t37 * t40, -t34 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
