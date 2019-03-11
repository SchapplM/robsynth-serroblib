% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(pkin(11));
t46 = cos(pkin(11));
t50 = sin(qJ(3));
t62 = cos(qJ(3));
t29 = t44 * t50 - t46 * t62;
t27 = t29 ^ 2;
t71 = -0.2e1 * t29;
t70 = 0.2e1 * t29;
t39 = -t46 * pkin(3) - pkin(4);
t52 = cos(qJ(5));
t34 = -t52 * pkin(5) + t39;
t69 = 0.2e1 * t34;
t68 = 0.2e1 * t50;
t67 = t29 * pkin(5);
t48 = sin(qJ(6));
t66 = t48 * pkin(5);
t51 = cos(qJ(6));
t65 = t51 * pkin(5);
t31 = t44 * t62 + t46 * t50;
t47 = cos(pkin(10));
t40 = -t47 * pkin(1) - pkin(2);
t35 = -t62 * pkin(3) + t40;
t10 = t29 * pkin(4) - t31 * pkin(8) + t35;
t49 = sin(qJ(5));
t45 = sin(pkin(10));
t38 = t45 * pkin(1) + pkin(7);
t57 = t62 * t38;
t24 = t62 * qJ(4) + t57;
t56 = (-qJ(4) - t38) * t50;
t15 = t46 * t24 + t44 * t56;
t59 = t52 * t15;
t5 = t59 + (-pkin(9) * t31 + t10) * t49;
t64 = t51 * t5;
t37 = t44 * pkin(3) + pkin(8);
t63 = pkin(9) + t37;
t32 = t48 * t49 - t51 * t52;
t18 = t29 * t32;
t33 = t48 * t52 + t51 * t49;
t19 = t33 * t29;
t22 = t49 * t29;
t61 = t49 * t31;
t60 = t49 * t52;
t58 = t52 * t31;
t6 = t52 * t10 - t49 * t15;
t4 = -pkin(9) * t58 + t6 + t67;
t1 = t51 * t4 - t48 * t5;
t13 = t44 * t24 - t46 * t56;
t55 = -t29 * t37 + t31 * t39;
t43 = t52 ^ 2;
t42 = t49 ^ 2;
t28 = t31 ^ 2;
t26 = t63 * t52;
t25 = t63 * t49;
t23 = t52 * t29;
t17 = -t48 * t25 + t51 * t26;
t16 = -t51 * t25 - t48 * t26;
t12 = -t48 * t61 + t51 * t58;
t11 = t33 * t31;
t8 = pkin(5) * t61 + t13;
t7 = t49 * t10 + t59;
t2 = t48 * t4 + t64;
t3 = [1, 0, 0 (t45 ^ 2 + t47 ^ 2) * pkin(1) ^ 2, t50 ^ 2, t62 * t68, 0, 0, 0, -0.2e1 * t40 * t62, t40 * t68, 0.2e1 * t13 * t31 - 0.2e1 * t15 * t29, t13 ^ 2 + t15 ^ 2 + t35 ^ 2, t43 * t28, -0.2e1 * t28 * t60, t58 * t70, t61 * t71, t27, 0.2e1 * t13 * t61 + 0.2e1 * t6 * t29, 0.2e1 * t13 * t58 - 0.2e1 * t7 * t29, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t70, t11 * t71, t27, 0.2e1 * t1 * t29 + 0.2e1 * t8 * t11, 0.2e1 * t8 * t12 - 0.2e1 * t2 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t29 + t15 * t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t28 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t50, t62, 0, -t50 * t38, -t57 (-t29 * t44 - t31 * t46) * pkin(3) (-t13 * t46 + t15 * t44) * pkin(3), t49 * t58 (-t42 + t43) * t31, t22, t23, 0, -t13 * t52 + t55 * t49, t13 * t49 + t55 * t52, t12 * t33, -t33 * t11 - t12 * t32, t19, -t18, 0, t34 * t11 + t16 * t29 + t8 * t32, t34 * t12 - t17 * t29 + t8 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t50, 0 (-t29 * t46 + t31 * t44) * pkin(3), 0, 0, 0, 0, 0, -t23, t22, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t44 ^ 2 + t46 ^ 2) * pkin(3) ^ 2, t42, 0.2e1 * t60, 0, 0, 0, -0.2e1 * t39 * t52, 0.2e1 * t39 * t49, t33 ^ 2, -0.2e1 * t33 * t32, 0, 0, 0, t32 * t69, t33 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, t23, -t22, 0, 0, 0, 0, 0, -t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t61, t29, t6, -t7, 0, 0, t12, -t11, t29, t29 * t65 + t1, -t64 + (-t4 - t67) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t58, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t52, 0, -t49 * t37, -t52 * t37, 0, 0, t33, -t32, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t49, 0, 0, 0, 0, 0, -t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
