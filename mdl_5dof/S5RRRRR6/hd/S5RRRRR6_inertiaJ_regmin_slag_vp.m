% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t46 = sin(qJ(4));
t47 = sin(qJ(3));
t50 = cos(qJ(4));
t51 = cos(qJ(3));
t28 = t46 * t47 - t50 * t51;
t40 = -t51 * pkin(3) - pkin(2);
t19 = t28 * pkin(4) + t40;
t58 = cos(qJ(2)) * pkin(1);
t18 = t19 - t58;
t68 = 0.2e1 * t18;
t67 = 0.2e1 * t19;
t31 = t40 - t58;
t66 = 0.2e1 * t31;
t65 = 0.2e1 * t40;
t64 = 0.2e1 * t51;
t29 = t46 * t51 + t50 * t47;
t63 = t29 * pkin(9);
t45 = sin(qJ(5));
t62 = t45 * pkin(4);
t61 = t46 * pkin(3);
t60 = sin(qJ(2)) * pkin(1);
t59 = t51 * pkin(7);
t39 = -pkin(2) - t58;
t57 = pkin(2) - t39;
t37 = pkin(7) + t60;
t56 = t51 * t37;
t55 = t18 + t19;
t54 = t31 + t40;
t49 = cos(qJ(5));
t53 = t49 * t61;
t25 = (-pkin(8) - t37) * t47;
t43 = t51 * pkin(8);
t26 = t43 + t56;
t10 = t50 * t25 - t46 * t26;
t32 = (-pkin(7) - pkin(8)) * t47;
t33 = t43 + t59;
t16 = t50 * t32 - t46 * t33;
t42 = t50 * pkin(3);
t38 = t42 + pkin(4);
t21 = t49 * t38 - t45 * t61;
t11 = -t46 * t25 - t50 * t26;
t17 = -t46 * t32 - t50 * t33;
t44 = t47 ^ 2;
t41 = t49 * pkin(4);
t35 = t47 * t64;
t27 = t29 ^ 2;
t24 = t28 * pkin(9);
t22 = -t45 * t38 - t53;
t15 = -0.2e1 * t29 * t28;
t14 = -t45 * t28 + t49 * t29;
t13 = t49 * t28 + t45 * t29;
t12 = t14 ^ 2;
t9 = -t17 - t24;
t8 = t16 - t63;
t7 = -t11 - t24;
t6 = t10 - t63;
t5 = -0.2e1 * t14 * t13;
t4 = -t45 * t8 - t49 * t9;
t3 = -t45 * t9 + t49 * t8;
t2 = -t45 * t6 - t49 * t7;
t1 = -t45 * t7 + t49 * t6;
t20 = [1, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t60, t44, t35, 0, 0, 0, -0.2e1 * t39 * t51, 0.2e1 * t39 * t47, t27, t15, 0, 0, 0, t28 * t66, t29 * t66, t12, t5, 0, 0, 0, t13 * t68, t14 * t68; 0, 0, 0, 1, t58, -t60, t44, t35, 0, 0, 0, t57 * t51, -t57 * t47, t27, t15, 0, 0, 0, t54 * t28, t54 * t29, t12, t5, 0, 0, 0, t55 * t13, t55 * t14; 0, 0, 0, 1, 0, 0, t44, t35, 0, 0, 0, pkin(2) * t64, -0.2e1 * pkin(2) * t47, t27, t15, 0, 0, 0, t28 * t65, t29 * t65, t12, t5, 0, 0, 0, t13 * t67, t14 * t67; 0, 0, 0, 0, 0, 0, 0, 0, t47, t51, 0, -t47 * t37, -t56, 0, 0, t29, -t28, 0, t10, t11, 0, 0, t14, -t13, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t47, t51, 0, -t47 * pkin(7), -t59, 0, 0, t29, -t28, 0, t16, t17, 0, 0, t14, -t13, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t61, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t10, t11, 0, 0, t14, -t13, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t16, t17, 0, 0, t14, -t13, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t42, -t61, 0, 0, 0, 0, 1, t21 + t41, -t53 + (-pkin(4) - t38) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t20;
