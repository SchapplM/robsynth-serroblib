% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR7
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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t68 = t56 ^ 2 + t57 ^ 2;
t69 = t68 * qJ(4);
t59 = sin(qJ(3));
t60 = sin(qJ(2));
t75 = cos(qJ(3));
t76 = cos(qJ(2));
t34 = t59 * t60 - t75 * t76;
t83 = -0.2e1 * t34;
t65 = t75 * pkin(2);
t49 = -t65 - pkin(3);
t79 = t57 * pkin(4);
t37 = t49 - t79;
t82 = 0.2e1 * t37;
t47 = -pkin(3) - t79;
t81 = 0.2e1 * t47;
t50 = -t76 * pkin(2) - pkin(1);
t80 = 0.2e1 * t50;
t78 = t59 * pkin(2);
t77 = pkin(3) - t49;
t42 = (-pkin(7) - pkin(6)) * t60;
t66 = t76 * pkin(6);
t43 = t76 * pkin(7) + t66;
t27 = -t75 * t42 + t59 * t43;
t74 = t27 * t57;
t35 = t59 * t76 + t75 * t60;
t73 = t56 * t35;
t72 = t57 * t35;
t19 = t34 * pkin(3) - t35 * qJ(4) + t50;
t28 = t59 * t42 + t75 * t43;
t8 = t56 * t19 + t57 * t28;
t71 = t37 + t47;
t46 = qJ(4) + t78;
t70 = t68 * t46;
t67 = 0.2e1 * t76;
t7 = t57 * t19 - t56 * t28;
t3 = -t7 * t56 + t8 * t57;
t64 = -pkin(3) * t35 - qJ(4) * t34;
t63 = -t34 * t46 + t35 * t49;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t33 = t61 * t56 + t58 * t57;
t32 = t58 * t56 - t61 * t57;
t53 = t57 * pkin(8);
t39 = t57 * qJ(4) + t53;
t38 = (-pkin(8) - qJ(4)) * t56;
t31 = t33 ^ 2;
t30 = t57 * t46 + t53;
t29 = (-pkin(8) - t46) * t56;
t26 = t33 * t34;
t25 = t32 * t34;
t24 = t58 * t38 + t61 * t39;
t23 = t61 * t38 - t58 * t39;
t22 = -0.2e1 * t33 * t32;
t21 = t27 * t56;
t18 = t58 * t29 + t61 * t30;
t17 = t61 * t29 - t58 * t30;
t14 = t32 * t35;
t13 = t33 * t35;
t12 = pkin(4) * t73 + t27;
t11 = t14 * t33;
t10 = t12 * t33;
t9 = t12 * t32;
t6 = -pkin(8) * t73 + t8;
t5 = t34 * pkin(4) - pkin(8) * t72 + t7;
t4 = -t33 * t13 + t14 * t32;
t2 = t58 * t5 + t61 * t6;
t1 = t61 * t5 - t58 * t6;
t15 = [1, 0, 0, t60 ^ 2, t60 * t67, 0, 0, 0, pkin(1) * t67, -0.2e1 * pkin(1) * t60, t35 ^ 2, t35 * t83, 0, 0, 0, t34 * t80, t35 * t80, 0.2e1 * t27 * t73 + 0.2e1 * t7 * t34, 0.2e1 * t27 * t72 - 0.2e1 * t8 * t34, 0.2e1 * (-t56 * t8 - t57 * t7) * t35, t27 ^ 2 + t7 ^ 2 + t8 ^ 2, t14 ^ 2, 0.2e1 * t14 * t13, t14 * t83, t13 * t83, t34 ^ 2, 0.2e1 * t1 * t34 + 0.2e1 * t12 * t13, -0.2e1 * t12 * t14 - 0.2e1 * t2 * t34; 0, 0, 0, 0, 0, t60, t76, 0, -t60 * pkin(6), -t66, 0, 0, t35, -t34, 0, -t27, -t28, t63 * t56 - t74, t63 * t57 + t21, t3, t27 * t49 + t3 * t46, -t11, t4, t26, -t25, 0, t37 * t13 + t17 * t34 + t9, -t37 * t14 - t18 * t34 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t78, -0.2e1 * t49 * t57, 0.2e1 * t49 * t56, 0.2e1 * t70, t68 * t46 ^ 2 + t49 ^ 2, t31, t22, 0, 0, 0, t32 * t82, t33 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, -t27, -t28, t64 * t56 - t74, t64 * t57 + t21, t3, -t27 * pkin(3) + t3 * qJ(4), -t11, t4, t26, -t25, 0, t47 * t13 + t23 * t34 + t9, -t47 * t14 - t24 * t34 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t78, t77 * t57, -t77 * t56, t69 + t70, -t49 * pkin(3) + t46 * t69, t31, t22, 0, 0, 0, t71 * t32, t71 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t57, -0.2e1 * pkin(3) * t56, 0.2e1 * t69, t68 * qJ(4) ^ 2 + pkin(3) ^ 2, t31, t22, 0, 0, 0, t32 * t81, t33 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, t27, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, 0, t49, 0, 0, 0, 0, 0, t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, 0, -pkin(3), 0, 0, 0, 0, 0, t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t13, t34, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
