% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = cos(qJ(4));
t49 = sin(qJ(4));
t64 = t49 * qJ(5);
t79 = pkin(4) + pkin(5);
t19 = t79 * t52 + pkin(3) + t64;
t84 = 0.2e1 * t19;
t57 = -t52 * pkin(4) - t64;
t27 = -pkin(3) + t57;
t83 = -0.2e1 * t27;
t50 = sin(qJ(3));
t82 = 0.2e1 * t50;
t53 = cos(qJ(3));
t81 = -0.2e1 * t53;
t80 = 0.2e1 * t53;
t78 = pkin(3) * t49;
t77 = pkin(3) * t52;
t76 = pkin(4) * t49;
t46 = sin(pkin(10));
t32 = t46 * pkin(1) + pkin(7);
t75 = t32 * t49;
t74 = t32 * t52;
t42 = t49 ^ 2;
t73 = t42 * t50;
t72 = t49 * t50;
t71 = t49 * t52;
t70 = t49 * t53;
t51 = cos(qJ(6));
t69 = t51 * t49;
t36 = t52 * t50;
t37 = t52 * t53;
t48 = sin(qJ(6));
t20 = t48 * t49 + t51 * t52;
t68 = t53 * t20;
t67 = t53 * t32;
t47 = cos(pkin(10));
t33 = -t47 * pkin(1) - pkin(2);
t18 = -t53 * pkin(3) - t50 * pkin(8) + t33;
t66 = -t52 * t18 + t49 * t67;
t9 = t49 * t18 + t52 * t67;
t44 = t52 ^ 2;
t65 = t42 + t44;
t63 = t52 * qJ(5);
t62 = t53 * qJ(5);
t61 = t50 * t80;
t41 = t53 * pkin(4);
t6 = t41 + t66;
t60 = t65 * pkin(8);
t39 = t49 * pkin(8);
t59 = -t49 * pkin(9) + t39;
t5 = -t62 + t9;
t3 = t53 * pkin(5) - pkin(9) * t36 + t6;
t4 = pkin(9) * t72 + t5;
t1 = t51 * t3 - t48 * t4;
t2 = t48 * t3 + t51 * t4;
t58 = t6 * t49 + t5 * t52;
t56 = t63 - t76;
t45 = t53 ^ 2;
t43 = t50 ^ 2;
t40 = t52 * pkin(8);
t35 = t44 * t50;
t34 = t44 * t43;
t31 = pkin(8) * t70;
t30 = t50 * t63;
t28 = -t52 * pkin(9) + t40;
t26 = t51 * qJ(5) - t48 * t79;
t25 = t48 * qJ(5) + t51 * t79;
t21 = -t48 * t52 + t69;
t17 = t21 * t53;
t14 = t20 * t50;
t13 = t48 * t36 - t50 * t69;
t12 = -t30 + (t32 + t76) * t50;
t11 = t51 * t28 + t48 * t59;
t10 = t48 * t28 - t51 * t59;
t7 = t30 + (-t79 * t49 - t32) * t50;
t8 = [1, 0, 0 (t46 ^ 2 + t47 ^ 2) * pkin(1) ^ 2, t43, t61, 0, 0, 0, t33 * t81, t33 * t82, t34, -0.2e1 * t43 * t71, -0.2e1 * t50 * t37, t49 * t61, t45, 0.2e1 * t43 * t75 + 0.2e1 * t53 * t66, 0.2e1 * t43 * t74 + 0.2e1 * t9 * t53, 0.2e1 * t12 * t72 + 0.2e1 * t6 * t53 (-t49 * t5 + t52 * t6) * t82, -0.2e1 * t12 * t36 - 0.2e1 * t5 * t53, t12 ^ 2 + t5 ^ 2 + t6 ^ 2, t14 ^ 2, -0.2e1 * t14 * t13, t14 * t80, t13 * t81, t45, 0.2e1 * t1 * t53 + 0.2e1 * t7 * t13, 0.2e1 * t7 * t14 - 0.2e1 * t2 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t53 + t58 * t50, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t43 + t34 + t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t50, t53, 0, -t50 * t32, -t67, t49 * t36, t35 - t73, -t70, -t37, 0, t31 + (-t74 - t78) * t50, pkin(8) * t37 + (t75 - t77) * t50, -t12 * t52 + t27 * t72 + t31, t58, -t12 * t49 + (-pkin(8) * t53 - t27 * t50) * t52, t58 * pkin(8) + t12 * t27, t14 * t21, -t21 * t13 - t14 * t20, t17, -t68, 0, -t10 * t53 + t19 * t13 + t7 * t20, -t11 * t53 + t19 * t14 + t7 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t50, 0, 0, 0, 0, 0, t37, -t70, t37, t35 + t73, t70, -t53 * t27 + t50 * t60, 0, 0, 0, 0, 0, t68, t17; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t42, 0.2e1 * t71, 0, 0, 0, 0.2e1 * t77, -0.2e1 * t78, t52 * t83, 0.2e1 * t60, t49 * t83, t65 * pkin(8) ^ 2 + t27 ^ 2, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t84, t21 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t72, -t53, -t66, -t9, -0.2e1 * t41 - t66, t57 * t50, -0.2e1 * t62 + t9, -t6 * pkin(4) + t5 * qJ(5), 0, 0, -t14, t13, -t53, -t25 * t53 - t1, -t26 * t53 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t36, -t72, 0, t36, -pkin(4) * t72 + t30, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t52, 0, -t39, -t40, -t39, t56, t40, t56 * pkin(8), 0, 0, -t21, t20, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t25, 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t36, 0, t6, 0, 0, 0, 0, 0, t51 * t53, -t48 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, -t51, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t53, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
