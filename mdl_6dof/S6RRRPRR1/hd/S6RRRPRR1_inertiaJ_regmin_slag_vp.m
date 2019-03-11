% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t62 = sin(qJ(3));
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t80 = cos(qJ(3));
t42 = t62 * t63 - t80 * t66;
t53 = -t66 * pkin(2) - pkin(1);
t34 = t42 * pkin(3) + t53;
t43 = t62 * t66 + t80 * t63;
t58 = sin(pkin(11));
t59 = cos(pkin(11));
t69 = t59 * t42 + t58 * t43;
t18 = t69 * pkin(4) + t34;
t88 = 0.2e1 * t18;
t87 = 0.2e1 * t53;
t86 = 0.2e1 * t66;
t85 = pkin(7) + pkin(8);
t60 = sin(qJ(6));
t84 = pkin(5) * t60;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t47 = t85 * t63;
t48 = t85 * t66;
t30 = -t80 * t47 - t62 * t48;
t26 = -t43 * qJ(4) + t30;
t31 = t62 * t47 - t80 * t48;
t27 = -t42 * qJ(4) - t31;
t10 = t59 * t26 - t58 * t27;
t29 = -t58 * t42 + t59 * t43;
t68 = -t29 * pkin(9) + t10;
t11 = t58 * t26 + t59 * t27;
t9 = -t69 * pkin(9) + t11;
t4 = t61 * t9 - t65 * t68;
t64 = cos(qJ(6));
t83 = t4 * t64;
t82 = t58 * pkin(3);
t81 = t62 * pkin(2);
t54 = t80 * pkin(2);
t52 = t54 + pkin(3);
t38 = t59 * t52 - t58 * t81;
t35 = pkin(4) + t38;
t32 = t65 * t35;
t40 = t58 * t52 + t59 * t81;
t24 = -t61 * t40 + t32;
t20 = -pkin(5) - t24;
t79 = t20 * t64;
t50 = t59 * pkin(3) + pkin(4);
t46 = t65 * t50;
t39 = -t61 * t82 + t46;
t36 = -pkin(5) - t39;
t78 = t36 * t64;
t16 = t61 * t29 + t65 * t69;
t13 = t60 * t16;
t17 = t65 * t29 - t61 * t69;
t77 = t60 * t17;
t76 = t60 * t64;
t75 = t64 * t17;
t74 = -0.2e1 * t17 * t16;
t73 = -t40 - t82;
t72 = -pkin(5) * t17 - pkin(10) * t16;
t25 = -t61 * t35 - t65 * t40;
t21 = pkin(10) - t25;
t71 = -t16 * t21 + t17 * t20;
t41 = -t61 * t50 - t65 * t82;
t37 = pkin(10) - t41;
t70 = -t16 * t37 + t17 * t36;
t57 = t64 ^ 2;
t56 = t60 ^ 2;
t55 = pkin(5) * t64;
t49 = 0.2e1 * t76;
t33 = t36 * t60;
t19 = t20 * t60;
t15 = t17 ^ 2;
t14 = t64 * t16;
t12 = t60 * t75;
t7 = (-t56 + t57) * t17;
t6 = t16 * pkin(5) - t17 * pkin(10) + t18;
t5 = t61 * t68 + t65 * t9;
t3 = t4 * t60;
t2 = t64 * t5 + t60 * t6;
t1 = -t60 * t5 + t64 * t6;
t8 = [1, 0, 0, t63 ^ 2, t63 * t86, 0, 0, 0, pkin(1) * t86, -0.2e1 * pkin(1) * t63, t43 ^ 2, -0.2e1 * t43 * t42, 0, 0, 0, t42 * t87, t43 * t87, -0.2e1 * t10 * t29 - 0.2e1 * t11 * t69, t10 ^ 2 + t11 ^ 2 + t34 ^ 2, t15, t74, 0, 0, 0, t16 * t88, t17 * t88, t57 * t15, -0.2e1 * t15 * t76, 0.2e1 * t16 * t75, t60 * t74, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t77, -0.2e1 * t2 * t16 + 0.2e1 * t4 * t75; 0, 0, 0, 0, 0, t63, t66, 0, -t63 * pkin(7), -t66 * pkin(7), 0, 0, t43, -t42, 0, t30, t31, -t38 * t29 - t40 * t69, t10 * t38 + t11 * t40, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t71 * t60 - t83, t71 * t64 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t81, 0, t38 ^ 2 + t40 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t24, 0.2e1 * t25, t56, t49, 0, 0, 0, -0.2e1 * t79, 0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42, 0, t30, t31 (-t59 * t29 - t58 * t69) * pkin(3) (t10 * t59 + t11 * t58) * pkin(3), 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t70 * t60 - t83, t70 * t64 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t54, -t81, 0 (t38 * t59 + t40 * t58) * pkin(3), 0, 0, 0, 0, 1, t73 * t61 + t32 + t46, t73 * t65 + (-t35 - t50) * t61, t56, t49, 0, 0, 0 (-t20 - t36) * t64, t33 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t58 ^ 2 + t59 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t39, 0.2e1 * t41, t56, t49, 0, 0, 0, -0.2e1 * t78, 0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t72 * t60 - t83, t72 * t64 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, t25, t56, t49, 0, 0, 0, t55 - t79, t19 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t39, t41, t56, t49, 0, 0, 0, t55 - t78, t33 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, t49, 0, 0, 0, 0.2e1 * t55, -0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t77, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * t21, -t64 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * t37, -t64 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * pkin(10), -t64 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
