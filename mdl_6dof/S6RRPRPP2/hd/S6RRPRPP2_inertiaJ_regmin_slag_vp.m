% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t47 = sin(qJ(4));
t41 = t47 * qJ(5);
t49 = cos(qJ(4));
t69 = t49 * pkin(4) + t41;
t50 = pkin(4) + pkin(5);
t79 = t49 * t50 + t41;
t68 = t49 * qJ(5);
t78 = t50 * t47 - t68;
t77 = 0.2e1 * t47;
t76 = -0.2e1 * t49;
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t48 = sin(qJ(2));
t75 = cos(qJ(2));
t30 = t45 * t48 - t46 * t75;
t26 = t30 * pkin(4);
t31 = t45 * t75 + t46 * t48;
t40 = -t75 * pkin(2) - pkin(1);
t11 = t30 * pkin(3) - t31 * pkin(8) + t40;
t64 = t75 * pkin(7);
t33 = t75 * qJ(3) + t64;
t63 = (-qJ(3) - pkin(7)) * t48;
t17 = t46 * t33 + t45 * t63;
t7 = t47 * t11 + t49 * t17;
t38 = t45 * pkin(2) + pkin(8);
t74 = t30 * t38;
t20 = t47 * t30;
t73 = t47 * t31;
t72 = t47 * t38;
t71 = t47 * t49;
t21 = t49 * t30;
t22 = t49 * t31;
t6 = t49 * t11 - t47 * t17;
t43 = t47 ^ 2;
t44 = t49 ^ 2;
t35 = t43 + t44;
t67 = t49 * qJ(6);
t66 = 0.2e1 * t75;
t24 = t30 * qJ(5);
t3 = t24 + t7;
t65 = 0.2e1 * t24 + t7;
t4 = -t26 - t6;
t39 = -t46 * pkin(2) - pkin(3);
t15 = t45 * t33 - t46 * t63;
t55 = t31 * t67 - t4;
t1 = -t30 * pkin(5) - t55;
t19 = qJ(6) * t73;
t2 = t19 + t3;
t62 = -t1 * t49 + t2 * t47;
t61 = t3 * t49 + t4 * t47;
t60 = t3 * t47 - t4 * t49;
t59 = pkin(4) * t47 - t68;
t25 = t39 - t69;
t58 = t25 * t31 - t74;
t27 = (-qJ(6) + t38) * t47;
t34 = t49 * t38;
t28 = t34 - t67;
t57 = -t27 * t49 + t28 * t47;
t56 = t31 * t39 - t74;
t53 = qJ(5) ^ 2;
t52 = 0.2e1 * qJ(5);
t29 = t31 ^ 2;
t18 = t49 * pkin(5) - t25;
t12 = t35 * t31;
t8 = t59 * t31 + t15;
t5 = -t78 * t31 - t15;
t9 = [1, 0, 0, t48 ^ 2, t48 * t66, 0, 0, 0, pkin(1) * t66, -0.2e1 * pkin(1) * t48, 0.2e1 * t15 * t31 - 0.2e1 * t17 * t30, t15 ^ 2 + t17 ^ 2 + t40 ^ 2, t44 * t29, -0.2e1 * t29 * t71, 0.2e1 * t30 * t22, -0.2e1 * t30 * t73, t30 ^ 2, 0.2e1 * t15 * t73 + 0.2e1 * t6 * t30, 0.2e1 * t15 * t22 - 0.2e1 * t7 * t30, -0.2e1 * t4 * t30 + 0.2e1 * t8 * t73, -0.2e1 * t60 * t31, -0.2e1 * t8 * t22 + 0.2e1 * t3 * t30, t3 ^ 2 + t4 ^ 2 + t8 ^ 2, -0.2e1 * t1 * t30 - 0.2e1 * t5 * t73, 0.2e1 * t2 * t30 + 0.2e1 * t5 * t22, 0.2e1 * t62 * t31, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t48, t75, 0, -t48 * pkin(7), -t64 (-t30 * t45 - t31 * t46) * pkin(2) (-t15 * t46 + t17 * t45) * pkin(2), t47 * t22 (-t43 + t44) * t31, t20, t21, 0, -t15 * t49 + t56 * t47, t15 * t47 + t56 * t49, t58 * t47 - t8 * t49, t61, -t8 * t47 - t58 * t49, t8 * t25 + t61 * t38, -t18 * t73 - t27 * t30 + t5 * t49, t18 * t22 + t28 * t30 + t5 * t47, -t1 * t47 - t2 * t49 + t57 * t31, t1 * t27 + t5 * t18 + t2 * t28; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t45 ^ 2 + t46 ^ 2) * pkin(2) ^ 2, t43, 0.2e1 * t71, 0, 0, 0, t39 * t76, t39 * t77, t25 * t76, 0.2e1 * t35 * t38, -0.2e1 * t25 * t47, t35 * t38 ^ 2 + t25 ^ 2, 0.2e1 * t18 * t49, t18 * t77, -0.2e1 * t27 * t47 - 0.2e1 * t28 * t49, t18 ^ 2 + t27 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, t21, -t20, t21, -t12, t20, t60, t21, t20, t12, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t73, t30, t6, -t7, -t4 + t26, -t69 * t31, t65, -t4 * pkin(4) + t3 * qJ(5) (pkin(5) + t50) * t30 + t55, t19 + t65, t79 * t31, t2 * qJ(5) - t1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t49, 0, -t72, -t34, -t72, -t59, t34, -t59 * t38, -t27, t28, t78, t28 * qJ(5) - t27 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t47, t49, 0, t47, t69, t49, t47, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t52, pkin(4) ^ 2 + t53, 0.2e1 * t50, t52, 0, t50 ^ 2 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t22, 0, t4, -t30, 0, -t22, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t72, 0, 0, -t47, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t22, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t47, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
