% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t32 = -t56 * pkin(2) - t53 * qJ(3) - pkin(1);
t24 = t56 * pkin(3) - t32;
t82 = 0.2e1 * t24;
t81 = -0.2e1 * t53;
t80 = 0.2e1 * t56;
t49 = sin(pkin(10));
t50 = cos(pkin(10));
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t42 = t56 * pkin(7);
t65 = -t56 * pkin(8) + t42;
t41 = t53 * pkin(7);
t66 = -t53 * pkin(8) + t41;
t19 = t52 * t65 - t55 * t66;
t28 = -t56 * t52 + t53 * t55;
t60 = -t28 * qJ(5) - t19;
t20 = t52 * t66 + t55 * t65;
t70 = t53 * t52 + t56 * t55;
t9 = -t70 * qJ(5) + t20;
t4 = t49 * t9 - t50 * t60;
t51 = sin(qJ(6));
t79 = t4 * t51;
t54 = cos(qJ(6));
t78 = t4 * t54;
t25 = t49 * t52 - t50 * t55;
t77 = t25 * t51;
t76 = t25 * t54;
t12 = t49 * t28 + t50 * t70;
t75 = t51 * t12;
t13 = t50 * t28 - t49 * t70;
t74 = t51 * t13;
t73 = t51 * t54;
t72 = t54 * t13;
t57 = -pkin(2) - pkin(3);
t30 = t52 * qJ(3) - t55 * t57;
t29 = -pkin(4) - t30;
t31 = t55 * qJ(3) + t52 * t57;
t16 = t50 * t29 - t49 * t31;
t14 = pkin(5) - t16;
t36 = -t50 * pkin(4) - pkin(5);
t71 = t14 - t36;
t17 = t49 * t29 + t50 * t31;
t46 = t53 ^ 2;
t69 = t56 ^ 2 + t46;
t68 = -0.2e1 * t73;
t67 = t51 * t72;
t64 = -t53 * pkin(2) + t56 * qJ(3);
t15 = -pkin(9) + t17;
t63 = -t12 * t15 + t13 * t14;
t27 = t49 * t55 + t50 * t52;
t62 = -t27 * t12 + t25 * t13;
t35 = t49 * pkin(4) + pkin(9);
t61 = -t12 * t35 + t13 * t36;
t18 = t70 * pkin(4) + t24;
t47 = t54 ^ 2;
t45 = t51 ^ 2;
t33 = 0.2e1 * t73;
t11 = t13 ^ 2;
t10 = t54 * t12;
t7 = (t45 - t47) * t13;
t6 = t49 * t60 + t50 * t9;
t3 = t12 * pkin(5) - t13 * pkin(9) + t18;
t2 = t51 * t3 + t54 * t6;
t1 = t54 * t3 - t51 * t6;
t5 = [1, 0, 0, t46, t53 * t80, 0, 0, 0, pkin(1) * t80, pkin(1) * t81, -0.2e1 * t32 * t56, 0.2e1 * t69 * pkin(7), t32 * t81, t69 * pkin(7) ^ 2 + t32 ^ 2, t28 ^ 2, -0.2e1 * t28 * t70, 0, 0, 0, t70 * t82, t28 * t82, -0.2e1 * t6 * t12 + 0.2e1 * t4 * t13, t18 ^ 2 + t4 ^ 2 + t6 ^ 2, t47 * t11, t11 * t68, 0.2e1 * t12 * t72, -0.2e1 * t12 * t74, t12 ^ 2, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t74, -0.2e1 * t2 * t12 + 0.2e1 * t4 * t72; 0, 0, 0, 0, 0, t53, t56, 0, -t41, -t42, -t41, t64, t42, t64 * pkin(7), 0, 0, -t28, t70, 0, t19, t20, -t17 * t12 - t16 * t13, -t4 * t16 + t6 * t17, -t67, t7, -t75, -t10, 0, t63 * t51 + t78, t63 * t54 - t79; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t30, 0.2e1 * t31, 0, t16 ^ 2 + t17 ^ 2, t45, t33, 0, 0, 0, 0.2e1 * t14 * t54, -0.2e1 * t14 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t41, 0, 0, 0, 0, 0, 0, 0, t62, t4 * t25 + t6 * t27, 0, 0, 0, 0, 0, t62 * t51, t62 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t55, t52, 0, -t16 * t25 + t17 * t27, 0, 0, 0, 0, 0, t76, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t25 ^ 2 + t27 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t70, 0, -t19, -t20 (-t12 * t49 - t13 * t50) * pkin(4) (-t4 * t50 + t49 * t6) * pkin(4), t67, -t7, t75, t10, 0, t61 * t51 - t78, t61 * t54 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t30, -t31, 0 (t16 * t50 + t17 * t49) * pkin(4), -t45, t68, 0, 0, 0, -t71 * t54, t71 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t52, 0 (-t25 * t50 + t27 * t49) * pkin(4), 0, 0, 0, 0, 0, -t76, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t49 ^ 2 + t50 ^ 2) * pkin(4) ^ 2, t45, t33, 0, 0, 0, -0.2e1 * t36 * t54, 0.2e1 * t36 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, t10, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t74, t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t54, 0, -t51 * t15, -t54 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t27, -t54 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t54, 0, -t51 * t35, -t54 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
