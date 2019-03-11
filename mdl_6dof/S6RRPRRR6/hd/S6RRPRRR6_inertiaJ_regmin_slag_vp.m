% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t60 = cos(qJ(6));
t89 = t60 ^ 2;
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t36 = -t63 * pkin(2) - t59 * qJ(3) - pkin(1);
t28 = t63 * pkin(3) - t36;
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t74 = t59 * t58 + t63 * t62;
t22 = t74 * pkin(4) + t28;
t88 = 0.2e1 * t22;
t87 = 0.2e1 * t28;
t86 = -0.2e1 * t59;
t85 = 0.2e1 * t60;
t84 = 0.2e1 * t63;
t56 = sin(qJ(6));
t83 = pkin(5) * t56;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t49 = t59 * pkin(7);
t37 = -t59 * pkin(8) + t49;
t50 = t63 * pkin(7);
t38 = -t63 * pkin(8) + t50;
t23 = -t62 * t37 + t58 * t38;
t30 = -t63 * t58 + t59 * t62;
t66 = -t30 * pkin(9) - t23;
t24 = t58 * t37 + t62 * t38;
t9 = -t74 * pkin(9) + t24;
t5 = t57 * t9 - t61 * t66;
t82 = t5 * t60;
t81 = t61 * pkin(4);
t45 = -pkin(5) - t81;
t80 = pkin(5) - t45;
t29 = t57 * t58 - t61 * t62;
t79 = t29 * t60;
t78 = t45 * t56;
t17 = t61 * t30 - t57 * t74;
t77 = t56 * t17;
t76 = t56 * t60;
t75 = t60 * t17;
t54 = t59 ^ 2;
t73 = t63 ^ 2 + t54;
t16 = t57 * t30 + t61 * t74;
t72 = -0.2e1 * t17 * t16;
t40 = -0.2e1 * t76;
t64 = -pkin(2) - pkin(3);
t34 = t58 * qJ(3) - t62 * t64;
t33 = -pkin(4) - t34;
t35 = t62 * qJ(3) + t58 * t64;
t20 = -t61 * t33 + t57 * t35;
t71 = -pkin(5) * t17 - pkin(10) * t16;
t70 = -t59 * pkin(2) + t63 * qJ(3);
t18 = pkin(5) + t20;
t21 = t57 * t33 + t61 * t35;
t19 = -pkin(10) + t21;
t69 = -t16 * t19 + t17 * t18;
t31 = t57 * t62 + t61 * t58;
t68 = -t16 * t31 + t17 * t29;
t48 = t57 * pkin(4);
t44 = t48 + pkin(10);
t67 = -t16 * t44 + t17 * t45;
t53 = t56 ^ 2;
t39 = 0.2e1 * t76;
t26 = t29 * t56;
t15 = t17 ^ 2;
t14 = t18 * t56;
t13 = t60 * t16;
t11 = t56 * t16;
t10 = t56 * t75;
t7 = (t53 - t89) * t17;
t6 = t57 * t66 + t61 * t9;
t4 = t5 * t56;
t3 = t16 * pkin(5) - t17 * pkin(10) + t22;
t2 = t56 * t3 + t60 * t6;
t1 = t60 * t3 - t56 * t6;
t8 = [1, 0, 0, t54, t59 * t84, 0, 0, 0, pkin(1) * t84, pkin(1) * t86, -0.2e1 * t36 * t63, 0.2e1 * t73 * pkin(7), t36 * t86, t73 * pkin(7) ^ 2 + t36 ^ 2, t30 ^ 2, -0.2e1 * t30 * t74, 0, 0, 0, t74 * t87, t30 * t87, t15, t72, 0, 0, 0, t16 * t88, t17 * t88, t89 * t15, t15 * t40, 0.2e1 * t16 * t75, t56 * t72, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t5 * t77, -0.2e1 * t2 * t16 + 0.2e1 * t5 * t75; 0, 0, 0, 0, 0, t59, t63, 0, -t49, -t50, -t49, t70, t50, t70 * pkin(7), 0, 0, -t30, t74, 0, t23, t24, 0, 0, -t17, t16, 0, t5, t6, -t10, t7, -t11, -t13, 0, t69 * t56 + t82, t69 * t60 - t4; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t34, 0.2e1 * t35, 0, 0, 0, 0, 1, 0.2e1 * t20, 0.2e1 * t21, t53, t39, 0, 0, 0, t18 * t85, -0.2e1 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t56, t68 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t62, t58, 0, 0, 0, 0, 0, t29, t31, 0, 0, 0, 0, 0, t79, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t74, 0, -t23, -t24, 0, 0, t17, -t16, 0, -t5, -t6, t10, -t7, t11, t13, 0, t67 * t56 - t82, t67 * t60 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t34, -t35, 0, 0, 0, 0, -1, -t20 - t81, -t21 + t48, -t53, t40, 0, 0, 0 (-t18 + t45) * t60, t14 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t58, 0, 0, 0, 0, 0, -t29, -t31, 0, 0, 0, 0, 0, -t79, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t81, -0.2e1 * t48, t53, t39, 0, 0, 0, -0.2e1 * t45 * t60, 0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t5, -t6, t10, -t7, t11, t13, 0, t71 * t56 - t82, t71 * t60 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t20, -t21, -t53, t40, 0, 0, 0 (-pkin(5) - t18) * t60, t14 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t31, 0, 0, 0, 0, 0, -t79, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t81, -t48, t53, t39, 0, 0, 0, t80 * t60, -t80 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t53, t39, 0, 0, 0, pkin(5) * t85, -0.2e1 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t77, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t60, 0, -t56 * t19, -t60 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56 * t31, -t60 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t60, 0, -t56 * t44, -t60 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t60, 0, -t56 * pkin(10), -t60 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
