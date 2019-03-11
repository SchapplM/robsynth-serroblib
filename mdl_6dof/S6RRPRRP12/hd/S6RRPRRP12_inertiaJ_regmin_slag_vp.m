% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t49 = sin(qJ(5));
t50 = sin(qJ(4));
t52 = cos(qJ(4));
t85 = cos(qJ(5));
t68 = t85 * t52;
t23 = t49 * t50 - t68;
t94 = t23 ^ 2;
t80 = t49 * t52;
t58 = -t85 * t50 - t80;
t93 = -0.2e1 * t58;
t51 = sin(qJ(2));
t92 = -0.2e1 * t51;
t91 = 0.2e1 * t51;
t53 = cos(qJ(2));
t90 = 0.2e1 * t53;
t89 = 2 * qJ(3);
t54 = -pkin(2) - pkin(8);
t65 = -t51 * qJ(3) - pkin(1);
t21 = t54 * t53 + t65;
t67 = pkin(9) * t53 - t21;
t42 = t51 * pkin(7);
t29 = t51 * pkin(3) + t42;
t79 = t50 * t29;
t10 = -t67 * t52 + t79;
t26 = t52 * t29;
t88 = t51 * pkin(4);
t8 = t67 * t50 + t26 + t88;
t4 = t85 * t10 + t49 * t8;
t87 = t51 * pkin(5);
t86 = -pkin(9) + t54;
t27 = t86 * t50;
t14 = t49 * t27 - t86 * t68;
t84 = t14 * t51;
t15 = t85 * t27 + t86 * t80;
t83 = t15 * t51;
t17 = t58 * t53;
t82 = t23 * t17;
t20 = t23 * t51;
t81 = t58 * t51;
t78 = t50 * t51;
t77 = t50 * t53;
t76 = t51 * t53;
t75 = t52 * t50;
t74 = t52 * t53;
t43 = t53 * pkin(7);
t30 = t53 * pkin(3) + t43;
t46 = t51 ^ 2;
t48 = t53 ^ 2;
t73 = t46 + t48;
t72 = t53 * qJ(3);
t34 = t50 * pkin(4) + qJ(3);
t71 = -0.2e1 * t76;
t39 = t51 * qJ(6);
t1 = t39 + t4;
t19 = pkin(4) * t74 + t30;
t70 = t85 * pkin(4);
t69 = t58 ^ 2 + t94;
t66 = t49 * t10 - t85 * t8;
t2 = t66 - t87;
t64 = -t1 * t58 + t2 * t23;
t63 = -t51 * pkin(2) + t72;
t62 = t23 * pkin(5) + qJ(6) * t58;
t61 = t14 * t23 - t15 * t58;
t40 = t49 * pkin(4);
t33 = t40 + qJ(6);
t37 = t70 + pkin(5);
t60 = t23 * t37 + t33 * t58;
t59 = t51 * t54 + t72;
t56 = 0.2e1 * pkin(5);
t55 = 0.2e1 * qJ(6);
t47 = t52 ^ 2;
t45 = t50 ^ 2;
t36 = t52 * t51;
t28 = -t53 * pkin(2) + t65;
t16 = t49 * t77 - t53 * t68;
t13 = t52 * t21 + t79;
t12 = -t50 * t21 + t26;
t11 = -pkin(5) * t58 + t23 * qJ(6) + t34;
t5 = -t16 * pkin(5) - t17 * qJ(6) + t19;
t3 = [1, 0, 0, t46, 0.2e1 * t76, 0, 0, 0, pkin(1) * t90, pkin(1) * t92, 0.2e1 * t73 * pkin(7), t28 * t90, t28 * t92, t73 * pkin(7) ^ 2 + t28 ^ 2, t45 * t48, 0.2e1 * t48 * t75, t50 * t71, t52 * t71, t46, 0.2e1 * t12 * t51 + 0.2e1 * t30 * t74, -0.2e1 * t13 * t51 - 0.2e1 * t30 * t77, t17 ^ 2, 0.2e1 * t17 * t16, t17 * t91, t16 * t91, t46, -0.2e1 * t19 * t16 - 0.2e1 * t51 * t66, 0.2e1 * t19 * t17 - 0.2e1 * t4 * t51, -0.2e1 * t5 * t16 - 0.2e1 * t2 * t51, 0.2e1 * t1 * t16 + 0.2e1 * t2 * t17, 0.2e1 * t1 * t51 - 0.2e1 * t5 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t51, t53, 0, -t42, -t43, t63, t42, t43, t63 * pkin(7), -t50 * t74 (t45 - t47) * t53, t36, -t78, 0, t30 * t50 + t59 * t52, t30 * t52 - t59 * t50, -t82, -t23 * t16 + t17 * t58, -t20, t81, 0, -t34 * t16 - t19 * t58 - t84, t34 * t17 - t19 * t23 - t83, -t11 * t16 - t5 * t58 - t84, t14 * t17 + t15 * t16 - t64, -t11 * t17 + t5 * t23 + t83, t1 * t15 + t5 * t11 + t2 * t14; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t89, pkin(2) ^ 2 + (qJ(3) ^ 2) t47, -0.2e1 * t75, 0, 0, 0, t50 * t89, t52 * t89, t94, t23 * t93, 0, 0, 0, t34 * t93, -0.2e1 * t34 * t23, t11 * t93, -0.2e1 * t61, 0.2e1 * t11 * t23, t11 ^ 2 + t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, t42, 0, 0, 0, 0, 0, t36, -t78, 0, 0, 0, 0, 0, -t20, t81, -t20, -t16 * t58 + t82, -t81, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t74, t51, t12, -t13, 0, 0, t17, t16, t51, t51 * t70 - t66, -t49 * t88 - t4 (pkin(5) + t37) * t51 - t66, t33 * t16 - t37 * t17, t33 * t51 + t1, t1 * t33 - t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, 0, t52 * t54, -t50 * t54, 0, 0, -t23, t58, 0, -t14, -t15, -t14, t60, t15, -t14 * t37 + t15 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, 0, 0, 0, 0, 0, -t23, t58, -t23, 0, -t58, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t40, 0.2e1 * t37, 0, 0.2e1 * t33, t33 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t16, t51, -t66, -t4, -t66 + 0.2e1 * t87, -pkin(5) * t17 + t16 * qJ(6), 0.2e1 * t39 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t58, 0, -t14, -t15, -t14, t62, t15, -t14 * pkin(5) + t15 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t58, -t23, 0, -t58, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t40, t56 + t70, 0, t55 + t40, t37 * pkin(5) + t33 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, 0, t55, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
