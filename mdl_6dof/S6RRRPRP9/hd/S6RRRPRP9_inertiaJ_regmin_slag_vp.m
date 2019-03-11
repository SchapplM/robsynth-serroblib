% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t93 = -pkin(3) * t58 - qJ(4) * t55;
t56 = sin(qJ(2));
t43 = t58 * t56;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t72 = t57 * t55;
t21 = t43 * t54 - t56 * t72;
t92 = -0.2e1 * t21;
t26 = t54 * t55 + t57 * t58;
t91 = 0.2e1 * t26;
t27 = -t54 * t58 + t72;
t90 = -0.2e1 * t27;
t33 = -pkin(2) + t93;
t89 = -0.2e1 * t33;
t88 = -0.2e1 * t56;
t59 = cos(qJ(2));
t87 = 0.2e1 * t59;
t60 = -pkin(3) - pkin(4);
t86 = pkin(2) * t55;
t85 = pkin(2) * t58;
t84 = pkin(3) * t55;
t83 = pkin(7) * t55;
t82 = pkin(7) * t58;
t81 = t55 * pkin(8);
t80 = t59 * pkin(5);
t34 = -pkin(2) * t59 - pkin(8) * t56 - pkin(1);
t71 = t58 * t59;
t19 = pkin(7) * t71 + t34 * t55;
t67 = t59 * qJ(4);
t16 = -t67 + t19;
t75 = t55 * t56;
t12 = pkin(9) * t75 + t16;
t49 = t59 * pkin(3);
t73 = t55 * t59;
t70 = pkin(7) * t73 - t34 * t58;
t17 = t49 + t70;
t9 = pkin(4) * t59 - pkin(9) * t43 + t17;
t79 = t12 * t54 - t57 * t9;
t4 = t12 * t57 + t54 * t9;
t46 = t58 * pkin(8);
t35 = -pkin(9) * t58 + t46;
t64 = (pkin(8) - pkin(9)) * t55;
t14 = t35 * t54 - t57 * t64;
t78 = t14 * t59;
t15 = t57 * t35 + t54 * t64;
t77 = t15 * t59;
t76 = t54 * t59;
t74 = t55 * t58;
t50 = t55 ^ 2;
t52 = t58 ^ 2;
t69 = t50 + t52;
t68 = t58 * qJ(4);
t66 = t59 * qJ(6);
t65 = t56 * t87;
t23 = pkin(4) * t58 - t33;
t31 = qJ(4) * t54 - t57 * t60;
t63 = t68 - t84;
t62 = t16 * t58 + t17 * t55;
t32 = qJ(4) * t57 + t54 * t60;
t37 = t56 * t68;
t13 = t37 + (t55 * t60 - pkin(7)) * t56;
t53 = t59 ^ 2;
t51 = t56 ^ 2;
t42 = t57 * t59;
t39 = pkin(8) * t73;
t30 = pkin(5) + t31;
t29 = -qJ(6) + t32;
t22 = t26 * t56;
t20 = -t37 + (pkin(7) + t84) * t56;
t6 = pkin(5) * t26 - qJ(6) * t27 + t23;
t5 = t21 * pkin(5) - t22 * qJ(6) + t13;
t2 = t79 - t80;
t1 = t66 + t4;
t3 = [1, 0, 0, t51, t65, 0, 0, 0, pkin(1) * t87, pkin(1) * t88, t52 * t51, -0.2e1 * t51 * t74, t71 * t88, t55 * t65, t53, 0.2e1 * t51 * t83 + 0.2e1 * t59 * t70, 0.2e1 * t19 * t59 + 0.2e1 * t51 * t82, 0.2e1 * t17 * t59 + 0.2e1 * t20 * t75, 0.2e1 * (-t16 * t55 + t17 * t58) * t56, -0.2e1 * t16 * t59 - 0.2e1 * t20 * t43, t16 ^ 2 + t17 ^ 2 + t20 ^ 2, t22 ^ 2, t22 * t92, t22 * t87, t59 * t92, t53, 0.2e1 * t13 * t21 - 0.2e1 * t59 * t79, 0.2e1 * t13 * t22 - 0.2e1 * t4 * t59, -0.2e1 * t2 * t59 + 0.2e1 * t21 * t5, -0.2e1 * t1 * t21 + 0.2e1 * t2 * t22, 0.2e1 * t1 * t59 - 0.2e1 * t22 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t56, t59, 0, -t56 * pkin(7), -t59 * pkin(7), t55 * t43 (-t50 + t52) * t56, -t73, -t71, 0, t39 + (-t82 - t86) * t56, pkin(8) * t71 + (t83 - t85) * t56, -t20 * t58 + t33 * t75 + t39, t62, -t20 * t55 + (-pkin(8) * t59 - t33 * t56) * t58, pkin(8) * t62 + t20 * t33, t22 * t27, -t21 * t27 - t22 * t26, t27 * t59, -t26 * t59, 0, t13 * t26 + t21 * t23 - t78, t13 * t27 + t22 * t23 - t77, t21 * t6 + t26 * t5 - t78, -t1 * t26 + t14 * t22 - t15 * t21 + t2 * t27, -t22 * t6 - t27 * t5 + t77, t1 * t15 + t14 * t2 + t5 * t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, 0.2e1 * t74, 0, 0, 0, 0.2e1 * t85, -0.2e1 * t86, t58 * t89, 0.2e1 * t69 * pkin(8), t55 * t89, pkin(8) ^ 2 * t69 + t33 ^ 2, t27 ^ 2, t26 * t90, 0, 0, 0, t23 * t91, 0.2e1 * t23 * t27, t6 * t91, 0.2e1 * t14 * t27 - 0.2e1 * t15 * t26, t6 * t90, t14 ^ 2 + t15 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t75, -t59, -t70, -t19, -0.2e1 * t49 - t70, t93 * t56, -0.2e1 * t67 + t19, -pkin(3) * t17 + qJ(4) * t16, 0, 0, -t22, t21, -t59, -t31 * t59 + t79, -t32 * t59 + t4 (-pkin(5) - t30) * t59 + t79, -t21 * t29 + t22 * t30 (-qJ(6) + t29) * t59 - t4, t1 * t29 + t2 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t81, -t46, -t81, t63, t46, t63 * pkin(8), 0, 0, -t27, t26, 0, t14, t15, t14, -t26 * t29 + t27 * t30, -t15, t14 * t30 + t15 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t31, 0.2e1 * t32, 0.2e1 * t30, 0, -0.2e1 * t29, t29 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t43, 0, t17, 0, 0, 0, 0, 0, t42, -t76, t42, -t21 * t54 - t22 * t57, t76, t1 * t54 - t2 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t81, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t54 - t27 * t57, 0, -t14 * t57 + t15 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t57, t54, -t57, 0, -t54, t29 * t54 - t30 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 ^ 2 + t57 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t59, -t79, -t4, -t79 + 0.2e1 * t80, -pkin(5) * t22 - qJ(6) * t21, 0.2e1 * t66 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, -t14, -t15, -t14, -pkin(5) * t27 - qJ(6) * t26, t15, -pkin(5) * t14 + qJ(6) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t31, -t32, -0.2e1 * pkin(5) - t31, 0, -0.2e1 * qJ(6) + t32, -pkin(5) * t30 + qJ(6) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54, t57, 0, t54, pkin(5) * t57 + qJ(6) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
