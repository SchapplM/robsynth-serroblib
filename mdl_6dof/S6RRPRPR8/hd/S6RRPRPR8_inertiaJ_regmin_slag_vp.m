% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t58 = sin(pkin(10));
t59 = cos(pkin(10));
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t90 = -t61 * t58 + t64 * t59;
t39 = t64 * t58 + t61 * t59;
t50 = -t59 * pkin(3) - pkin(2);
t68 = t39 * qJ(5) - t50;
t84 = -pkin(4) - pkin(5);
t13 = -t84 * t90 + t68;
t89 = 0.2e1 * t13;
t88 = -0.2e1 * t39;
t87 = 0.2e1 * t50;
t65 = cos(qJ(2));
t86 = -0.2e1 * t65;
t85 = 0.2e1 * t65;
t83 = pkin(7) * t58;
t62 = sin(qJ(2));
t52 = t62 * pkin(7);
t82 = t65 * pkin(7);
t75 = pkin(8) + qJ(3);
t44 = t75 * t58;
t45 = t75 * t59;
t27 = t64 * t44 + t61 * t45;
t81 = t27 * t65;
t28 = -t61 * t44 + t64 * t45;
t80 = t28 * t65;
t79 = t58 * t62;
t78 = t59 * t62;
t43 = -t65 * pkin(2) - t62 * qJ(3) - pkin(1);
t36 = t59 * t43;
t23 = -pkin(8) * t78 + t36 + (-pkin(3) - t83) * t65;
t30 = t58 * t43 + t59 * t82;
t26 = -pkin(8) * t79 + t30;
t74 = -t64 * t23 + t61 * t26;
t11 = t61 * t23 + t64 * t26;
t40 = pkin(3) * t79 + t52;
t73 = t58 ^ 2 + t59 ^ 2;
t72 = t65 * qJ(5);
t53 = t65 * pkin(4);
t8 = t53 + t74;
t7 = -t72 + t11;
t33 = t90 * t62;
t71 = t33 * qJ(5) - t40;
t3 = t65 * pkin(5) - t33 * pkin(9) + t8;
t32 = t39 * t62;
t4 = t32 * pkin(9) + t7;
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t1 = t63 * t3 - t60 * t4;
t2 = t60 * t3 + t63 * t4;
t70 = -pkin(2) * t62 + qJ(3) * t65;
t29 = -t58 * t82 + t36;
t69 = -t29 * t58 + t30 * t59;
t67 = -t39 * pkin(9) + t27;
t57 = t65 ^ 2;
t56 = t62 ^ 2;
t42 = t63 * qJ(5) + t60 * t84;
t41 = t60 * qJ(5) - t63 * t84;
t22 = t63 * t39 - t60 * t90;
t21 = t60 * t39 + t63 * t90;
t20 = -pkin(4) * t90 - t68;
t17 = -pkin(9) * t90 + t28;
t16 = t60 * t32 + t63 * t33;
t15 = -t63 * t32 + t60 * t33;
t12 = t32 * pkin(4) - t71;
t9 = t84 * t32 + t71;
t6 = t63 * t17 + t60 * t67;
t5 = t60 * t17 - t63 * t67;
t10 = [1, 0, 0, t56, t62 * t85, 0, 0, 0, pkin(1) * t85, -0.2e1 * pkin(1) * t62, -0.2e1 * t29 * t65 + 0.2e1 * t56 * t83, 0.2e1 * t56 * pkin(7) * t59 + 0.2e1 * t30 * t65, 0.2e1 * (-t29 * t59 - t30 * t58) * t62, t56 * pkin(7) ^ 2 + t29 ^ 2 + t30 ^ 2, t33 ^ 2, -0.2e1 * t33 * t32, t33 * t86, t32 * t85, t57, 0.2e1 * t40 * t32 + 0.2e1 * t65 * t74, 0.2e1 * t11 * t65 + 0.2e1 * t40 * t33, 0.2e1 * t12 * t32 + 0.2e1 * t8 * t65, -0.2e1 * t7 * t32 + 0.2e1 * t8 * t33, -0.2e1 * t12 * t33 - 0.2e1 * t7 * t65, t12 ^ 2 + t7 ^ 2 + t8 ^ 2, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t85, t15 * t86, t57, 0.2e1 * t1 * t65 + 0.2e1 * t9 * t15, 0.2e1 * t9 * t16 - 0.2e1 * t2 * t65; 0, 0, 0, 0, 0, t62, t65, 0, -t52, -t82, -pkin(7) * t78 + t70 * t58, pkin(7) * t79 + t70 * t59, t69, -pkin(2) * t52 + t69 * qJ(3), t33 * t39, -t39 * t32 + t33 * t90, -t39 * t65, -t90 * t65, 0, t50 * t32 - t40 * t90 + t81, t50 * t33 + t40 * t39 + t80, -t12 * t90 + t20 * t32 + t81, t27 * t33 - t28 * t32 + t8 * t39 + t7 * t90, -t12 * t39 - t20 * t33 - t80, t12 * t20 + t8 * t27 + t7 * t28, t16 * t22, -t22 * t15 - t16 * t21, t22 * t65, -t21 * t65, 0, t13 * t15 + t9 * t21 - t5 * t65, t13 * t16 + t9 * t22 - t6 * t65; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t59, -0.2e1 * pkin(2) * t58, 0.2e1 * t73 * qJ(3), t73 * qJ(3) ^ 2 + pkin(2) ^ 2, t39 ^ 2, -t90 * t88, 0, 0, 0, -t90 * t87, t39 * t87, -0.2e1 * t20 * t90, 0.2e1 * t27 * t39 + 0.2e1 * t28 * t90, t20 * t88, t20 ^ 2 + t27 ^ 2 + t28 ^ 2, t22 ^ 2, -0.2e1 * t22 * t21, 0, 0, 0, t21 * t89, t22 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t78, 0, t52, 0, 0, 0, 0, 0, t32, t33, t32, 0, -t33, t12, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, -pkin(2), 0, 0, 0, 0, 0, -t90, t39, -t90, 0, -t39, t20, 0, 0, 0, 0, 0, -t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, -t65, -t74, -t11, -0.2e1 * t53 - t74, -pkin(4) * t33 - t32 * qJ(5), -0.2e1 * t72 + t11, -t8 * pkin(4) + t7 * qJ(5), 0, 0, -t16, t15, -t65, -t41 * t65 - t1, -t42 * t65 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t90, 0, -t27, -t28, -t27, -pkin(4) * t39 + qJ(5) * t90, t28, -t27 * pkin(4) + t28 * qJ(5), 0, 0, -t22, t21, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t41, 0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t33, 0, t8, 0, 0, 0, 0, 0, t63 * t65, -t60 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, -t63, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t65, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t41, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
