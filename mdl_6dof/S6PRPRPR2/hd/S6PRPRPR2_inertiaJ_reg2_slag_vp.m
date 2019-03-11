% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t55 = sin(pkin(12));
t58 = cos(pkin(12));
t61 = sin(qJ(6));
t64 = cos(qJ(6));
t103 = -t61 * t55 + t64 * t58;
t56 = sin(pkin(11));
t57 = sin(pkin(6));
t59 = cos(pkin(11));
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t20 = (t56 * t66 + t59 * t63) * t57;
t60 = cos(pkin(6));
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t13 = t20 * t62 - t60 * t65;
t12 = t13 ^ 2;
t86 = t57 * t66;
t87 = t57 * t63;
t18 = t56 * t87 - t59 * t86;
t102 = t18 ^ 2;
t48 = -t58 * pkin(5) - pkin(4);
t101 = 0.2e1 * t48;
t100 = 0.2e1 * t62;
t99 = -0.2e1 * t65;
t98 = t56 * pkin(2);
t97 = t59 * pkin(2);
t96 = t65 * pkin(4);
t95 = t13 * t65;
t32 = t64 * t55 + t61 * t58;
t23 = t32 * t62;
t94 = t32 * t23;
t93 = t32 * t65;
t49 = t55 ^ 2;
t92 = t49 * t62;
t46 = pkin(8) + t98;
t53 = t62 ^ 2;
t91 = t53 * t46;
t90 = t55 * t58;
t89 = t55 * t62;
t88 = t55 * t65;
t85 = t58 * t62;
t84 = t58 * t65;
t37 = t62 * t46;
t82 = t62 * t65;
t80 = t65 * t103;
t79 = t65 * t46;
t78 = pkin(9) + qJ(5);
t47 = -pkin(3) - t97;
t29 = -t62 * qJ(5) + t47 - t96;
t11 = t55 * t29 + t58 * t79;
t51 = t58 ^ 2;
t77 = t49 + t51;
t54 = t65 ^ 2;
t76 = t53 + t54;
t75 = 0.2e1 * t82;
t74 = t55 * t85;
t73 = t77 * qJ(5);
t15 = t20 * t65 + t60 * t62;
t5 = -t15 * t55 + t18 * t58;
t6 = t15 * t58 + t18 * t55;
t72 = -t5 * t55 + t6 * t58;
t71 = -pkin(4) * t62 + qJ(5) * t65;
t27 = t58 * t29;
t10 = -t55 * t79 + t27;
t70 = -t10 * t55 + t11 * t58;
t69 = t13 * t62 + t15 * t65;
t52 = t60 ^ 2;
t43 = t51 * t62;
t42 = t51 * t53;
t41 = t49 * t53;
t40 = t46 ^ 2;
t36 = t53 * t40;
t35 = t78 * t58;
t34 = t78 * t55;
t28 = pkin(5) * t89 + t37;
t25 = t103 * t62;
t22 = t25 ^ 2;
t21 = t23 ^ 2;
t17 = -t61 * t34 + t64 * t35;
t16 = -t64 * t34 - t61 * t35;
t9 = t25 * t103;
t8 = -pkin(9) * t89 + t11;
t7 = -pkin(9) * t85 + t27 + (-t46 * t55 - pkin(5)) * t65;
t4 = t61 * t7 + t64 * t8;
t3 = -t61 * t8 + t64 * t7;
t2 = t61 * t5 + t64 * t6;
t1 = t64 * t5 - t61 * t6;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 + (t63 ^ 2 + t66 ^ 2) * t57 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 ^ 2 + t102 + t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t102 + t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t87, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20, 0 (-t18 * t59 + t20 * t56) * pkin(2), 0, 0, 0, 0, 0, 0, -t18 * t65, t18 * t62, t69, t18 * t47 + t69 * t46, 0, 0, 0, 0, 0, 0, t13 * t89 - t5 * t65, t13 * t85 + t6 * t65 (-t5 * t58 - t55 * t6) * t62, t5 * t10 + t6 * t11 + t13 * t37, 0, 0, 0, 0, 0, 0, -t1 * t65 + t13 * t23, t13 * t25 + t2 * t65, -t1 * t25 - t2 * t23, t1 * t3 + t13 * t28 + t2 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t97, -0.2e1 * t98, 0 (t56 ^ 2 + t59 ^ 2) * pkin(2) ^ 2, t53, t75, 0, t54, 0, 0, t47 * t99, t47 * t100, 0.2e1 * t76 * t46, t54 * t40 + t47 ^ 2 + t36, t42, -0.2e1 * t53 * t90, -0.2e1 * t58 * t82, t41, t55 * t75, t54, -0.2e1 * t10 * t65 + 0.2e1 * t55 * t91, 0.2e1 * t11 * t65 + 0.2e1 * t58 * t91 (-t10 * t58 - t11 * t55) * t100, t10 ^ 2 + t11 ^ 2 + t36, t22, -0.2e1 * t25 * t23, t25 * t99, t21, -t23 * t99, t54, 0.2e1 * t28 * t23 - 0.2e1 * t3 * t65, 0.2e1 * t28 * t25 + 0.2e1 * t4 * t65, -0.2e1 * t4 * t23 - 0.2e1 * t3 * t25, t28 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t62 - t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t72 - t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t23 + t2 * t25 - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t70 - t79) * t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t23 + t4 * t25 - t28 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 + t41 + t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 + t21 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t58, t13 * t55, t72, -t13 * pkin(4) + qJ(5) * t72, 0, 0, 0, 0, 0, 0, -t13 * t103, t13 * t32, -t1 * t32 + t103 * t2, t1 * t16 + t13 * t48 + t2 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t65, 0, -t37, -t79, 0, 0, t74, t43 - t92, -t88, -t74, -t84, 0, -t58 * t37 + t55 * t71, t55 * t37 + t58 * t71, t70, -pkin(4) * t37 + qJ(5) * t70, t25 * t32, t9 - t94, -t93, -t23 * t103, -t80, 0, -t103 * t28 - t16 * t65 + t48 * t23, t17 * t65 + t48 * t25 + t28 * t32, t103 * t4 - t16 * t25 - t17 * t23 - t3 * t32, t3 * t16 + t4 * t17 + t28 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t88, t43 + t92, t62 * t73 + t96, 0, 0, 0, 0, 0, 0, t80, -t93, t9 + t94, -t23 * t16 + t25 * t17 - t65 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, 0.2e1 * t90, 0, t51, 0, 0, 0.2e1 * pkin(4) * t58, -0.2e1 * pkin(4) * t55, 0.2e1 * t73, t77 * qJ(5) ^ 2 + pkin(4) ^ 2, t32 ^ 2, 0.2e1 * t32 * t103, 0, t103 ^ 2, 0, 0, -t103 * t101, t32 * t101, 0.2e1 * t103 * t17 - 0.2e1 * t16 * t32, t16 ^ 2 + t17 ^ 2 + t48 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t85, 0, t37, 0, 0, 0, 0, 0, 0, t23, t25, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t55, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t103, t32, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, -t65, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t103, 0, t16, -t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t14;
