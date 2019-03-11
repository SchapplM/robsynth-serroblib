% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t68 = sin(pkin(10));
t64 = t68 ^ 2;
t70 = cos(pkin(10));
t66 = t70 ^ 2;
t57 = t64 + t66;
t69 = sin(pkin(9));
t71 = cos(pkin(9));
t73 = sin(qJ(3));
t97 = cos(qJ(3));
t51 = t97 * t69 + t73 * t71;
t107 = -0.2e1 * t51;
t72 = sin(qJ(6));
t74 = cos(qJ(6));
t50 = t74 * t68 - t72 * t70;
t91 = pkin(7) + qJ(2);
t56 = t91 * t71;
t85 = t91 * t69;
t24 = t73 * t56 + t97 * t85;
t106 = t24 ^ 2;
t48 = t73 * t69 - t97 * t71;
t105 = t48 ^ 2;
t84 = t68 * qJ(5) + pkin(3);
t98 = pkin(4) + pkin(5);
t38 = t98 * t70 + t84;
t104 = 0.2e1 * t38;
t103 = -0.2e1 * t48;
t102 = 0.2e1 * t48;
t60 = -t71 * pkin(2) - pkin(1);
t101 = 0.2e1 * t60;
t100 = -0.2e1 * t68;
t99 = 0.2e1 * t71;
t46 = t72 * t68 + t74 * t70;
t15 = t46 * t51;
t96 = t15 * t46;
t95 = t50 * t48;
t33 = t68 * t48;
t34 = t68 * t51;
t94 = t68 * t70;
t36 = t70 * t48;
t37 = t70 * t51;
t17 = t48 * pkin(3) - t51 * qJ(4) + t60;
t27 = t97 * t56 - t73 * t85;
t10 = t68 * t17 + t70 * t27;
t90 = t57 * qJ(4) ^ 2;
t65 = t69 ^ 2;
t67 = t71 ^ 2;
t89 = t65 + t67;
t88 = qJ(4) * t48;
t87 = t48 * t34;
t44 = t51 ^ 2;
t86 = t44 * t94;
t5 = t48 * qJ(5) + t10;
t20 = t68 * t27;
t9 = t70 * t17 - t20;
t83 = -qJ(5) * t37 + t24;
t6 = -t48 * pkin(4) - t9;
t82 = t5 * t70 + t6 * t68;
t81 = t5 * t68 - t6 * t70;
t80 = t10 * t70 - t9 * t68;
t79 = t10 * t68 + t9 * t70;
t78 = -pkin(3) * t51 - t88;
t53 = -t70 * pkin(4) - t84;
t77 = -t51 * t53 + t88;
t62 = t68 * qJ(4);
t55 = (-pkin(8) + qJ(4)) * t70;
t54 = -t68 * pkin(8) + t62;
t52 = 0.2e1 * t57 * qJ(4);
t43 = t50 ^ 2;
t40 = t46 ^ 2;
t35 = t66 * t44;
t32 = t64 * t44;
t29 = t68 * t37;
t28 = t46 * t48;
t26 = t72 * t54 + t74 * t55;
t23 = t74 * t54 - t72 * t55;
t22 = t37 * t102;
t19 = t57 * t51;
t18 = (t64 - t66) * t51;
t13 = t50 * t51;
t12 = t50 * t13;
t11 = pkin(4) * t34 + t83;
t7 = t98 * t34 + t83;
t4 = pkin(8) * t34 + t5;
t3 = t20 + (-pkin(8) * t51 - t17) * t70 - t98 * t48;
t2 = t72 * t3 + t74 * t4;
t1 = t74 * t3 - t72 * t4;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t65, t69 * t99, 0, t67, 0, 0, pkin(1) * t99, -0.2e1 * pkin(1) * t69, 0.2e1 * t89 * qJ(2), t89 * qJ(2) ^ 2 + pkin(1) ^ 2, t44, t51 * t103, 0, t105, 0, 0, t48 * t101, t51 * t101, 0.2e1 * t24 * t51 - 0.2e1 * t27 * t48, t27 ^ 2 + t60 ^ 2 + t106, t35, -0.2e1 * t86, t22, t32, -0.2e1 * t87, t105, 0.2e1 * t24 * t34 + 0.2e1 * t9 * t48, -0.2e1 * t10 * t48 + 0.2e1 * t24 * t37, t79 * t107, t10 ^ 2 + t9 ^ 2 + t106, t35, t22, 0.2e1 * t86, t105, 0.2e1 * t87, t32, 0.2e1 * t11 * t34 - 0.2e1 * t6 * t48, t81 * t107, -0.2e1 * t11 * t37 + 0.2e1 * t5 * t48, t11 ^ 2 + t5 ^ 2 + t6 ^ 2, t15 ^ 2, 0.2e1 * t15 * t13, t15 * t103, t13 ^ 2, -t13 * t102, t105, -0.2e1 * t1 * t48 + 0.2e1 * t7 * t13, -0.2e1 * t7 * t15 + 0.2e1 * t2 * t48, -0.2e1 * t1 * t15 + 0.2e1 * t2 * t13, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t48, t51, 0, t60, 0, 0, 0, 0, 0, 0, t36, -t33, -t19, t79, 0, 0, 0, 0, 0, 0, t36, -t19, t33, t81, 0, 0, 0, 0, 0, 0, t28, t95, t12 + t96, -t1 * t46 + t2 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t48, 0, -t24, -t27, 0, 0, t29, -t18, t33, -t29, t36, 0, -t24 * t70 + t78 * t68, t24 * t68 + t78 * t70, t80, -t24 * pkin(3) + t80 * qJ(4), t29, t33, t18, 0, -t36, -t29, -t11 * t70 - t68 * t77, t82, -t11 * t68 + t70 * t77, t82 * qJ(4) + t11 * t53, t15 * t50, t12 - t96, -t95, -t13 * t46, t28, 0, -t38 * t13 - t23 * t48 - t7 * t46, t38 * t15 + t26 * t48 - t7 * t50, -t1 * t50 + t26 * t13 - t23 * t15 - t2 * t46, t1 * t23 + t2 * t26 - t7 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t23 + t50 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t64, 0.2e1 * t94, 0, t66, 0, 0, 0.2e1 * pkin(3) * t70, pkin(3) * t100, t52, pkin(3) ^ 2 + t90, t64, 0, -0.2e1 * t94, 0, 0, t66, -0.2e1 * t53 * t70, t52, t53 * t100, t53 ^ 2 + t90, t43, -0.2e1 * t50 * t46, 0, t40, 0, 0, t46 * t104, t50 * t104, -0.2e1 * t23 * t50 - 0.2e1 * t26 * t46, t23 ^ 2 + t26 ^ 2 + t38 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t37, 0, t24, 0, 0, 0, 0, 0, 0, t34, 0, -t37, t11, 0, 0, 0, 0, 0, 0, t13, -t15, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t68, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t70, 0, -t68, t53, 0, 0, 0, 0, 0, 0, -t46, -t50, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t37, 0, t6, 0, 0, 0, 0, 0, 0, -t74 * t48, t72 * t48, t72 * t13 - t74 * t15, t1 * t74 + t2 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t74 + t50 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t46 - t74 * t50, t23 * t74 + t26 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 ^ 2 + t74 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, t13, -t48, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t46, 0, t23, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t8;
