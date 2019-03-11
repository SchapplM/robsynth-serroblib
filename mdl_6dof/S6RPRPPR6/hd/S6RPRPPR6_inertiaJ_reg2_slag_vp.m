% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t62 = sin(qJ(3));
t63 = cos(qJ(3));
t34 = t58 * t63 + t60 * t62;
t57 = sin(pkin(10));
t59 = cos(pkin(10));
t61 = sin(qJ(6));
t91 = cos(qJ(6));
t67 = -t61 * t57 + t91 * t59;
t104 = t67 * t34;
t38 = -t58 * t62 + t60 * t63;
t13 = t67 * t38;
t103 = -0.2e1 * t38;
t102 = (t34 * t58 + t38 * t60) * pkin(3);
t29 = t34 ^ 2;
t99 = t38 ^ 2;
t101 = t29 + t99;
t64 = -pkin(1) - pkin(7);
t79 = -qJ(4) + t64;
t42 = t79 * t62;
t73 = t79 * t63;
t20 = t58 * t42 - t60 * t73;
t100 = t20 ^ 2;
t98 = 0.2e1 * t34;
t93 = t60 * pkin(3);
t49 = -pkin(4) - t93;
t43 = -t59 * pkin(5) + t49;
t97 = 0.2e1 * t43;
t50 = t62 * pkin(3) + qJ(2);
t96 = 0.2e1 * t50;
t95 = 0.2e1 * qJ(2);
t94 = t58 * pkin(3);
t46 = qJ(5) + t94;
t92 = pkin(8) + t46;
t90 = t13 * t67;
t89 = t20 * t38;
t40 = t91 * t57 + t61 * t59;
t88 = t38 * t40;
t87 = t38 * t49;
t86 = t40 * t34;
t85 = t57 * t34;
t84 = t57 * t38;
t83 = t57 * t59;
t82 = t59 * t38;
t18 = t34 * pkin(4) - t38 * qJ(5) + t50;
t22 = t60 * t42 + t58 * t73;
t6 = t57 * t18 + t59 * t22;
t52 = t57 ^ 2;
t53 = t59 ^ 2;
t80 = t52 + t53;
t54 = t62 ^ 2;
t55 = t63 ^ 2;
t44 = t54 + t55;
t78 = t34 * t103;
t77 = t57 * t82;
t74 = t80 * t46;
t5 = t59 * t18 - t57 * t22;
t72 = t5 * t59 + t6 * t57;
t71 = -t5 * t57 + t6 * t59;
t70 = -t22 * t34 + t89;
t69 = -t34 * t46 + t87;
t65 = qJ(2) ^ 2;
t41 = t44 * t64;
t33 = t40 ^ 2;
t30 = t67 ^ 2;
t27 = t92 * t59;
t26 = t92 * t57;
t25 = t59 * t34;
t17 = -t61 * t26 + t91 * t27;
t16 = -t91 * t26 - t61 * t27;
t8 = pkin(5) * t84 + t20;
t7 = t40 * t88;
t4 = -pkin(8) * t84 + t6;
t3 = t34 * pkin(5) - pkin(8) * t82 + t5;
t2 = t61 * t3 + t91 * t4;
t1 = t91 * t3 - t61 * t4;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t95 (pkin(1) ^ 2) + t65, t55, -0.2e1 * t63 * t62, 0, t54, 0, 0, t62 * t95, t63 * t95, -0.2e1 * t41, t44 * t64 ^ 2 + t65, t99, t78, 0, t29, 0, 0, t34 * t96, t38 * t96, 0.2e1 * t70, t22 ^ 2 + t50 ^ 2 + t100, t53 * t99, -0.2e1 * t99 * t83, t82 * t98, t52 * t99, t57 * t78, t29, 0.2e1 * t20 * t84 + 0.2e1 * t5 * t34, 0.2e1 * t20 * t82 - 0.2e1 * t6 * t34, t72 * t103, t5 ^ 2 + t6 ^ 2 + t100, t13 ^ 2, -0.2e1 * t13 * t88, t13 * t98, t88 ^ 2, -t88 * t98, t29, 0.2e1 * t1 * t34 + 0.2e1 * t8 * t88, 0.2e1 * t8 * t13 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t13 - 0.2e1 * t2 * t88, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t44, t41, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t70, 0, 0, 0, 0, 0, 0, -t101 * t57, -t101 * t59, 0, t71 * t34 - t89, 0, 0, 0, 0, 0, 0, -t34 * t86 - t38 * t88, -t104 * t34 - t38 * t13, -t104 * t88 + t13 * t86, -t1 * t86 + t104 * t2 - t8 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t29 + t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 ^ 2 + t86 ^ 2 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, -t62, 0, t63 * t64, -t62 * t64, 0, 0, 0, 0, t38, 0, -t34, 0, -t20, -t22, -t102 (-t20 * t60 + t22 * t58) * pkin(3), t77 (-t52 + t53) * t38, t85, -t77, t25, 0, -t20 * t59 + t57 * t69, t20 * t57 + t59 * t69, t71, t20 * t49 + t71 * t46, t13 * t40, -t7 + t90, t86, -t88 * t67, t104, 0, t16 * t34 + t43 * t88 - t67 * t8, t43 * t13 - t17 * t34 + t8 * t40, -t1 * t40 - t16 * t13 - t17 * t88 + t2 * t67, t1 * t16 + t2 * t17 + t8 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t62, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t34, 0, t102, 0, 0, 0, 0, 0, 0, t82, -t84, t80 * t34, t34 * t74 - t87, 0, 0, 0, 0, 0, 0, t13, -t88, t104 * t67 + t40 * t86, t104 * t17 - t16 * t86 - t38 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t93, -0.2e1 * t94, 0 (t58 ^ 2 + t60 ^ 2) * pkin(3) ^ 2, t52, 0.2e1 * t83, 0, t53, 0, 0, -0.2e1 * t49 * t59, 0.2e1 * t49 * t57, 0.2e1 * t74, t80 * t46 ^ 2 + t49 ^ 2, t33, 0.2e1 * t40 * t67, 0, t30, 0, 0, -t67 * t97, t40 * t97, -0.2e1 * t16 * t40 + 0.2e1 * t17 * t67, t16 ^ 2 + t17 ^ 2 + t43 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t38, 0, t50, 0, 0, 0, 0, 0, 0, t25, -t85, -t80 * t38, t72, 0, 0, 0, 0, 0, 0, t104, -t86, -t7 - t90, t1 * t67 + t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * t40 - t67 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t67 + t17 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t82, 0, t20, 0, 0, 0, 0, 0, 0, t88, t13, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t57, 0, t49, 0, 0, 0, 0, 0, 0, -t67, t40, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t88, t34, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t104, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t67, 0, t16, -t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t9;
