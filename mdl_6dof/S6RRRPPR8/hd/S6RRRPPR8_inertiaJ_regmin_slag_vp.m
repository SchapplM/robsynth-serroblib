% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:16:27
% EndTime: 2019-05-07 06:16:32
% DurationCPUTime: 1.09s
% Computational Cost: add. (860->155), mult. (1942->275), div. (0->0), fcn. (2101->8), ass. (0->89)
t62 = cos(pkin(6));
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t61 = sin(pkin(6));
t66 = sin(qJ(2));
t97 = t61 * t66;
t32 = t62 * t65 + t68 * t97;
t108 = -0.2e1 * t32;
t107 = -0.2e1 * t65;
t106 = -0.2e1 * t68;
t105 = 0.2e1 * t68;
t104 = -pkin(4) - pkin(10);
t103 = pkin(1) * t66;
t102 = t68 * pkin(9);
t31 = -t62 * t68 + t65 * t97;
t69 = cos(qJ(2));
t47 = t61 * t69;
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t17 = t31 * t67 + t64 * t47;
t101 = t17 * t64;
t100 = t32 * t64;
t26 = t32 * t65;
t99 = t32 * t67;
t55 = t61 ^ 2;
t98 = t55 * t69;
t96 = t62 * t66;
t56 = -pkin(3) + t104;
t95 = t64 * t56;
t94 = t64 * t65;
t93 = t64 * t67;
t92 = t64 * t68;
t91 = t67 * t56;
t90 = t67 * t65;
t89 = t67 * t68;
t63 = qJ(4) + pkin(5);
t82 = pkin(8) * t47;
t24 = t82 + (pkin(9) + t103) * t62;
t25 = (-pkin(2) * t69 - pkin(9) * t66 - pkin(1)) * t61;
t88 = t65 * t24 - t68 * t25;
t13 = t68 * t24 + t65 * t25;
t33 = t62 * t69 * pkin(1) - pkin(8) * t97;
t58 = t65 ^ 2;
t60 = t68 ^ 2;
t87 = t58 + t60;
t86 = t31 * qJ(4);
t85 = t68 * qJ(4);
t84 = 0.2e1 * t47;
t83 = t65 * t105;
t36 = -t68 * pkin(3) - t65 * qJ(4) - pkin(2);
t81 = t65 * t47;
t80 = t68 * t47;
t45 = pkin(3) * t47;
t11 = t45 + t88;
t78 = qJ(4) * t47;
t79 = -0.2e1 * t78 + t13;
t23 = -t62 * pkin(2) - t33;
t77 = pkin(9) * t80;
t35 = t68 * pkin(4) - t36;
t9 = t31 * pkin(3) - t32 * qJ(4) + t23;
t76 = -pkin(3) * t65 + t85;
t10 = -t78 + t13;
t75 = t10 * t68 + t11 * t65;
t74 = -t56 * t65 - t63 * t68;
t6 = pkin(4) * t47 - t32 * qJ(5) + t11;
t72 = qJ(4) ^ 2;
t71 = 0.2e1 * qJ(4);
t70 = -pkin(3) - pkin(4);
t59 = t67 ^ 2;
t57 = t64 ^ 2;
t52 = t65 * pkin(9);
t39 = -t68 * qJ(5) + t102;
t38 = -t65 * qJ(5) + t52;
t37 = pkin(9) * t81;
t34 = pkin(1) * t96 + t82;
t30 = t32 ^ 2;
t27 = t31 * qJ(5);
t18 = t65 * pkin(5) + t68 * pkin(10) + t35;
t16 = t31 * t64 - t67 * t47;
t15 = t64 * t18 + t67 * t38;
t14 = t67 * t18 - t64 * t38;
t8 = -t10 - t27;
t7 = -t31 * pkin(4) - t9;
t5 = -t63 * t47 + t13 + t27;
t4 = pkin(10) * t47 + t6;
t3 = t32 * pkin(5) + t104 * t31 - t9;
t2 = t64 * t3 + t67 * t4;
t1 = t67 * t3 - t64 * t4;
t12 = [1, 0, 0, t55 * t66 ^ 2, 0.2e1 * t66 * t98, 0.2e1 * t61 * t96, t62 * t84, t62 ^ 2, 0.2e1 * pkin(1) * t98 + 0.2e1 * t33 * t62, -0.2e1 * t55 * t103 - 0.2e1 * t34 * t62, t30, t31 * t108, t47 * t108, t31 * t84, t55 * t69 ^ 2, 0.2e1 * t23 * t31 + 0.2e1 * t47 * t88, 0.2e1 * t13 * t47 + 0.2e1 * t23 * t32, 0.2e1 * t11 * t47 + 0.2e1 * t9 * t31, -0.2e1 * t10 * t31 + 0.2e1 * t11 * t32, -0.2e1 * t10 * t47 - 0.2e1 * t9 * t32, t10 ^ 2 + t11 ^ 2 + t9 ^ 2, 0.2e1 * t7 * t32 + 0.2e1 * t8 * t47, 0.2e1 * t7 * t31 - 0.2e1 * t6 * t47, -0.2e1 * t8 * t31 - 0.2e1 * t6 * t32, t6 ^ 2 + t7 ^ 2 + t8 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t32, t16 * t108, t30, 0.2e1 * t1 * t32 + 0.2e1 * t5 * t16, 0.2e1 * t5 * t17 - 0.2e1 * t2 * t32; 0, 0, 0, 0, 0, t97, t47, t62, t33, -t34, t26, -t65 * t31 + t32 * t68, -t81, -t80, 0, -pkin(2) * t31 - t23 * t68 + t37, -pkin(2) * t32 + t23 * t65 + t77, t36 * t31 - t9 * t68 + t37 (-t31 * t68 + t26) * pkin(9) + t75, -t36 * t32 - t9 * t65 - t77, pkin(9) * t75 + t9 * t36, t35 * t32 - t39 * t47 + t7 * t65, t35 * t31 - t38 * t47 - t7 * t68, t39 * t31 - t38 * t32 - t6 * t65 + t8 * t68, t7 * t35 + t6 * t38 - t8 * t39, -t17 * t89 (t16 * t67 + t101) * t68, t17 * t65 - t32 * t89, -t16 * t65 + t32 * t92, t26, t1 * t65 + t14 * t32 + t39 * t16 - t5 * t92, -t15 * t32 + t39 * t17 - t2 * t65 - t5 * t89; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, t83, 0, 0, 0, pkin(2) * t105, pkin(2) * t107, t36 * t106, 0.2e1 * t87 * pkin(9), t36 * t107, t87 * pkin(9) ^ 2 + t36 ^ 2, 0.2e1 * t35 * t65, t35 * t106, -0.2e1 * t38 * t65 - 0.2e1 * t39 * t68, t35 ^ 2 + t38 ^ 2 + t39 ^ 2, t59 * t60, -0.2e1 * t60 * t93, t89 * t107, t64 * t83, t58, 0.2e1 * t14 * t65 - 0.2e1 * t39 * t92, -0.2e1 * t15 * t65 - 0.2e1 * t39 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, -t47, -t88, -t13, -0.2e1 * t45 - t88, -t32 * pkin(3) - t86, t79, -t11 * pkin(3) + t10 * qJ(4), t27 + t79, -t70 * t47 + t6, -t70 * t32 + t86, -t8 * qJ(4) + t6 * t70, -t101, t64 * t16 - t17 * t67, -t100, -t99, 0, t63 * t16 - t32 * t95 + t5 * t67, t63 * t17 - t32 * t91 - t5 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t68, 0, -t52, -t102, -t52, t76, t102, t76 * pkin(9), t39, t38, -t70 * t65 - t85, t39 * qJ(4) + t38 * t70, t64 * t89 (-t57 + t59) * t68, -t94, -t90, 0, t39 * t67 + t64 * t74, -t39 * t64 + t67 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t71, pkin(3) ^ 2 + t72, t71, 0.2e1 * t70, 0, t70 ^ 2 + t72, t57, 0.2e1 * t93, 0, 0, 0, 0.2e1 * t63 * t67, -0.2e1 * t63 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t32, 0, t11, 0, -t47, -t32, t6, 0, 0, 0, 0, 0, -t100, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t52, 0, 0, -t65, t38, 0, 0, 0, 0, 0, -t94, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 1, 0, t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t31, 0, t7, 0, 0, 0, 0, 0, t99, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t68, 0, t35, 0, 0, 0, 0, 0, t90, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t32, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t92, t65, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t67, 0, -t95, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
