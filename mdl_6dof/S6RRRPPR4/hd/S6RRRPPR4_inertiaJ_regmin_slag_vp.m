% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:46:09
% EndTime: 2019-05-07 04:46:12
% DurationCPUTime: 0.70s
% Computational Cost: add. (971->131), mult. (1935->228), div. (0->0), fcn. (2127->8), ass. (0->76)
t64 = sin(pkin(10));
t65 = cos(pkin(10));
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t41 = t64 * t67 - t65 * t70;
t42 = t64 * t70 + t65 * t67;
t57 = -t70 * pkin(3) - pkin(2);
t74 = t42 * qJ(5) - t57;
t92 = -pkin(4) - pkin(5);
t13 = t92 * t41 + t74;
t96 = 0.2e1 * t13;
t68 = sin(qJ(2));
t38 = t42 * t68;
t86 = t70 * t68;
t89 = t67 * t68;
t39 = -t64 * t89 + t65 * t86;
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t16 = -t69 * t38 + t66 * t39;
t95 = -0.2e1 * t16;
t94 = -0.2e1 * t68;
t71 = cos(qJ(2));
t93 = 0.2e1 * t71;
t91 = pkin(2) * t70;
t90 = pkin(7) * t67;
t88 = t67 * t70;
t87 = t67 * t71;
t85 = t70 * t71;
t84 = -qJ(4) - pkin(8);
t49 = -t71 * pkin(2) - t68 * pkin(8) - pkin(1);
t43 = t70 * t49;
t83 = qJ(4) * t68;
t25 = -t70 * t83 + t43 + (-pkin(3) - t90) * t71;
t80 = pkin(7) * t85;
t28 = t80 + (t49 - t83) * t67;
t10 = t65 * t25 - t64 * t28;
t11 = t64 * t25 + t65 * t28;
t50 = t84 * t70;
t79 = t84 * t67;
t32 = -t65 * t50 + t64 * t79;
t58 = t68 * pkin(7);
t47 = pkin(3) * t89 + t58;
t82 = t68 * t93;
t30 = -t64 * t50 - t65 * t79;
t81 = t30 ^ 2 + t32 ^ 2;
t9 = t71 * pkin(4) - t10;
t55 = t65 * pkin(3) + pkin(4);
t78 = t30 * t39 - t32 * t38;
t77 = -pkin(5) - t55;
t7 = -t71 * qJ(5) + t11;
t76 = t39 * qJ(5) - t47;
t3 = t71 * pkin(5) - t39 * pkin(9) + t9;
t4 = t38 * pkin(9) + t7;
t1 = t69 * t3 - t66 * t4;
t2 = t66 * t3 + t69 * t4;
t75 = 0.2e1 * t30 * t42 - 0.2e1 * t32 * t41;
t73 = -t42 * pkin(9) + t30;
t63 = t71 ^ 2;
t62 = t70 ^ 2;
t61 = t68 ^ 2;
t60 = t67 ^ 2;
t53 = t64 * pkin(3) + qJ(5);
t37 = t69 * t53 + t66 * t77;
t36 = t66 * t53 - t69 * t77;
t35 = t67 * t49 + t80;
t34 = -pkin(7) * t87 + t43;
t24 = t66 * t41 + t69 * t42;
t23 = -t69 * t41 + t66 * t42;
t20 = t41 * pkin(4) - t74;
t18 = t41 * pkin(9) + t32;
t17 = t66 * t38 + t69 * t39;
t12 = t38 * pkin(4) - t76;
t8 = t92 * t38 + t76;
t6 = t69 * t18 + t66 * t73;
t5 = t66 * t18 - t69 * t73;
t14 = [1, 0, 0, t61, t82, 0, 0, 0, pkin(1) * t93, pkin(1) * t94, t62 * t61, -0.2e1 * t61 * t88, t85 * t94, t67 * t82, t63, -0.2e1 * t34 * t71 + 0.2e1 * t61 * t90, 0.2e1 * t61 * pkin(7) * t70 + 0.2e1 * t35 * t71, -0.2e1 * t10 * t39 - 0.2e1 * t11 * t38, t10 ^ 2 + t11 ^ 2 + t47 ^ 2, 0.2e1 * t12 * t38 + 0.2e1 * t9 * t71, -0.2e1 * t7 * t38 + 0.2e1 * t9 * t39, -0.2e1 * t12 * t39 - 0.2e1 * t7 * t71, t12 ^ 2 + t7 ^ 2 + t9 ^ 2, t17 ^ 2, t17 * t95, t17 * t93, t71 * t95, t63, 0.2e1 * t1 * t71 + 0.2e1 * t8 * t16, 0.2e1 * t8 * t17 - 0.2e1 * t2 * t71; 0, 0, 0, 0, 0, t68, t71, 0, -t58, -t71 * pkin(7), t67 * t86 (-t60 + t62) * t68, -t87, -t85, 0, -pkin(7) * t86 + (-pkin(2) * t68 + pkin(8) * t71) * t67, pkin(8) * t85 + (t90 - t91) * t68, -t10 * t42 - t11 * t41 + t78, -t10 * t30 + t11 * t32 + t47 * t57, t12 * t41 + t20 * t38 + t30 * t71, -t7 * t41 + t9 * t42 + t78, -t12 * t42 - t20 * t39 - t32 * t71, t12 * t20 + t9 * t30 + t7 * t32, t17 * t24, -t24 * t16 - t17 * t23, t24 * t71, -t23 * t71, 0, t13 * t16 + t8 * t23 - t5 * t71, t13 * t17 + t8 * t24 - t6 * t71; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t60, 0.2e1 * t88, 0, 0, 0, 0.2e1 * t91, -0.2e1 * pkin(2) * t67, t75, t57 ^ 2 + t81, 0.2e1 * t20 * t41, t75, -0.2e1 * t20 * t42, t20 ^ 2 + t81, t24 ^ 2, -0.2e1 * t24 * t23, 0, 0, 0, t23 * t96, t24 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t89, -t71, t34, -t35 (-t38 * t64 - t39 * t65) * pkin(3) (t10 * t65 + t11 * t64) * pkin(3), -t55 * t71 - t9, -t53 * t38 - t55 * t39 (-qJ(5) - t53) * t71 + t11, t7 * t53 - t9 * t55, 0, 0, -t17, t16, -t71, -t36 * t71 - t1, -t37 * t71 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t67 * pkin(8), -t70 * pkin(8) (-t41 * t64 - t42 * t65) * pkin(3) (-t30 * t65 + t32 * t64) * pkin(3), -t30, -t53 * t41 - t55 * t42, t32, -t30 * t55 + t32 * t53, 0, 0, -t24, t23, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t64 ^ 2 + t65 ^ 2) * pkin(3) ^ 2, 0.2e1 * t55, 0, 0.2e1 * t53, t53 ^ 2 + t55 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t36, 0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t38, 0, -t39, t12, 0, 0, 0, 0, 0, -t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t41, 0, -t42, t20, 0, 0, 0, 0, 0, -t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t39, 0, t9, 0, 0, 0, 0, 0, t69 * t71, -t66 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t55, 0, 0, 0, 0, 0, -t69, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t71, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t14;
