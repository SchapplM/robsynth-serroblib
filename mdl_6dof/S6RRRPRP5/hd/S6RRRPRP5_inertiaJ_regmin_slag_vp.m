% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:55:18
% EndTime: 2019-05-07 07:55:20
% DurationCPUTime: 0.88s
% Computational Cost: add. (1624->142), mult. (3228->247), div. (0->0), fcn. (3592->8), ass. (0->78)
t56 = sin(pkin(10));
t57 = cos(pkin(10));
t60 = sin(qJ(2));
t61 = cos(qJ(3));
t79 = t61 * t60;
t59 = sin(qJ(3));
t82 = t59 * t60;
t35 = -t56 * t82 + t57 * t79;
t58 = sin(qJ(5));
t40 = t56 * t61 + t57 * t59;
t67 = t40 * t60;
t85 = cos(qJ(5));
t18 = t35 * t85 - t58 * t67;
t94 = -0.2e1 * t18;
t68 = t56 * t59 - t57 * t61;
t23 = t40 * t85 - t58 * t68;
t93 = -0.2e1 * t23;
t50 = -t61 * pkin(3) - pkin(2);
t33 = pkin(4) * t68 + t50;
t92 = 0.2e1 * t33;
t91 = -0.2e1 * t60;
t62 = cos(qJ(2));
t90 = 0.2e1 * t62;
t89 = pkin(2) * t61;
t88 = pkin(7) * t59;
t87 = t56 * pkin(3);
t86 = t62 * pkin(5);
t45 = -pkin(2) * t62 - pkin(8) * t60 - pkin(1);
t41 = t61 * t45;
t76 = qJ(4) * t60;
t24 = -t61 * t76 + t41 + (-pkin(3) - t88) * t62;
t78 = t61 * t62;
t73 = pkin(7) * t78;
t27 = t73 + (t45 - t76) * t59;
t15 = t24 * t56 + t27 * t57;
t13 = -pkin(9) * t67 + t15;
t14 = t24 * t57 - t27 * t56;
t8 = -pkin(4) * t62 - pkin(9) * t35 + t14;
t4 = t13 * t85 + t58 * t8;
t77 = -qJ(4) - pkin(8);
t70 = t77 * t59;
t71 = t77 * t61;
t29 = t56 * t70 - t57 * t71;
t19 = -pkin(9) * t68 + t29;
t28 = t56 * t71 + t57 * t70;
t66 = -t40 * pkin(9) + t28;
t11 = t19 * t58 - t66 * t85;
t84 = t11 * t62;
t12 = t19 * t85 + t58 * t66;
t83 = t12 * t62;
t81 = t59 * t61;
t80 = t59 * t62;
t51 = t60 * pkin(7);
t44 = pkin(3) * t82 + t51;
t75 = t62 * qJ(6);
t74 = t60 * t90;
t72 = t58 * t13 - t8 * t85;
t49 = pkin(3) * t57 + pkin(4);
t69 = -t49 * t85 + t58 * t87;
t38 = t49 * t58 + t85 * t87;
t25 = pkin(4) * t67 + t44;
t64 = 0.2e1 * pkin(5);
t63 = 0.2e1 * qJ(6);
t55 = t62 ^ 2;
t54 = t61 ^ 2;
t53 = t60 ^ 2;
t52 = t59 ^ 2;
t36 = -pkin(5) + t69;
t34 = qJ(6) + t38;
t32 = t45 * t59 + t73;
t31 = -pkin(7) * t80 + t41;
t22 = t40 * t58 + t68 * t85;
t17 = t35 * t58 + t67 * t85;
t9 = t22 * pkin(5) - t23 * qJ(6) + t33;
t5 = t17 * pkin(5) - t18 * qJ(6) + t25;
t2 = t72 + t86;
t1 = -t75 + t4;
t3 = [1, 0, 0, t53, t74, 0, 0, 0, pkin(1) * t90, pkin(1) * t91, t54 * t53, -0.2e1 * t53 * t81, t78 * t91, t59 * t74, t55, -0.2e1 * t31 * t62 + 0.2e1 * t53 * t88, 0.2e1 * pkin(7) * t53 * t61 + 0.2e1 * t32 * t62, -0.2e1 * t14 * t35 - 0.2e1 * t15 * t67, t14 ^ 2 + t15 ^ 2 + t44 ^ 2, t18 ^ 2, t17 * t94, t62 * t94, t17 * t90, t55, 0.2e1 * t17 * t25 + 0.2e1 * t62 * t72, 0.2e1 * t18 * t25 + 0.2e1 * t4 * t62, 0.2e1 * t17 * t5 + 0.2e1 * t2 * t62, -0.2e1 * t1 * t17 + 0.2e1 * t18 * t2, -0.2e1 * t1 * t62 - 0.2e1 * t18 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t60, t62, 0, -t51, -t62 * pkin(7), t59 * t79 (-t52 + t54) * t60, -t80, -t78, 0, -pkin(7) * t79 + (-pkin(2) * t60 + pkin(8) * t62) * t59, pkin(8) * t78 + (t88 - t89) * t60, -t14 * t40 - t15 * t68 - t28 * t35 - t29 * t67, t14 * t28 + t15 * t29 + t44 * t50, t18 * t23, -t17 * t23 - t18 * t22, -t23 * t62, t22 * t62, 0, t17 * t33 + t22 * t25 + t84, t18 * t33 + t23 * t25 + t83, t17 * t9 + t22 * t5 + t84, -t1 * t22 + t11 * t18 - t12 * t17 + t2 * t23, -t18 * t9 - t23 * t5 - t83, t1 * t12 + t11 * t2 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, 0.2e1 * t81, 0, 0, 0, 0.2e1 * t89, -0.2e1 * pkin(2) * t59, -0.2e1 * t28 * t40 - 0.2e1 * t29 * t68, t28 ^ 2 + t29 ^ 2 + t50 ^ 2, t23 ^ 2, t22 * t93, 0, 0, 0, t22 * t92, t23 * t92, 0.2e1 * t9 * t22, 0.2e1 * t11 * t23 - 0.2e1 * t12 * t22, t9 * t93, t11 ^ 2 + t12 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t82, -t62, t31, -t32 (-t57 * t35 - t56 * t67) * pkin(3) (t14 * t57 + t15 * t56) * pkin(3), 0, 0, t18, -t17, -t62, t62 * t69 - t72, t38 * t62 - t4 (-pkin(5) + t36) * t62 - t72, -t17 * t34 + t18 * t36 (-qJ(6) - t34) * t62 + t4, t1 * t34 + t2 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t61, 0, -t59 * pkin(8), -t61 * pkin(8) (-t57 * t40 - t56 * t68) * pkin(3) (t28 * t57 + t29 * t56) * pkin(3), 0, 0, t23, -t22, 0, -t11, -t12, -t11, -t22 * t34 + t23 * t36, t12, t11 * t36 + t12 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t56 ^ 2 + t57 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, -0.2e1 * t69, -0.2e1 * t38, -0.2e1 * t36, 0, 0.2e1 * t34, t34 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, t17, t18, t17, 0, -t18, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, t22, t23, t22, 0, -t23, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t62, -t72, -t4, -t72 - 0.2e1 * t86, -pkin(5) * t18 - qJ(6) * t17, -0.2e1 * t75 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t11, -t12, -t11, -pkin(5) * t23 - qJ(6) * t22, t12, -pkin(5) * t11 + qJ(6) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t69, -t38, t64 - t69, 0, t63 + t38, -pkin(5) * t36 + qJ(6) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t64, 0, t63, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
