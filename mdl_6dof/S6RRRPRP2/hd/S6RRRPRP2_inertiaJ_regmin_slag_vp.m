% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP2
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:32:56
% EndTime: 2019-05-07 07:32:58
% DurationCPUTime: 0.88s
% Computational Cost: add. (1480->102), mult. (2720->194), div. (0->0), fcn. (3250->8), ass. (0->84)
t56 = sin(pkin(10));
t49 = t56 * pkin(3) + pkin(9);
t58 = sin(qJ(5));
t54 = t58 ^ 2;
t61 = cos(qJ(5));
t55 = t61 ^ 2;
t77 = t54 + t55;
t78 = t77 * t49;
t63 = cos(qJ(2));
t52 = -t63 * pkin(2) - pkin(1);
t101 = 0.2e1 * t52;
t100 = -0.2e1 * t58;
t99 = 0.2e1 * t58;
t98 = -0.2e1 * t61;
t97 = 0.2e1 * t63;
t96 = pkin(7) + pkin(8);
t59 = sin(qJ(3));
t60 = sin(qJ(2));
t62 = cos(qJ(3));
t38 = t59 * t63 + t62 * t60;
t57 = cos(pkin(10));
t67 = t59 * t60 - t62 * t63;
t23 = t56 * t38 + t57 * t67;
t24 = t57 * t38 - t56 * t67;
t30 = t67 * pkin(3) + t52;
t11 = t23 * pkin(4) - t24 * pkin(9) + t30;
t43 = t96 * t60;
t44 = t96 * t63;
t27 = t59 * t43 - t62 * t44;
t17 = -t67 * qJ(4) - t27;
t26 = -t62 * t43 - t59 * t44;
t65 = -t38 * qJ(4) + t26;
t14 = t57 * t17 + t56 * t65;
t5 = t58 * t11 + t61 * t14;
t95 = t23 * pkin(5);
t94 = t57 * pkin(3);
t93 = t59 * pkin(2);
t12 = t56 * t17 - t57 * t65;
t72 = pkin(5) * t58 - t61 * qJ(6);
t6 = t72 * t24 + t12;
t92 = t6 * t58;
t91 = t6 * t61;
t90 = t12 * t61;
t53 = t62 * pkin(2);
t51 = t53 + pkin(3);
t35 = t56 * t51 + t57 * t93;
t33 = pkin(9) + t35;
t89 = t23 * t33;
t88 = t23 * t49;
t19 = t58 * t23;
t87 = t58 * t24;
t86 = t58 * t33;
t85 = t58 * t49;
t84 = t58 * t61;
t20 = t61 * t23;
t21 = t61 * t24;
t83 = t61 * t33;
t82 = t61 * t49;
t34 = t57 * t51 - t56 * t93;
t73 = t61 * pkin(5) + t58 * qJ(6);
t66 = -pkin(4) - t73;
t25 = -t34 + t66;
t37 = t66 - t94;
t81 = -t25 - t37;
t80 = t77 * t33;
t32 = -pkin(4) - t34;
t50 = -pkin(4) - t94;
t79 = t32 + t50;
t76 = t23 * qJ(6);
t75 = -t61 * t11 + t58 * t14;
t2 = t76 + t5;
t3 = t75 - t95;
t1 = t2 * t61 + t3 * t58;
t74 = t2 * t58 - t3 * t61;
t71 = -t24 * t25 + t89;
t70 = t24 * t32 - t89;
t69 = -t24 * t37 + t88;
t68 = t24 * t50 - t88;
t48 = 0.2e1 * t84;
t22 = t24 ^ 2;
t18 = t58 * t21;
t15 = (-t54 + t55) * t24;
t10 = t12 * t58;
t4 = [1, 0, 0, t60 ^ 2, t60 * t97, 0, 0, 0, pkin(1) * t97, -0.2e1 * pkin(1) * t60, t38 ^ 2, -0.2e1 * t38 * t67, 0, 0, 0, t67 * t101, t38 * t101, 0.2e1 * t12 * t24 - 0.2e1 * t14 * t23, t12 ^ 2 + t14 ^ 2 + t30 ^ 2, t55 * t22, -0.2e1 * t22 * t84, 0.2e1 * t23 * t21, -0.2e1 * t23 * t87, t23 ^ 2, 0.2e1 * t12 * t87 - 0.2e1 * t23 * t75, 0.2e1 * t12 * t21 - 0.2e1 * t5 * t23, -0.2e1 * t3 * t23 + 0.2e1 * t6 * t87, -0.2e1 * t74 * t24, 0.2e1 * t2 * t23 - 0.2e1 * t6 * t21, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t60, t63, 0, -t60 * pkin(7), -t63 * pkin(7), 0, 0, t38, -t67, 0, t26, t27, -t35 * t23 - t34 * t24, -t12 * t34 + t14 * t35, t18, t15, t19, t20, 0, t70 * t58 - t90, t70 * t61 + t10, -t71 * t58 - t91, t1, t71 * t61 - t92, t1 * t33 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t93, 0, t34 ^ 2 + t35 ^ 2, t54, t48, 0, 0, 0, t32 * t98, t32 * t99, t25 * t98, 0.2e1 * t80, t25 * t100, t77 * t33 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t67, 0, t26, t27 (-t23 * t56 - t24 * t57) * pkin(3) (-t12 * t57 + t14 * t56) * pkin(3), t18, t15, t19, t20, 0, t68 * t58 - t90, t68 * t61 + t10, -t69 * t58 - t91, t1, t69 * t61 - t92, t1 * t49 + t6 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t93, 0 (t34 * t57 + t35 * t56) * pkin(3), t54, t48, 0, 0, 0, -t79 * t61, t79 * t58, t81 * t61, t78 + t80, t81 * t58, t25 * t37 + t33 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t56 ^ 2 + t57 ^ 2) * pkin(3) ^ 2, t54, t48, 0, 0, 0, t50 * t98, t50 * t99, t37 * t98, 0.2e1 * t78, t37 * t100, t77 * t49 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, t20, -t19, t20, -t77 * t24, t19, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t87, t23, -t75, -t5, -t75 + 0.2e1 * t95, -t73 * t24, 0.2e1 * t76 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t61, 0, -t86, -t83, -t86, -t72, t83, -t72 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t61, 0, -t85, -t82, -t85, -t72, t82, -t72 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t58, t61, 0, t58, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
