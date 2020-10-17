% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:38:47
% EndTime: 2019-05-06 12:38:50
% DurationCPUTime: 0.71s
% Computational Cost: add. (858->106), mult. (1437->183), div. (0->0), fcn. (1498->6), ass. (0->70)
t57 = sin(pkin(9));
t58 = cos(pkin(9));
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t29 = -t57 * t59 + t58 * t61;
t67 = t57 * t61 + t58 * t59;
t99 = (t29 * t58 + t57 * t67) * pkin(4);
t80 = t29 ^ 2 + t67 ^ 2;
t93 = pkin(2) + pkin(8);
t79 = -qJ(5) - t93;
t34 = t79 * t59;
t73 = t79 * t61;
t20 = t57 * t34 - t58 * t73;
t22 = t58 * t34 + t57 * t73;
t76 = -t20 * t29 + t22 * t67;
t46 = t59 * pkin(4) + qJ(3);
t15 = pkin(5) * t67 - t29 * qJ(6) + t46;
t97 = -0.2e1 * t15;
t60 = sin(qJ(2));
t96 = -0.2e1 * t60;
t62 = cos(qJ(2));
t95 = 0.2e1 * t62;
t94 = 0.2e1 * qJ(3);
t50 = t60 * pkin(7);
t36 = t60 * pkin(3) + t50;
t33 = t61 * t36;
t77 = -t60 * qJ(3) - pkin(1);
t27 = -t93 * t62 + t77;
t74 = qJ(5) * t62 - t27;
t12 = t60 * pkin(4) + t74 * t59 + t33;
t90 = t59 * t36;
t14 = -t74 * t61 + t90;
t4 = t57 * t12 + t58 * t14;
t89 = t59 * t60;
t88 = t59 * t62;
t87 = t60 * t62;
t86 = t61 * t59;
t85 = t61 * t62;
t51 = t62 * pkin(7);
t37 = t62 * pkin(3) + t51;
t54 = t60 ^ 2;
t56 = t62 ^ 2;
t84 = t54 + t56;
t83 = t62 * qJ(3);
t82 = -0.2e1 * t87;
t81 = t20 ^ 2 + t22 ^ 2;
t1 = t60 * qJ(6) + t4;
t25 = pkin(4) * t85 + t37;
t3 = t58 * t12 - t57 * t14;
t23 = t57 * t88 - t58 * t85;
t24 = t67 * t62;
t78 = -t20 * t24 + t22 * t23;
t75 = t23 * t67 + t29 * t24;
t2 = -t60 * pkin(5) - t3;
t72 = t1 * t67 - t2 * t29;
t71 = t3 * t29 + t4 * t67;
t70 = -t60 * pkin(2) + t83;
t40 = t57 * pkin(4) + qJ(6);
t44 = t58 * pkin(4) + pkin(5);
t69 = t44 * t29 + t40 * t67;
t66 = 0.2e1 * t76;
t65 = -t93 * t60 + t83;
t55 = t61 ^ 2;
t53 = t59 ^ 2;
t47 = t61 * t60;
t35 = -t62 * pkin(2) + t77;
t18 = t61 * t27 + t90;
t17 = -t59 * t27 + t33;
t5 = -t23 * pkin(5) + t24 * qJ(6) + t25;
t6 = [1, 0, 0, t54, 0.2e1 * t87, 0, 0, 0, pkin(1) * t95, pkin(1) * t96, 0.2e1 * t84 * pkin(7), t35 * t95, t35 * t96, t84 * pkin(7) ^ 2 + t35 ^ 2, t53 * t56, 0.2e1 * t56 * t86, t59 * t82, t61 * t82, t54, 0.2e1 * t17 * t60 + 0.2e1 * t37 * t85, -0.2e1 * t18 * t60 - 0.2e1 * t37 * t88, 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24, t25 ^ 2 + t3 ^ 2 + t4 ^ 2, -0.2e1 * t2 * t60 - 0.2e1 * t5 * t23, 0.2e1 * t1 * t23 - 0.2e1 * t2 * t24, 0.2e1 * t1 * t60 + 0.2e1 * t5 * t24, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t60, t62, 0, -t50, -t51, t70, t50, t51, t70 * pkin(7), -t59 * t85 (t53 - t55) * t62, t47, -t89, 0, t37 * t59 + t65 * t61, t37 * t61 - t65 * t59, -t71 + t78, -t3 * t20 + t4 * t22 + t25 * t46, -t15 * t23 - t20 * t60 + t5 * t67, -t72 + t78, t15 * t24 + t22 * t60 - t5 * t29, t1 * t22 + t5 * t15 + t2 * t20; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t94, pkin(2) ^ 2 + qJ(3) ^ 2, t55, -0.2e1 * t86, 0, 0, 0, t59 * t94, t61 * t94, -t66, t46 ^ 2 + t81, -t67 * t97, -t66, t29 * t97, t15 ^ 2 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, t50, 0, 0, 0, 0, 0, t47, -t89, t75, t71, t29 * t60, t75, t67 * t60, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, -t80, t76, 0, -t80, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t85, t60, t17, -t18 (t23 * t57 + t24 * t58) * pkin(4) (t3 * t58 + t4 * t57) * pkin(4) (pkin(5) + t44) * t60 + t3, t40 * t23 + t44 * t24, t40 * t60 + t1, t1 * t40 - t2 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t59, 0, -t61 * t93, t59 * t93, -t99 (-t20 * t58 + t22 * t57) * pkin(4), -t20, -t69, t22, -t20 * t44 + t22 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t59, 0, t99, t29, 0, t67, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t57 ^ 2 + t58 ^ 2) * pkin(4) ^ 2, 0.2e1 * t44, 0, 0.2e1 * t40, t40 ^ 2 + t44 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23, 0, t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t67, 0, -t29, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
