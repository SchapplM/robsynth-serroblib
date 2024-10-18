% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:16
% DurationCPUTime: 0.48s
% Computational Cost: add. (625->99), mult. (1510->133), div. (0->0), fcn. (1533->10), ass. (0->86)
t72 = sin(pkin(5));
t75 = sin(qJ(4));
t62 = t72 * t75;
t79 = cos(qJ(4));
t63 = t72 * t79;
t74 = sin(qJ(5));
t78 = cos(qJ(5));
t40 = t74 * t62 - t78 * t63;
t104 = -0.2e1 * t40;
t73 = cos(pkin(5));
t103 = 0.2e1 * t73;
t102 = pkin(4) * t79;
t101 = pkin(10) * t72;
t67 = t73 * pkin(4);
t100 = t74 * pkin(4);
t76 = sin(qJ(3));
t99 = t76 * pkin(2);
t98 = sin(qJ(2)) * pkin(1);
t97 = t78 * pkin(4);
t69 = cos(qJ(2)) * pkin(1);
t65 = t69 + pkin(2);
t80 = cos(qJ(3));
t82 = t80 * t98;
t45 = -t76 * t65 - t82;
t66 = t72 * pkin(9);
t37 = -t45 + t66;
t44 = t80 * t65 - t76 * t98;
t42 = pkin(3) + t44;
t91 = t73 * t75;
t15 = t79 * t37 + t42 * t91;
t58 = pkin(10) * t63;
t12 = t15 + t58;
t90 = t73 * t79;
t33 = t42 * t90;
t8 = t33 + t67 + (-t37 - t101) * t75;
t2 = -t74 * t12 + t78 * t8;
t28 = (-t42 - t102) * t72;
t96 = t2 * t73 + t28 * t40;
t68 = t80 * pkin(2);
t64 = t68 + pkin(3);
t43 = (-t64 - t102) * t72;
t50 = t64 * t90;
t54 = t66 + t99;
t20 = t50 + t67 + (-t54 - t101) * t75;
t31 = t79 * t54 + t64 * t91;
t23 = t31 + t58;
t5 = t78 * t20 - t74 * t23;
t95 = t43 * t40 + t5 * t73;
t60 = pkin(3) * t90;
t29 = t60 + t67 + (-pkin(9) - pkin(10)) * t62;
t47 = pkin(3) * t91 + pkin(9) * t63;
t36 = t47 + t58;
t10 = t78 * t29 - t74 * t36;
t48 = (-pkin(3) - t102) * t72;
t94 = t10 * t73 + t48 * t40;
t70 = t72 ^ 2;
t93 = t70 * t75;
t92 = t70 * t79;
t89 = t78 * t12;
t88 = t78 * t23;
t87 = t78 * t36;
t14 = -t75 * t37 + t33;
t86 = t14 * t73 + t42 * t92;
t30 = -t75 * t54 + t50;
t85 = t30 * t73 + t64 * t92;
t46 = -pkin(9) * t62 + t60;
t84 = pkin(3) * t92 + t46 * t73;
t83 = t72 * t103;
t71 = t73 ^ 2;
t61 = t70 * t75 ^ 2;
t59 = t73 * t97;
t53 = 0.2e1 * t75 * t92;
t52 = t79 * t83;
t51 = t75 * t83;
t41 = (t74 * t79 + t75 * t78) * t72;
t38 = t41 ^ 2;
t35 = t41 * t103;
t34 = t73 * t104;
t25 = t48 * t41;
t22 = t43 * t41;
t19 = t41 * t104;
t17 = t28 * t41;
t11 = t74 * t29 + t87;
t6 = t74 * t20 + t88;
t3 = t74 * t8 + t89;
t1 = [1, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t98, 1, 0.2e1 * t44, 0.2e1 * t45, t61, t53, t51, t52, t71, 0.2e1 * t86, -0.2e1 * t15 * t73 - 0.2e1 * t42 * t93, t38, t19, t35, t34, t71, 0.2e1 * t96, -0.2e1 * t3 * t73 + 0.2e1 * t17; 0, 0, 0, 1, t69, -t98, 1, t44 + t68, -t82 + (-pkin(2) - t65) * t76, t61, t53, t51, t52, t71, t85 + t86, (-t15 - t31) * t73 + (-t42 - t64) * t93, t38, t19, t35, t34, t71, t95 + t96, t17 + t22 + (-t3 - t6) * t73; 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t99, t61, t53, t51, t52, t71, 0.2e1 * t85, -0.2e1 * t31 * t73 - 0.2e1 * t64 * t93, t38, t19, t35, t34, t71, 0.2e1 * t95, -0.2e1 * t6 * t73 + 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 1, t44, t45, t61, t53, t51, t52, t71, t84 + t86, (-t15 - t47) * t73 + (-pkin(3) - t42) * t93, t38, t19, t35, t34, t71, t94 + t96, t17 + t25 + (-t11 - t3) * t73; 0, 0, 0, 0, 0, 0, 1, t68, -t99, t61, t53, t51, t52, t71, t84 + t85, (-t31 - t47) * t73 + (-pkin(3) - t64) * t93, t38, t19, t35, t34, t71, t94 + t95, t22 + t25 + (-t11 - t6) * t73; 0, 0, 0, 0, 0, 0, 1, 0, 0, t61, t53, t51, t52, t71, 0.2e1 * t84, -0.2e1 * pkin(3) * t93 - 0.2e1 * t47 * t73, t38, t19, t35, t34, t71, 0.2e1 * t94, -0.2e1 * t11 * t73 + 0.2e1 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t63, t73, t14, -t15, 0, 0, t41, -t40, t73, t2 + t59, -t89 + (-t8 - t67) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t63, t73, t30, -t31, 0, 0, t41, -t40, t73, t5 + t59, -t88 + (-t20 - t67) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t63, t73, t46, -t47, 0, 0, t41, -t40, t73, t10 + t59, -t87 + (-t29 - t67) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t97, -0.2e1 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, t73, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, t73, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, t73, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t97, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
