% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:25:02
% EndTime: 2019-05-06 17:25:05
% DurationCPUTime: 0.85s
% Computational Cost: add. (1403->98), mult. (2608->181), div. (0->0), fcn. (3146->8), ass. (0->80)
t54 = cos(pkin(10));
t47 = t54 * pkin(2) + pkin(3);
t56 = sin(qJ(4));
t84 = cos(qJ(4));
t53 = sin(pkin(10));
t90 = t53 * pkin(2);
t35 = -t56 * t47 - t84 * t90;
t33 = pkin(9) - t35;
t55 = sin(qJ(5));
t51 = t55 ^ 2;
t58 = cos(qJ(5));
t52 = t58 ^ 2;
t73 = t51 + t52;
t75 = t73 * t33;
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t37 = t53 * t59 + t54 * t57;
t62 = t53 * t57 - t54 * t59;
t24 = t84 * t37 - t56 * t62;
t97 = -0.2e1 * t24;
t48 = -t59 * pkin(2) - pkin(1);
t30 = t62 * pkin(3) + t48;
t96 = 0.2e1 * t30;
t95 = -0.2e1 * t55;
t94 = -0.2e1 * t58;
t93 = 0.2e1 * t59;
t23 = t56 * t37 + t84 * t62;
t11 = t23 * pkin(4) - t24 * pkin(9) + t30;
t77 = -qJ(3) - pkin(7);
t42 = t77 * t57;
t44 = t77 * t59;
t26 = t54 * t42 + t53 * t44;
t16 = -t37 * pkin(8) + t26;
t27 = t53 * t42 - t54 * t44;
t17 = -t62 * pkin(8) + t27;
t13 = t56 * t16 + t84 * t17;
t5 = t55 * t11 + t58 * t13;
t92 = pkin(9) * t23;
t91 = t23 * pkin(5);
t89 = t55 * pkin(9);
t88 = t58 * pkin(9);
t12 = -t84 * t16 + t56 * t17;
t65 = pkin(5) * t55 - t58 * qJ(6);
t6 = t65 * t24 + t12;
t87 = t6 * t55;
t86 = t6 * t58;
t34 = t84 * t47 - t56 * t90;
t32 = -pkin(4) - t34;
t85 = pkin(4) - t32;
t83 = t12 * t58;
t82 = t23 * t33;
t19 = t55 * t23;
t81 = t55 * t24;
t80 = t55 * t33;
t79 = t55 * t58;
t20 = t58 * t23;
t21 = t58 * t24;
t78 = t58 * t33;
t66 = t58 * pkin(5) + t55 * qJ(6);
t41 = -pkin(4) - t66;
t25 = t41 - t34;
t76 = -t25 - t41;
t74 = t73 * pkin(9);
t72 = t23 * qJ(6);
t71 = t23 * t97;
t70 = -t58 * t11 + t55 * t13;
t69 = -pkin(4) * t24 - t92;
t2 = t72 + t5;
t3 = t70 - t91;
t1 = t2 * t58 + t3 * t55;
t68 = t2 * t55 - t3 * t58;
t67 = -t24 * t41 + t92;
t64 = -t24 * t25 + t82;
t63 = t24 * t32 - t82;
t46 = 0.2e1 * t79;
t22 = t24 ^ 2;
t18 = t55 * t21;
t14 = (-t51 + t52) * t24;
t10 = t12 * t55;
t4 = [1, 0, 0, t57 ^ 2, t57 * t93, 0, 0, 0, pkin(1) * t93, -0.2e1 * pkin(1) * t57, -0.2e1 * t26 * t37 - 0.2e1 * t27 * t62, t26 ^ 2 + t27 ^ 2 + t48 ^ 2, t22, t71, 0, 0, 0, t23 * t96, t24 * t96, t52 * t22, -0.2e1 * t22 * t79, 0.2e1 * t23 * t21, t55 * t71, t23 ^ 2, 0.2e1 * t12 * t81 - 0.2e1 * t23 * t70, 0.2e1 * t12 * t21 - 0.2e1 * t5 * t23, -0.2e1 * t3 * t23 + 0.2e1 * t6 * t81, t68 * t97, 0.2e1 * t2 * t23 - 0.2e1 * t6 * t21, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t57, t59, 0, -t57 * pkin(7), -t59 * pkin(7) (-t54 * t37 - t53 * t62) * pkin(2) (t26 * t54 + t27 * t53) * pkin(2), 0, 0, t24, -t23, 0, -t12, -t13, t18, t14, t19, t20, 0, t63 * t55 - t83, t63 * t58 + t10, -t55 * t64 - t86, t1, t58 * t64 - t87, t1 * t33 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t53 ^ 2 + t54 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t34, 0.2e1 * t35, t51, t46, 0, 0, 0, t32 * t94, 0.2e1 * t32 * t55, t25 * t94, 0.2e1 * t75, t25 * t95, t73 * t33 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t20, -t19, t20, -t73 * t24, t19, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t12, -t13, t18, t14, t19, t20, 0, t69 * t55 - t83, t58 * t69 + t10, -t55 * t67 - t86, t1, t58 * t67 - t87, pkin(9) * t1 + t6 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, t35, t51, t46, 0, 0, 0, t85 * t58, -t85 * t55, t76 * t58, t74 + t75, t76 * t55, pkin(9) * t75 + t25 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t51, t46, 0, 0, 0, 0.2e1 * pkin(4) * t58, pkin(4) * t95, t41 * t94, 0.2e1 * t74, t41 * t95, pkin(9) ^ 2 * t73 + t41 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t81, t23, -t70, -t5, -t70 + 0.2e1 * t91, -t66 * t24, 0.2e1 * t72 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t80, -t78, -t80, -t65, t78, -t65 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55, t58, 0, t55, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t89, -t88, -t89, -t65, t88, -t65 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
