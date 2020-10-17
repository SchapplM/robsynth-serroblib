% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:09:17
% EndTime: 2019-05-07 04:09:20
% DurationCPUTime: 0.85s
% Computational Cost: add. (1497->111), mult. (2818->203), div. (0->0), fcn. (3427->10), ass. (0->88)
t73 = sin(pkin(10));
t63 = pkin(3) * t73 + qJ(5);
t72 = sin(pkin(11));
t74 = cos(pkin(11));
t88 = t72 ^ 2 + t74 ^ 2;
t89 = t88 * t63;
t78 = sin(qJ(2));
t98 = pkin(7) + pkin(8);
t58 = t98 * t78;
t81 = cos(qJ(2));
t59 = t98 * t81;
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t37 = t77 * t58 - t80 * t59;
t84 = t77 * t78 - t80 * t81;
t28 = -qJ(4) * t84 - t37;
t75 = cos(pkin(10));
t36 = -t58 * t80 - t59 * t77;
t53 = t77 * t81 + t78 * t80;
t83 = -qJ(4) * t53 + t36;
t18 = t28 * t73 - t75 * t83;
t106 = t18 ^ 2;
t32 = t53 * t75 - t73 * t84;
t76 = sin(qJ(6));
t79 = cos(qJ(6));
t52 = t72 * t79 + t74 * t76;
t21 = t52 * t32;
t105 = -0.2e1 * t21;
t69 = t80 * pkin(2);
t66 = t69 + pkin(3);
t96 = t77 * pkin(2);
t45 = t66 * t75 - t73 * t96;
t44 = -pkin(4) - t45;
t97 = t74 * pkin(5);
t38 = t44 - t97;
t104 = 0.2e1 * t38;
t65 = -pkin(3) * t75 - pkin(4);
t54 = t65 - t97;
t103 = 0.2e1 * t54;
t67 = -t81 * pkin(2) - pkin(1);
t102 = 0.2e1 * t67;
t101 = 0.2e1 * t72;
t100 = -0.2e1 * t74;
t99 = 0.2e1 * t81;
t95 = t18 * t74;
t31 = t53 * t73 + t75 * t84;
t26 = t52 * t31;
t94 = t72 * t32;
t93 = t74 * t32;
t41 = pkin(3) * t84 + t67;
t17 = t31 * pkin(4) - t32 * qJ(5) + t41;
t20 = t28 * t75 + t73 * t83;
t8 = t17 * t72 + t20 * t74;
t92 = t38 + t54;
t46 = t66 * t73 + t75 * t96;
t43 = qJ(5) + t46;
t91 = t88 * t43;
t90 = t44 + t65;
t7 = t17 * t74 - t20 * t72;
t87 = t7 * t74 + t72 * t8;
t3 = -t7 * t72 + t74 * t8;
t86 = -t31 * t43 + t32 * t44;
t85 = -t31 * t63 + t32 * t65;
t51 = t72 * t76 - t74 * t79;
t68 = t74 * pkin(9);
t50 = t52 ^ 2;
t49 = t63 * t74 + t68;
t48 = (-pkin(9) - t63) * t72;
t35 = t43 * t74 + t68;
t34 = (-pkin(9) - t43) * t72;
t33 = -0.2e1 * t52 * t51;
t30 = t48 * t76 + t49 * t79;
t29 = t48 * t79 - t49 * t76;
t25 = t51 * t31;
t24 = t34 * t76 + t35 * t79;
t23 = t34 * t79 - t35 * t76;
t22 = t51 * t32;
t16 = t18 * t72;
t12 = t22 * t52;
t11 = pkin(5) * t94 + t18;
t10 = t11 * t52;
t9 = t11 * t51;
t6 = -t21 * t52 + t22 * t51;
t5 = -pkin(9) * t94 + t8;
t4 = pkin(5) * t31 - pkin(9) * t93 + t7;
t2 = t4 * t76 + t5 * t79;
t1 = t4 * t79 - t5 * t76;
t13 = [1, 0, 0, t78 ^ 2, t78 * t99, 0, 0, 0, pkin(1) * t99, -0.2e1 * pkin(1) * t78, t53 ^ 2, -0.2e1 * t53 * t84, 0, 0, 0, t84 * t102, t53 * t102, 0.2e1 * t18 * t32 - 0.2e1 * t20 * t31, t20 ^ 2 + t41 ^ 2 + t106, 0.2e1 * t18 * t94 + 0.2e1 * t31 * t7, 0.2e1 * t18 * t93 - 0.2e1 * t31 * t8, -0.2e1 * t87 * t32, t7 ^ 2 + t8 ^ 2 + t106, t22 ^ 2, -t22 * t105, -0.2e1 * t22 * t31, t31 * t105, t31 ^ 2, 0.2e1 * t1 * t31 + 0.2e1 * t11 * t21, -0.2e1 * t11 * t22 - 0.2e1 * t2 * t31; 0, 0, 0, 0, 0, t78, t81, 0, -t78 * pkin(7), -t81 * pkin(7), 0, 0, t53, -t84, 0, t36, t37, -t31 * t46 - t32 * t45, -t18 * t45 + t20 * t46, t72 * t86 - t95, t74 * t86 + t16, t3, t18 * t44 + t3 * t43, -t12, t6, t26, -t25, 0, t21 * t38 + t23 * t31 + t9, -t22 * t38 - t24 * t31 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t96, 0, t45 ^ 2 + t46 ^ 2, t44 * t100, t44 * t101, 0.2e1 * t91, t43 ^ 2 * t88 + t44 ^ 2, t50, t33, 0, 0, 0, t51 * t104, t52 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t84, 0, t36, t37 (-t31 * t73 - t32 * t75) * pkin(3) (-t18 * t75 + t20 * t73) * pkin(3), t72 * t85 - t95, t74 * t85 + t16, t3, t18 * t65 + t3 * t63, -t12, t6, t26, -t25, 0, t21 * t54 + t29 * t31 + t9, -t22 * t54 - t30 * t31 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t96, 0 (t45 * t75 + t46 * t73) * pkin(3), -t90 * t74, t90 * t72, t89 + t91, t43 * t89 + t44 * t65, t50, t33, 0, 0, 0, t92 * t51, t92 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t73 ^ 2 + t75 ^ 2) * pkin(3) ^ 2, t65 * t100, t65 * t101, 0.2e1 * t89, t63 ^ 2 * t88 + t65 ^ 2, t50, t33, 0, 0, 0, t51 * t103, t52 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t74 * t31, -t72 * t31, -t88 * t32, t87, 0, 0, 0, 0, 0, -t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t93, 0, t18, 0, 0, 0, 0, 0, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t72, 0, t44, 0, 0, 0, 0, 0, t51, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t72, 0, t65, 0, 0, 0, 0, 0, t51, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21, t31, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t51, 0, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t51, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
