% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:05:21
% EndTime: 2019-05-05 08:05:25
% DurationCPUTime: 0.89s
% Computational Cost: add. (1237->140), mult. (3101->297), div. (0->0), fcn. (3749->14), ass. (0->92)
t65 = cos(pkin(7));
t68 = sin(qJ(4));
t72 = cos(qJ(4));
t62 = sin(pkin(7));
t69 = sin(qJ(3));
t96 = t62 * t69;
t42 = -t72 * t65 + t68 * t96;
t43 = t68 * t65 + t72 * t96;
t61 = sin(pkin(13));
t64 = cos(pkin(13));
t27 = -t61 * t42 + t64 * t43;
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t73 = cos(qJ(3));
t95 = t62 * t73;
t21 = t67 * t27 + t71 * t95;
t103 = -0.2e1 * t21;
t102 = -0.2e1 * t43;
t101 = 0.2e1 * t72;
t100 = pkin(2) * t69;
t99 = pkin(2) * t73;
t22 = t71 * t27 - t67 * t95;
t98 = t22 * t67;
t58 = t62 ^ 2;
t97 = t58 * t73;
t63 = sin(pkin(6));
t70 = sin(qJ(2));
t94 = t63 * t70;
t74 = cos(qJ(2));
t93 = t63 * t74;
t92 = t65 * t69;
t91 = t65 * t74;
t26 = t64 * t42 + t61 * t43;
t90 = t67 * t26;
t47 = t61 * t68 - t64 * t72;
t89 = t67 * t47;
t48 = t61 * t72 + t64 * t68;
t88 = t67 * t48;
t55 = t61 * pkin(4) + pkin(11);
t87 = t67 * t55;
t86 = t67 * t71;
t85 = t71 * t48;
t84 = t71 * t55;
t83 = -qJ(5) - pkin(10);
t81 = pkin(9) * t95;
t37 = t81 + (pkin(10) + t100) * t65;
t38 = (-pkin(3) * t73 - pkin(10) * t69 - pkin(2)) * t62;
t23 = -t68 * t37 + t72 * t38;
t14 = -pkin(4) * t95 - t43 * qJ(5) + t23;
t24 = t72 * t37 + t68 * t38;
t18 = -t42 * qJ(5) + t24;
t8 = t61 * t14 + t64 * t18;
t82 = 0.2e1 * t95;
t80 = t68 * t95;
t79 = t72 * t95;
t57 = -t72 * pkin(4) - pkin(3);
t78 = t83 * t68;
t7 = t64 * t14 - t61 * t18;
t66 = cos(pkin(6));
t31 = t66 * t96 + (t69 * t91 + t70 * t73) * t63;
t40 = -t62 * t93 + t66 * t65;
t77 = t31 * t68 - t40 * t72;
t56 = -t64 * pkin(4) - pkin(5);
t76 = -t47 * t55 + t48 * t56;
t53 = pkin(9) * t96;
t36 = t53 + (-pkin(3) - t99) * t65;
t28 = t42 * pkin(4) + t36;
t60 = t71 ^ 2;
t59 = t67 ^ 2;
t51 = t83 * t72;
t46 = t48 ^ 2;
t45 = pkin(2) * t92 + t81;
t44 = t65 * t99 - t53;
t41 = t71 * t47;
t34 = -t64 * t51 + t61 * t78;
t32 = -t61 * t51 - t64 * t78;
t30 = -t63 * t73 * t91 - t66 * t95 + t69 * t94;
t29 = t47 * pkin(5) - t48 * pkin(11) + t57;
t25 = t71 * t26;
t20 = t31 * t72 + t40 * t68;
t16 = t67 * t29 + t71 * t34;
t15 = t71 * t29 - t67 * t34;
t12 = t64 * t20 - t61 * t77;
t10 = t61 * t20 + t64 * t77;
t9 = t26 * pkin(5) - t27 * pkin(11) + t28;
t6 = -pkin(11) * t95 + t8;
t5 = pkin(5) * t95 - t7;
t4 = t71 * t12 + t30 * t67;
t3 = -t67 * t12 + t30 * t71;
t2 = t71 * t6 + t67 * t9;
t1 = -t67 * t6 + t71 * t9;
t11 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t12 ^ 2 + t30 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t93, -t94, 0, 0, 0, 0, 0, -t30 * t65 - t40 * t95, -t31 * t65 + t40 * t96, 0, 0, 0, 0, 0, t30 * t42 + t77 * t95, t20 * t95 + t30 * t43, t10 * t27 - t12 * t26, -t10 * t7 + t12 * t8 + t30 * t28, 0, 0, 0, 0, 0, t10 * t21 + t3 * t26, t10 * t22 - t4 * t26; 0, 1, 0, 0, t58 * t69 ^ 2, 0.2e1 * t69 * t97, 0.2e1 * t62 * t92, t65 * t82, t65 ^ 2, 0.2e1 * pkin(2) * t97 + 0.2e1 * t44 * t65, -0.2e1 * t58 * t100 - 0.2e1 * t45 * t65, t43 ^ 2, t42 * t102, t95 * t102, t42 * t82, t58 * t73 ^ 2, -0.2e1 * t23 * t95 + 0.2e1 * t36 * t42, 0.2e1 * t24 * t95 + 0.2e1 * t36 * t43, -0.2e1 * t8 * t26 - 0.2e1 * t7 * t27, t28 ^ 2 + t7 ^ 2 + t8 ^ 2, t22 ^ 2, t22 * t103, 0.2e1 * t22 * t26, t26 * t103, t26 ^ 2, 0.2e1 * t1 * t26 + 0.2e1 * t5 * t21, -0.2e1 * t2 * t26 + 0.2e1 * t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t31, 0, 0, 0, 0, 0, -t30 * t72, t30 * t68, t10 * t48 - t12 * t47, t10 * t32 + t12 * t34 + t30 * t57, 0, 0, 0, 0, 0, t10 * t88 + t3 * t47, t10 * t85 - t4 * t47; 0, 0, 0, 0, 0, 0, t96, t95, t65, t44, -t45, t43 * t68, -t68 * t42 + t43 * t72, -t80, -t79, 0, -pkin(3) * t42 + pkin(10) * t80 - t36 * t72, -pkin(3) * t43 + pkin(10) * t79 + t36 * t68, -t34 * t26 + t32 * t27 - t8 * t47 - t7 * t48, t28 * t57 - t7 * t32 + t8 * t34, t22 * t85 (-t21 * t71 - t98) * t48, t22 * t47 + t26 * t85, -t21 * t47 - t26 * t88, t26 * t47, t1 * t47 + t15 * t26 + t32 * t21 + t5 * t88, -t16 * t26 - t2 * t47 + t32 * t22 + t5 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t68 ^ 2, t68 * t101, 0, 0, 0, pkin(3) * t101, -0.2e1 * pkin(3) * t68, 0.2e1 * t32 * t48 - 0.2e1 * t34 * t47, t32 ^ 2 + t34 ^ 2 + t57 ^ 2, t60 * t46, -0.2e1 * t46 * t86, 0.2e1 * t47 * t85, -0.2e1 * t47 * t88, t47 ^ 2, 0.2e1 * t15 * t47 + 0.2e1 * t32 * t88, -0.2e1 * t16 * t47 + 0.2e1 * t32 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t20, 0 (-t10 * t64 + t12 * t61) * pkin(4), 0, 0, 0, 0, 0, -t10 * t71, t10 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42, -t95, t23, -t24 (-t26 * t61 - t27 * t64) * pkin(4) (t61 * t8 + t64 * t7) * pkin(4), t98, -t67 * t21 + t22 * t71, t90, t25, 0, t56 * t21 - t26 * t87 - t5 * t71, t56 * t22 - t26 * t84 + t5 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t72, 0, -t68 * pkin(10), -t72 * pkin(10) (-t47 * t61 - t48 * t64) * pkin(4) (-t32 * t64 + t34 * t61) * pkin(4), t67 * t85 (-t59 + t60) * t48, t89, t41, 0, -t32 * t71 + t76 * t67, t32 * t67 + t76 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t61 ^ 2 + t64 ^ 2) * pkin(4) ^ 2, t59, 0.2e1 * t86, 0, 0, 0, -0.2e1 * t56 * t71, 0.2e1 * t56 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, t25, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, t41, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t88, t47, t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t71, 0, -t87, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
