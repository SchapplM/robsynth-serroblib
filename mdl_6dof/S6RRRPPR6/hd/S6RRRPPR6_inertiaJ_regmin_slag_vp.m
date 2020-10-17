% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:33:31
% EndTime: 2019-05-07 05:33:34
% DurationCPUTime: 0.96s
% Computational Cost: add. (1514->146), mult. (3431->283), div. (0->0), fcn. (3946->10), ass. (0->93)
t64 = cos(pkin(6));
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t62 = sin(pkin(6));
t67 = sin(qJ(2));
t94 = t62 * t67;
t40 = -t64 * t69 + t66 * t94;
t41 = t64 * t66 + t69 * t94;
t61 = sin(pkin(11));
t63 = cos(pkin(11));
t27 = -t40 * t61 + t41 * t63;
t105 = 0.2e1 * t27;
t45 = t61 * t66 - t63 * t69;
t46 = t61 * t69 + t63 * t66;
t57 = -t69 * pkin(3) - pkin(2);
t74 = -t46 * qJ(5) + t57;
t29 = pkin(4) * t45 + t74;
t104 = -0.2e1 * t29;
t103 = -0.2e1 * t41;
t52 = pkin(3) * t61 + qJ(5);
t102 = 0.2e1 * t52;
t101 = 0.2e1 * t69;
t100 = pkin(4) + pkin(10);
t99 = pkin(1) * t67;
t70 = cos(qJ(2));
t98 = pkin(1) * t70;
t26 = t40 * t63 + t41 * t61;
t65 = sin(qJ(6));
t68 = cos(qJ(6));
t93 = t62 * t70;
t19 = t26 * t65 - t68 * t93;
t97 = t19 * t68;
t96 = t52 * t45;
t58 = t62 ^ 2;
t95 = t58 * t70;
t92 = t64 * t67;
t91 = t65 * t27;
t90 = t65 * t45;
t89 = t65 * t46;
t56 = -pkin(3) * t63 - pkin(4);
t51 = -pkin(10) + t56;
t88 = t65 * t51;
t23 = t68 * t27;
t87 = t68 * t45;
t39 = t68 * t46;
t86 = t68 * t51;
t85 = t68 * t65;
t84 = -qJ(4) - pkin(9);
t81 = pkin(8) * t93;
t36 = t81 + (pkin(9) + t99) * t64;
t37 = (-pkin(2) * t70 - pkin(9) * t67 - pkin(1)) * t62;
t21 = -t36 * t66 + t37 * t69;
t14 = -pkin(3) * t93 - qJ(4) * t41 + t21;
t22 = t36 * t69 + t37 * t66;
t17 = -qJ(4) * t40 + t22;
t9 = t14 * t61 + t17 * t63;
t83 = 0.2e1 * t45 * t46;
t82 = 0.2e1 * t93;
t80 = t66 * t93;
t79 = t69 * t93;
t48 = t84 * t69;
t77 = t84 * t66;
t30 = -t48 * t61 - t63 * t77;
t32 = -t63 * t48 + t61 * t77;
t78 = t30 ^ 2 + t32 ^ 2;
t8 = t14 * t63 - t17 * t61;
t7 = pkin(4) * t93 - t8;
t76 = -t26 * t32 + t27 * t30;
t75 = -t46 * t51 + t96;
t6 = qJ(5) * t93 - t9;
t49 = pkin(8) * t94;
t35 = t49 + (-pkin(2) - t98) * t64;
t73 = 0.2e1 * t30 * t46 - 0.2e1 * t32 * t45;
t28 = t40 * pkin(3) + t35;
t72 = -t27 * qJ(5) + t28;
t60 = t68 ^ 2;
t59 = t65 ^ 2;
t44 = t45 ^ 2;
t43 = pkin(1) * t92 + t81;
t42 = t64 * t98 - t49;
t25 = -t45 * pkin(5) + t32;
t24 = pkin(5) * t46 + t30;
t20 = t100 * t45 + t74;
t18 = t26 * t68 + t65 * t93;
t12 = t20 * t68 + t24 * t65;
t11 = -t20 * t65 + t24 * t68;
t10 = t26 * pkin(4) + t72;
t5 = t100 * t26 + t72;
t4 = -pkin(5) * t26 - t6;
t3 = pkin(5) * t27 + pkin(10) * t93 + t7;
t2 = t3 * t65 + t5 * t68;
t1 = t3 * t68 - t5 * t65;
t13 = [1, 0, 0, t58 * t67 ^ 2, 0.2e1 * t67 * t95, 0.2e1 * t62 * t92, t64 * t82, t64 ^ 2, 0.2e1 * pkin(1) * t95 + 0.2e1 * t42 * t64, -0.2e1 * t43 * t64 - 0.2e1 * t58 * t99, t41 ^ 2, t40 * t103, t93 * t103, t40 * t82, t58 * t70 ^ 2, -0.2e1 * t21 * t93 + 0.2e1 * t35 * t40, 0.2e1 * t22 * t93 + 0.2e1 * t35 * t41, -0.2e1 * t26 * t9 - 0.2e1 * t27 * t8, t28 ^ 2 + t8 ^ 2 + t9 ^ 2, 0.2e1 * t26 * t6 + 0.2e1 * t27 * t7, -0.2e1 * t10 * t26 - 0.2e1 * t7 * t93, -0.2e1 * t10 * t27 + 0.2e1 * t6 * t93, t10 ^ 2 + t6 ^ 2 + t7 ^ 2, t19 ^ 2, 0.2e1 * t19 * t18, t19 * t105, t18 * t105, t27 ^ 2, 0.2e1 * t1 * t27 - 0.2e1 * t18 * t4, 0.2e1 * t19 * t4 - 0.2e1 * t2 * t27; 0, 0, 0, 0, 0, t94, t93, t64, t42, -t43, t41 * t66, -t40 * t66 + t41 * t69, -t80, -t79, 0, -pkin(2) * t40 + pkin(9) * t80 - t35 * t69, -pkin(2) * t41 + pkin(9) * t79 + t35 * t66, -t45 * t9 - t46 * t8 + t76, t28 * t57 - t30 * t8 + t32 * t9, t45 * t6 + t46 * t7 + t76, -t10 * t45 - t26 * t29 - t30 * t93, -t10 * t46 - t27 * t29 - t32 * t93, t10 * t29 + t30 * t7 - t32 * t6, t19 * t90 (t18 * t65 + t97) * t45, t19 * t46 + t27 * t90, t18 * t46 + t27 * t87, t27 * t46, t1 * t46 + t11 * t27 - t18 * t25 - t4 * t87, -t12 * t27 + t19 * t25 - t2 * t46 + t4 * t90; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t66 ^ 2, t66 * t101, 0, 0, 0, pkin(2) * t101, -0.2e1 * pkin(2) * t66, t73, t57 ^ 2 + t78, t73, t45 * t104, t46 * t104, t29 ^ 2 + t78, t59 * t44, 0.2e1 * t44 * t85, t65 * t83, t68 * t83, t46 ^ 2, 0.2e1 * t11 * t46 - 0.2e1 * t25 * t87, -0.2e1 * t12 * t46 + 0.2e1 * t25 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, -t93, t21, -t22 (-t26 * t61 - t27 * t63) * pkin(3) (t61 * t9 + t63 * t8) * pkin(3), -t26 * t52 + t27 * t56, -t56 * t93 + t7 (-qJ(5) - t52) * t93 + t9, -t52 * t6 + t56 * t7, t97, t18 * t68 - t19 * t65, t23, -t91, 0, -t18 * t52 + t27 * t86 + t4 * t65, t19 * t52 - t27 * t88 + t4 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t69, 0, -t66 * pkin(9), -t69 * pkin(9) (-t45 * t61 - t46 * t63) * pkin(3) (-t30 * t63 + t32 * t61) * pkin(3), t46 * t56 - t96, t30, t32, t30 * t56 + t32 * t52, t45 * t85 (-t59 + t60) * t45, t39, -t89, 0, t25 * t65 - t68 * t75, t25 * t68 + t65 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t61 ^ 2 + t63 ^ 2) * pkin(3) ^ 2, 0, 0.2e1 * t56, t102, t52 ^ 2 + t56 ^ 2, t60, -0.2e1 * t85, 0, 0, 0, t65 * t102, t68 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, -t27, t10, 0, 0, 0, 0, 0, -t91, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t45, -t46, t29, 0, 0, 0, 0, 0, -t89, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t93, 0, t7, 0, 0, 0, 0, 0, t23, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, t30, 0, 0, 0, 0, 0, t39, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18, t27, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t87, t46, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t65, 0, t86, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
