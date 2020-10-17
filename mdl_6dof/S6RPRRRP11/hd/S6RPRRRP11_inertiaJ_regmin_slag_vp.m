% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:08:30
% EndTime: 2019-05-06 02:08:33
% DurationCPUTime: 1.18s
% Computational Cost: add. (2459->161), mult. (6647->324), div. (0->0), fcn. (7731->12), ass. (0->95)
t63 = sin(pkin(6));
t62 = sin(pkin(7));
t78 = cos(pkin(6));
t73 = t78 * t62;
t64 = cos(pkin(12));
t65 = cos(pkin(7));
t91 = t64 * t65;
t108 = t63 * t91 + t73;
t61 = sin(pkin(12));
t74 = pkin(1) * t78;
t80 = qJ(2) * t63;
t41 = t61 * t74 + t64 * t80;
t25 = pkin(9) * t108 + t41;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t53 = t64 * t74;
t95 = t61 * t63;
t29 = t78 * pkin(2) + t53 + (-pkin(9) * t65 - qJ(2)) * t95;
t34 = (-pkin(9) * t61 * t62 - pkin(2) * t64 - pkin(1)) * t63;
t72 = t29 * t65 + t34 * t62;
t16 = -t25 * t68 + t71 * t72;
t28 = t68 * t73 + (t61 * t71 + t68 * t91) * t63;
t92 = t63 * t64;
t38 = t62 * t92 - t65 * t78;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t21 = t28 * t67 + t38 * t70;
t107 = -0.2e1 * t21;
t27 = -t108 * t71 + t68 * t95;
t106 = 0.2e1 * t27;
t105 = -0.2e1 * t28;
t104 = -0.2e1 * t67;
t103 = 0.2e1 * t70;
t57 = t63 ^ 2;
t102 = pkin(1) * t57;
t69 = cos(qJ(5));
t101 = pkin(4) * t69;
t66 = sin(qJ(5));
t100 = pkin(10) * t66;
t20 = -t29 * t62 + t34 * t65;
t12 = pkin(3) * t27 - pkin(10) * t28 + t20;
t17 = t25 * t71 + t68 * t72;
t15 = -t38 * pkin(10) + t17;
t8 = t12 * t70 - t15 * t67;
t6 = -pkin(4) * t27 - t8;
t99 = t6 * t66;
t98 = t6 * t69;
t97 = t66 * pkin(5);
t22 = t28 * t70 - t38 * t67;
t19 = t22 * t69 + t27 * t66;
t96 = t19 * t66;
t94 = t62 * t68;
t93 = t62 * t71;
t90 = t66 * t21;
t89 = t66 * t67;
t88 = t66 * t69;
t87 = t66 * t70;
t86 = t67 * t27;
t85 = t69 * t21;
t84 = t69 * t67;
t83 = t69 * t70;
t82 = t70 * t27;
t81 = -qJ(6) - pkin(11);
t79 = qJ(6) * t67;
t77 = t67 * t103;
t76 = pkin(10) * t83;
t14 = t38 * pkin(3) - t16;
t11 = t21 * pkin(4) - t22 * pkin(11) + t14;
t9 = t12 * t67 + t15 * t70;
t7 = pkin(11) * t27 + t9;
t3 = t11 * t69 - t66 * t7;
t4 = t11 * t66 + t69 * t7;
t60 = t69 ^ 2;
t59 = t67 ^ 2;
t58 = t66 ^ 2;
t56 = -pkin(5) * t69 - pkin(4);
t49 = t81 * t69;
t48 = t81 * t66;
t47 = -pkin(4) * t70 - pkin(11) * t67 - pkin(3);
t46 = (pkin(10) + t97) * t67;
t44 = t69 * t47;
t43 = t65 * t67 + t70 * t94;
t42 = -t65 * t70 + t67 * t94;
t40 = -t61 * t80 + t53;
t36 = t47 * t66 + t76;
t35 = -pkin(10) * t87 + t44;
t32 = t76 + (t47 - t79) * t66;
t31 = t43 * t69 - t66 * t93;
t30 = -t43 * t66 - t69 * t93;
t26 = -t69 * t79 + t44 + (-pkin(5) - t100) * t70;
t18 = t22 * t66 - t27 * t69;
t5 = pkin(5) * t18 + t6;
t2 = -qJ(6) * t18 + t4;
t1 = pkin(5) * t21 - qJ(6) * t19 + t3;
t10 = [1, 0, 0, 0.2e1 * t102 * t64 + 0.2e1 * t40 * t78, -0.2e1 * t102 * t61 - 0.2e1 * t41 * t78, 0.2e1 * (-t40 * t61 + t41 * t64) * t63, pkin(1) ^ 2 * t57 + t40 ^ 2 + t41 ^ 2, t28 ^ 2, t27 * t105, t38 * t105, t38 * t106, t38 ^ 2, -0.2e1 * t16 * t38 + 0.2e1 * t20 * t27, 0.2e1 * t17 * t38 + 0.2e1 * t20 * t28, t22 ^ 2, t22 * t107, t22 * t106, t27 * t107, t27 ^ 2, 0.2e1 * t14 * t21 + 0.2e1 * t27 * t8, 0.2e1 * t14 * t22 - 0.2e1 * t27 * t9, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t21, t18 * t107, t21 ^ 2, 0.2e1 * t18 * t6 + 0.2e1 * t21 * t3, 0.2e1 * t19 * t6 - 0.2e1 * t21 * t4, -0.2e1 * t1 * t19 - 0.2e1 * t18 * t2, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t92, t95, 0, -t63 * pkin(1), 0, 0, 0, 0, 0, t27 * t65 - t38 * t93, t28 * t65 + t38 * t94, 0, 0, 0, 0, 0, -t21 * t93 - t27 * t42, -t22 * t93 - t27 * t43, 0, 0, 0, 0, 0, t18 * t42 + t21 * t30, t19 * t42 - t21 * t31, -t18 * t31 - t19 * t30, t1 * t30 + t2 * t31 + t42 * t5; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 ^ 2 + t31 ^ 2 + t42 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, -t38, t16, -t17, t22 * t67, -t21 * t67 + t22 * t70, t86, t82, 0, -pkin(3) * t21 - pkin(10) * t86 - t14 * t70, -pkin(3) * t22 - pkin(10) * t82 + t14 * t67, t19 * t84 (-t18 * t69 - t96) * t67, -t19 * t70 + t21 * t84, t18 * t70 - t21 * t89, -t21 * t70, t21 * t35 - t3 * t70 + (pkin(10) * t18 + t99) * t67, -t21 * t36 + t4 * t70 + (pkin(10) * t19 + t98) * t67, -t18 * t32 - t19 * t26 + (-t1 * t69 - t2 * t66) * t67, t1 * t26 + t2 * t32 + t46 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t94, 0, 0, 0, 0, 0, t70 * t93, -t67 * t93, 0, 0, 0, 0, 0, -t30 * t70 + t42 * t89, t31 * t70 + t42 * t84 (-t30 * t69 - t31 * t66) * t67, t26 * t30 + t31 * t32 + t42 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, t77, 0, 0, 0, pkin(3) * t103, pkin(3) * t104, t60 * t59, -0.2e1 * t59 * t88, t83 * t104, t66 * t77, t70 ^ 2, 0.2e1 * t100 * t59 - 0.2e1 * t35 * t70, 0.2e1 * pkin(10) * t59 * t69 + 0.2e1 * t36 * t70, 0.2e1 * (-t26 * t69 - t32 * t66) * t67, t26 ^ 2 + t32 ^ 2 + t46 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t27, t8, -t9, t96, -t18 * t66 + t19 * t69, t90, t85, 0, -pkin(4) * t18 - pkin(11) * t90 - t98, -pkin(4) * t19 - pkin(11) * t85 + t99, -t1 * t66 + t18 * t49 - t19 * t48 + t2 * t69, t1 * t48 - t2 * t49 + t5 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t43, 0, 0, 0, 0, 0, -t42 * t69, t42 * t66, -t30 * t66 + t31 * t69, t30 * t48 - t31 * t49 + t42 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t67 * pkin(10), -t70 * pkin(10), t66 * t84 (-t58 + t60) * t67, -t87, -t83, 0, -pkin(10) * t84 + (-pkin(4) * t67 + pkin(11) * t70) * t66, pkin(11) * t83 + (t100 - t101) * t67 (-t48 * t67 + t32) * t69 + (t49 * t67 - t26) * t66, t26 * t48 - t32 * t49 + t46 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, 0.2e1 * t88, 0, 0, 0, 0.2e1 * t101, -0.2e1 * pkin(4) * t66, -0.2e1 * t48 * t66 - 0.2e1 * t49 * t69, t48 ^ 2 + t49 ^ 2 + t56 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t21, t3, -t4, -pkin(5) * t19, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, 0, t30 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t89, -t70, t35, -t36, -pkin(5) * t84, t26 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t69, 0, -t66 * pkin(11), -t69 * pkin(11), -t97, t48 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
