% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:21:15
% EndTime: 2019-05-06 11:21:18
% DurationCPUTime: 0.80s
% Computational Cost: add. (675->118), mult. (1410->208), div. (0->0), fcn. (1579->8), ass. (0->75)
t64 = sin(pkin(10));
t65 = cos(pkin(10));
t93 = t64 ^ 2 + t65 ^ 2;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t37 = t70 * t64 - t67 * t65;
t41 = -t65 * pkin(3) - t64 * qJ(4) - pkin(2);
t33 = t65 * pkin(4) - t41;
t36 = t67 * t64 + t70 * t65;
t20 = t36 * pkin(5) + t33;
t92 = 0.2e1 * t20;
t91 = 0.2e1 * t33;
t90 = -0.2e1 * t41;
t68 = sin(qJ(2));
t89 = 0.2e1 * t68;
t71 = cos(qJ(2));
t88 = -0.2e1 * t71;
t87 = 0.2e1 * t71;
t86 = pkin(2) * t64;
t85 = pkin(7) * t65;
t66 = sin(qJ(6));
t84 = t66 * pkin(5);
t58 = t68 * pkin(7);
t69 = cos(qJ(6));
t83 = t69 * pkin(5);
t30 = t37 * t68;
t42 = -t71 * pkin(2) - t68 * qJ(3) - pkin(1);
t80 = t71 * pkin(7);
t48 = t64 * t80;
t59 = t71 * pkin(3);
t16 = t71 * pkin(4) + t48 + t59 + (-pkin(8) * t68 - t42) * t65;
t27 = t64 * t42 + t65 * t80;
t24 = -t71 * qJ(4) + t27;
t50 = t64 * t68;
t19 = pkin(8) * t50 + t24;
t9 = t67 * t16 + t70 * t19;
t5 = t30 * pkin(9) + t9;
t82 = t69 * t5;
t81 = t71 * pkin(5);
t51 = t65 * t68;
t77 = qJ(4) * t51 - t58;
t76 = t93 * qJ(3) ^ 2;
t75 = qJ(3) * t71;
t54 = t64 * qJ(3);
t31 = t36 * t68;
t8 = t70 * t16 - t67 * t19;
t4 = -t31 * pkin(9) + t8 + t81;
t1 = t69 * t4 - t66 * t5;
t43 = -t64 * pkin(8) + t54;
t44 = (-pkin(8) + qJ(3)) * t65;
t21 = t70 * t43 - t67 * t44;
t26 = t65 * t42 - t48;
t25 = -t26 + t59;
t74 = t24 * t65 + t25 * t64;
t73 = -t26 * t64 + t27 * t65;
t22 = t67 * t43 + t70 * t44;
t23 = (-pkin(3) - pkin(4)) * t50 + t77;
t63 = t71 ^ 2;
t62 = t68 ^ 2;
t46 = t71 * t54;
t40 = t66 * t70 + t69 * t67;
t39 = -t66 * t67 + t69 * t70;
t38 = 0.2e1 * t93 * qJ(3);
t29 = pkin(3) * t50 - t77;
t18 = -t66 * t36 + t69 * t37;
t17 = t69 * t36 + t66 * t37;
t14 = -t36 * pkin(9) + t22;
t13 = -t37 * pkin(9) + t21;
t12 = t66 * t30 + t69 * t31;
t11 = -t69 * t30 + t66 * t31;
t10 = -t30 * pkin(5) + t23;
t7 = t66 * t13 + t69 * t14;
t6 = t69 * t13 - t66 * t14;
t2 = t66 * t4 + t82;
t3 = [1, 0, 0, t62, t68 * t87, 0, 0, 0, pkin(1) * t87, -0.2e1 * pkin(1) * t68, 0.2e1 * t62 * pkin(7) * t64 - 0.2e1 * t26 * t71, 0.2e1 * t27 * t71 + 0.2e1 * t62 * t85 (-t26 * t65 - t27 * t64) * t89, t62 * pkin(7) ^ 2 + t26 ^ 2 + t27 ^ 2, 0.2e1 * t25 * t71 + 0.2e1 * t29 * t50 (-t24 * t64 + t25 * t65) * t89, -0.2e1 * t24 * t71 - 0.2e1 * t29 * t51, t24 ^ 2 + t25 ^ 2 + t29 ^ 2, t31 ^ 2, 0.2e1 * t31 * t30, t31 * t87, -t30 * t88, t63, -0.2e1 * t23 * t30 + 0.2e1 * t8 * t71, 0.2e1 * t23 * t31 - 0.2e1 * t9 * t71, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t87, t11 * t88, t63, 0.2e1 * t1 * t71 + 0.2e1 * t10 * t11, 0.2e1 * t10 * t12 - 0.2e1 * t2 * t71; 0, 0, 0, 0, 0, t68, t71, 0, -t58, -t80, t46 + (-t85 - t86) * t68, pkin(7) * t50 + (-pkin(2) * t68 + t75) * t65, t73, -pkin(2) * t58 + t73 * qJ(3), -t29 * t65 + t41 * t50 + t46, t74, -t29 * t64 + (-t41 * t68 - t75) * t65, t74 * qJ(3) + t29 * t41, t31 * t37, t37 * t30 - t31 * t36, t37 * t71, -t36 * t71, 0, t21 * t71 + t23 * t36 - t33 * t30, -t22 * t71 + t23 * t37 + t33 * t31, t12 * t18, -t18 * t11 - t12 * t17, t18 * t71, -t17 * t71, 0, t10 * t17 + t20 * t11 + t6 * t71, t10 * t18 + t20 * t12 - t7 * t71; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t65, -0.2e1 * t86, t38, pkin(2) ^ 2 + t76, t65 * t90, t38, t64 * t90, t41 ^ 2 + t76, t37 ^ 2, -0.2e1 * t37 * t36, 0, 0, 0, t36 * t91, t37 * t91, t18 ^ 2, -0.2e1 * t18 * t17, 0, 0, 0, t17 * t92, t18 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, 0, t58, t50, 0, -t51, t29, 0, 0, 0, 0, 0, t30, -t31, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t64, 0, -pkin(2), -t65, 0, -t64, t41, 0, 0, 0, 0, 0, -t36, -t37, 0, 0, 0, 0, 0, -t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t51, 0, t25, 0, 0, 0, 0, 0, t70 * t71, -t67 * t71, 0, 0, 0, 0, 0, t39 * t71, -t40 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t30, t71, t8, -t9, 0, 0, t12, -t11, t71, t69 * t81 + t1, -t82 + (-t4 - t81) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, t21, -t22, 0, 0, t18, -t17, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t67, 0, 0, 0, 0, 0, t39, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t83, -0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t71, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t83, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
