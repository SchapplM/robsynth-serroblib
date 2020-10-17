% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:50:39
% EndTime: 2019-05-05 22:50:41
% DurationCPUTime: 0.64s
% Computational Cost: add. (1240->106), mult. (2457->205), div. (0->0), fcn. (2996->10), ass. (0->75)
t64 = sin(pkin(11));
t66 = cos(pkin(11));
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t45 = -t64 * t69 + t66 * t71;
t59 = -t71 * pkin(4) - pkin(3);
t35 = -t45 * pkin(5) + t59;
t92 = 0.2e1 * t35;
t65 = sin(pkin(10));
t67 = cos(pkin(10));
t70 = sin(qJ(3));
t87 = cos(qJ(3));
t46 = t70 * t65 - t67 * t87;
t91 = -0.2e1 * t46;
t90 = 0.2e1 * t46;
t58 = -t67 * pkin(2) - pkin(1);
t89 = 0.2e1 * t58;
t88 = pkin(4) * t64;
t48 = t65 * t87 + t70 * t67;
t28 = t46 * pkin(3) - t48 * pkin(8) + t58;
t79 = pkin(7) + qJ(2);
t52 = t79 * t65;
t53 = t79 * t67;
t32 = -t70 * t52 + t53 * t87;
t18 = t71 * t28 - t69 * t32;
t76 = qJ(5) * t48;
t10 = t46 * pkin(4) - t71 * t76 + t18;
t81 = t71 * t32;
t14 = t81 + (t28 - t76) * t69;
t7 = t64 * t10 + t66 * t14;
t86 = cos(qJ(6));
t47 = t64 * t71 + t66 * t69;
t68 = sin(qJ(6));
t30 = t68 * t45 + t47 * t86;
t85 = t30 * t46;
t84 = t69 * t46;
t83 = t69 * t48;
t82 = t69 * t71;
t80 = t71 * t48;
t78 = -qJ(5) - pkin(8);
t54 = t78 * t69;
t55 = t78 * t71;
t34 = t64 * t54 - t66 * t55;
t77 = t65 ^ 2 + t67 ^ 2;
t75 = t48 * t91;
t26 = -t64 * t83 + t66 * t80;
t6 = t66 * t10 - t64 * t14;
t4 = t46 * pkin(5) - t26 * pkin(9) + t6;
t25 = t47 * t48;
t5 = -t25 * pkin(9) + t7;
t1 = t86 * t4 - t68 * t5;
t33 = t66 * t54 + t64 * t55;
t74 = -pkin(3) * t48 - pkin(8) * t46;
t31 = t52 * t87 + t70 * t53;
t21 = pkin(4) * t83 + t31;
t2 = t68 * t4 + t5 * t86;
t63 = t71 ^ 2;
t62 = t69 ^ 2;
t57 = t66 * pkin(4) + pkin(5);
t43 = t48 ^ 2;
t42 = t46 ^ 2;
t41 = t71 * t46;
t39 = t68 * t57 + t86 * t88;
t38 = t57 * t86 - t68 * t88;
t29 = -t45 * t86 + t68 * t47;
t24 = t45 * pkin(9) + t34;
t23 = -t47 * pkin(9) + t33;
t20 = t29 * t46;
t19 = t69 * t28 + t81;
t17 = t25 * pkin(5) + t21;
t16 = -t68 * t25 + t26 * t86;
t15 = t25 * t86 + t68 * t26;
t13 = t68 * t23 + t24 * t86;
t12 = t23 * t86 - t68 * t24;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t67, -0.2e1 * pkin(1) * t65, 0.2e1 * t77 * qJ(2), qJ(2) ^ 2 * t77 + pkin(1) ^ 2, t43, t75, 0, 0, 0, t46 * t89, t48 * t89, t63 * t43, -0.2e1 * t43 * t82, t80 * t90, t69 * t75, t42, 0.2e1 * t18 * t46 + 0.2e1 * t31 * t83, -0.2e1 * t19 * t46 + 0.2e1 * t31 * t80, -0.2e1 * t7 * t25 - 0.2e1 * t6 * t26, t21 ^ 2 + t6 ^ 2 + t7 ^ 2, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t90, t15 * t91, t42, 0.2e1 * t1 * t46 + 0.2e1 * t17 * t15, 0.2e1 * t17 * t16 - 0.2e1 * t2 * t46; 0, 0, 0, -t67, t65, 0, -pkin(1), 0, 0, 0, 0, 0, t46, t48, 0, 0, 0, 0, 0, t41, -t84, -t47 * t25 - t45 * t26, t6 * t45 + t7 * t47, 0, 0, 0, 0, 0, -t20, -t85; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 ^ 2 + t47 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t46, 0, -t31, -t32, t69 * t80 (-t62 + t63) * t48, t84, t41, 0, -t31 * t71 + t69 * t74, t31 * t69 + t71 * t74, -t34 * t25 - t33 * t26 + t7 * t45 - t6 * t47, t21 * t59 + t6 * t33 + t7 * t34, t16 * t30, -t30 * t15 - t16 * t29, t85, -t20, 0, t12 * t46 + t35 * t15 + t17 * t29, -t13 * t46 + t35 * t16 + t17 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t33 + t47 * t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, 0.2e1 * t82, 0, 0, 0, 0.2e1 * pkin(3) * t71, -0.2e1 * pkin(3) * t69, -0.2e1 * t33 * t47 + 0.2e1 * t34 * t45, t33 ^ 2 + t34 ^ 2 + t59 ^ 2, t30 ^ 2, -0.2e1 * t30 * t29, 0, 0, 0, t29 * t92, t30 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t83, t46, t18, -t19 (-t25 * t64 - t26 * t66) * pkin(4) (t6 * t66 + t64 * t7) * pkin(4), 0, 0, t16, -t15, t46, t38 * t46 + t1, -t39 * t46 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t69, 0 (t45 * t66 + t47 * t64) * pkin(4), 0, 0, 0, 0, 0, -t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t71, 0, -t69 * pkin(8), -t71 * pkin(8) (t45 * t64 - t47 * t66) * pkin(4) (t33 * t66 + t34 * t64) * pkin(4), 0, 0, t30, -t29, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t64 ^ 2 + t66 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, t29, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t46, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
