% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:35:50
% EndTime: 2019-05-07 08:35:53
% DurationCPUTime: 0.76s
% Computational Cost: add. (661->121), mult. (1253->221), div. (0->0), fcn. (1295->6), ass. (0->70)
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t86 = -t59 * pkin(3) - t56 * qJ(4);
t57 = sin(qJ(2));
t43 = t59 * t57;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t71 = t58 * t56;
t22 = t55 * t43 - t57 * t71;
t85 = -0.2e1 * t22;
t33 = -pkin(2) + t86;
t24 = t59 * pkin(4) - t33;
t84 = 0.2e1 * t24;
t83 = -0.2e1 * t33;
t82 = -0.2e1 * t57;
t60 = cos(qJ(2));
t81 = 0.2e1 * t60;
t80 = -pkin(3) - pkin(4);
t79 = pkin(2) * t56;
t78 = pkin(2) * t59;
t77 = pkin(3) * t56;
t76 = pkin(7) * t56;
t75 = pkin(7) * t59;
t50 = t60 * pkin(3);
t34 = -t60 * pkin(2) - t57 * pkin(8) - pkin(1);
t72 = t56 * t60;
t69 = pkin(7) * t72 - t59 * t34;
t18 = t50 + t69;
t10 = t60 * pkin(4) - pkin(9) * t43 + t18;
t70 = t59 * t60;
t20 = pkin(7) * t70 + t56 * t34;
t66 = t60 * qJ(4);
t17 = -t66 + t20;
t74 = t56 * t57;
t12 = pkin(9) * t74 + t17;
t4 = t55 * t10 + t58 * t12;
t73 = t56 * t59;
t51 = t56 ^ 2;
t53 = t59 ^ 2;
t68 = t51 + t53;
t67 = t59 * qJ(4);
t65 = t57 * t81;
t64 = -t58 * t10 + t55 * t12;
t46 = t56 * pkin(8);
t35 = -t56 * pkin(9) + t46;
t47 = t59 * pkin(8);
t36 = -t59 * pkin(9) + t47;
t15 = -t58 * t35 + t55 * t36;
t31 = t55 * qJ(4) - t58 * t80;
t63 = t67 - t77;
t62 = t17 * t59 + t18 * t56;
t16 = t55 * t35 + t58 * t36;
t27 = t55 * t56 + t58 * t59;
t38 = t57 * t67;
t14 = t38 + (t80 * t56 - pkin(7)) * t57;
t54 = t60 ^ 2;
t52 = t57 ^ 2;
t40 = pkin(8) * t72;
t32 = t58 * qJ(4) + t55 * t80;
t30 = -pkin(5) - t31;
t28 = -t55 * t59 + t71;
t23 = t27 * t57;
t21 = -t38 + (pkin(7) + t77) * t57;
t13 = t27 * pkin(5) + t24;
t7 = -t27 * qJ(6) + t16;
t6 = -t28 * qJ(6) - t15;
t5 = t22 * pkin(5) + t14;
t2 = -t22 * qJ(6) + t4;
t1 = t60 * pkin(5) - t23 * qJ(6) - t64;
t3 = [1, 0, 0, t52, t65, 0, 0, 0, pkin(1) * t81, pkin(1) * t82, t53 * t52, -0.2e1 * t52 * t73, t70 * t82, t56 * t65, t54, 0.2e1 * t52 * t76 + 0.2e1 * t60 * t69, 0.2e1 * t20 * t60 + 0.2e1 * t52 * t75, 0.2e1 * t18 * t60 + 0.2e1 * t21 * t74, 0.2e1 * (-t17 * t56 + t18 * t59) * t57, -0.2e1 * t17 * t60 - 0.2e1 * t21 * t43, t17 ^ 2 + t18 ^ 2 + t21 ^ 2, t23 ^ 2, t23 * t85, t23 * t81, t60 * t85, t54, 0.2e1 * t14 * t22 - 0.2e1 * t60 * t64, 0.2e1 * t14 * t23 - 0.2e1 * t4 * t60, -0.2e1 * t1 * t23 - 0.2e1 * t2 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t57, t60, 0, -t57 * pkin(7), -t60 * pkin(7), t56 * t43 (-t51 + t53) * t57, -t72, -t70, 0, t40 + (-t75 - t79) * t57, pkin(8) * t70 + (t76 - t78) * t57, -t21 * t59 + t33 * t74 + t40, t62, -t21 * t56 + (-pkin(8) * t60 - t33 * t57) * t59, pkin(8) * t62 + t21 * t33, t23 * t28, -t28 * t22 - t23 * t27, t60 * t28, -t60 * t27, 0, t14 * t27 - t15 * t60 + t24 * t22, t14 * t28 - t16 * t60 + t24 * t23, -t1 * t28 - t2 * t27 - t7 * t22 - t6 * t23, t1 * t6 + t5 * t13 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t51, 0.2e1 * t73, 0, 0, 0, 0.2e1 * t78, -0.2e1 * t79, t59 * t83, 0.2e1 * t68 * pkin(8), t56 * t83, pkin(8) ^ 2 * t68 + t33 ^ 2, t28 ^ 2, -0.2e1 * t28 * t27, 0, 0, 0, t27 * t84, t28 * t84, -0.2e1 * t7 * t27 - 0.2e1 * t6 * t28, t13 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t74, -t60, -t69, -t20, -0.2e1 * t50 - t69, t86 * t57, -0.2e1 * t66 + t20, -t18 * pkin(3) + t17 * qJ(4), 0, 0, -t23, t22, -t60, -t31 * t60 + t64, -t32 * t60 + t4, -t32 * t22 - t30 * t23, t1 * t30 + t2 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t59, 0, -t46, -t47, -t46, t63, t47, t63 * pkin(8), 0, 0, -t28, t27, 0, t15, t16, -t32 * t27 - t30 * t28, t6 * t30 + t7 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t31, 0.2e1 * t32, 0, t30 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t43, 0, t18, 0, 0, 0, 0, 0, t58 * t60, -t55 * t60, -t55 * t22 - t58 * t23, t1 * t58 + t2 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, t46, 0, 0, 0, 0, 0, 0, 0, -t55 * t27 - t58 * t28, t7 * t55 + t6 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t58, t55, 0, t30 * t58 + t32 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t55 ^ 2 + t58 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t60, -t64, -t4, -pkin(5) * t23, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t15, -t16, -pkin(5) * t28, t6 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t31, -t32, 0, t30 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55, 0, t58 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
