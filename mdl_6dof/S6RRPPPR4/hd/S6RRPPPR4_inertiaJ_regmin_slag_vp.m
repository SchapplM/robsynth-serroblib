% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:44:15
% EndTime: 2019-05-06 08:44:17
% DurationCPUTime: 0.60s
% Computational Cost: add. (422->93), mult. (741->157), div. (0->0), fcn. (743->6), ass. (0->62)
t56 = sin(qJ(2));
t81 = -0.2e1 * t56;
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t36 = t53 ^ 2 + t54 ^ 2;
t75 = sin(qJ(6));
t76 = cos(qJ(6));
t25 = -t76 * t53 + t75 * t54;
t82 = 0.2e1 * t25;
t57 = cos(qJ(2));
t80 = 0.2e1 * t57;
t79 = 2 * qJ(3);
t78 = -pkin(4) - pkin(5);
t55 = -pkin(2) - qJ(4);
t77 = pkin(8) + t55;
t60 = t75 * t53 + t76 * t54;
t74 = t60 * t56;
t73 = t25 * t56;
t72 = t53 * t56;
t71 = t53 * t57;
t41 = t54 * t57;
t70 = t55 * t56;
t63 = -t56 * qJ(3) - pkin(1);
t22 = t55 * t57 + t63;
t45 = t56 * pkin(7);
t32 = t56 * pkin(3) + t45;
t10 = t54 * t22 + t53 * t32;
t69 = t36 * t55 ^ 2;
t46 = t57 * pkin(7);
t33 = t57 * pkin(3) + t46;
t51 = t56 ^ 2;
t68 = t57 ^ 2 + t51;
t67 = qJ(5) * t53;
t66 = t57 * qJ(3);
t7 = t56 * qJ(5) + t10;
t9 = -t53 * t22 + t54 * t32;
t62 = t54 * qJ(5) - qJ(3);
t8 = -t56 * pkin(4) - t9;
t3 = t7 * t53 - t8 * t54;
t4 = t10 * t53 + t9 * t54;
t61 = -t56 * pkin(2) + t66;
t58 = qJ(3) ^ 2;
t40 = t54 * t56;
t34 = t54 * t70;
t31 = -t57 * pkin(2) + t63;
t30 = t53 * pkin(4) - t62;
t29 = t77 * t54;
t28 = t77 * t53;
t23 = t36 * t55;
t21 = t78 * t53 + t62;
t20 = 0.2e1 * t23;
t16 = t25 * t57;
t15 = t60 * t57;
t14 = (pkin(4) * t54 + t67) * t57 + t33;
t13 = t76 * t28 - t75 * t29;
t12 = -t75 * t28 - t76 * t29;
t11 = (t78 * t54 - t67) * t57 - t33;
t6 = pkin(8) * t41 + t7;
t5 = pkin(8) * t71 + t78 * t56 - t9;
t2 = t75 * t5 + t76 * t6;
t1 = t76 * t5 - t75 * t6;
t17 = [1, 0, 0, t51, t56 * t80, 0, 0, 0, pkin(1) * t80, pkin(1) * t81, 0.2e1 * t68 * pkin(7), t31 * t80, t31 * t81, t68 * pkin(7) ^ 2 + t31 ^ 2, 0.2e1 * t33 * t41 + 0.2e1 * t9 * t56, -0.2e1 * t10 * t56 - 0.2e1 * t33 * t71 (-t10 * t54 + t53 * t9) * t80, t10 ^ 2 + t33 ^ 2 + t9 ^ 2, 0.2e1 * t14 * t41 - 0.2e1 * t8 * t56 (-t53 * t8 - t54 * t7) * t80, 0.2e1 * t14 * t71 + 0.2e1 * t7 * t56, t14 ^ 2 + t7 ^ 2 + t8 ^ 2, t16 ^ 2, 0.2e1 * t16 * t15, t16 * t81, t15 * t81, t51, -0.2e1 * t1 * t56 - 0.2e1 * t11 * t15, 0.2e1 * t11 * t16 + 0.2e1 * t2 * t56; 0, 0, 0, 0, 0, t56, t57, 0, -t45, -t46, t61, t45, t46, t61 * pkin(7), t33 * t53 + t54 * t66 + t34, t33 * t54 + (-t66 - t70) * t53, -t4, t33 * qJ(3) + t4 * t55, t14 * t53 + t30 * t41 + t34, -t3, -t14 * t54 + (t30 * t57 + t70) * t53, t14 * t30 + t3 * t55, t16 * t60, t15 * t60 - t16 * t25, -t74, t73, 0, t11 * t25 - t12 * t56 - t21 * t15, t11 * t60 + t13 * t56 + t21 * t16; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t79, pkin(2) ^ 2 + t58, t53 * t79, t54 * t79, -t20, t58 + t69, 0.2e1 * t30 * t53, -t20, -0.2e1 * t30 * t54, t30 ^ 2 + t69, t60 ^ 2, -t60 * t82, 0, 0, 0, t21 * t82, 0.2e1 * t21 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, t45, t40, -t72, 0, t4, t40, 0, t72, t3, 0, 0, 0, 0, 0, t74, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -t36, t23, 0, -t36, 0, t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t36, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t71, 0, t33, t41, 0, t71, t14, 0, 0, 0, 0, 0, t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t54, 0, qJ(3), t53, 0, -t54, t30, 0, 0, 0, 0, 0, -t25, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t71, 0, t8, 0, 0, 0, 0, 0, -t76 * t56, t75 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t54 * t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, -t56, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t25, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t17;
