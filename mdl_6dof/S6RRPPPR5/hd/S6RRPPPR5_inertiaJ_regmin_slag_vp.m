% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:54:16
% EndTime: 2019-05-06 08:54:19
% DurationCPUTime: 0.68s
% Computational Cost: add. (442->113), mult. (855->187), div. (0->0), fcn. (856->6), ass. (0->66)
t50 = sin(pkin(9));
t51 = cos(pkin(9));
t81 = t50 ^ 2 + t51 ^ 2;
t65 = pkin(3) + qJ(5);
t80 = t65 * t51;
t53 = sin(qJ(2));
t37 = t51 * t53;
t55 = cos(qJ(2));
t79 = pkin(4) * t37 + t55 * qJ(5);
t52 = sin(qJ(6));
t54 = cos(qJ(6));
t25 = t52 * t50 - t54 * t51;
t78 = -t65 * t50 - pkin(7);
t24 = t54 * t50 + t52 * t51;
t19 = t24 * t53;
t77 = -0.2e1 * t19;
t76 = 0.2e1 * t25;
t75 = -0.2e1 * t50;
t74 = 0.2e1 * t51;
t73 = 0.2e1 * t53;
t72 = 0.2e1 * t55;
t49 = t53 ^ 2;
t71 = pkin(7) * t49;
t70 = t53 * pkin(7);
t69 = t55 * pkin(7);
t36 = t50 * t53;
t66 = t54 * t55;
t64 = -pkin(5) - qJ(4);
t28 = -t55 * pkin(2) - t53 * qJ(3) - pkin(1);
t15 = t50 * t28 + t51 * t69;
t63 = t81 * qJ(3) ^ 2;
t40 = t50 * qJ(3);
t29 = t50 * pkin(4) + t40;
t30 = (pkin(4) + qJ(3)) * t51;
t62 = qJ(3) * t55;
t61 = t50 * qJ(4) + pkin(2);
t33 = t50 * t69;
t14 = t51 * t28 - t33;
t46 = t55 * pkin(3);
t13 = -t14 + t46;
t12 = t55 * qJ(4) - t15;
t60 = -pkin(2) * t53 + t62;
t59 = -t12 * t51 + t13 * t50;
t58 = -t14 * t50 + t15 * t51;
t27 = -t51 * pkin(3) - t61;
t57 = -t27 * t53 - t62;
t38 = t52 * t55;
t32 = qJ(4) * t37;
t26 = 0.2e1 * t81 * qJ(3);
t22 = t51 * pkin(8) + t30;
t21 = t50 * pkin(8) + t29;
t20 = t61 + t80;
t18 = t25 * t53;
t17 = -t32 + (pkin(3) * t50 + pkin(7)) * t53;
t16 = t64 * t50 - pkin(2) - t80;
t10 = -t53 * t78 - t32;
t9 = -pkin(4) * t36 - t12;
t8 = -t32 + (-pkin(5) * t51 - t78) * t53;
t7 = t54 * t21 + t52 * t22;
t6 = -t52 * t21 + t54 * t22;
t5 = t13 + t79;
t4 = t64 * t55 + (-pkin(4) - pkin(8)) * t36 + t15;
t3 = t33 + t46 + (pkin(8) * t53 - t28) * t51 + t79;
t2 = t54 * t3 + t52 * t4;
t1 = -t52 * t3 + t54 * t4;
t11 = [1, 0, 0, t49, t53 * t72, 0, 0, 0, pkin(1) * t72, -0.2e1 * pkin(1) * t53, -0.2e1 * t14 * t55 + 0.2e1 * t50 * t71, 0.2e1 * t15 * t55 + 0.2e1 * t51 * t71 (-t14 * t51 - t15 * t50) * t73, t49 * pkin(7) ^ 2 + t14 ^ 2 + t15 ^ 2 (t12 * t50 + t13 * t51) * t73, -0.2e1 * t13 * t55 - 0.2e1 * t17 * t36, 0.2e1 * t12 * t55 - 0.2e1 * t17 * t37, t12 ^ 2 + t13 ^ 2 + t17 ^ 2, -0.2e1 * t10 * t37 - 0.2e1 * t9 * t55 (-t5 * t51 + t50 * t9) * t73, 0.2e1 * t10 * t36 + 0.2e1 * t5 * t55, t10 ^ 2 + t5 ^ 2 + t9 ^ 2, t19 ^ 2, t18 * t77, t55 * t77, t18 * t72, t55 ^ 2, -0.2e1 * t1 * t55 + 0.2e1 * t8 * t18, 0.2e1 * t8 * t19 + 0.2e1 * t2 * t55; 0, 0, 0, 0, 0, t53, t55, 0, -t70, -t69, -pkin(7) * t37 + t60 * t50, pkin(7) * t36 + t60 * t51, t58, -pkin(2) * t70 + t58 * qJ(3), t59, t17 * t51 + t57 * t50, -t17 * t50 + t57 * t51, t59 * qJ(3) + t17 * t27, -t10 * t50 + t20 * t37 - t30 * t55 (-t29 * t53 - t9) * t51 + (t30 * t53 - t5) * t50, -t10 * t51 - t20 * t36 + t29 * t55, -t10 * t20 + t5 * t29 + t9 * t30, t19 * t25, -t25 * t18 + t19 * t24, -t25 * t55, -t24 * t55, 0, t16 * t18 - t8 * t24 - t6 * t55, t16 * t19 + t8 * t25 + t7 * t55; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, pkin(2) * t74, pkin(2) * t75, t26, pkin(2) ^ 2 + t63, t26, t27 * t74, t27 * t75, t27 ^ 2 + t63, 0.2e1 * t20 * t50, -0.2e1 * t29 * t50 - 0.2e1 * t30 * t51, t20 * t74, t20 ^ 2 + t29 ^ 2 + t30 ^ 2, t25 ^ 2, t24 * t76, 0, 0, 0, -0.2e1 * t16 * t24, t16 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t37, 0, t70, 0, -t36, -t37, t17, -t37, 0, t36, t10, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, -pkin(2), 0, t51, -t50, t27, -t50, 0, -t51, -t20, 0, 0, 0, 0, 0, -t24, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t55, 0, t13, 0, -t37, t55, t5, 0, 0, 0, 0, 0, t38, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, t40, 0, -t50, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t36, 0, t9, 0, 0, 0, 0, 0, -t66, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, -t55, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t24, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
