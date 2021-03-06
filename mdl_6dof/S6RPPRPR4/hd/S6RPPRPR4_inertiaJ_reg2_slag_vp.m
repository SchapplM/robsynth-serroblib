% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:16:01
% EndTime: 2019-05-05 14:16:03
% DurationCPUTime: 1.05s
% Computational Cost: add. (838->104), mult. (1314->183), div. (0->0), fcn. (1496->8), ass. (0->65)
t47 = sin(pkin(9));
t51 = sin(qJ(4));
t43 = t51 ^ 2;
t53 = cos(qJ(4));
t45 = t53 ^ 2;
t67 = t43 + t45;
t83 = t67 * t47;
t49 = cos(pkin(9));
t54 = -pkin(1) - pkin(2);
t30 = t49 * qJ(2) + t47 * t54;
t27 = -pkin(7) + t30;
t66 = qJ(5) - t27;
t10 = t66 * t53;
t46 = sin(pkin(10));
t48 = cos(pkin(10));
t61 = t66 * t51;
t3 = -t46 * t10 - t48 * t61;
t82 = t3 ^ 2;
t23 = t46 * t53 + t48 * t51;
t11 = t23 * t47;
t81 = t11 ^ 2;
t21 = t46 * t51 - t48 * t53;
t80 = t21 ^ 2;
t79 = t23 ^ 2;
t78 = -0.2e1 * t23;
t77 = 0.2e1 * t53;
t76 = t3 * t11;
t75 = t3 * t21;
t74 = t46 * pkin(4);
t73 = t48 * pkin(4);
t72 = t11 * t21;
t50 = sin(qJ(6));
t16 = t21 * t50;
t52 = cos(qJ(6));
t15 = t21 * t52;
t71 = t50 * t23;
t70 = t50 * t52;
t69 = t52 * t23;
t28 = t47 * qJ(2) - t49 * t54;
t42 = t50 ^ 2;
t44 = t52 ^ 2;
t68 = t42 + t44;
t65 = t21 * t78;
t26 = pkin(3) + t28;
t64 = t50 * t69;
t63 = t68 * t23;
t33 = pkin(8) + t74;
t62 = t68 * t33;
t14 = t53 * pkin(4) + t26;
t5 = -t48 * t10 + t46 * t61;
t6 = -t21 * pkin(5) + t23 * pkin(8) + t14;
t1 = -t50 * t5 + t52 * t6;
t2 = t52 * t5 + t50 * t6;
t60 = t1 * t52 + t2 * t50;
t59 = -t1 * t50 + t2 * t52;
t13 = t21 * t47;
t7 = t50 * t13 - t52 * t49;
t8 = -t52 * t13 - t50 * t49;
t58 = t8 * t50 + t7 * t52;
t57 = -t7 * t50 + t8 * t52;
t34 = -pkin(5) - t73;
t56 = t21 * t33 - t23 * t34;
t41 = t49 ^ 2;
t40 = t47 ^ 2;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, 0.2e1 * t30, 0, t28 ^ 2 + t30 ^ 2, t43, t51 * t77, 0, t45, 0, 0, t26 * t77, -0.2e1 * t26 * t51, -0.2e1 * t67 * t27, t67 * t27 ^ 2 + t26 ^ 2, t79, t65, 0, t80, 0, 0, -0.2e1 * t14 * t21, t14 * t78, 0.2e1 * t5 * t21 - 0.2e1 * t3 * t23, t14 ^ 2 + t5 ^ 2 + t82, t44 * t79, -0.2e1 * t79 * t70, 0.2e1 * t21 * t69, t42 * t79, t50 * t65, t80, -0.2e1 * t1 * t21 - 0.2e1 * t3 * t71, 0.2e1 * t2 * t21 - 0.2e1 * t3 * t69, 0.2e1 * t60 * t23, t1 ^ 2 + t2 ^ 2 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t49, t47, 0, -t28 * t49 + t30 * t47, 0, 0, 0, 0, 0, 0, -t49 * t53, t51 * t49, -t83, -t26 * t49 + t27 * t83, 0, 0, 0, 0, 0, 0, t49 * t21, t49 * t23, -t11 * t23 - t13 * t21, -t5 * t13 - t14 * t49 + t76, 0, 0, 0, 0, 0, 0, -t11 * t71 - t7 * t21, -t11 * t69 + t8 * t21, t58 * t23, t1 * t7 + t2 * t8 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 + t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t40 + t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t41 + t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t23 + t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t23 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t23 + t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t23 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 + t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t79 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, -t53, 0, -t51 * t27, -t53 * t27, 0, 0, 0, 0, -t23, 0, t21, 0, -t3, -t5 (t21 * t46 + t23 * t48) * pkin(4) (-t3 * t48 + t46 * t5) * pkin(4), -t64 (t42 - t44) * t23, -t16, t64, -t15, 0, -t3 * t52 + t56 * t50, t3 * t50 + t56 * t52, t59, t3 * t34 + t59 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t47, -t53 * t47, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t13, 0 (-t11 * t48 - t13 * t46) * pkin(4), 0, 0, 0, 0, 0, 0, -t11 * t52, t11 * t50, t57, t11 * t34 + t57 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t51, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, 0 (-t21 * t48 + t23 * t46) * pkin(4), 0, 0, 0, 0, 0, 0, -t15, t16, t63, t21 * t34 + t23 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t74, 0 (t46 ^ 2 + t48 ^ 2) * pkin(4) ^ 2, t42, 0.2e1 * t70, 0, t44, 0, 0, -0.2e1 * t34 * t52, 0.2e1 * t34 * t50, 0.2e1 * t62, t68 * t33 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, 0, t14, 0, 0, 0, 0, 0, 0, -t15, t16, t63, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, t71, -t21, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t69, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t52, 0, -t50 * t33, -t52 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t4;
