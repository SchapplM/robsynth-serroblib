% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:00:32
% EndTime: 2019-05-05 03:00:35
% DurationCPUTime: 1.07s
% Computational Cost: add. (375->95), mult. (723->149), div. (0->0), fcn. (787->8), ass. (0->70)
t46 = cos(pkin(6));
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t45 = sin(pkin(6));
t50 = sin(qJ(2));
t76 = t45 * t50;
t12 = -t46 * t52 + t49 * t76;
t14 = t46 * t49 + t52 * t76;
t63 = t12 * t49 + t14 * t52;
t42 = t49 ^ 2;
t44 = t52 ^ 2;
t86 = t42 + t44;
t11 = t14 ^ 2;
t39 = t45 ^ 2;
t53 = cos(qJ(2));
t27 = t39 * t53 ^ 2;
t85 = t12 ^ 2 + t11 + t27;
t80 = t52 * pkin(8);
t21 = -t52 * qJ(5) + t80;
t84 = t21 ^ 2;
t83 = -0.2e1 * t49;
t82 = -0.2e1 * t52;
t81 = 0.2e1 * t52;
t54 = -pkin(3) - pkin(4);
t78 = t14 * t21;
t48 = sin(qJ(6));
t41 = t48 ^ 2;
t77 = t41 * t52;
t28 = t45 * t53;
t75 = t48 * t49;
t51 = cos(qJ(6));
t74 = t48 * t51;
t73 = t48 * t52;
t72 = t49 * t52;
t71 = t51 * t49;
t70 = t51 * t52;
t68 = t86 * pkin(8) ^ 2;
t43 = t51 ^ 2;
t25 = t43 + t41;
t67 = t14 * qJ(4);
t66 = t52 * qJ(4);
t65 = -0.2e1 * t72;
t26 = 0.2e1 * t72;
t19 = -t52 * pkin(3) - t49 * qJ(4) - pkin(2);
t64 = t48 * t70;
t16 = t52 * pkin(4) - t19;
t62 = t63 * pkin(8);
t34 = t49 * pkin(8);
t20 = -t49 * qJ(5) + t34;
t7 = t49 * pkin(5) + t52 * pkin(9) + t16;
t3 = -t48 * t20 + t51 * t7;
t4 = t51 * t20 + t48 * t7;
t61 = t3 * t51 + t4 * t48;
t1 = -t3 * t48 + t4 * t51;
t5 = -t12 * t48 + t51 * t28;
t6 = t12 * t51 + t48 * t28;
t60 = t6 * t48 + t5 * t51;
t2 = -t5 * t48 + t6 * t51;
t59 = -t49 * pkin(3) + t66;
t40 = -pkin(9) + t54;
t47 = qJ(4) + pkin(5);
t58 = -t40 * t49 - t47 * t52;
t56 = qJ(4) ^ 2;
t55 = 0.2e1 * qJ(4);
t30 = t43 * t52;
t24 = t52 * t28;
t23 = t49 * t28;
t18 = 0.2e1 * t86 * pkin(8);
t17 = t25 * t40;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t50 ^ 2 + t46 ^ 2 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t76, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, t63, pkin(2) * t28 + t62, 0, 0, 0, 0, 0, 0, t24, t63, t23, -t19 * t28 + t62, 0, 0, 0, 0, 0, 0, t23, -t24, -t63, t12 * t20 + t16 * t28 + t78, 0, 0, 0, 0, 0, 0, -t14 * t73 + t5 * t49, -t14 * t70 - t6 * t49, t60 * t52, t5 * t3 + t6 * t4 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t42, t26, 0, t44, 0, 0, pkin(2) * t81, pkin(2) * t83, t18, pkin(2) ^ 2 + t68, t42, 0, t65, 0, 0, t44, t19 * t82, t18, t19 * t83, t19 ^ 2 + t68, t44, t26, 0, t42, 0, 0, 0.2e1 * t16 * t49, t16 * t82, -0.2e1 * t20 * t49 - 0.2e1 * t21 * t52, t16 ^ 2 + t20 ^ 2 + t84, t43 * t44, -0.2e1 * t44 * t74, t51 * t65, t41 * t44, t48 * t26, t42, -0.2e1 * t21 * t73 + 0.2e1 * t3 * t49, -0.2e1 * t21 * t70 - 0.2e1 * t4 * t49, t61 * t81, t3 ^ 2 + t4 ^ 2 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, t14, -t12 * pkin(3) + t67, 0, 0, 0, 0, 0, 0, t14, t12, 0, t12 * t54 + t67, 0, 0, 0, 0, 0, 0, t14 * t51, -t14 * t48, -t2, t14 * t47 + t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t52, 0, -t34, -t80, 0, 0, 0, t49, 0, 0, -t52, 0, -t34, t59, t80, t59 * pkin(8), 0, 0, t52, 0, t49, 0, t21, t20, -t54 * t49 - t66, t21 * qJ(4) + t20 * t54, t64, t30 - t77, -t75, -t64, -t71, 0, t21 * t51 + t58 * t48, -t21 * t48 + t58 * t51, -t1, t1 * t40 + t21 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t55, pkin(3) ^ 2 + t56, 0, 0, 0, 0, 0, 1, t55, 0.2e1 * t54, 0, t54 ^ 2 + t56, t41, 0.2e1 * t74, 0, t43, 0, 0, 0.2e1 * t47 * t51, -0.2e1 * t47 * t48, -0.2e1 * t17, t25 * t40 ^ 2 + t47 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t20, 0, 0, 0, 0, 0, 0, -t75, -t71, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 1, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t52, 0, t16, 0, 0, 0, 0, 0, 0, t71, -t75, t30 + t77, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, t73, t49, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, -t51, 0, -t48 * t40, -t51 * t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t8;
