% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:21:49
% EndTime: 2019-05-07 04:21:51
% DurationCPUTime: 0.58s
% Computational Cost: add. (840->83), mult. (1568->153), div. (0->0), fcn. (1866->8), ass. (0->71)
t51 = sin(qJ(3));
t52 = sin(qJ(2));
t54 = cos(qJ(3));
t55 = cos(qJ(2));
t31 = t51 * t52 - t54 * t55;
t32 = t51 * t55 + t54 * t52;
t48 = sin(pkin(10));
t49 = cos(pkin(10));
t19 = t49 * t31 + t48 * t32;
t20 = -t48 * t31 + t49 * t32;
t44 = -t55 * pkin(2) - pkin(1);
t23 = t31 * pkin(3) + t44;
t57 = -t20 * qJ(5) + t23;
t8 = t19 * pkin(4) + t57;
t79 = -0.2e1 * t8;
t78 = 0.2e1 * t44;
t50 = sin(qJ(6));
t77 = 0.2e1 * t50;
t53 = cos(qJ(6));
t76 = 0.2e1 * t53;
t75 = 0.2e1 * t55;
t74 = pkin(7) + pkin(8);
t73 = t51 * pkin(2);
t45 = t54 * pkin(2);
t43 = t45 + pkin(3);
t72 = pkin(3) + t43;
t63 = t49 * t73;
t29 = t48 * t43 + t63;
t25 = qJ(5) + t29;
t71 = t25 * t19;
t39 = t48 * pkin(3) + qJ(5);
t70 = t39 * t19;
t69 = t50 * t19;
t68 = t50 * t20;
t67 = t53 * t19;
t17 = t53 * t20;
t66 = t53 * t50;
t65 = t25 + t39;
t64 = 0.2e1 * t19 * t20;
t34 = t74 * t52;
t35 = t74 * t55;
t22 = t51 * t34 - t54 * t35;
t15 = -t31 * qJ(4) - t22;
t21 = -t54 * t34 - t51 * t35;
t59 = -t32 * qJ(4) + t21;
t11 = t49 * t15 + t48 * t59;
t9 = t48 * t15 - t49 * t59;
t62 = t11 ^ 2 + t9 ^ 2;
t41 = -t49 * pkin(3) - pkin(4);
t36 = t48 * t73;
t28 = t49 * t43 - t36;
t27 = -pkin(4) - t28;
t24 = -pkin(9) + t27;
t61 = -t20 * t24 + t71;
t37 = -pkin(9) + t41;
t60 = -t20 * t37 + t70;
t58 = -0.2e1 * t11 * t19 + 0.2e1 * t9 * t20;
t47 = t53 ^ 2;
t46 = t50 ^ 2;
t38 = -0.2e1 * t66;
t18 = t19 ^ 2;
t16 = t19 * t66;
t13 = (-t46 + t47) * t19;
t7 = -t19 * pkin(5) + t11;
t6 = t20 * pkin(5) + t9;
t5 = (pkin(4) + pkin(9)) * t19 + t57;
t4 = t7 * t53;
t3 = t7 * t50;
t2 = t53 * t5 + t50 * t6;
t1 = -t50 * t5 + t53 * t6;
t10 = [1, 0, 0, t52 ^ 2, t52 * t75, 0, 0, 0, pkin(1) * t75, -0.2e1 * pkin(1) * t52, t32 ^ 2, -0.2e1 * t32 * t31, 0, 0, 0, t31 * t78, t32 * t78, t58, t23 ^ 2 + t62, t58, t19 * t79, t20 * t79, t8 ^ 2 + t62, t46 * t18, 0.2e1 * t18 * t66, t50 * t64, t53 * t64, t20 ^ 2, 0.2e1 * t1 * t20 - 0.2e1 * t7 * t67, -0.2e1 * t2 * t20 + 0.2e1 * t7 * t69; 0, 0, 0, 0, 0, t52, t55, 0, -t52 * pkin(7), -t55 * pkin(7), 0, 0, t32, -t31, 0, t21, t22, -t29 * t19 - t28 * t20, t11 * t29 - t9 * t28, t27 * t20 - t71, t9, t11, t11 * t25 + t9 * t27, t16, t13, t17, -t68, 0, -t61 * t53 + t3, t61 * t50 + t4; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t73, 0, t28 ^ 2 + t29 ^ 2, 0, 0.2e1 * t27, 0.2e1 * t25, t25 ^ 2 + t27 ^ 2, t47, t38, 0, 0, 0, t25 * t77, t25 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t21, t22 (-t19 * t48 - t20 * t49) * pkin(3) (t11 * t48 - t49 * t9) * pkin(3), t41 * t20 - t70, t9, t11, t11 * t39 + t9 * t41, t16, t13, t17, -t68, 0, -t60 * t53 + t3, t60 * t50 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, -t73, 0 (t28 * t49 + t29 * t48) * pkin(3), 0, -t72 * t49 - 0.2e1 * pkin(4) + t36, t72 * t48 + 0.2e1 * qJ(5) + t63, t25 * t39 + t27 * t41, t47, t38, 0, 0, 0, t65 * t50, t65 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t48 ^ 2 + t49 ^ 2) * pkin(3) ^ 2, 0, 0.2e1 * t41, 0.2e1 * t39, t39 ^ 2 + t41 ^ 2, t47, t38, 0, 0, 0, t39 * t77, t39 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t19, -t20, t8, 0, 0, 0, 0, 0, -t68, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, t9, 0, 0, 0, 0, 0, t17, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t67, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t50, 0, t53 * t24, -t50 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t50, 0, t53 * t37, -t50 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
