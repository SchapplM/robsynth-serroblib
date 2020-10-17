% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:54:47
% EndTime: 2019-05-06 10:54:48
% DurationCPUTime: 0.58s
% Computational Cost: add. (606->88), mult. (1086->154), div. (0->0), fcn. (1311->8), ass. (0->69)
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t34 = -pkin(2) * t54 - qJ(3) * t51 - pkin(1);
t22 = pkin(3) * t54 - t34;
t47 = sin(pkin(10));
t48 = cos(pkin(10));
t24 = t47 * t51 + t48 * t54;
t17 = pkin(4) * t24 + t22;
t81 = 0.2e1 * t17;
t80 = 0.2e1 * t22;
t49 = sin(qJ(6));
t79 = -0.2e1 * t49;
t78 = -0.2e1 * t51;
t52 = cos(qJ(6));
t77 = 0.2e1 * t52;
t76 = 0.2e1 * t54;
t75 = -pkin(2) - pkin(3);
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t40 = t54 * pkin(7);
t35 = -qJ(4) * t54 + t40;
t60 = (pkin(7) - qJ(4)) * t51;
t18 = t35 * t47 - t48 * t60;
t25 = -t47 * t54 + t48 * t51;
t7 = -pkin(8) * t25 - t18;
t20 = t35 * t48 + t47 * t60;
t8 = -pkin(8) * t24 + t20;
t4 = t50 * t8 - t53 * t7;
t74 = t4 * t49;
t73 = t4 * t52;
t72 = t51 * pkin(7);
t31 = qJ(3) * t47 - t48 * t75;
t30 = -pkin(4) - t31;
t33 = qJ(3) * t48 + t47 * t75;
t15 = -t30 * t53 + t33 * t50;
t13 = pkin(5) + t15;
t71 = pkin(5) + t13;
t23 = t47 * t50 - t48 * t53;
t70 = t23 * t49;
t69 = t23 * t52;
t11 = t24 * t53 + t25 * t50;
t68 = t49 * t11;
t12 = -t24 * t50 + t25 * t53;
t67 = t49 * t12;
t66 = t49 * t52;
t65 = t52 * t12;
t44 = t51 ^ 2;
t64 = t54 ^ 2 + t44;
t63 = -0.2e1 * t12 * t11;
t62 = -0.2e1 * t66;
t61 = t49 * t65;
t59 = -pkin(5) * t12 - pkin(9) * t11;
t58 = -pkin(2) * t51 + qJ(3) * t54;
t16 = t30 * t50 + t33 * t53;
t14 = -pkin(9) + t16;
t57 = -t11 * t14 + t12 * t13;
t26 = t47 * t53 + t48 * t50;
t56 = -t11 * t26 + t12 * t23;
t45 = t52 ^ 2;
t43 = t49 ^ 2;
t36 = 0.2e1 * t66;
t10 = t12 ^ 2;
t9 = t52 * t11;
t6 = (t43 - t45) * t12;
t5 = t50 * t7 + t53 * t8;
t3 = pkin(5) * t11 - pkin(9) * t12 + t17;
t2 = t3 * t49 + t5 * t52;
t1 = t3 * t52 - t49 * t5;
t19 = [1, 0, 0, t44, t51 * t76, 0, 0, 0, pkin(1) * t76, pkin(1) * t78, -0.2e1 * t34 * t54, 0.2e1 * t64 * pkin(7), t34 * t78, pkin(7) ^ 2 * t64 + t34 ^ 2, t24 * t80, t25 * t80, 0.2e1 * t18 * t25 - 0.2e1 * t20 * t24, t18 ^ 2 + t20 ^ 2 + t22 ^ 2, t10, t63, 0, 0, 0, t11 * t81, t12 * t81, t45 * t10, t10 * t62, 0.2e1 * t11 * t65, t49 * t63, t11 ^ 2, 0.2e1 * t1 * t11 + 0.2e1 * t4 * t67, -0.2e1 * t11 * t2 + 0.2e1 * t4 * t65; 0, 0, 0, 0, 0, t51, t54, 0, -t72, -t40, -t72, t58, t40, t58 * pkin(7), t18, t20, -t24 * t33 + t25 * t31, t18 * t31 + t20 * t33, 0, 0, -t12, t11, 0, t4, t5, -t61, t6, -t68, -t9, 0, t49 * t57 + t73, t52 * t57 - t74; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0.2e1 * t31, 0.2e1 * t33, 0, t31 ^ 2 + t33 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t15, 0.2e1 * t16, t43, t36, 0, 0, 0, t13 * t77, t13 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t72, 0, 0, -t24 * t47 - t25 * t48, -t18 * t48 + t20 * t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t49, t56 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), -t48, t47, 0, -t31 * t48 + t33 * t47, 0, 0, 0, 0, 0, t23, t26, 0, 0, 0, 0, 0, t69, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t47 ^ 2 + t48 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, 0, t22, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t9, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, -t4, -t5, t61, -t6, t68, t9, 0, t49 * t59 - t73, t52 * t59 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t15, -t16, -t43, t62, 0, 0, 0, -t71 * t52, t71 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t26, 0, 0, 0, 0, 0, -t69, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t43, t36, 0, 0, 0, pkin(5) * t77, pkin(5) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t67, t11, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t52, 0, -t49 * t14, -t52 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49 * t26, -t52 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t52, 0, -t49 * pkin(9), -t52 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t19;
