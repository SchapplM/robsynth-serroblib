% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:27:53
% EndTime: 2019-05-06 09:27:55
% DurationCPUTime: 0.66s
% Computational Cost: add. (740->108), mult. (1292->176), div. (0->0), fcn. (1367->6), ass. (0->66)
t49 = sin(pkin(9));
t52 = sin(qJ(5));
t50 = cos(pkin(9));
t78 = cos(qJ(5));
t65 = t78 * t50;
t26 = t52 * t49 - t65;
t86 = t26 ^ 2;
t70 = t52 * t50;
t57 = -t78 * t49 - t70;
t85 = -0.2e1 * t57;
t53 = sin(qJ(2));
t84 = -0.2e1 * t53;
t83 = 0.2e1 * t53;
t54 = cos(qJ(2));
t82 = 0.2e1 * t54;
t81 = 2 * qJ(3);
t80 = t53 * pkin(5);
t51 = -pkin(2) - qJ(4);
t79 = -pkin(8) + t51;
t63 = -t53 * qJ(3) - pkin(1);
t24 = t51 * t54 + t63;
t41 = t53 * pkin(7);
t33 = t53 * pkin(3) + t41;
t14 = t50 * t24 + t49 * t33;
t72 = t50 * t54;
t11 = -pkin(8) * t72 + t14;
t30 = t50 * t33;
t9 = t53 * pkin(4) + t30 + (pkin(8) * t54 - t24) * t49;
t4 = t78 * t11 + t52 * t9;
t31 = t79 * t49;
t15 = t52 * t31 - t79 * t65;
t77 = t15 * t53;
t16 = t78 * t31 + t79 * t70;
t76 = t16 * t53;
t18 = t57 * t54;
t75 = t26 * t18;
t21 = t26 * t53;
t74 = t57 * t53;
t73 = t49 * t54;
t42 = t54 * pkin(7);
t34 = t54 * pkin(3) + t42;
t36 = t49 ^ 2 + t50 ^ 2;
t47 = t53 ^ 2;
t69 = t54 ^ 2 + t47;
t68 = t53 * qJ(6);
t67 = t54 * qJ(3);
t38 = t49 * pkin(4) + qJ(3);
t22 = pkin(4) * t72 + t34;
t66 = t57 ^ 2 + t86;
t64 = t52 * t11 - t78 * t9;
t1 = t68 + t4;
t2 = t64 - t80;
t62 = -t1 * t57 + t2 * t26;
t61 = -t53 * pkin(2) + t67;
t60 = t26 * pkin(5) + qJ(6) * t57;
t13 = -t49 * t24 + t30;
t5 = t13 * t50 + t14 * t49;
t59 = t15 * t26 - t16 * t57;
t58 = t51 * t53 + t67;
t55 = qJ(3) ^ 2;
t32 = -t54 * pkin(2) + t63;
t25 = t36 * t51;
t17 = t26 * t54;
t12 = -pkin(5) * t57 + t26 * qJ(6) + t38;
t6 = -t17 * pkin(5) - t18 * qJ(6) + t22;
t3 = [1, 0, 0, t47, t53 * t82, 0, 0, 0, pkin(1) * t82, pkin(1) * t84, 0.2e1 * t69 * pkin(7), t32 * t82, t32 * t84, t69 * pkin(7) ^ 2 + t32 ^ 2, 0.2e1 * t13 * t53 + 0.2e1 * t34 * t72, -0.2e1 * t14 * t53 - 0.2e1 * t34 * t73 (t13 * t49 - t14 * t50) * t82, t13 ^ 2 + t14 ^ 2 + t34 ^ 2, t18 ^ 2, 0.2e1 * t18 * t17, t18 * t83, t17 * t83, t47, -0.2e1 * t22 * t17 - 0.2e1 * t53 * t64, 0.2e1 * t22 * t18 - 0.2e1 * t4 * t53, -0.2e1 * t6 * t17 - 0.2e1 * t2 * t53, 0.2e1 * t1 * t17 + 0.2e1 * t2 * t18, 0.2e1 * t1 * t53 - 0.2e1 * t6 * t18, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t53, t54, 0, -t41, -t42, t61, t41, t42, t61 * pkin(7), t34 * t49 + t58 * t50, t34 * t50 - t58 * t49, -t5, t34 * qJ(3) + t5 * t51, -t75, -t26 * t17 + t18 * t57, -t21, t74, 0, -t38 * t17 - t22 * t57 - t77, t38 * t18 - t22 * t26 - t76, -t12 * t17 - t57 * t6 - t77, t15 * t18 + t16 * t17 - t62, -t12 * t18 + t6 * t26 + t76, t1 * t16 + t6 * t12 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t81, pkin(2) ^ 2 + t55, t49 * t81, t50 * t81, -0.2e1 * t25, t36 * t51 ^ 2 + t55, t86, t26 * t85, 0, 0, 0, t38 * t85, -0.2e1 * t38 * t26, t12 * t85, -0.2e1 * t59, 0.2e1 * t12 * t26, t12 ^ 2 + t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, t41, t50 * t53, -t49 * t53, 0, t5, 0, 0, 0, 0, 0, -t21, t74, -t21, -t17 * t57 + t75, -t74, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -t36, t25, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t73, 0, t34, 0, 0, 0, 0, 0, -t17, t18, -t17, 0, -t18, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t50, 0, qJ(3), 0, 0, 0, 0, 0, -t57, -t26, -t57, 0, t26, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, t53, -t64, -t4, -t64 + 0.2e1 * t80, -pkin(5) * t18 + t17 * qJ(6), 0.2e1 * t68 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t57, 0, -t15, -t16, -t15, t60, t16, -t15 * pkin(5) + t16 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t57, -t26, 0, -t57, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
