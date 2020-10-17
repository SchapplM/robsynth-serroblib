% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:04:53
% EndTime: 2019-05-06 09:04:55
% DurationCPUTime: 0.70s
% Computational Cost: add. (1238->106), mult. (2327->197), div. (0->0), fcn. (2739->8), ass. (0->63)
t52 = cos(pkin(10));
t50 = sin(pkin(10));
t54 = sin(qJ(5));
t68 = t54 * t50;
t74 = cos(qJ(5));
t82 = t74 * t52 - t68;
t55 = sin(qJ(2));
t40 = (-qJ(3) - pkin(7)) * t55;
t75 = cos(qJ(2));
t64 = t75 * pkin(7);
t41 = t75 * qJ(3) + t64;
t51 = sin(pkin(9));
t53 = cos(pkin(9));
t26 = -t53 * t40 + t51 * t41;
t81 = t26 ^ 2;
t36 = t51 * t75 + t53 * t55;
t63 = t74 * t50;
t37 = t54 * t52 + t63;
t16 = t37 * t36;
t80 = -0.2e1 * t16;
t79 = -0.2e1 * t37;
t46 = -t53 * pkin(2) - pkin(3);
t39 = -t52 * pkin(4) + t46;
t78 = 0.2e1 * t39;
t34 = t51 * t55 - t53 * t75;
t47 = -t75 * pkin(2) - pkin(1);
t22 = t34 * pkin(3) - t36 * qJ(4) + t47;
t28 = t51 * t40 + t53 * t41;
t12 = t50 * t22 + t52 * t28;
t70 = t50 * t36;
t10 = -pkin(8) * t70 + t12;
t11 = t52 * t22 - t50 * t28;
t69 = t52 * t36;
t8 = t34 * pkin(4) - pkin(8) * t69 + t11;
t4 = t74 * t10 + t54 * t8;
t77 = t34 * pkin(5);
t43 = t51 * pkin(2) + qJ(4);
t76 = pkin(8) + t43;
t32 = t76 * t52;
t20 = t54 * t32 + t76 * t63;
t73 = t20 * t34;
t21 = t74 * t32 - t76 * t68;
t72 = t21 * t34;
t17 = t82 * t36;
t71 = t82 * t17;
t24 = t82 * t34;
t25 = t37 * t34;
t67 = t50 ^ 2 + t52 ^ 2;
t66 = t34 * qJ(6);
t65 = 0.2e1 * t75;
t61 = t54 * t10 - t74 * t8;
t14 = pkin(4) * t70 + t26;
t60 = pkin(5) * t82 + t37 * qJ(6);
t59 = t11 * t52 + t12 * t50;
t58 = -t11 * t50 + t12 * t52;
t57 = -t34 * t43 + t36 * t46;
t33 = t37 ^ 2;
t15 = t39 - t60;
t13 = t37 * t16;
t5 = t16 * pkin(5) - t17 * qJ(6) + t14;
t2 = t61 - t77;
t1 = t66 + t4;
t3 = [1, 0, 0, t55 ^ 2, t55 * t65, 0, 0, 0, pkin(1) * t65, -0.2e1 * pkin(1) * t55, 0.2e1 * t26 * t36 - 0.2e1 * t28 * t34, t28 ^ 2 + t47 ^ 2 + t81, 0.2e1 * t11 * t34 + 0.2e1 * t26 * t70, -0.2e1 * t12 * t34 + 0.2e1 * t26 * t69, -0.2e1 * t59 * t36, t11 ^ 2 + t12 ^ 2 + t81, t17 ^ 2, t17 * t80, 0.2e1 * t17 * t34, t34 * t80, t34 ^ 2, 0.2e1 * t14 * t16 - 0.2e1 * t34 * t61, 0.2e1 * t14 * t17 - 0.2e1 * t4 * t34, 0.2e1 * t5 * t16 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t16 + 0.2e1 * t2 * t17, 0.2e1 * t1 * t34 - 0.2e1 * t5 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t55, t75, 0, -t55 * pkin(7), -t64 (-t34 * t51 - t36 * t53) * pkin(2) (-t26 * t53 + t28 * t51) * pkin(2), -t26 * t52 + t50 * t57, t26 * t50 + t52 * t57, t58, t26 * t46 + t43 * t58, t17 * t37, -t13 + t71, t25, t24, 0, -t14 * t82 + t39 * t16 - t73, t14 * t37 + t39 * t17 - t72, t15 * t16 - t5 * t82 - t73, t1 * t82 - t21 * t16 + t20 * t17 + t2 * t37, -t15 * t17 - t5 * t37 + t72, t1 * t21 + t5 * t15 + t2 * t20; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t51 ^ 2 + t53 ^ 2) * pkin(2) ^ 2, -0.2e1 * t46 * t52, 0.2e1 * t46 * t50, 0.2e1 * t67 * t43, t43 ^ 2 * t67 + t46 ^ 2, t33, -t82 * t79, 0, 0, 0, -t82 * t78, t37 * t78, -0.2e1 * t15 * t82, 0.2e1 * t20 * t37 + 0.2e1 * t21 * t82, t15 * t79, t15 ^ 2 + t20 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t52 * t34, -t50 * t34, -t67 * t36, t59, 0, 0, 0, 0, 0, t24, -t25, t24, -t13 - t71, t25, t1 * t37 - t2 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t82 + t21 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 ^ 2 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, 0, t26, 0, 0, 0, 0, 0, t16, t17, t16, 0, -t17, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t50, 0, t46, 0, 0, 0, 0, 0, -t82, t37, -t82, 0, -t37, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t34, -t61, -t4, -t61 + 0.2e1 * t77, -pkin(5) * t17 - t16 * qJ(6), 0.2e1 * t66 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t82, 0, -t20, -t21, -t20, -pkin(5) * t37 + qJ(6) * t82, t21, -t20 * pkin(5) + t21 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t37, t82, 0, t37, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
