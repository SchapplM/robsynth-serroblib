% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:55:47
% EndTime: 2019-05-05 09:55:50
% DurationCPUTime: 0.76s
% Computational Cost: add. (792->140), mult. (1728->227), div. (0->0), fcn. (1990->10), ass. (0->75)
t48 = sin(qJ(5));
t50 = sin(qJ(3));
t52 = cos(qJ(4));
t77 = cos(qJ(5));
t60 = t77 * t52;
t49 = sin(qJ(4));
t70 = t49 * t50;
t23 = -t48 * t70 + t50 * t60;
t86 = -0.2e1 * t23;
t61 = t77 * t49;
t28 = t48 * t52 + t61;
t85 = -0.2e1 * t28;
t39 = -t52 * pkin(4) - pkin(3);
t84 = 0.2e1 * t39;
t83 = -0.2e1 * t50;
t53 = cos(qJ(3));
t82 = 0.2e1 * t53;
t81 = -pkin(10) - pkin(9);
t80 = pkin(3) * t52;
t79 = pkin(8) * t49;
t40 = t48 * pkin(4);
t78 = t53 * pkin(5);
t31 = -t53 * pkin(3) - t50 * pkin(9) - pkin(2);
t26 = t52 * t31;
t67 = t52 * t50;
t13 = -pkin(10) * t67 + t26 + (-pkin(4) - t79) * t53;
t66 = t52 * t53;
t63 = pkin(8) * t66;
t17 = t63 + (-pkin(10) * t50 + t31) * t49;
t4 = t48 * t13 + t77 * t17;
t32 = t81 * t52;
t18 = -t48 * t32 - t81 * t61;
t76 = t18 * t53;
t71 = t48 * t49;
t19 = -t77 * t32 + t81 * t71;
t75 = t19 * t53;
t47 = cos(pkin(6));
t46 = sin(pkin(6));
t73 = t46 * sin(qJ(2));
t24 = -t47 * t53 + t50 * t73;
t74 = t24 * t28;
t72 = t46 * cos(qJ(2));
t69 = t49 * t52;
t68 = t49 * t53;
t41 = t50 * pkin(8);
t30 = pkin(4) * t70 + t41;
t65 = t53 * qJ(6);
t64 = t50 * t82;
t62 = t77 * pkin(4);
t22 = t28 * t50;
t25 = t47 * t50 + t53 * t73;
t16 = t25 * t52 - t49 * t72;
t57 = t25 * t49 + t52 * t72;
t5 = t48 * t16 + t77 * t57;
t59 = t24 * t22 + t5 * t53;
t3 = t77 * t13 - t48 * t17;
t6 = t77 * t16 - t48 * t57;
t58 = t24 * t23 + t6 * t53;
t56 = 0.2e1 * pkin(5);
t55 = 0.2e1 * qJ(6);
t45 = t53 ^ 2;
t44 = t52 ^ 2;
t43 = t50 ^ 2;
t42 = t49 ^ 2;
t37 = t62 + pkin(5);
t35 = t40 + qJ(6);
t27 = -t60 + t71;
t21 = t49 * t31 + t63;
t20 = -pkin(8) * t68 + t26;
t15 = t24 * t27;
t10 = t27 * pkin(5) - t28 * qJ(6) + t39;
t7 = t22 * pkin(5) - t23 * qJ(6) + t30;
t2 = -t3 + t78;
t1 = -t65 + t4;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, t72, -t73, 0, 0, 0, 0, 0, t53 * t72, -t50 * t72, 0, 0, 0, 0, 0, t24 * t70 + t57 * t53, t16 * t53 + t24 * t67, 0, 0, 0, 0, 0, t59, t58, t59, -t6 * t22 + t5 * t23, -t58, t6 * t1 + t5 * t2 + t24 * t7; 0, 1, 0, 0, t43, t64, 0, 0, 0, pkin(2) * t82, pkin(2) * t83, t44 * t43, -0.2e1 * t43 * t69, t66 * t83, t49 * t64, t45, -0.2e1 * t20 * t53 + 0.2e1 * t43 * t79, 0.2e1 * t43 * pkin(8) * t52 + 0.2e1 * t21 * t53, t23 ^ 2, t22 * t86, t53 * t86, t22 * t82, t45, 0.2e1 * t30 * t22 - 0.2e1 * t3 * t53, 0.2e1 * t30 * t23 + 0.2e1 * t4 * t53, 0.2e1 * t2 * t53 + 0.2e1 * t7 * t22, -0.2e1 * t1 * t22 + 0.2e1 * t2 * t23, -0.2e1 * t1 * t53 - 0.2e1 * t7 * t23, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, 0, 0, 0, 0, -t24 * t52, t24 * t49, 0, 0, 0, 0, 0, t15, t74, t15, -t6 * t27 + t5 * t28, -t74, t24 * t10 + t5 * t18 + t6 * t19; 0, 0, 0, 0, 0, 0, t50, t53, 0, -t41, -t53 * pkin(8), t49 * t67 (-t42 + t44) * t50, -t68, -t66, 0, -pkin(8) * t67 + (-pkin(3) * t50 + pkin(9) * t53) * t49, pkin(9) * t66 + (t79 - t80) * t50, t23 * t28, -t28 * t22 - t23 * t27, -t28 * t53, t27 * t53, 0, t39 * t22 + t30 * t27 + t76, t39 * t23 + t30 * t28 + t75, t10 * t22 + t7 * t27 + t76, -t1 * t27 + t18 * t23 - t19 * t22 + t2 * t28, -t10 * t23 - t7 * t28 - t75, t1 * t19 + t7 * t10 + t2 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t42, 0.2e1 * t69, 0, 0, 0, 0.2e1 * t80, -0.2e1 * pkin(3) * t49, t28 ^ 2, t27 * t85, 0, 0, 0, t27 * t84, t28 * t84, 0.2e1 * t10 * t27, 0.2e1 * t18 * t28 - 0.2e1 * t19 * t27, t10 * t85, t10 ^ 2 + t18 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t16, 0, 0, 0, 0, 0, -t5, -t6, -t5, 0, t6, t6 * t35 - t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t70, -t53, t20, -t21, 0, 0, t23, -t22, -t53, -t53 * t62 + t3, t53 * t40 - t4 (-pkin(5) - t37) * t53 + t3, -t35 * t22 - t37 * t23 (-qJ(6) - t35) * t53 + t4, t1 * t35 - t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t52, 0, -t49 * pkin(9), -t52 * pkin(9), 0, 0, t28, -t27, 0, -t18, -t19, -t18, -t35 * t27 - t37 * t28, t19, -t18 * t37 + t19 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t40, 0.2e1 * t37, 0, 0.2e1 * t35, t35 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -t5, 0, t6, -t5 * pkin(5) + t6 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, -t53, t3, -t4, t3 - 0.2e1 * t78, -pkin(5) * t23 - t22 * qJ(6), -0.2e1 * t65 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t18, -t19, -t18, -pkin(5) * t28 - t27 * qJ(6), t19, -t18 * pkin(5) + t19 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t62, -t40, t56 + t62, 0, t55 + t40, t37 * pkin(5) + t35 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, 0, t55, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t23, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
