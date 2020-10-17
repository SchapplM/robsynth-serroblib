% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:25:35
% EndTime: 2019-05-06 01:25:36
% DurationCPUTime: 0.61s
% Computational Cost: add. (1009->89), mult. (1904->158), div. (0->0), fcn. (2331->8), ass. (0->74)
t46 = sin(pkin(10));
t47 = cos(pkin(10));
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t28 = t50 * t46 - t52 * t47;
t29 = t52 * t46 + t50 * t47;
t49 = sin(qJ(4));
t72 = cos(qJ(4));
t21 = -t49 * t28 + t72 * t29;
t81 = -0.2e1 * t21;
t37 = -t47 * pkin(2) - pkin(1);
t22 = t28 * pkin(3) + t37;
t80 = 0.2e1 * t22;
t79 = 0.2e1 * t37;
t48 = sin(qJ(5));
t78 = t48 * pkin(5);
t77 = t49 * pkin(3);
t51 = cos(qJ(5));
t76 = t51 * pkin(5);
t75 = t51 * pkin(9);
t66 = pkin(7) + qJ(2);
t31 = t66 * t46;
t32 = t66 * t47;
t60 = -t52 * t31 - t50 * t32;
t14 = -t29 * pkin(8) + t60;
t55 = t50 * t31 - t52 * t32;
t15 = -t28 * pkin(8) - t55;
t9 = -t72 * t14 + t49 * t15;
t74 = t9 * t51;
t62 = t72 * pkin(3);
t39 = -t62 - pkin(4);
t73 = pkin(4) - t39;
t20 = t72 * t28 + t49 * t29;
t17 = t48 * t20;
t71 = t48 * t21;
t70 = t48 * t51;
t10 = t49 * t14 + t72 * t15;
t69 = t51 * t10;
t68 = t51 * t21;
t38 = pkin(9) + t77;
t67 = t51 * t38;
t65 = t46 ^ 2 + t47 ^ 2;
t44 = t48 ^ 2;
t45 = t51 ^ 2;
t64 = t44 + t45;
t41 = t51 * qJ(6);
t63 = t20 * t81;
t11 = t20 * pkin(4) - t21 * pkin(9) + t22;
t4 = -t48 * t10 + t51 * t11;
t1 = t20 * pkin(5) - t21 * t41 + t4;
t3 = t69 + (-qJ(6) * t21 + t11) * t48;
t61 = -t1 * t48 + t3 * t51;
t59 = -pkin(4) * t21 - pkin(9) * t20;
t58 = t1 * t51 + t3 * t48;
t57 = -t20 * t38 + t21 * t39;
t25 = (-qJ(6) - t38) * t48;
t26 = t41 + t67;
t56 = t51 * t25 + t48 * t26;
t34 = (-qJ(6) - pkin(9)) * t48;
t35 = t41 + t75;
t54 = t51 * t34 + t48 * t35;
t40 = -pkin(4) - t76;
t36 = 0.2e1 * t70;
t33 = t39 - t76;
t30 = t35 * t51;
t23 = t26 * t51;
t19 = t21 ^ 2;
t18 = t51 * t20;
t16 = t48 * t68;
t12 = (-t44 + t45) * t21;
t7 = t9 * t48;
t6 = pkin(5) * t71 + t9;
t5 = t48 * t11 + t69;
t2 = [1, 0, 0, 0.2e1 * pkin(1) * t47, -0.2e1 * pkin(1) * t46, 0.2e1 * t65 * qJ(2), t65 * qJ(2) ^ 2 + pkin(1) ^ 2, t29 ^ 2, -0.2e1 * t28 * t29, 0, 0, 0, t28 * t79, t29 * t79, t19, t63, 0, 0, 0, t20 * t80, t21 * t80, t45 * t19, -0.2e1 * t19 * t70, 0.2e1 * t20 * t68, t48 * t63, t20 ^ 2, 0.2e1 * t4 * t20 + 0.2e1 * t9 * t71, -0.2e1 * t5 * t20 + 0.2e1 * t9 * t68, t58 * t81, t1 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, -t47, t46, 0, -pkin(1), 0, 0, 0, 0, 0, t28, t29, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, t18, -t17, -t64 * t21, t58; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t60, t55, 0, 0, t21, -t20, 0, -t9, -t10, t16, t12, t17, t18, 0, t57 * t48 - t74, t57 * t51 + t7, -t56 * t21 + t61, t1 * t25 + t3 * t26 + t6 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t77, t44, t36, 0, 0, 0, -0.2e1 * t39 * t51, 0.2e1 * t39 * t48, -0.2e1 * t25 * t48 + 0.2e1 * t23, t25 ^ 2 + t26 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t9, -t10, t16, t12, t17, t18, 0, t59 * t48 - t74, t59 * t51 + t7, -t54 * t21 + t61, t1 * t34 + t3 * t35 + t6 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t62, -t77, t44, t36, 0, 0, 0, t73 * t51, -t73 * t48, t23 + t30 + (-t25 - t34) * t48, t25 * t34 + t26 * t35 + t33 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t44, t36, 0, 0, 0, 0.2e1 * pkin(4) * t51, -0.2e1 * pkin(4) * t48, -0.2e1 * t34 * t48 + 0.2e1 * t30, t34 ^ 2 + t35 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t71, t20, t4, -t5, -pkin(5) * t68, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * t38, -t67, -t78, t25 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * pkin(9), -t75, -t78, t34 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
