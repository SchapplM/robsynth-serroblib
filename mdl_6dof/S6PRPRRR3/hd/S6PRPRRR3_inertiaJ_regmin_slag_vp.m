% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:42:23
% EndTime: 2019-05-05 00:42:25
% DurationCPUTime: 0.52s
% Computational Cost: add. (546->79), mult. (1182->142), div. (0->0), fcn. (1546->12), ass. (0->67)
t42 = sin(pkin(12));
t44 = cos(pkin(12));
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t29 = t42 * t48 - t44 * t51;
t35 = -pkin(3) * t44 - pkin(2);
t24 = pkin(4) * t29 + t35;
t75 = 0.2e1 * t24;
t74 = 0.2e1 * t35;
t47 = sin(qJ(5));
t73 = t47 * pkin(4);
t50 = cos(qJ(6));
t63 = pkin(8) + qJ(3);
t31 = t63 * t42;
t32 = t63 * t44;
t55 = t31 * t48 - t32 * t51;
t15 = -pkin(9) * t29 - t55;
t30 = t42 * t51 + t44 * t48;
t59 = -t31 * t51 - t32 * t48;
t54 = -pkin(9) * t30 + t59;
t70 = cos(qJ(5));
t6 = t15 * t47 - t54 * t70;
t72 = t6 * t50;
t60 = t70 * pkin(4);
t37 = -t60 - pkin(5);
t71 = pkin(5) - t37;
t45 = cos(pkin(6));
t43 = sin(pkin(6));
t68 = t43 * sin(qJ(2));
t25 = -t42 * t68 + t44 * t45;
t26 = t42 * t45 + t44 * t68;
t16 = t25 * t51 - t26 * t48;
t17 = t25 * t48 + t26 * t51;
t10 = -t16 * t70 + t17 * t47;
t69 = t10 * t50;
t52 = cos(qJ(2));
t67 = t43 * t52;
t22 = t29 * t70 + t30 * t47;
t46 = sin(qJ(6));
t19 = t46 * t22;
t23 = -t29 * t47 + t30 * t70;
t66 = t46 * t23;
t65 = t46 * t50;
t64 = t50 * t23;
t62 = t42 ^ 2 + t44 ^ 2;
t61 = -0.2e1 * t23 * t22;
t58 = -pkin(5) * t23 - pkin(10) * t22;
t36 = pkin(10) + t73;
t57 = -t22 * t36 + t23 * t37;
t56 = -t25 * t42 + t26 * t44;
t41 = t50 ^ 2;
t40 = t46 ^ 2;
t33 = 0.2e1 * t65;
t21 = t23 ^ 2;
t20 = t50 * t22;
t18 = t46 * t64;
t12 = (-t40 + t41) * t23;
t11 = t16 * t47 + t17 * t70;
t9 = t10 * t46;
t8 = pkin(5) * t22 - pkin(10) * t23 + t24;
t7 = t15 * t70 + t47 * t54;
t5 = t6 * t46;
t4 = t11 * t50 - t46 * t67;
t3 = -t11 * t46 - t50 * t67;
t2 = t46 * t8 + t50 * t7;
t1 = -t46 * t7 + t50 * t8;
t13 = [1, 0, 0, 0, 0, 0, 0, t43 ^ 2 * t52 ^ 2 + t25 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t67, -t68, t44 * t67, -t42 * t67, t56, pkin(2) * t67 + qJ(3) * t56, 0, 0, 0, 0, 0, -t29 * t67, -t30 * t67, 0, 0, 0, 0, 0, -t22 * t67, -t23 * t67, 0, 0, 0, 0, 0, t10 * t66 + t22 * t3, t10 * t64 - t22 * t4; 0, 1, 0, 0, 0.2e1 * pkin(2) * t44, -0.2e1 * pkin(2) * t42, 0.2e1 * t62 * qJ(3), qJ(3) ^ 2 * t62 + pkin(2) ^ 2, t30 ^ 2, -0.2e1 * t30 * t29, 0, 0, 0, t29 * t74, t30 * t74, t21, t61, 0, 0, 0, t22 * t75, t23 * t75, t41 * t21, -0.2e1 * t21 * t65, 0.2e1 * t22 * t64, t46 * t61, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t6 * t66, -0.2e1 * t2 * t22 + 0.2e1 * t6 * t64; 0, 0, 0, 0, 0, 0, 0, -t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t44, t42, 0, -pkin(2), 0, 0, 0, 0, 0, t29, t30, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, t20, -t19; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t69, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, t59, t55, 0, 0, t23, -t22, 0, -t6, -t7, t18, t12, t19, t20, 0, t46 * t57 - t72, t50 * t57 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t73, t40, t33, 0, 0, 0, -0.2e1 * t37 * t50, 0.2e1 * t37 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t69, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t6, -t7, t18, t12, t19, t20, 0, t46 * t58 - t72, t50 * t58 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t60, -t73, t40, t33, 0, 0, 0, t71 * t50, -t71 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t40, t33, 0, 0, 0, 0.2e1 * pkin(5) * t50, -0.2e1 * pkin(5) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t66, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t50, 0, -t46 * t36, -t50 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t50, 0, -t46 * pkin(10), -t50 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
