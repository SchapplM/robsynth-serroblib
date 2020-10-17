% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:04:46
% EndTime: 2019-05-05 06:04:47
% DurationCPUTime: 0.56s
% Computational Cost: add. (311->91), mult. (696->165), div. (0->0), fcn. (828->10), ass. (0->75)
t45 = sin(qJ(5));
t32 = pkin(5) * t45 + qJ(4);
t81 = 0.2e1 * t32;
t46 = sin(qJ(3));
t80 = -0.2e1 * t46;
t79 = 0.2e1 * t46;
t50 = cos(qJ(3));
t78 = 0.2e1 * t50;
t77 = 0.2e1 * qJ(4);
t52 = -pkin(3) - pkin(9);
t44 = sin(qJ(6));
t76 = t44 * pkin(5);
t75 = t46 * pkin(5);
t48 = cos(qJ(6));
t74 = t48 * pkin(5);
t49 = cos(qJ(5));
t57 = -t46 * qJ(4) - pkin(2);
t22 = t50 * t52 + t57;
t58 = pkin(10) * t50 - t22;
t35 = t46 * pkin(8);
t29 = pkin(4) * t46 + t35;
t69 = t45 * t29;
t7 = -t49 * t58 + t69;
t73 = t48 * t7;
t23 = t44 * t49 + t45 * t48;
t72 = t23 * t46;
t42 = sin(pkin(6));
t71 = t42 * sin(qJ(2));
t51 = cos(qJ(2));
t70 = t42 * t51;
t68 = t45 * t46;
t67 = t45 * t50;
t66 = t46 * t50;
t65 = t49 * t45;
t64 = t49 * t50;
t36 = t50 * pkin(8);
t30 = pkin(4) * t50 + t36;
t39 = t46 ^ 2;
t41 = t50 ^ 2;
t63 = t39 + t41;
t62 = qJ(4) * t50;
t61 = -0.2e1 * t66;
t60 = t46 * t70;
t59 = t50 * t70;
t25 = t49 * t29;
t6 = t45 * t58 + t25 + t75;
t1 = -t44 * t7 + t48 * t6;
t56 = -pkin(3) * t46 + t62;
t43 = cos(pkin(6));
t18 = -t43 * t50 + t46 * t71;
t19 = t43 * t46 + t50 * t71;
t55 = t18 * t46 + t19 * t50;
t54 = t46 * t52 + t62;
t40 = t49 ^ 2;
t38 = t45 ^ 2;
t34 = t49 * t52;
t33 = t49 * t46;
t28 = -pkin(3) * t50 + t57;
t27 = -pkin(10) * t49 + t34;
t26 = (-pkin(10) + t52) * t45;
t24 = -t44 * t45 + t48 * t49;
t21 = t24 * t46;
t17 = pkin(5) * t64 + t30;
t16 = t23 * t50;
t15 = t44 * t67 - t48 * t64;
t13 = t26 * t48 + t27 * t44;
t12 = -t26 * t44 + t27 * t48;
t11 = -t18 * t45 + t49 * t70;
t10 = t18 * t49 + t45 * t70;
t9 = t22 * t49 + t69;
t8 = -t22 * t45 + t25;
t4 = t10 * t44 - t11 * t48;
t3 = t10 * t48 + t11 * t44;
t2 = t44 * t6 + t73;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 ^ 2 * t51 ^ 2 + t18 ^ 2 + t19 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t70, -t71, 0, 0, 0, 0, 0, t59, -t60, t55, -t59, t60, pkin(8) * t55 - t28 * t70, 0, 0, 0, 0, 0, t10 * t46 + t19 * t64, t11 * t46 - t19 * t67, 0, 0, 0, 0, 0, -t15 * t19 + t3 * t46, -t16 * t19 - t4 * t46; 0, 1, 0, 0, t39, 0.2e1 * t66, 0, 0, 0, pkin(2) * t78, pkin(2) * t80, 0.2e1 * t63 * pkin(8), t28 * t78, t28 * t80, pkin(8) ^ 2 * t63 + t28 ^ 2, t38 * t41, 0.2e1 * t41 * t65, t45 * t61, t49 * t61, t39, 0.2e1 * t30 * t64 + 0.2e1 * t46 * t8, -0.2e1 * t30 * t67 - 0.2e1 * t46 * t9, t16 ^ 2, -0.2e1 * t16 * t15, -t16 * t79, t15 * t79, t39, 0.2e1 * t1 * t46 - 0.2e1 * t15 * t17, -0.2e1 * t16 * t17 - 0.2e1 * t2 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, t18, t19, -pkin(3) * t18 + qJ(4) * t19, 0, 0, 0, 0, 0, t19 * t45, t19 * t49, 0, 0, 0, 0, 0, t19 * t23, t19 * t24; 0, 0, 0, 0, 0, 0, t46, t50, 0, -t35, -t36, t56, t35, t36, t56 * pkin(8), -t45 * t64 (t38 - t40) * t50, t33, -t68, 0, t30 * t45 + t49 * t54, t30 * t49 - t45 * t54, -t16 * t24, t15 * t24 + t16 * t23, t21, -t72, 0, t12 * t46 - t15 * t32 + t17 * t23, -t13 * t46 - t16 * t32 + t17 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t77, pkin(3) ^ 2 + qJ(4) ^ 2, t40, -0.2e1 * t65, 0, 0, 0, t45 * t77, t49 * t77, t24 ^ 2, -0.2e1 * t24 * t23, 0, 0, 0, t23 * t81, t24 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, t35, 0, 0, 0, 0, 0, t33, -t68, 0, 0, 0, 0, 0, t21, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t64, t46, t8, -t9, 0, 0, -t16, t15, t46, t46 * t74 + t1, -t73 + (-t6 - t75) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t45, 0, t34, -t45 * t52, 0, 0, t24, -t23, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t45, 0, 0, 0, 0, 0, t24, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t74, -0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, t46, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t74, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
