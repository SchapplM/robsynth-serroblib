% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:23
% EndTime: 2021-01-16 03:28:26
% DurationCPUTime: 0.61s
% Computational Cost: add. (625->88), mult. (1347->163), div. (0->0), fcn. (1726->12), ass. (0->69)
t44 = sin(pkin(12));
t46 = cos(pkin(12));
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t34 = t44 * t50 - t46 * t53;
t41 = -pkin(3) * t53 - pkin(2);
t26 = pkin(4) * t34 + t41;
t74 = 0.2e1 * t26;
t73 = 0.2e1 * t41;
t72 = 0.2e1 * t53;
t71 = t44 * pkin(3);
t70 = t46 * pkin(3);
t52 = cos(qJ(6));
t60 = qJ(4) + pkin(8);
t37 = t60 * t50;
t38 = t60 * t53;
t25 = -t37 * t44 + t38 * t46;
t15 = -pkin(9) * t34 + t25;
t49 = sin(qJ(5));
t24 = -t46 * t37 - t38 * t44;
t35 = t44 * t53 + t46 * t50;
t56 = -pkin(9) * t35 + t24;
t67 = cos(qJ(5));
t7 = t49 * t15 - t56 * t67;
t69 = t7 * t52;
t40 = pkin(4) + t70;
t29 = t40 * t67 - t49 * t71;
t27 = -pkin(5) - t29;
t68 = pkin(5) - t27;
t47 = cos(pkin(6));
t45 = sin(pkin(6));
t65 = t45 * sin(qJ(2));
t32 = t47 * t53 - t50 * t65;
t33 = t47 * t50 + t53 * t65;
t16 = t32 * t46 - t33 * t44;
t17 = t32 * t44 + t33 * t46;
t10 = -t16 * t67 + t49 * t17;
t66 = t10 * t52;
t54 = cos(qJ(2));
t64 = t45 * t54;
t22 = t34 * t67 + t49 * t35;
t48 = sin(qJ(6));
t19 = t48 * t22;
t23 = -t49 * t34 + t35 * t67;
t63 = t48 * t23;
t62 = t48 * t52;
t61 = t52 * t23;
t59 = -0.2e1 * t23 * t22;
t58 = -pkin(5) * t23 - pkin(10) * t22;
t30 = -t49 * t40 - t67 * t71;
t28 = pkin(10) - t30;
t57 = -t22 * t28 + t23 * t27;
t43 = t52 ^ 2;
t42 = t48 ^ 2;
t39 = 0.2e1 * t62;
t21 = t23 ^ 2;
t20 = t52 * t22;
t18 = t48 * t61;
t12 = (-t42 + t43) * t23;
t11 = t49 * t16 + t17 * t67;
t9 = t10 * t48;
t8 = t15 * t67 + t49 * t56;
t6 = pkin(5) * t22 - pkin(10) * t23 + t26;
t5 = t7 * t48;
t4 = t11 * t52 - t48 * t64;
t3 = -t11 * t48 - t52 * t64;
t2 = t48 * t6 + t52 * t8;
t1 = -t48 * t8 + t52 * t6;
t13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 ^ 2 * t54 ^ 2 + t16 ^ 2 + t17 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t64, -t65, 0, 0, 0, 0, 0, t53 * t64, -t50 * t64, -t34 * t64, -t35 * t64, -t16 * t35 - t17 * t34, t16 * t24 + t17 * t25 - t41 * t64, 0, 0, 0, 0, 0, -t22 * t64, -t23 * t64, 0, 0, 0, 0, 0, t10 * t63 + t22 * t3, t10 * t61 - t22 * t4; 0, 1, 0, 0, t50 ^ 2, t50 * t72, 0, 0, 0, pkin(2) * t72, -0.2e1 * pkin(2) * t50, t34 * t73, t35 * t73, -0.2e1 * t24 * t35 - 0.2e1 * t25 * t34, t24 ^ 2 + t25 ^ 2 + t41 ^ 2, t21, t59, 0, 0, 0, t22 * t74, t23 * t74, t43 * t21, -0.2e1 * t21 * t62, 0.2e1 * t22 * t61, t48 * t59, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t63 * t7, -0.2e1 * t2 * t22 + 0.2e1 * t61 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33, t16, -t17, 0, (t16 * t46 + t17 * t44) * pkin(3), 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t66, t9; 0, 0, 0, 0, 0, 0, t50, t53, 0, -t50 * pkin(8), -t53 * pkin(8), t24, -t25, (-t34 * t44 - t35 * t46) * pkin(3), (t24 * t46 + t25 * t44) * pkin(3), 0, 0, t23, -t22, 0, -t7, -t8, t18, t12, t19, t20, 0, t48 * t57 - t69, t52 * t57 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t70, -0.2e1 * t71, 0, (t44 ^ 2 + t46 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t29, 0.2e1 * t30, t42, t39, 0, 0, 0, -0.2e1 * t27 * t52, 0.2e1 * t27 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, t41, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t66, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t7, -t8, t18, t12, t19, t20, 0, t48 * t58 - t69, t52 * t58 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, t30, t42, t39, 0, 0, 0, t68 * t52, -t68 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t42, t39, 0, 0, 0, 0.2e1 * pkin(5) * t52, -0.2e1 * pkin(5) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t63, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t52, 0, -t48 * t28, -t52 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t52, 0, -t48 * pkin(10), -t52 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t13;
