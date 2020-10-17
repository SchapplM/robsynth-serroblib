% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:52
% EndTime: 2019-12-31 20:33:54
% DurationCPUTime: 0.46s
% Computational Cost: add. (519->86), mult. (1142->167), div. (0->0), fcn. (1311->8), ass. (0->60)
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t32 = t47 * t50 - t53 * t48;
t41 = -pkin(3) * t48 - pkin(2);
t25 = pkin(4) * t32 + t41;
t71 = 0.2e1 * t25;
t70 = 0.2e1 * t41;
t54 = cos(qJ(2));
t69 = -0.2e1 * t54;
t68 = 0.2e1 * t54;
t67 = pkin(6) * t47;
t49 = sin(qJ(5));
t66 = t49 * pkin(4);
t51 = sin(qJ(2));
t42 = t51 * pkin(6);
t52 = cos(qJ(5));
t65 = t52 * pkin(4);
t33 = t47 * t53 + t48 * t50;
t26 = t33 * t51;
t35 = -pkin(2) * t54 - qJ(3) * t51 - pkin(1);
t30 = t48 * t35;
t60 = t48 * t51;
t17 = -pkin(7) * t60 + t30 + (-pkin(3) - t67) * t54;
t62 = t54 * pkin(6);
t23 = t47 * t35 + t48 * t62;
t61 = t47 * t51;
t19 = -pkin(7) * t61 + t23;
t9 = t17 * t50 + t19 * t53;
t7 = -pkin(8) * t26 + t9;
t64 = t52 * t7;
t63 = t54 * pkin(4);
t59 = pkin(7) + qJ(3);
t34 = pkin(3) * t61 + t42;
t58 = t47 ^ 2 + t48 ^ 2;
t27 = t32 * t51;
t8 = t53 * t17 - t19 * t50;
t4 = pkin(8) * t27 - t63 + t8;
t1 = t52 * t4 - t49 * t7;
t36 = t59 * t47;
t37 = t59 * t48;
t20 = -t53 * t36 - t37 * t50;
t57 = -pkin(2) * t51 + qJ(3) * t54;
t22 = -t47 * t62 + t30;
t56 = -t22 * t47 + t23 * t48;
t21 = -t36 * t50 + t37 * t53;
t46 = t54 ^ 2;
t45 = t51 ^ 2;
t18 = pkin(4) * t26 + t34;
t16 = -t32 * t49 + t33 * t52;
t15 = t52 * t32 + t33 * t49;
t13 = -pkin(8) * t32 + t21;
t12 = -pkin(8) * t33 + t20;
t11 = -t26 * t49 - t27 * t52;
t10 = t52 * t26 - t27 * t49;
t6 = t12 * t49 + t13 * t52;
t5 = t12 * t52 - t13 * t49;
t2 = t4 * t49 + t64;
t3 = [1, 0, 0, t45, t51 * t68, 0, 0, 0, pkin(1) * t68, -0.2e1 * pkin(1) * t51, -0.2e1 * t22 * t54 + 0.2e1 * t45 * t67, 0.2e1 * pkin(6) * t45 * t48 + 0.2e1 * t23 * t54, 0.2e1 * (-t22 * t48 - t23 * t47) * t51, pkin(6) ^ 2 * t45 + t22 ^ 2 + t23 ^ 2, t27 ^ 2, 0.2e1 * t27 * t26, -t27 * t69, t26 * t68, t46, 0.2e1 * t26 * t34 - 0.2e1 * t54 * t8, -0.2e1 * t27 * t34 + 0.2e1 * t54 * t9, t11 ^ 2, -0.2e1 * t11 * t10, t11 * t69, t10 * t68, t46, -0.2e1 * t1 * t54 + 0.2e1 * t10 * t18, 0.2e1 * t11 * t18 + 0.2e1 * t2 * t54; 0, 0, 0, 0, 0, t51, t54, 0, -t42, -t62, -pkin(6) * t60 + t47 * t57, pkin(6) * t61 + t48 * t57, t56, -pkin(2) * t42 + qJ(3) * t56, -t27 * t33, -t26 * t33 + t27 * t32, -t33 * t54, t32 * t54, 0, -t20 * t54 + t26 * t41 + t32 * t34, t21 * t54 - t27 * t41 + t33 * t34, t11 * t16, -t10 * t16 - t11 * t15, -t16 * t54, t15 * t54, 0, t10 * t25 + t15 * t18 - t5 * t54, t11 * t25 + t16 * t18 + t54 * t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t48, -0.2e1 * pkin(2) * t47, 0.2e1 * t58 * qJ(3), qJ(3) ^ 2 * t58 + pkin(2) ^ 2, t33 ^ 2, -0.2e1 * t33 * t32, 0, 0, 0, t32 * t70, t33 * t70, t16 ^ 2, -0.2e1 * t16 * t15, 0, 0, 0, t15 * t71, t16 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, 0, t42, 0, 0, 0, 0, 0, t26, -t27, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, -pkin(2), 0, 0, 0, 0, 0, t32, t33, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, -t54, t8, -t9, 0, 0, t11, -t10, -t54, -t52 * t63 + t1, -t64 + (-t4 + t63) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, t20, -t21, 0, 0, t16, -t15, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, -t54, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
