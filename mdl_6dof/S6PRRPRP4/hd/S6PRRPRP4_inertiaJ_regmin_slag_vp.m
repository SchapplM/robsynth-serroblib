% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:48
% EndTime: 2021-01-16 03:12:51
% DurationCPUTime: 0.57s
% Computational Cost: add. (343->90), mult. (702->162), div. (0->0), fcn. (776->8), ass. (0->66)
t39 = cos(pkin(6));
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t38 = sin(pkin(6));
t66 = t38 * sin(qJ(2));
t17 = t39 * t41 + t44 * t66;
t73 = t17 ^ 2;
t40 = sin(qJ(5));
t28 = pkin(5) * t40 + qJ(4);
t72 = 0.2e1 * t28;
t71 = -0.2e1 * t41;
t70 = 0.2e1 * t44;
t69 = 0.2e1 * qJ(4);
t46 = -pkin(3) - pkin(9);
t68 = t41 * pkin(5);
t43 = cos(qJ(5));
t67 = t43 * pkin(5);
t45 = cos(qJ(2));
t65 = t38 * t45;
t31 = t41 * pkin(8);
t24 = pkin(4) * t41 + t31;
t64 = t40 * t24;
t63 = t40 * t41;
t62 = t40 * t44;
t61 = t41 * t44;
t60 = t43 * t40;
t59 = t43 * t44;
t32 = t44 * pkin(8);
t25 = t44 * pkin(4) + t32;
t35 = t41 ^ 2;
t37 = t44 ^ 2;
t58 = t35 + t37;
t57 = qJ(4) * t44;
t56 = qJ(6) * t44;
t55 = -0.2e1 * t61;
t54 = t41 * t65;
t53 = t44 * t65;
t52 = -qJ(4) * t41 - pkin(2);
t19 = t44 * t46 + t52;
t6 = -t19 * t40 + t43 * t24;
t51 = t40 * t56 + t6;
t50 = -pkin(3) * t41 + t57;
t16 = -t39 * t44 + t41 * t66;
t49 = t16 * t41 + t17 * t44;
t48 = t41 * t46 + t57;
t36 = t43 ^ 2;
t34 = t40 ^ 2;
t30 = t43 * t46;
t29 = t43 * t41;
t26 = t34 + t36;
t23 = -pkin(3) * t44 + t52;
t22 = -qJ(6) * t43 + t30;
t21 = (-qJ(6) + t46) * t40;
t15 = pkin(5) * t59 + t25;
t14 = t17 * t43;
t13 = t17 * t40;
t10 = -t16 * t40 + t43 * t65;
t9 = t16 * t43 + t40 * t65;
t8 = t21 * t40 + t22 * t43;
t7 = t19 * t43 + t64;
t5 = t64 + (t19 - t56) * t43;
t4 = t51 + t68;
t3 = t10 * t41 - t17 * t62;
t2 = t17 * t59 + t41 * t9;
t1 = -t10 * t40 + t43 * t9;
t11 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 ^ 2 * t45 ^ 2 + t16 ^ 2 + t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t9 ^ 2 + t73; 0, 0, t65, -t66, 0, 0, 0, 0, 0, t53, -t54, t49, -t53, t54, pkin(8) * t49 - t23 * t65, 0, 0, 0, 0, 0, t2, t3, t2, t3, (t10 * t43 + t40 * t9) * t44, -t10 * t5 + t15 * t17 + t4 * t9; 0, 1, 0, 0, t35, 0.2e1 * t61, 0, 0, 0, pkin(2) * t70, pkin(2) * t71, 0.2e1 * t58 * pkin(8), t23 * t70, t23 * t71, pkin(8) ^ 2 * t58 + t23 ^ 2, t34 * t37, 0.2e1 * t37 * t60, t40 * t55, t43 * t55, t35, 0.2e1 * t25 * t59 + 0.2e1 * t41 * t6, -0.2e1 * t25 * t62 - 0.2e1 * t41 * t7, 0.2e1 * t15 * t59 + 0.2e1 * t4 * t41, -0.2e1 * t15 * t62 - 0.2e1 * t41 * t5, (t4 * t40 - t43 * t5) * t70, t15 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, t16, t17, -pkin(3) * t16 + qJ(4) * t17, 0, 0, 0, 0, 0, t13, t14, t13, t14, -t1, -t10 * t21 + t17 * t28 + t22 * t9; 0, 0, 0, 0, 0, 0, t41, t44, 0, -t31, -t32, t50, t31, t32, t50 * pkin(8), -t40 * t59, (t34 - t36) * t44, t29, -t63, 0, t25 * t40 + t43 * t48, t25 * t43 - t40 * t48, t15 * t40 + t22 * t41 + t28 * t59, t15 * t43 - t21 * t41 - t28 * t62, (-t21 * t44 - t4) * t43 + (t22 * t44 - t5) * t40, t15 * t28 + t21 * t5 + t22 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t69, pkin(3) ^ 2 + qJ(4) ^ 2, t36, -0.2e1 * t60, 0, 0, 0, t40 * t69, t43 * t69, t40 * t72, t43 * t72, -0.2e1 * t8, t21 ^ 2 + t22 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, t31, 0, 0, 0, 0, 0, t29, -t63, t29, -t63, 0, t4 * t43 + t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, t9, t10, 0, t9 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t59, t41, t6, -t7, t51 + 0.2e1 * t68, -t5, pkin(5) * t62, t4 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, 0, t30, -t40 * t46, t22, -t21, -t67, t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, t43, -t40, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t62, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t11;
