% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:04
% EndTime: 2021-01-16 02:06:07
% DurationCPUTime: 0.62s
% Computational Cost: add. (694->106), mult. (1507->200), div. (0->0), fcn. (1866->12), ass. (0->69)
t52 = sin(qJ(3));
t48 = sin(pkin(6));
t63 = t48 * sin(qJ(2));
t66 = cos(pkin(6));
t73 = cos(qJ(3));
t28 = t66 * t52 + t73 * t63;
t47 = sin(pkin(11));
t50 = cos(pkin(11));
t56 = -t52 * t63 + t66 * t73;
t14 = t47 * t28 - t50 * t56;
t81 = t14 ^ 2;
t64 = t73 * pkin(8);
t37 = t73 * qJ(4) + t64;
t62 = (-pkin(8) - qJ(4)) * t52;
t24 = t47 * t37 - t50 * t62;
t80 = t24 ^ 2;
t33 = t47 * t73 + t50 * t52;
t46 = sin(pkin(12));
t49 = cos(pkin(12));
t51 = sin(qJ(6));
t53 = cos(qJ(6));
t34 = t53 * t46 + t51 * t49;
t12 = t34 * t33;
t79 = -0.2e1 * t12;
t75 = t50 * pkin(3);
t42 = -pkin(4) - t75;
t36 = -t49 * pkin(5) + t42;
t78 = 0.2e1 * t36;
t43 = -t73 * pkin(3) - pkin(2);
t77 = 0.2e1 * t43;
t76 = t47 * pkin(3);
t39 = qJ(5) + t76;
t74 = pkin(9) + t39;
t72 = t14 * t24;
t31 = t47 * t52 - t50 * t73;
t71 = t34 * t31;
t70 = t46 * t33;
t54 = cos(qJ(2));
t69 = t48 * t54;
t68 = t49 * t33;
t21 = t31 * pkin(4) - t33 * qJ(5) + t43;
t26 = t50 * t37 + t47 * t62;
t8 = t46 * t21 + t49 * t26;
t67 = t46 ^ 2 + t49 ^ 2;
t65 = 0.2e1 * t73;
t7 = t49 * t21 - t46 * t26;
t61 = t8 * t46 + t7 * t49;
t60 = -t7 * t46 + t8 * t49;
t16 = t50 * t28 + t47 * t56;
t10 = t49 * t16 - t46 * t69;
t9 = -t46 * t16 - t49 * t69;
t59 = t10 * t49 - t9 * t46;
t58 = t10 * t46 + t9 * t49;
t57 = -t31 * t39 + t33 * t42;
t32 = t51 * t46 - t53 * t49;
t30 = t74 * t49;
t29 = t74 * t46;
t23 = t32 * t31;
t20 = -t51 * t29 + t53 * t30;
t19 = -t53 * t29 - t51 * t30;
t13 = t32 * t33;
t11 = pkin(5) * t70 + t24;
t6 = -pkin(9) * t70 + t8;
t5 = t31 * pkin(5) - pkin(9) * t68 + t7;
t4 = t53 * t10 + t51 * t9;
t3 = -t51 * t10 + t53 * t9;
t2 = t51 * t5 + t53 * t6;
t1 = t53 * t5 - t51 * t6;
t15 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 ^ 2 * t54 ^ 2 + t16 ^ 2 + t81, 0, 0, 0, t10 ^ 2 + t9 ^ 2 + t81, 0, 0, 0, 0, 0, 0, 0; 0, 0, t69, -t63, 0, 0, 0, 0, 0, t73 * t69, -t52 * t69, -t31 * t69, -t33 * t69, t14 * t33 - t16 * t31, t16 * t26 - t43 * t69 + t72, t14 * t70 + t9 * t31, -t10 * t31 + t14 * t68, -t58 * t33, t10 * t8 + t9 * t7 + t72, 0, 0, 0, 0, 0, t14 * t12 + t3 * t31, -t14 * t13 - t4 * t31; 0, 1, 0, 0, t52 ^ 2, t52 * t65, 0, 0, 0, pkin(2) * t65, -0.2e1 * pkin(2) * t52, t31 * t77, t33 * t77, 0.2e1 * t24 * t33 - 0.2e1 * t26 * t31, t26 ^ 2 + t43 ^ 2 + t80, 0.2e1 * t24 * t70 + 0.2e1 * t7 * t31, 0.2e1 * t24 * t68 - 0.2e1 * t8 * t31, -0.2e1 * t61 * t33, t7 ^ 2 + t8 ^ 2 + t80, t13 ^ 2, -t13 * t79, -0.2e1 * t13 * t31, t31 * t79, t31 ^ 2, 0.2e1 * t1 * t31 + 0.2e1 * t11 * t12, -0.2e1 * t11 * t13 - 0.2e1 * t2 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t28, -t14, -t16, 0, (-t14 * t50 + t16 * t47) * pkin(3), -t14 * t49, t14 * t46, t59, t14 * t42 + t59 * t39, 0, 0, 0, 0, 0, t14 * t32, t14 * t34; 0, 0, 0, 0, 0, 0, t52, t73, 0, -t52 * pkin(8), -t64, -t24, -t26, (-t31 * t47 - t33 * t50) * pkin(3), (-t24 * t50 + t26 * t47) * pkin(3), -t24 * t49 + t57 * t46, t24 * t46 + t57 * t49, t60, t24 * t42 + t60 * t39, -t13 * t34, -t34 * t12 + t13 * t32, t71, -t23, 0, t11 * t32 + t36 * t12 + t19 * t31, t11 * t34 - t36 * t13 - t20 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t75, -0.2e1 * t76, 0, (t47 ^ 2 + t50 ^ 2) * pkin(3) ^ 2, -0.2e1 * t42 * t49, 0.2e1 * t42 * t46, 0.2e1 * t67 * t39, t67 * t39 ^ 2 + t42 ^ 2, t34 ^ 2, -0.2e1 * t34 * t32, 0, 0, 0, t32 * t78, t34 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, t43, t49 * t31, -t46 * t31, -t67 * t33, t61, 0, 0, 0, 0, 0, -t23, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t68, 0, t24, 0, 0, 0, 0, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t46, 0, t42, 0, 0, 0, 0, 0, t32, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, t31, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
