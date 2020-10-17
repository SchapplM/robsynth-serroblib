% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR16_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:26
% EndTime: 2019-12-31 20:47:28
% DurationCPUTime: 0.56s
% Computational Cost: add. (399->92), mult. (1015->206), div. (0->0), fcn. (1109->8), ass. (0->76)
t79 = -2 * pkin(2);
t39 = cos(pkin(5));
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t38 = sin(pkin(5));
t45 = cos(qJ(2));
t68 = t38 * t45;
t18 = t39 * t41 + t44 * t68;
t78 = -0.2e1 * t18;
t77 = 0.2e1 * t38;
t43 = cos(qJ(5));
t76 = 0.2e1 * t43;
t75 = 2 * qJ(3);
t46 = -pkin(2) - pkin(8);
t42 = sin(qJ(2));
t74 = pkin(1) * t42;
t73 = pkin(1) * t45;
t30 = t38 * t42;
t48 = -qJ(3) * t42 - pkin(1);
t12 = (t46 * t45 + t48) * t38;
t25 = pkin(7) * t30;
t49 = -pkin(2) - t73;
t8 = pkin(3) * t30 + t25 + (-pkin(8) + t49) * t39;
t5 = -t41 * t12 + t44 * t8;
t3 = -pkin(4) * t30 - t5;
t40 = sin(qJ(5));
t72 = t3 * t40;
t71 = t3 * t43;
t19 = t39 * t44 - t41 * t68;
t10 = t19 * t43 + t40 * t30;
t70 = t10 * t40;
t33 = t38 ^ 2;
t69 = t33 * t45;
t67 = t39 * t45;
t66 = t40 * t18;
t65 = t40 * t41;
t64 = t40 * t43;
t63 = t40 * t44;
t62 = t41 * t46;
t61 = t43 * t18;
t60 = t43 * t41;
t31 = t43 * t44;
t59 = t43 * t46;
t58 = t44 * t10;
t57 = t44 * t18;
t56 = t44 * t41;
t55 = t44 * t46;
t21 = pkin(7) * t68 + t39 * t74;
t35 = t41 ^ 2;
t37 = t44 ^ 2;
t54 = -t35 - t37;
t53 = 0.2e1 * t30;
t52 = -0.2e1 * t56;
t51 = t41 * t30;
t50 = t46 * t30;
t32 = t39 * qJ(3);
t15 = -t32 - t21;
t11 = pkin(3) * t68 - t15;
t47 = -pkin(4) * t44 - pkin(9) * t41;
t6 = t44 * t12 + t41 * t8;
t36 = t43 ^ 2;
t34 = t40 ^ 2;
t29 = t33 * t42 ^ 2;
t24 = t44 * t30;
t22 = t41 * pkin(4) - t44 * pkin(9) + qJ(3);
t20 = pkin(1) * t67 - t25;
t17 = t49 * t39 + t25;
t16 = (-pkin(2) * t45 + t48) * t38;
t14 = t40 * t22 + t41 * t59;
t13 = t43 * t22 - t40 * t62;
t9 = t19 * t40 - t43 * t30;
t7 = t18 * pkin(4) - t19 * pkin(9) + t11;
t4 = pkin(9) * t30 + t6;
t2 = t43 * t4 + t40 * t7;
t1 = -t40 * t4 + t43 * t7;
t23 = [1, 0, 0, t29, 0.2e1 * t42 * t69, t39 * t53, t67 * t77, t39 ^ 2, 0.2e1 * pkin(1) * t69 + 0.2e1 * t20 * t39, -0.2e1 * t21 * t39 - 0.2e1 * t33 * t74, (-t15 * t45 + t17 * t42) * t77, 0.2e1 * t16 * t68 + 0.2e1 * t17 * t39, -0.2e1 * t15 * t39 - 0.2e1 * t16 * t30, t15 ^ 2 + t16 ^ 2 + t17 ^ 2, t19 ^ 2, t19 * t78, t19 * t53, t30 * t78, t29, 0.2e1 * t11 * t18 + 0.2e1 * t5 * t30, 0.2e1 * t11 * t19 - 0.2e1 * t6 * t30, t10 ^ 2, -0.2e1 * t10 * t9, 0.2e1 * t10 * t18, t9 * t78, t18 ^ 2, 0.2e1 * t1 * t18 + 0.2e1 * t3 * t9, 0.2e1 * t3 * t10 - 0.2e1 * t2 * t18; 0, 0, 0, 0, 0, t30, t68, t39, t20, -t21, (-pkin(2) * t42 + qJ(3) * t45) * t38, t25 + (t79 - t73) * t39, 0.2e1 * t32 + t21, -t17 * pkin(2) - t15 * qJ(3), t19 * t44, -t19 * t41 - t57, t24, -t51, 0, qJ(3) * t18 + t11 * t41 + t44 * t50, qJ(3) * t19 + t11 * t44 - t41 * t50, t43 * t58, (-t43 * t9 - t70) * t44, t10 * t41 + t43 * t57, -t40 * t57 - t9 * t41, t18 * t41, t1 * t41 + t13 * t18 + (-t46 * t9 + t72) * t44, -t14 * t18 - t2 * t41 + (-t10 * t46 + t71) * t44; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t79, t75, pkin(2) ^ 2 + qJ(3) ^ 2, t37, t52, 0, 0, 0, t41 * t75, t44 * t75, t36 * t37, -0.2e1 * t37 * t64, t56 * t76, t40 * t52, t35, -0.2e1 * t37 * t46 * t40 + 0.2e1 * t13 * t41, -0.2e1 * t14 * t41 - 0.2e1 * t37 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t39, 0, t17, 0, 0, 0, 0, 0, t24, -t51, 0, 0, 0, 0, 0, -t18 * t65 - t44 * t9, -t18 * t60 - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t40, t54 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t30, t5, -t6, t70, t10 * t43 - t40 * t9, t66, t61, 0, -pkin(4) * t9 - pkin(9) * t66 - t71, -pkin(4) * t10 - pkin(9) * t61 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, t55, -t62, t40 * t31, (-t34 + t36) * t44, t65, t60, 0, t47 * t40 + t43 * t55, -t40 * t55 + t47 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, 0, 0, 0, 0, t31, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t34, 0.2e1 * t64, 0, 0, 0, pkin(4) * t76, -0.2e1 * pkin(4) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t63, t41, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(9), -t43 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t23;
